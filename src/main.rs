#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use isatty::stdout_isatty;
use std::process::ExitCode;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::PathBuf,
    str::FromStr,
};
use structopt::clap::arg_enum;

use anyhow::{bail, ensure, Context, Result};
use env_logger::Env;
use flate2::bufread::MultiGzDecoder;
use flume::{bounded, Receiver, Sender};
use fqgrep_lib::matcher::{validate_fixed_pattern, Matcher, MatcherFactory, MatcherOpts};
use fqgrep_lib::{is_fastq_path, is_gzip_path};
use gzp::BUFSIZE;
use itertools::{self, izip, Itertools};
use lazy_static::lazy_static;
use parking_lot::Mutex;
use proglog::{CountFormatterKind, ProgLog, ProgLogBuilder};
use rayon::prelude::*;
use seq_io::fastq::{self, OwnedRecord, Record};
use structopt::{clap::AppSettings, StructOpt};

/// The number of reads in a chunk
const CHUNKSIZE: usize = 5000; // * num_cpus::get();
/// The number of chunks allowed in a channel
const READER_CHANNEL_SIZE: usize = 100;
const WRITER_CHANNEL_SIZE: usize = 2000;

lazy_static! {
    /// Return the number of cpus as a String
    pub static ref NUM_CPU: String = num_cpus::get().to_string();
}

pub mod built_info {
    use lazy_static::lazy_static;
    include!(concat!(env!("OUT_DIR"), "/built.rs"));

    /// Get a software version string including
    ///   - Git commit hash
    ///   - Git dirty info (whether the repo had uncommitted changes)
    ///   - Cargo package version if no git info found
    fn get_software_version() -> String {
        let prefix = if let Some(s) = GIT_COMMIT_HASH {
            format!("{}-{}", PKG_VERSION, s[0..8].to_owned())
        } else {
            // This shouldn't happen
            PKG_VERSION.to_string()
        };
        let suffix = match GIT_DIRTY {
            Some(true) => "-dirty",
            _ => "",
        };
        format!("{}{}", prefix, suffix)
    }

    lazy_static! {
        /// Version of the software with git hash
        pub static ref VERSION: String = get_software_version();
    }
}
#[derive(Debug)]
struct FastqWriter {
    tx: Sender<Vec<OwnedRecord>>,
    lock: Mutex<()>,
}

impl FastqWriter {
    fn new(count_tx: Sender<usize>, count: bool, paired: bool, output: Option<PathBuf>) -> Self {
        // TODO: try making this unbounded
        let (tx, rx): (Sender<Vec<OwnedRecord>>, Receiver<Vec<OwnedRecord>>) =
            bounded(WRITER_CHANNEL_SIZE);

        std::thread::spawn(move || {
            let mut maybe_writer: Option<Box<dyn Write>> = {
                if count {
                    None
                } else if let Some(file_path) = output {
                    Some(Box::new(BufWriter::with_capacity(
                        BUFSIZE,
                        File::create(file_path).unwrap(),
                    )))
                } else {
                    Some(Box::new(BufWriter::with_capacity(
                        BUFSIZE,
                        std::io::stdout(),
                    )))
                }
            };

            let mut num_matches = 0;
            while let Ok(reads) = rx.recv() {
                num_matches += reads.len();
                if let Some(ref mut writer) = maybe_writer {
                    for read in reads {
                        fastq::write_to(&mut *writer, &read.head, &read.seq, &read.qual)
                            .expect("failed writing read");
                    }
                };
            }
            if paired {
                num_matches /= 2;
            }

            if count {
                std::io::stdout()
                    .write_all(format!("{}\n", num_matches).as_bytes())
                    .unwrap();
            }

            if let Some(mut writer) = maybe_writer {
                writer.flush().expect("Error flushing writer");
            };
            count_tx
                .send(num_matches)
                .expect("failed sending final count");
        });
        let lock = Mutex::new(());
        Self { tx, lock }
    }
}

fn spawn_reader(file: PathBuf, decompress: bool) -> Receiver<Vec<OwnedRecord>> {
    let (tx, rx) = bounded(READER_CHANNEL_SIZE);
    std::thread::spawn(move || {
        // Open the file or standad input
        let raw_handle = if file.as_os_str() == "-" {
            Box::new(std::io::stdin()) as Box<dyn Read>
        } else {
            let handle = File::open(&file)
                .with_context(|| format!("Error opening input: {}", file.display()))
                .unwrap();
            Box::new(handle) as Box<dyn Read>
        };
        // Wrap it in a buffer
        let buf_handle = BufReader::with_capacity(BUFSIZE, raw_handle);
        // Maybe wrap it in a decompressor
        let maybe_decoder_handle = {
            let is_gzip = is_gzip_path(&file) || (!is_fastq_path(&file) && decompress);
            if is_gzip {
                Box::new(MultiGzDecoder::new(buf_handle)) as Box<dyn Read>
            } else {
                Box::new(buf_handle) as Box<dyn Read>
            }
        };
        // Open a FASTQ reader, get an iterator over the records, and chunk them
        let fastq_reader = fastq::Reader::with_capacity(maybe_decoder_handle, BUFSIZE)
            .into_records()
            .chunks(CHUNKSIZE * num_cpus::get());
        // Iterate over the chunks
        for chunk in &fastq_reader {
            tx.send(chunk.map(|r| r.expect("Error reading")).collect())
                .expect("Error sending");
        }
    });
    rx
}

arg_enum! {
    #[derive(PartialEq, Debug)]
    enum Color {
    Never,
    Always,
    Auto,
}
}

/// The fqgrep utility searches any given input FASTQ files, selecting records whose bases match
/// one or more patterns.  By default, a pattern matches the bases in a FASTQ record if the
/// regular expression (RE) in the pattern matches the bases.  An empty expression matches every
/// line.  Each FASTQ record that matches at least one of the patterns is written to the standard
/// output.
///
/// INPUT COMPRESSION
///
/// By default, the input files are assumed to be uncompressed with the following exceptions: (1)
/// If the input files are real files and end with `.gz` or `.bgz`, they are assumed to be GZIP
/// compressed, or (2) if they end with `.fastq` or `.fq`, they are assumed to be uncompressed, or
/// (3) if the `-Z/--decompress` option is specified then any unrecongized inputs (including
/// standard input) are assumed to be GZIP compressed.
///
/// THREADS
///
/// The `--threads` option controls the number of threads used to _search_ the reads.
/// Independently, for single end reads or interleaved paired end reads, a single thread will be
/// used to read each input FASTQ.  For paired end reads across pairs of FASTQs, two threads will
/// be used to read the FASTQs for each end of a pair.  Finally, a single thread will be created
/// for the writer.
///
/// EXIT STATUS
///
/// The fqgrep utility exits with one of the following values:
/// 0 if one or more lines were selected, 1 if no lines were selected, and >1 if an error occurred.
///
#[derive(StructOpt, Debug)]
#[structopt(
    name = "fqgrep", 
    global_setting(AppSettings::ColoredHelp),
    global_setting(AppSettings::DeriveDisplayOrder),
    version = built_info::VERSION.as_str())
 ]
#[allow(clippy::struct_excessive_bools)]
struct Opts {
    /// The number of threads to use for matching reads against pattern.  See the full usage for
    /// threads specific to reading and writing.
    #[structopt(long, short = "t", default_value = NUM_CPU.as_str())]
    threads: usize,

    // TODO: support GREP_COLOR(S) (or FQGREP_COLOR(S)) environment variables
    /// Mark up the matching text.  The possible values of when are “never”, “always” and “auto”.
    #[structopt(long = "color", default_value = "never")]
    color: Color,

    /// Only a count of selected lines is written to standard output.
    #[structopt(long, short = "c")]
    count: bool,

    /// Specify a pattern used during the search of the input: an input line is selected if it
    /// matches any of the specified patterns.  This option is most useful when multiple `-e`
    /// options are used to specify multiple patterns.
    #[structopt(long, short = "e")]
    regexp: Vec<String>,

    /// Interpret pattern as a set of fixed strings
    #[structopt(long, short = "F")]
    fixed_strings: bool,

    /// Read one or more newline separated patterns from file.  Empty pattern lines match every
    /// input line.  Newlines are not considered part of a pattern.  If file is empty, nothing
    /// is matched.
    #[structopt(long, short = "f")]
    file: Option<PathBuf>,

    /// Selected lines are those not matching any of the specified patterns
    #[structopt(short = "v")]
    invert_match: bool,

    /// Assume all unrecognized inputs are GZIP compressed.
    #[structopt(short = "Z", long)]
    decompress: bool,

    /// Treat the input files as paired.  The number of input files must be a multiple of two,
    /// with the first file being R1, second R2, third R1, fourth R2, and so on.  If the pattern
    /// matches either R1 or R2, then both R1 and R2 will be output (interleaved).  If the input
    /// is standard input, then treat the input as interlaved paired end reads.
    #[structopt(long)]
    paired: bool,

    /// Search the reverse complement for matches.
    #[structopt(long)]
    reverse_complement: bool,

    /// Write progress information
    #[structopt(long)]
    progress: bool,

    /// The first argument is the pattern to match, with the remaining arguments containing the
    /// files to match.  If `-e` is given, then all the arguments are files to match.
    /// Use standard input if either no files are given or `-` is given.
    ///
    args: Vec<String>,

    /// Hidden option to capture stdout for testing
    ///
    #[structopt(long, hidden = true)]
    output: Option<PathBuf>,
}

fn read_patterns(file: &PathBuf) -> Result<Vec<String>> {
    let f = File::open(&file).expect("error in opening file");
    let r = BufReader::new(f);
    let mut v = Vec::new();
    for result in r.lines() {
        v.push(result?);
    }
    Ok(v)
}

fn main() -> ExitCode {
    // Set the exit code:
    // - exit code 0 if there were matches
    // - exit code 1 if there were no matches
    // - exit code 2 if fqgrep returned an error
    // - exit code 101 if fqgrep panicked
    let outer = std::panic::catch_unwind(fqgrep);
    match outer {
        Err(_) => {
            eprintln!("Error: fqgrep panicked.  Please report this as a bug!");
            ExitCode::from(101)
        }
        Ok(inner) => match inner {
            Ok(0) => ExitCode::from(1),
            Ok(_) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("Error: {}", e);
                ExitCode::from(2)
            }
        },
    }
}

fn fqgrep() -> Result<usize> {
    // define instance of opts from setup() -- relies on env var
    let mut opts = setup();
    fqgrep_from_opts(&mut opts)
}

#[allow(clippy::too_many_lines)]
fn fqgrep_from_opts(opts: &mut Opts) -> Result<usize> {

    // Add patterns from a file if given
    if let Some(file) = &opts.file {
        for pattern in read_patterns(file)? {
            opts.regexp.push(pattern);
        }
    }

    // Inspect the positional arguments to extract a fixed pattern
    let (pattern, mut files): (Option<String>, Vec<PathBuf>) = {
        let (pattern, file_strings): (Option<String>, Vec<String>) = if opts.regexp.is_empty() {
            // No patterns given by -e, so assume the first positional argument is the pattern and
            // the rest are files
            ensure!(
                !opts.args.is_empty(),
                "Pattern must be given with -e or as the first positional argument "
            );
            let files = opts.args.iter().skip(1).cloned().collect();
            (Some(opts.args[0].clone()), files)
        } else {
            // Patterns given by -e, so assume all positional arguments are files
            (None, opts.args.clone())
        };

        // Convert file strings into paths
        let files = file_strings
            .iter()
            .map(|s| {
                PathBuf::from_str(s)
                    .with_context(|| format!("Cannot create a path from: {}", s))
                    .unwrap()
            })
            .collect();
        (pattern, files)
    };

    // Ensure that if multiple files are given, its a multiple of two.
    if opts.paired {
        ensure!(
            files.len() <= 1 || files.len() % 2 == 0,
            "Input files must be a multiple of two, or either a single file or standard input (assume interleaved) with --paired"
        );
    }

    // Validate the fixed string pattern, if fixed-strings are specified
    if opts.fixed_strings {
        if let Some(pattern) = &pattern {
            validate_fixed_pattern(pattern)?;
        } else if !opts.regexp.is_empty() {
            for pattern in &opts.regexp {
                validate_fixed_pattern(pattern)?;
            }
        } else {
            bail!("A pattern must be given as a positional argument or with -e/--regexp")
        }
    }

    // Set up a progress logger if desired
    let progress_logger = if opts.progress {
        Some(
            ProgLogBuilder::new()
                .name("fqgrep-progress")
                .noun("reads")
                .verb("Searched")
                .unit(50_000_000)
                .count_formatter(CountFormatterKind::Comma)
                .build(),
        )
    } else {
        None
    };

    // Build the common pattern matching options
    let match_opts = MatcherOpts {
        invert_match: opts.invert_match,
        reverse_complement: opts.reverse_complement,
        color: opts.color == Color::Always || (opts.color == Color::Auto && stdout_isatty()),
    };

    // The matcher used in the primary search
    let matcher: Box<dyn Matcher + Sync + Send> =
        MatcherFactory::new_matcher(&pattern, opts.fixed_strings, &opts.regexp, match_opts);

    // The thread pool from which threads are spanwed.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.threads)
        .build()
        .unwrap();

    // Sender and receiver channels for the final count of matching records
    let (count_tx, count_rx): (Sender<usize>, Receiver<usize>) = bounded(1);

    // The writer of final counts or matching records
    let writer = FastqWriter::new(count_tx, opts.count, opts.paired, opts.output.clone());

    // The main loop
    pool.install(|| {
        // If no files, use "-" to signify standard input.
        if files.is_empty() {
            // read from standard input
            files.push(PathBuf::from_str("-").unwrap());
        }
        if opts.paired {
            // Either an interleaved paired end FASTQ, or pairs of FASTQs
            if files.len() == 1 {
                // Interleaved paired end FASTQ
                // The channel FASTQ record chunks are received after being read in
                let rx = spawn_reader(files[0].clone(), opts.decompress);
                for reads in izip!(rx.iter()) {
                    let paired_reads = reads
                        .into_iter()
                        .tuples::<(OwnedRecord, OwnedRecord)>()
                        .collect_vec();
                    process_paired_reads(paired_reads, &matcher, &writer, &progress_logger);
                }
            } else {
                // Pairs of FASTQ files
                for file_pairs in files.chunks_exact(2) {
                    // The channels for R1 and R2 with FASTQ record chunks that are received after being read in
                    let rx1 = spawn_reader(file_pairs[0].clone(), opts.decompress);
                    let rx2 = spawn_reader(file_pairs[1].clone(), opts.decompress);
                    for (reads1, reads2) in izip!(rx1.iter(), rx2.iter()) {
                        let paired_reads = reads1.into_iter().zip(reads2.into_iter()).collect_vec();
                        process_paired_reads(paired_reads, &matcher, &writer, &progress_logger);
                    }
                }
            }
        } else {
            // Process one FASTQ at a time
            for file in files {
                // The channel FASTQ record chunks are received after being read in
                let rx = spawn_reader(file.clone(), opts.decompress);
                for reads in rx.iter() {
                    // Get the matched reads
                    let matched_reads: Vec<OwnedRecord> = reads
                        .into_par_iter()
                        .map(|mut read| -> Option<OwnedRecord> {
                            if let Some(progress) = &progress_logger {
                                progress.record();
                            }
                            if matcher.read_match(&mut read) {
                                Some(read)
                            } else {
                                None
                            }
                        })
                        .flatten()
                        .collect();

                    let _lock = writer.lock.lock();
                    writer.tx.send(matched_reads).expect("Failed to send read");
                }
            }
        }
    });

    drop(writer); // so count_tx.send will execute
                  // Get the final count of records matched

    match count_rx.recv() {
        Ok(count) => Ok(count),
        Err(error) => Err(error).with_context(|| "failed receive final match counts"),
    }
}

/// Process a chunk of paired end records in parallel.
#[allow(clippy::borrowed_box)] // FIXME: remove me later and solve
fn process_paired_reads(
    reads: Vec<(OwnedRecord, OwnedRecord)>,
    matcher: &Box<dyn Matcher + Sync + Send>,
    writer: &FastqWriter,
    progress_logger: &Option<ProgLog>,
) {
    reads
        .into_par_iter()
        .map(|(mut read1, mut read2)| {
            if let Some(progress) = progress_logger {
                progress.record();
                progress.record();
            }
            assert!(
                read1.head() == read2.head(),
                "Mismatching read pair!  R1: {} R2: {}",
                std::str::from_utf8(read1.head()).unwrap(),
                std::str::from_utf8(read2.head()).unwrap()
            );
            // NB: if the output is to be colored, always call read_match on read2, regardless of
            // whether or not read1 had a match, so that read2 is always colored.  If the output
            // isn't to be colored, only search for a match in read2 if read1 does not have a match
            let match1 = matcher.read_match(&mut read1);
            let match2 = (!matcher.opts().color && match1) || matcher.read_match(&mut read2);
            if match1 || match2 {
                Some((read1, read2))
            } else {
                None
            }
        })
        .flatten()
        .fold(std::vec::Vec::new, |mut matched_reads, (r1, r2)| {
            matched_reads.push(r1);
            matched_reads.push(r2);
            matched_reads
        })
        .for_each(|matched_read| {
            let _lock = writer.lock.lock();
            writer.tx.send(matched_read).expect("Failed to send read");
        });
}

/// Parse args and set up logging / tracing
fn setup() -> Opts {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    Opts::from_args()
}

// Tests
#[cfg(test)]
pub mod tests {
    use crate::*;
    use fgoxide::io::Io;
    use tempfile::TempDir;

    /// Returns the path(s)  of type `Vec<String>` to the fastq(s) written from the provided sequence
    ///
    /// # Arguments
    ///
    /// * `temp_dir` - A temp directory that must be created in the actual test fucntion
    /// * `sequences` - A &Vec<Vec<&str>> where the outer Vec contains the the reads for a given FASTQ file
    ///
    /// # Examples
    /// if - sequences = vec![vec!["AAGTCTGAATCCATGGAAAGCTATTG", "GGGTCTGAATCCATGGAAAGCTATTG"], vec!["AAGTCTGAATCCATGGAAAGCTATTG", "GGGTCTGAATCCATGGAAAGCTATTG"]]
    /// write_fastq() will return Vec</temp/path/to/first.fa, /temp/path/to/second.fq> where 'paths' are Strings.
    fn write_fastq(temp_dir: &TempDir, sequences_per_fastq: &Vec<Vec<&str>>) -> Vec<String> {
        // Set io (from fgoxide())
        let io = Io::default();

        // Initialize a Vec<String>
        let mut fastq_paths = Vec::new();

        // Iterate through each FASTQ file we are to build
        for (fastq_index, fastq_sequences) in sequences.iter().enumerate() {
            let name = format!("sample_{i}.fq");
            let fastq_path = temp_dir.path().join(name);
            let mut lines = Vec::new();

            // Second loop through &str in &Vec<Vec<&str>>
            for (num, seq) in s.iter().enumerate() {
                // Seq ID, Sequence, Sep, and Qual
                lines.push(format!("@{num}"));
                lines.push(seq.to_string());
                lines.push("+".to_string());
                lines.push((String::from("-")).repeat(seq.len()));

                io.write_lines(&f1, &lines).unwrap();
            }
            // Convert PathBuf to String - Opts expects Vec<String>
            // as_string() is not a method of PathBuf
            let path_as_string = f1.as_path().display().to_string();
            path_vec.push(path_as_string);
        }
        path_vec
    }

    /// Returns a path (PathBuf) to a file with patterns to read from
    ///
    /// # Arguments
    ///
    /// * `temp_dir` - A temp directory that must be created in the actual test fucntion
    /// * `pattern` - One or more patterns
    ///
    /// # Examples
    /// if - pattern = &Vec<String::from("^A"), String::from("^G")>
    /// write_pattern() will return '/temp/path/to/pattern/txt' where pattern.txt has two lines ^A/n^G
    fn write_pattern(temp_dir: &TempDir, pattern: &Vec<String>) -> PathBuf {
        // Set io
        let io = Io::default();
        // File name and path
        let name = String::from("pattern.txt");
        let f1 = temp_dir.path().join(name);
        // Simply pass pattern to io.write_lines()
        io.write_lines(&f1, &*pattern).unwrap();
        // Return
        f1
    }


    /// Builds the command line options used for testing.
    /// Builds an instance of StructOpts to be passed to fqgrep_from_opts().  This will also write the
    /// FASTQ(s) (via write_fastq) and optionally writes the search pattern to a file (via write_pattern) if pattern_from_file is true.  For the latter, opts.regex and opts.file will set appropriately.
    ///
    /// # Arguments
    ///
    /// * `temp_dir` - A temp directory that must be created in the actual test fucntion
    /// * `seqs` - A &Vec<Vec<&str>> in which the number of outer Vec indicate the number of fastq files, the inner Vecs indicate the number of records in each fastq (passed to write_fastq())
    /// * `regexp` - One or more patterns to search for
    /// * `pattern_from_file` - True to write the pattern(s) to file and set opts.file, otherwise use opts.regexp
    /// * `output_path` - Optionally the path to the output, or None to output to standard output
    ///
    /// # Examples
    /// let seqs = vec![vec!["AAGTCTGAATCCATGGAAAGCTATTG", "GGGTCTGAATCCATGGAAAGCTATTG"], vec!["AAGTCTGAATCCATGGAAAGCTATTG", "GGGTCTGAATCCATGGAAAGCTATTG"]];
    /// let pattern = Vec<String::from("AGTG")>;
    /// let dir = = TempDir::new().unwrap();
    /// let ex_opts = call_opts(dir, seqs, pattern, true)
    /// ex_opts will be an instance of StructOpts where ex_opts.file is a PathBuf for a file with one line 'AGTG' and ex_opts.args is Vec<String> with two fastq paths
    fn build_opts(
        dir: &TempDir,
        seqs: &Vec<Vec<&str>>,
        regexp: &Vec<String>,
        pattern_from_file: bool,
        output: Option<PathBuf>,
    ) -> Opts {
        // write fastq
        let fq_path = write_fastq(&dir, &seqs);

        // define pattern as from file or from string (opts.file or opts.regexp)
        if pattern_from_file {
            // Write pattern.txt from provided patterns
            let pattern_path: PathBuf = write_pattern(&dir, &regexp);
            let regex: Vec<String> = vec![];
            // Write instance of opts
            let file_test_opts = write_opts(Some(pattern_path), &fq_path, &regex, output);
            file_test_opts
        } else {
            // Call write_opts
            let nofile_test_opts: Opts = write_opts(None, &fq_path, &regexp, output);
            nofile_test_opts
        }
    }

    /// Returns an instance of StructOpts within call_opts().
    ///
    /// # Arguments
    ///
    /// * `pattern_file_path` - Either None or PathBuf to 'pattern.txt' created in write_pattern()
    /// * `fastq_file_path` - &Vec<String> created in write_fastq()
    /// * `regexp_str` - Either an empty Vec<String> or Vec<String> where String = reg expressions
    /// * `output_path` - Either None or PathBuf to 'output.fq'
    ///
    fn write_opts(
        pattern_file_path: Option<PathBuf>,
        fastq_file_path: &Vec<String>,
        regex_str: &Vec<String>,
        output_path: Option<PathBuf>,
    ) -> Opts {
        // Takes a PathBuf as Opts.file, &Vec<String> for Opts.args, and &Vec<String> for Opts.regexp
        let return_opts = Opts {
            // todo count must be false when in combination with an output path
            threads: 4,
            color: Color::Never,
            count: false,
            regexp: regex_str.to_vec(),
            fixed_strings: false,
            file: pattern_file_path,
            invert_match: false,
            decompress: false,
            paired: false,
            reverse_complement: false,
            progress: true,
            args: fastq_file_path.to_vec(),
            output: output_path,
        };
        return_opts
    }

    /// Tests a single fq file for a seq starting with either A or G
    ///
    /// Workflow
    /// Define a temp dir
    /// Define the sequences to search in
    /// Define the pattern to seach for
    /// Call call_opts() and pass resulting StructOpts instance to fqgrep_from_opts()
    /// Assert that fqgrep() returns two matches
    ///  
    #[test]
    fn test_single_fq() {
        let dir = TempDir::new().unwrap();
        let seqs = vec![vec![
            "AAGTCTGAATCCATGGAAAGCTATTG",
            "GGGTCTGAATCCATGGAAAGCTATTG",
        ]];
        let test_pattern = vec![String::from("^A"), String::from("^G")];
        let mut opts_testcase = call_opts(&dir, &seqs, &test_pattern, false, None);
        // TODO: Automatically set count as t/f depending on opts.output
        opts_testcase.count = true;
        let result = fqgrep_from_opts(&mut opts_testcase);
        assert_eq!(result.unwrap(), 2)
    }

    /// Tests two files (paired reads) for seqs starting with either A or G
    ///
    /// Workflow
    /// Define a temp dir
    /// Define the sequences to search in
    /// Define the pattern to seach for
    /// Call call_opts() and pass resulting StructOpts instance to fqgrep_from_opts()
    /// Set opts.testcase.paired to true
    /// Assert that fqgrep() returns two matches
    ///
    #[test]
    fn test_multiple_fq() {
        let dir = TempDir::new().unwrap();
        let seqs = vec![
            vec!["AAGTCTGAATCCATGGAAAGCTATTG", "GGGTCTGAATCCATGGAAAGCTATTG"],
            vec!["AAGTCTGAATCCATGGAAAGCTATTG", "GGGTCTGAATCCATGGAAAGCTATTG"],
        ];
        let test_pattern = vec![String::from("^A"), String::from("^G")];
        let mut opts_testcase = call_opts(&dir, &seqs, &test_pattern, true, None);
        opts_testcase.count = true;
        opts_testcase.paired = true;
        let result = fqgrep_from_opts(&mut opts_testcase);
        assert_eq!(result.unwrap(), 2)
    }

    /// Tests two fastqs for a seq that starts with A and check the output
    ///
    /// Workflow
    /// Set io as io from fgoxide::io::Io
    /// Define a temp dir and add an output path to the dir
    /// Clone the output path so it can be read from
    ///
    /// Define the sequences to search in
    /// Define the pattern to seach for
    /// Call call_opts() and pass resulting StructOpts instance to fqgrep_from_opts()
    /// Assert that fqgrep() returns two matches
    ///
    #[test]
    fn test_get_output() {
        let io = Io::default();
        let dir: TempDir = TempDir::new().unwrap();
        let out_path = dir.path().join(String::from("output.fq"));
        // Clone path to read from below
        // TODO: refactor
        let x = &out_path.clone();
        let seqs = vec![
            vec!["AAGTCTGAATCCATGGAAAGCTATTG"],
            vec!["GGGTCTGAATCCATGGAAAGCTATTG"],
        ];
        let test_pattern = vec![String::from("^A")];
        let mut opts_testcase = call_opts(&dir, &seqs, &test_pattern, true, Some(out_path));
        // run fqgrep
        let _result = fqgrep_from_opts(&mut opts_testcase);

        // read seqs back from the output path
        let mut lines = io.read_lines(&x).unwrap();
        // Get only the seqs (every 4th line skipping the first line)
        let mut count = -2;
        lines.retain(|_| {
            count += 1;
            return count % 4 == 0;
        });
        // make expected
        let expected = vec!["AAGTCTGAATCCATGGAAAGCTATTG"];
        // compare
        assert_eq!(lines, expected);
    }
}
