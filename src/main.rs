#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use anyhow::{Context, Result, bail, ensure};
use clap::builder::styling;
use clap::{ColorChoice, Parser as ClapParser, ValueEnum};
use env_logger::Env;
use flate2::bufread::MultiGzDecoder;
use fqgrep_lib::matcher::{Matcher, MatcherFactory, MatcherOpts, validate_fixed_pattern};
use fqgrep_lib::seq_io::{
    InterleavedFastqReader, PairedFastqReader, parallel_interleaved_fastq, parallel_paired_fastq,
};
use fqgrep_lib::{is_fastq_path, is_gzip_path};
use gzp::BUFSIZE;
use isatty::stdout_isatty;
use lazy_static::lazy_static;
use proglog::{CountFormatterKind, ProgLog, ProgLogBuilder};
use seq_io::fastq::{self, Record, RefRecord};
use seq_io::parallel::parallel_fastq;
use std::process::ExitCode;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::PathBuf,
    str::FromStr,
};

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

fn spawn_reader(file: PathBuf, decompress: bool) -> fastq::Reader<Box<dyn std::io::Read + Send>> {
    // Open the file or standad input
    let raw_handle = if file.as_os_str() == "-" {
        Box::new(std::io::stdin()) as Box<dyn Read + Send>
    } else {
        let handle = File::open(&file)
            .with_context(|| format!("Error opening input: {}", file.display()))
            .unwrap();
        Box::new(handle) as Box<dyn Read + Send>
    };
    // Wrap it in a buffer
    let buf_handle: BufReader<Box<dyn std::io::Read + Send>> =
        BufReader::with_capacity(BUFSIZE, raw_handle);
    // Maybe wrap it in a decompressor
    let maybe_decoder_handle: Box<dyn std::io::Read + Send> = {
        let is_gzip = is_gzip_path(&file) || (!is_fastq_path(&file) && decompress);
        if is_gzip {
            Box::new(MultiGzDecoder::new(buf_handle)) as Box<dyn std::io::Read + Send>
        } else {
            Box::new(buf_handle) as Box<dyn std::io::Read + Send>
        }
    };
    // Open a FASTQ reader
    fastq::Reader::with_capacity(maybe_decoder_handle, BUFSIZE)
}

#[derive(ValueEnum, PartialEq, Debug, Clone)]
enum Color {
    Never,
    Always,
    Auto,
}

/// The style for the usage
const STYLES: styling::Styles = styling::Styles::styled()
    .header(styling::AnsiColor::Yellow.on_default().bold())
    .usage(styling::AnsiColor::Yellow.on_default().bold())
    .literal(styling::AnsiColor::Blue.on_default().bold())
    .placeholder(styling::AnsiColor::Cyan.on_default());

/// # OVERVIEW
///
/// The fqgrep utility searches any given input FASTQ files, selecting records whose bases match
/// one or more patterns.  By default, a pattern matches the bases in a FASTQ record if the
/// regular expression (RE) in the pattern matches the bases.  An empty expression matches every
/// line.  Each FASTQ record that matches at least one of the patterns is written to the standard
/// output.
///
/// # INPUT COMPRESSION
///
/// By default, the input files are assumed to be uncompressed with the following exceptions: (1)
/// If the input files are real files and end with `.gz` or `.bgz`, they are assumed to be GZIP
/// compressed, or (2) if they end with `.fastq` or `.fq`, they are assumed to be uncompressed, or
/// (3) if the `-Z/--decompress` option is specified then any unrecongized inputs (including
/// standard input) are assumed to be GZIP compressed.
///
/// # THREADS
///
/// The `--threads` option controls the number of threads used to _search_ the reads.
/// Independently, for single end reads or interleaved paired end reads, a single thread will be
/// used to read each input FASTQ.  For paired end reads across pairs of FASTQs, two threads will
/// be used to read the FASTQs for each end of a pair.  Finally, a single thread will be created
/// for the writer.
///
/// # EXIT STATUS
///
/// The fqgrep utility exits with one of the following values:
/// 0 if one or more lines were selected, 1 if no lines were selected, and >1 if an error occurred.
///
#[derive(ClapParser, Debug, Clone)]
#[clap(
    name = "fqgrep",
    color = ColorChoice::Auto,
    styles = STYLES,
    version = built_info::VERSION.as_str())
 ]
#[allow(clippy::struct_excessive_bools)]
struct Opts {
    /// The number of threads to use for matching reads against pattern.  See the full usage for
    /// threads specific to reading and writing.
    #[clap(long, short = 't', default_value = NUM_CPU.as_str())]
    threads: usize,

    // TODO: support GREP_COLOR(S) (or FQGREP_COLOR(S)) environment variables
    /// Mark up the matching text.  The possible values of when are “never”, “always” and “auto”.
    #[clap(long = "color", default_value = "never")]
    color: Color,

    /// Only a count of selected lines is written to standard output.
    #[clap(long, short = 'c')]
    count: bool,

    /// Specify a pattern used during the search of the input: an input line is selected if it
    /// matches any of the specified patterns.  This option is most useful when multiple `-e`
    /// options are used to specify multiple patterns.
    #[clap(long, short = 'e')]
    regexp: Vec<String>,

    /// Interpret pattern as a set of fixed strings
    #[clap(long, short = 'F')]
    fixed_strings: bool,

    /// Read one or more newline separated patterns from file.  Empty pattern lines match every
    /// input line.  Newlines are not considered part of a pattern.  If file is empty, nothing
    /// is matched.
    #[clap(long, short = 'f')]
    file: Option<PathBuf>,

    /// Selected lines are those not matching any of the specified patterns
    #[clap(short = 'v')]
    invert_match: bool,

    /// Assume all unrecognized inputs are GZIP compressed.
    #[clap(short = 'Z', long)]
    decompress: bool,

    /// Treat the input files as paired.  The number of input files must be a multiple of two,
    /// with the first file being R1, second R2, third R1, fourth R2, and so on.  If the pattern
    /// matches either R1 or R2, then both R1 and R2 will be output (interleaved).  If the input
    /// is standard input, then treat the input as interlaved paired end reads.
    #[clap(long)]
    paired: bool,

    /// Search the reverse complement for matches.
    #[clap(long)]
    reverse_complement: bool,

    /// Write progress information
    #[clap(long)]
    progress: bool,

    /// The first argument is the pattern to match, with the remaining arguments containing the
    /// files to match.  If `-e` is given, then all the arguments are files to match.
    /// Use standard input if either no files are given or `-` is given.
    ///
    args: Vec<String>,

    /// Hidden option to capture stdout for testing
    ///
    #[clap(long, hide_short_help = true, hide_long_help = true)]
    output: Option<PathBuf>,
}

fn read_patterns(file: &PathBuf) -> Result<Vec<String>> {
    let f = File::open(file).expect("error in opening file");
    let r = BufReader::new(f);
    let mut v = Vec::new();
    for result in r.lines() {
        v.push(result?);
    }
    Ok(v)
}

/// Runs fqgrep and converts Option<u8> to ExitCode
///
/// Set the exit code:
/// - exit code SUCCESS if there were matches
/// - exit code 1 if there were no matches
/// - exit code 2 if fqgrep returned an error
/// - exit code 101 if panic
fn main() -> ExitCode {
    // Receives u8 from
    if let Some(fqgrep_output) = fqgrep(&setup()) {
        ExitCode::from(fqgrep_output)
    } else {
        ExitCode::SUCCESS
    }
}

/// Runs fqgrep_from_opts and returns None upon success and an error number if error or zero matches
///
/// - None if there were matches
/// - 1 if there were no matches
/// - 2 if fqgrep_from_opts returned an error
/// - 101 if fqgrep_from_opts panicked
fn fqgrep(opts: &Opts) -> Option<u8> {
    let outer = std::panic::catch_unwind(|| fqgrep_from_opts(opts));
    match outer {
        Err(_) => {
            eprintln!("Error: fqgrep panicked.  Please report this as a bug!");
            Some(101)
        }
        Ok(inner) => match inner {
            Ok(0) => Some(1),
            Ok(_) => None,
            Err(e) => {
                eprintln!("Error: {}", e);
                Some(2)
            }
        },
    }
}

#[allow(clippy::too_many_lines)]
fn fqgrep_from_opts(opts: &Opts) -> Result<usize> {
    let mut opts = opts.clone();
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
    };

    let color = opts.color == Color::Always || (opts.color == Color::Auto && stdout_isatty());

    let mut maybe_writer: Option<Box<dyn Write>> = {
        if opts.count {
            None
        } else if let Some(file_path) = opts.output {
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
    let mut num_matches = 0usize;

    // The matcher used in the primary search
    let matcher: Box<dyn Matcher + Sync + Send> =
        MatcherFactory::new_matcher(&pattern, opts.fixed_strings, &opts.regexp, match_opts);

    // The main loop
    // If no files, use "-" to signify standard input.
    if files.is_empty() {
        // read from standard input
        files.push(PathBuf::from_str("-").unwrap());
    }

    if opts.paired {
        // Either an interleaved paired end FASTQ, or pairs of FASTQs
        if files.len() == 1 {
            // Interleaved paired end FASTQ
            parallel_interleaved_fastq(
                InterleavedFastqReader {
                    reader: spawn_reader(files.first().unwrap().clone(), opts.decompress),
                },
                opts.threads as u32,
                opts.threads,
                |(read1, read2), found| {
                    *found = process_paired_reads(&read1, &read2, &matcher, &progress_logger);
                },
                |(read1, read2), found| {
                    if *found > 0 {
                        num_matches += 1;
                        if let Some(writer) = &mut maybe_writer {
                            if color {
                                let found1 = *found == 1;
                                let found2 = *found == 2 || matcher.read_match(&read2);
                                let mut read1 = read1.to_owned_record();
                                let mut read2 = read2.to_owned_record();
                                matcher.color(&mut read1, found1);
                                matcher.color(&mut read2, found2);
                                fastq::write_to(
                                    &mut *writer,
                                    read1.head(),
                                    read1.seq(),
                                    read1.qual(),
                                )
                                .unwrap();
                                fastq::write_to(
                                    &mut *writer,
                                    read2.head(),
                                    read2.seq(),
                                    read2.qual(),
                                )
                                .unwrap();
                            } else {
                                fastq::write_to(
                                    &mut *writer,
                                    read1.head(),
                                    read1.seq(),
                                    read1.qual(),
                                )
                                .unwrap();
                                fastq::write_to(
                                    &mut *writer,
                                    read2.head(),
                                    read2.seq(),
                                    read2.qual(),
                                )
                                .unwrap();
                            }
                        }
                    }
                    None::<()>
                },
            )
            .unwrap();
        } else {
            // // Pairs of FASTQ files
            for file_pairs in files.chunks_exact(2) {
                let reader1: fastq::Reader<Box<dyn Read + Send>> =
                    spawn_reader(file_pairs[0].clone(), opts.decompress);
                let reader2: fastq::Reader<Box<dyn Read + Send>> =
                    spawn_reader(file_pairs[1].clone(), opts.decompress);

                parallel_paired_fastq(
                    PairedFastqReader { reader1, reader2 },
                    opts.threads as u32,
                    opts.threads,
                    |(read1, read2), found| {
                        *found = process_paired_reads(&read1, &read2, &matcher, &progress_logger);
                    },
                    |(read1, read2), found| {
                        if *found > 0 {
                            num_matches += 1;
                            if let Some(writer) = &mut maybe_writer {
                                if color {
                                    let found1 = *found == 1;
                                    let found2 = *found == 2 || matcher.read_match(&read2);
                                    let mut read1 = read1.to_owned_record();
                                    let mut read2 = read2.to_owned_record();
                                    matcher.color(&mut read1, found1);
                                    matcher.color(&mut read2, found2);
                                    fastq::write_to(
                                        &mut *writer,
                                        read1.head(),
                                        read1.seq(),
                                        read1.qual(),
                                    )
                                    .unwrap();
                                    fastq::write_to(
                                        &mut *writer,
                                        read2.head(),
                                        read2.seq(),
                                        read2.qual(),
                                    )
                                    .unwrap();
                                } else {
                                    fastq::write_to(
                                        &mut *writer,
                                        read1.head(),
                                        read1.seq(),
                                        read1.qual(),
                                    )
                                    .unwrap();
                                    fastq::write_to(
                                        &mut *writer,
                                        read2.head(),
                                        read2.seq(),
                                        read2.qual(),
                                    )
                                    .unwrap();
                                }
                            }
                        }
                        None::<()>
                    },
                )
                .unwrap();
            }
        }
    } else {
        // Process one FASTQ at a time
        for file in files {
            let reader: fastq::Reader<Box<dyn Read + Send>> =
                spawn_reader(file.clone(), opts.decompress);
            parallel_fastq(
                reader,
                opts.threads as u32,
                opts.threads,
                |record, found| {
                    if let Some(progress) = &progress_logger {
                        progress.record();
                    }
                    *found = matcher.read_match(&record);
                },
                |record, found| {
                    if *found {
                        num_matches += 1;
                        if let Some(writer) = &mut maybe_writer {
                            if color {
                                let mut record = record.to_owned_record();
                                matcher.color(&mut record, *found);
                                fastq::write_to(writer, record.head(), record.seq(), record.qual())
                                    .unwrap();
                            } else {
                                fastq::write_to(writer, record.head(), record.seq(), record.qual())
                                    .unwrap();
                            }
                        }
                    }
                    None::<()>
                },
            )
            .unwrap();
        }
    }

    if opts.count {
        std::io::stdout()
            .write_all(format!("{num_matches}\n").as_bytes())
            .unwrap();
    }

    // Get the final count of records matched
    Ok(num_matches)
}

#[allow(clippy::ref_option)] // FIXME: remove me later and solve
#[allow(clippy::borrowed_box)] // FIXME: remove me later and solve
fn process_paired_reads(
    read1: &RefRecord,
    read2: &RefRecord,
    matcher: &Box<dyn Matcher + Sync + Send>, //&(dyn Matcher + Sync + Send),
    progress_logger: &Option<ProgLog>,
) -> u32 {
    if let Some(progress) = progress_logger {
        progress.record();
        progress.record();
    }
    assert_eq!(
        read1.id_bytes(),
        read2.id_bytes(),
        "Mismatching read pair!  R1: {} R2: {}",
        read1.id().unwrap(),
        read2.id().unwrap()
    );
    // NB: only search for a match in read2 if read1 does not have a match
    if matcher.read_match(read1) {
        1
    } else if matcher.read_match(read2) {
        2
    } else {
        0
    }
}

/// Parse args and set up logging / tracing
fn setup() -> Opts {
    if std::env::var("RUST_LOG").is_err() {
        unsafe {
            std::env::set_var("RUST_LOG", "info");
        }
    }
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    Opts::parse()
}

// Tests
#[cfg(test)]
pub mod tests {
    use crate::*;
    use fgoxide::io::Io;
    use rstest::rstest;
    use seq_io::fastq::{OwnedRecord, Record};
    use tempfile::TempDir;

    /// Returns the path(s)  of type `Vec<String>` to the fastq(s) written from the provided sequence
    ///
    /// # Arguments
    ///
    /// * `temp_dir` - A temp directory that must be created in the actual test fucntion
    /// * `sequences` - A &Vec<Vec<&str>> where the outer Vec contains the the reads for a given FASTQ file
    /// * `file_extension` - The desired file extension i.e. .fq, .fq.gz etc.
    ///
    /// # Examples
    /// if - sequences = vec![vec!["AAGTCTGAATCCATGGAAAGCTATTG", "GGGTCTGAATCCATGGAAAGCTATTG"], vec!["AAGTCTGAATCCATGGAAAGCTATTG", "GGGTCTGAATCCATGGAAAGCTATTG"]]
    /// write_fastq() will return Vec</temp/path/to/first.fa, /temp/path/to/second.fq> where 'paths' are Strings.
    fn write_fastq(
        temp_dir: &TempDir,
        sequences_per_fastq: &Vec<Vec<&str>>,
        file_extension: String,
    ) -> Vec<String> {
        // Initialize a Vec<String>
        let mut fastq_paths = Vec::new();

        // Iterate through each FASTQ file we are to build
        for (fastq_index, fastq_sequences) in sequences_per_fastq.iter().enumerate() {
            let name = format!("sample_{fastq_index}{file_extension}");
            let fastq_path = temp_dir.path().join(name);
            let io = Io::default();
            let mut writer = io.new_writer(&fastq_path).unwrap();

            // Second loop through &str in &Vec<Vec<&str>>
            for (num, seq) in fastq_sequences.iter().enumerate() {
                let read = OwnedRecord {
                    head: format!("@{num}").as_bytes().to_vec(),
                    seq: seq.as_bytes().to_vec(),
                    qual: vec![b'X'; seq.len()],
                };
                fastq::write_to(&mut writer, &read.head, &read.seq, &read.qual)
                    .expect("failed writing read");
            }
            // Convert PathBuf to String - Opts expects Vec<String>
            // as_string() is not a method of PathBuf
            let path_as_string = fastq_path.as_path().display().to_string();
            fastq_paths.push(path_as_string);
        }
        fastq_paths
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
        let pattern_path = temp_dir.path().join(name);
        // Simply pass pattern to io.write_lines()
        io.write_lines(&pattern_path, pattern).unwrap();
        // Return
        pattern_path
    }

    /// Builds the command line options used for testing.
    /// Builds an instance of StructOpts to be passed to fqgrep_from_opts().  This will also write the
    /// FASTQ(s) (via write_fastq) and optionally writes the search pattern to a file (via write_pattern) if pattern_from_file is true.  For the latter, opts.regex and opts.file will set appropriately.
    ///
    /// # Arguments
    ///
    /// * `temp_dir` - A temp directory
    /// * `seqs` - Test sequences to be written to fastq(s)
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
        compression: String,
    ) -> Opts {
        let fq_path = write_fastq(dir, seqs, compression);

        let (pattern_string, pattern_file) = {
            if pattern_from_file {
                (vec![], Some(write_pattern(dir, regexp)))
            } else {
                (regexp.clone(), None)
            }
        };

        Opts {
            threads: 4,
            color: Color::Never,
            count: output.is_none(),
            regexp: pattern_string,
            fixed_strings: false,
            file: pattern_file,
            invert_match: false,
            decompress: false,
            paired: false,
            reverse_complement: false,
            progress: true,
            args: fq_path.clone(),
            output,
        }
    }

    /// Returns sequences from fastq
    ///
    /// # Arguments
    /// * result_path" path to fastq of matches
    ///
    /// # Example
    /// let dir: TempDir = TempDir::new().unwrap();
    /// let out_path = dir.path().join(String::from("output.fq"));
    /// let return_seq = slurp_output(out_path);
    fn slurp_output(result_path: PathBuf) -> Vec<String> {
        let mut return_seqs = vec![];
        let mut reader = fastq::Reader::from_path(result_path).unwrap();
        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            // convert bytes to str
            let seq = std::str::from_utf8(record.seq()).unwrap();
            return_seqs.push(seq.to_owned());
        }
        return_seqs
    }

    // ############################################################################################
    // Tests match with unpaired and paired reads when count is true
    // ############################################################################################
    #[rstest]
    // Unpaired reads with fixed strings
    #[case(false, vec!["AA"], 1)] // unpaired: fixed string with one match
    #[case(false, vec!["CCCG"], 0)] // unpaired: fixed string with zero matches
    #[case(false, vec!["T"], 3)] // unpaired: fixed string with multiple matches
    // Unpaired reads with regex
    #[case(false, vec!["^A"], 1)] // unpaired: regex with one match
    #[case(false, vec!["T.....G"], 0)] // unpaired: regex with zero matches
    #[case(false, vec!["G"], 4)] // unpaired: regex with multiple matches
    // Unpaired reads with mixed patterns
    #[case(false, vec!["GGCC", "G..C"], 1)] // unpaired: mixed set with one match
    #[case(false, vec!["Z", "A.....G"], 0)] // unpaired: mixed set with zero matches
    #[case(false,vec!["^T", "AA"], 3)] // unpaired: mixed set with multiple matches
    // Paired reads with fixed strings
    #[case(true, vec!["AA"], 1)] // paired: fixed string with one match
    #[case(true, vec!["CCCG"], 0)] // paired: fixed string with zero matches
    #[case(true, vec!["CC"], 2)] // paired: fixed string with multiple matches
    // Paired reads with regex
    #[case(true, vec!["^A"], 1)] // paired: regex with one match
    #[case(true, vec!["T.....G"], 0)] // paired: regex with zero matches
    #[case(true, vec!["G"], 3)] // paired: regex with multiple matches
    // Paired reads with mixed patterns
    #[case(true, vec!["GGCC", "G..C"], 1)] // paired: mixed set with one match
    #[case(true, vec!["Z", "A.....G"], 0)] // paired: mixed set with zero matches
    #[case(true, vec!["^T", "AA"], 3)] // paired: mixed set with multiple matches
    fn test_reads_when_count_true(
        #[case] paired: bool,
        #[case] pattern: Vec<&str>,
        #[case] expected: usize,
    ) {
        let dir = TempDir::new().unwrap();
        let seqs = vec![
            vec!["AAAA", "TTTT"],
            vec!["GGGG", "CCCC"],
            vec!["TTCT", "CGCG"],
            vec!["GGTT", "GGCC"],
        ];
        let pattern = pattern.iter().map(|&s| s.to_owned()).collect::<Vec<_>>();
        let mut opts = build_opts(&dir, &seqs, &pattern, true, None, String::from(".fq"));
        opts.paired = paired;
        let result = fqgrep_from_opts(&opts);
        assert_eq!(result.unwrap(), expected);
    }

    // ############################################################################################
    //Tests match with unpaired and paired reads when count is false
    // ############################################################################################
    #[rstest]
    // Unpaired reads with fixed strings
    #[case(false, vec!["A"], vec!["AAAA", "ATAT", "TATA", "AAAT", "TTTA"], true)] // unpaired: fixed string with one match
    #[case(false, vec!["A", "G"], vec!["AAAA", "ATAT", "TATA", "AAAT", "TTTA", "GGGG", "CGCG", "GCGC", "CCCG", "GGGC"], true)] // unpaired: fixed string set with two matches
    // Unpaired reads with regex
    #[case(false, vec!["^A"], vec!["AAAA", "ATAT", "AAAT"], true)] // unpaired: regex with one match
    #[case(false, vec!["^A", "^G"], vec!["AAAA", "ATAT", "AAAT", "GGGG", "GCGC", "GGGC"], true)] // unpaired: regex set with two matches
    // Paired reads with fixed string sets
    #[case(true, vec!["AAAA"], vec!["AAAA", "CCCC"], true)] // paired: fixed string with one match
    #[case(true, vec!["A", "G"], vec!["AAAA", "CCCC", "TTTT", "GGGG", "ATAT", "CGCG", "TATA", "GCGC", "AAAT", "CCCG", "TTTA", "GGGC"], true)]
    // paired: fixed string set with two matches in correct interleave order
    #[case(true, vec!["A", "G"], vec!["AAAA", "GGGG", "TTTT", "CCCC", "CGCG", "ATAT", "GCGC", "TATA", "CCCG", "AAAT", "GGGC", "TTTA"], false)]
    // paired: fixed string set with two matches in incorrect interleave order
    // Paired reads with regex sets
    #[case(true, vec!["^AAAA"], vec!["AAAA", "CCCC"], true)] // paired: regex with one match
    #[case(true, vec!["^AA", "^GG"], vec!["AAAA", "CCCC", "TTTT", "GGGG", "AAAT", "CCCG", "TTTA", "GGGC"], true)] // paired: regex set with two matches in correct interleave order
    #[case(true, vec!["^AA", "^GG"], vec!["AAAA", "GGGG", "TTTT", "CCCC", "TTTA", "CCCG", "AAAT", "GGGC"], false)] // paired: regex set with two matches in incorrect interleave order
    fn test_reads_when_count_false(
        #[case] paired: bool,
        #[case] pattern: Vec<&str>,
        #[case] expected_seq: Vec<&str>,
        #[case] expected_bool: bool,
    ) {
        let dir = TempDir::new().unwrap();
        let seqs = vec![
            vec!["AAAA", "TTTT", "ATAT", "TATA", "AAAT", "TTTA"],
            vec!["CCCC", "GGGG", "CGCG", "GCGC", "CCCG", "GGGC"],
        ];
        let out_path = dir.path().join(String::from("output.fq"));
        let result_path = &out_path.clone();
        let pattern = pattern.iter().map(|&s| s.to_owned()).collect::<Vec<_>>();
        let mut opts = build_opts(
            &dir,
            &seqs,
            &pattern,
            true,
            Some(out_path),
            String::from(".fq"),
        );

        opts.paired = paired;
        let _result = fqgrep_from_opts(&opts);
        let return_sequences = slurp_output(result_path.to_path_buf());
        let sequence_match = expected_seq == return_sequences;
        assert_eq!(sequence_match, expected_bool);
    }

    // ############################################################################################
    // Tests two fastqs for 'TGGATTCAGACTT' which is only found once in the reverse complement
    // Tests inverse_match
    // Tests both paired and unpaired reads
    // ############################################################################################
    #[rstest]
    // Unpaired reads
    #[case(false, false, false, 0)] // unpaired: zero matches when invert_match and reverse_complement are false
    #[case(false, false, true, 1)] //  unpaired: one match when invert_match is false and reverse_complement is true
    #[case(false, true, false, 4)] //  unpaired: four matches when invert_match is true and reverse_complement is false
    #[case(false, true, true, 3)] // unpaired: three matches when invert_match and reverse_complement are true
    // Paired reads
    #[case(true, false, false, 0)] // paired: zero matches when invert_match and reverse_complement are false
    #[case(true, false, true, 1)] //  paired: one match when invert_match is false and reverse_complement is true
    #[case(true, true, false, 2)] //  paired: two matches when invert_match is true and reverse_complement is false
    #[case(true, true, true, 2)] // paired: two matches when invert_match and reverse_complement are true
    fn test_reverse_complement_and_invert_match(
        #[case] paired: bool,
        #[case] invert_match: bool,
        #[case] reverse_complement: bool,
        #[case] expected: usize,
    ) {
        let dir = TempDir::new().unwrap();
        let seqs = vec![vec!["GGGG", "GGGG"], vec!["AAAA", "CCCC"]];
        let pattern = vec![String::from("TTTT")];
        let mut opts = build_opts(&dir, &seqs, &pattern, false, None, String::from(".fq"));

        opts.paired = paired;
        opts.invert_match = invert_match;
        opts.reverse_complement = reverse_complement;

        let result = fqgrep_from_opts(&opts);
        assert_eq!(result.unwrap(), expected);
    }

    // ############################################################################################
    // Tests that an error is returned when fixed_strings is true and regex is present
    // ############################################################################################
    #[rstest]
    #[should_panic(
        expected = "called `Result::unwrap()` on an `Err` value: Fixed pattern must contain only DNA bases:  .. [^] .. G"
    )]
    #[case(true, vec![String::from("^G")])] // panic with single regex
    #[should_panic(
        expected = "called `Result::unwrap()` on an `Err` value: Fixed pattern must contain only DNA bases:  .. [^] .. A"
    )]
    #[case(true, vec![String::from("^A"),String::from("AA")])] // panic with combination of regex and fixed string
    fn test_regexp_from_fixed_string_fails_with_regex(
        #[case] fixed_strings: bool,
        #[case] pattern: Vec<String>,
    ) {
        let dir = TempDir::new().unwrap();
        let seqs = vec![vec!["GGGG", "TTTT"], vec!["AAAA", "CCCC"]];

        let mut opts = build_opts(&dir, &seqs, &pattern, true, None, String::from(".fq"));

        opts.fixed_strings = fixed_strings;
        let _result = fqgrep_from_opts(&opts);
        _result.unwrap();
    }
    // ############################################################################################
    // Tests error is returned from main when three records are defined as paired
    // ############################################################################################
    #[test]
    #[should_panic(
        expected = "Input files must be a multiple of two, or either a single file or standard input (assume interleaved) with --paired"
    )]
    fn test_paired_should_panic_with_three_records() {
        let dir = TempDir::new().unwrap();
        let seqs = vec![
            vec!["GTCAGCTCGAGCATCAGCTACGCACT"],
            vec!["AGTGCGTAGCTGATGCTCGAGCTGAC"],
            vec!["GGGTCTGAATCCATGGAAAGCTATTG"],
        ];

        let test_pattern = vec![String::from("A")];
        let mut opts_test = build_opts(&dir, &seqs, &test_pattern, true, None, String::from(".fq"));

        opts_test.paired = true;
        let _num_matches = fqgrep_from_opts(&opts_test);
        _num_matches.unwrap();
    }

    // ############################################################################################
    // Tests that correct match count is returned regardless of pattern provided via file or
    // string
    // ############################################################################################
    #[test]
    fn test_regexp_matches_from_file_and_string() {
        let dir = TempDir::new().unwrap();
        let seqs = vec![
            vec!["GTCAGCTCGAGCATCAGCTACGCACT"],
            vec!["AGTGCGTAGCTGATGCTCGAGCTGAC"],
            vec!["GGGTCTGAATCCATGGAAAGCTATTG"],
        ];

        let test_pattern = vec![String::from("^G")];
        let mut opts_test = build_opts(&dir, &seqs, &test_pattern, true, None, String::from(".fq"));

        // Test pattern from file
        let result = fqgrep_from_opts(&opts_test);
        assert_eq!(result.unwrap(), 2);

        // Test pattern from string
        opts_test.regexp = test_pattern;
        let result = fqgrep_from_opts(&opts_test);
        assert_eq!(result.unwrap(), 2);
    }

    // ############################################################################################
    // Tests that correct match count is returned from .fq, .fastq, .fq.gz, and .fq.bgz
    // ############################################################################################
    // Three matches are found from all file types
    #[rstest]
    #[case(String::from(".fq"), 3)]
    #[case(String::from(".fastq"), 3)]
    #[case(String::from(".fq.gz"), 3)]
    #[case(String::from(".fq.bgz"), 3)]
    fn test_fastq_compression(#[case] extension: String, #[case] expected: usize) {
        let dir = TempDir::new().unwrap();
        let seqs = vec![
            vec!["GTCAG", "ACGT", "GGTG"],
            vec!["AGTGCGTAGCTGATGCTCGAGCTGAC"],
            vec!["GGGTCTGAATCCATGGAAAGCTATTG"],
        ];

        let test_pattern = vec![String::from("^G")];

        let opts = build_opts(&dir, &seqs, &test_pattern, true, None, extension);
        let result = fqgrep_from_opts(&opts);
        assert_eq!(result.unwrap(), expected);
    }

    // ############################################################################################
    // Tests ExitCode status 101, 1, None, and 2
    // ############################################################################################

    #[rstest]
    #[case(false, vec![String::from("*")], Some(101))] // invalid regex will panic - ExitCode(101)
    #[case(false, vec![String::from("^T")], Some(1))] // zero matches - ExitCode(1)
    #[case(true, vec![String::from("GTCAGC")], None)] // one match - ExitCode(SUCCESS) None returned from fqgrep()
    #[case(true, vec![String::from("^T")], Some(2))] // returns inner error when regex is declared fixed_string - ExitCode(2)
    fn test_exit_code_catching(
        #[case] fixed_strings: bool,
        #[case] pattern: Vec<String>,
        #[case] expected: Option<u8>,
    ) {
        let dir = TempDir::new().unwrap();
        let seqs = vec![vec!["GTCAGC"], vec!["AGTGCG"], vec!["GGGTCTG"]];
        let mut opts = build_opts(&dir, &seqs, &pattern, true, None, String::from(".fq"));
        opts.fixed_strings = fixed_strings;
        assert_eq!(fqgrep(&opts), expected);
    }
}
