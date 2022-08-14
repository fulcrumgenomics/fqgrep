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
use rayon::{prelude::*, ThreadPool};
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

struct FastqWriter {
    tx: Sender<Vec<OwnedRecord>>,
    lock: Mutex<()>,
}

impl FastqWriter {
    fn new(count_tx: Sender<usize>, pool: &ThreadPool, count: bool, paired: bool) -> Self {
        // TODO: try making this unbounded
        let (tx, rx): (Sender<Vec<OwnedRecord>>, Receiver<Vec<OwnedRecord>>) =
            bounded(WRITER_CHANNEL_SIZE);
        if count {
            pool.spawn(move || {
                let mut num_matches = 0;
                while let Ok(reads) = rx.recv() {
                    num_matches += reads.len();
                }
                if paired {
                    num_matches /= 2;
                }
                std::io::stdout()
                    .write_all(format!("{}\n", num_matches).as_bytes())
                    .unwrap();

                count_tx
                    .send(num_matches)
                    .expect("failed sending final count");
            });
        } else {
            pool.spawn(move || {
                let mut num_matches = 0;
                let mut writer = BufWriter::with_capacity(BUFSIZE, std::io::stdout());
                while let Ok(reads) = rx.recv() {
                    num_matches += reads.len();
                    for read in reads {
                        fastq::write_to(&mut writer, &read.head, &read.seq, &read.qual)
                            .expect("failed writing read");
                    }
                }
                if paired {
                    num_matches /= 2;
                }
                writer.flush().expect("Error flushing writer");
                count_tx
                    .send(num_matches)
                    .expect("failed sending final count");
            });
        }
        let lock = Mutex::new(());
        Self { tx, lock }
    }
}

fn spawn_reader(file: PathBuf, decompress: bool) -> Receiver<Vec<OwnedRecord>> {
    let (tx, rx) = bounded(READER_CHANNEL_SIZE);
    rayon::spawn(move || {
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
/// one or more patterns.   By default, a pattern matches the bases in a FASTQ record if the
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
    /// The number of threads to use for matching reads against pattern
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
    /// Input files must be gzip compressed unless `--plain` is given
    args: Vec<String>,
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

#[allow(clippy::too_many_lines)]
fn fqgrep() -> Result<usize> {
    let mut opts = setup();

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
    let writer = FastqWriter::new(count_tx, &pool, opts.count, opts.paired);

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
            // Process one FSATQ at a time
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
