#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::PathBuf,
    str::FromStr,
    thread::JoinHandle,
};

use anyhow::{bail, ensure, Context, Result};
use env_logger::Env;
use flate2::bufread::MultiGzDecoder;
use flume::{bounded, Receiver, Sender};
use fqgrep_lib::matcher::{
    FixedStringMatcher, FixedStringSetMatcher, Matcher, MatcherOpts, RegexMatcher, RegexSetMatcher,
};
use gzp::BUFSIZE;
use itertools::{self, izip, Itertools};
use lazy_static::lazy_static;
use parking_lot::Mutex;
use proglog::{CountFormatterKind, ProgLog, ProgLogBuilder};
use rayon::{prelude::*, Scope, ThreadPool};
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
    fn new(pool: &ThreadPool, count: bool, paired: bool) -> Self {
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
            });
        } else {
            pool.spawn(move || {
                let mut writer = BufWriter::with_capacity(BUFSIZE, std::io::stdout());
                while let Ok(reads) = rx.recv() {
                    for read in reads {
                        fastq::write_to(&mut writer, &read.head, &read.seq, &read.qual)
                            .expect("failed writing read");
                    }
                }
                writer.flush().expect("Error flushing writer");
            });
        }
        let lock = Mutex::new(());
        Self { tx, lock }
    }
}

struct ThreadReader {
    handle: JoinHandle<()>,
    rx: Receiver<Vec<OwnedRecord>>,
}

impl ThreadReader {
    fn new(file: PathBuf, plain: bool) -> Self {
        let (tx, rx) = bounded(READER_CHANNEL_SIZE);
        let handle = std::thread::spawn(move || {
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
            let maybe_decoder_handle = if plain {
                Box::new(buf_handle) as Box<dyn Read>
            } else {
                Box::new(MultiGzDecoder::new(buf_handle)) as Box<dyn Read>
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

        Self { handle, rx }
    }
}

/// A small utility to "grep" a pair of gzipped FASTQ files, outputting the read pair if the sequence
/// matches either the given ref sequence, the given alt sequence, or both.
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

    /// Treat the input files as paired.  The number of input files must be a multiple of two,
    /// with the first file being R1, second R2, third R1, fourth R2, and so on.  If the pattern
    /// matches either R1 or R2, then both R1 and R2 will be output (interleaved).  If the input
    /// is standard input, then treat the input as interlaved paired end reads.
    #[structopt(long)]
    paired: bool,

    /// Only a count of selected lines is written to standard output.
    #[structopt(long, short = "c")]
    count: bool,

    /// Specify a pattern used during the search of the input: an input line is selected if it
    /// matches any of the specified patterns.  This option is most useful when multiple `-e`
    /// options are used to specify multiple patterns.
    #[structopt(long, short = "e")]
    regexp: Vec<String>,

    // Interpret pattern as a set of fixed strings
    #[structopt(long, short = "F")]
    fixed_strings: bool,

    /// Read one or more newline separated patterns from file.  Empty pattern lines match every
    /// input line.  Newlines are not considered part of a pattern.  If file is empty, nothing
    /// is matched.
    #[structopt(long, short = "f")]
    file: Option<PathBuf>,

    /// Search the reverse complement for matches.
    #[structopt(long)]
    reverse_complement: bool,

    /// Selected lines are those not matching any of the specified patterns
    #[structopt(short = "v")]
    invert_match: bool,

    /// Write progress information
    #[structopt(long)]
    progress: bool,

    /// The input FASTQ(s) is not compressed
    #[structopt(long)]
    plain: bool,

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

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
    let mut opts = setup();

    if let Some(file) = &opts.file {
        for pattern in read_patterns(file)? {
            opts.regexp.push(pattern);
        }
    }

    let (pattern, mut files): (Option<String>, Vec<PathBuf>) = {
        let (pattern, file_strings): (Option<String>, Vec<String>) = if opts.regexp.is_empty() {
            ensure!(
                !opts.args.is_empty(),
                "Pattern must be given with -e or as the first positional argument "
            );
            let files = opts.args.iter().skip(1).cloned().collect();
            (Some(opts.args[0].clone()), files)
        } else {
            (None, opts.args.clone())
        };

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

    if opts.paired {
        ensure!(
            files.len() <= 1 || files.len() % 2 == 0,
            "Input files must be a multiple of two, or either a single file or standard input (assume interleaved) with --paired"
        );
    }

    if opts.fixed_strings {
        if let Some(pattern) = &pattern {
            FixedStringMatcher::validate(pattern)?;
        } else if !opts.regexp.is_empty() {
            for pattern in &opts.regexp {
                FixedStringMatcher::validate(pattern)?;
            }
        } else {
            bail!("A pattern must be given as a positional argument or with -e/--regexp")
        }
    }

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

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.threads)
        .build()
        .unwrap();

    let writer = FastqWriter::new(&pool, opts.count, opts.paired);

    pool.install(|| {
        let match_opts = MatcherOpts {
            invert_match: opts.invert_match,
            reverse_complement: opts.reverse_complement,
        };
        let matcher: Box<dyn Matcher + Sync> = match (opts.fixed_strings, &pattern) {
            (true, Some(pattern)) => Box::new(FixedStringMatcher::new(pattern, match_opts)),
            (false, Some(pattern)) => Box::new(RegexMatcher::new(pattern, match_opts)),
            (true, None) => Box::new(FixedStringSetMatcher::new(&opts.regexp, match_opts)),
            (false, None) => Box::new(RegexSetMatcher::new(&opts.regexp, match_opts)),
        };

        if opts.paired {
            if files.is_empty() {
                // read from standad input
                files.push(PathBuf::from_str("-").unwrap());
            }
            if files.len() == 1 {
                // interleaved paired end FASTQ
                let reader = ThreadReader::new(files[0].clone(), opts.plain);
                for reads in izip!(reader.rx.iter()) {
                    let paired_reads = reads
                        .into_iter()
                        .tuples::<(OwnedRecord, OwnedRecord)>()
                        .collect_vec();
                    process_paired_reads(paired_reads, &matcher, &writer, &progress_logger);
                }
                reader.handle.join().expect("Failed to join reader");
            } else {
                // pairs of FASTQ files
                for file_pairs in files.chunks_exact(2) {
                    let reader1 = ThreadReader::new(file_pairs[0].clone(), opts.plain);
                    let reader2 = ThreadReader::new(file_pairs[1].clone(), opts.plain);
                    for (reads1, reads2) in izip!(reader1.rx.iter(), reader2.rx.iter()) {
                        let paired_reads = reads1.into_iter().zip(reads2.into_iter()).collect_vec();
                        process_paired_reads(paired_reads, &matcher, &writer, &progress_logger);
                    }
                    reader1.handle.join().expect("Failed to join reader");
                    reader2.handle.join().expect("Failed to join reader");
                }
            }
        } else {
            if files.is_empty() {
                // read from standad input
                files.push(PathBuf::from_str("-").unwrap());
            }
            for file in files {
                let reader = ThreadReader::new(file.clone(), opts.plain);
                for reads in reader.rx.iter() {
                    // Get the matched reads
                    let matched_reads: Vec<OwnedRecord> = reads
                        .into_par_iter()
                        .map(|read| -> Option<OwnedRecord> {
                            if let Some(progress) = &progress_logger {
                                progress.record();
                            }
                            if matcher.read_match(&read) {
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

                reader.handle.join().expect("Failed to join reader");
            }
        }
    });

    while !writer.tx.is_empty() {}

    Ok(())
}

#[allow(clippy::borrowed_box)] // FIXME: remove me later and solve
fn process_paired_reads(
    reads: Vec<(OwnedRecord, OwnedRecord)>,
    matcher: &Box<dyn Matcher + Sync>,
    writer: &FastqWriter,
    progress_logger: &Option<ProgLog>,
) {
    reads
        .into_par_iter()
        .map(|(read1, read2)| {
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
            if matcher.read_match(&read1) || matcher.read_match(&read2) {
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
