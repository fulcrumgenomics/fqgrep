#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::{
    borrow::Borrow,
    fs,
    fs::File,
    io::{BufReader, BufWriter},
    path::{Path, PathBuf},
    thread::JoinHandle,
};

use anyhow::{Error, Result};
use bstr::ByteSlice;
use env_logger::Env;
use flate2::bufread::MultiGzDecoder;
use flume::{bounded, Receiver, Sender};
use gzp::{deflate::Bgzf, Compression, ZBuilder, ZWriter, BUFSIZE};
use itertools::{self, izip, Itertools};
use lazy_static::lazy_static;
use log::info;
use parking_lot::Mutex;
use rayon::prelude::*;
use regex::{Regex, RegexBuilder};
use seq_io::fastq::{self, OwnedRecord};
use serde::{Deserialize, Serialize};
use structopt::{clap::AppSettings::ColoredHelp, StructOpt};

/// The number of reads in a chunk
const CHUNKSIZE: usize = 5000; // * num_cpus::get();
/// The number of chunks allowed in a channel
const READER_CHANNEL_SIZE: usize = 100;
const WRITER_CHANNEL_SIZE: usize = 2000;

const REF_PREFIX: &str = "REF";
const ALT_PREFIX: &str = "ALT";
const BOTH_PREFIX: &str = "BOTH";

/// Helper type to represent the accumlated reads and what match type they are.
///
/// `p` are ref reads
/// `1` are alt reads
/// `2` are reads that matched both
type MatchedReads = (
    (Vec<OwnedRecord>, Vec<OwnedRecord>),
    (Vec<OwnedRecord>, Vec<OwnedRecord>),
    (Vec<OwnedRecord>, Vec<OwnedRecord>),
);

lazy_static! {
    /// Return the number of cpus as a String
    pub static ref NUM_CPU: String = num_cpus::get().to_string();
}

lazy_static! {
    static ref R1_REGEX: Regex = RegexBuilder::new(r"r1")
        .case_insensitive(true)
        .build()
        .expect("Failed to compile R1 regex");
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

lazy_static! {
    static ref COMPLEMENT: [u8; 256] = {
        let mut comp = [0; 256];
        for (v, a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in b"AGCTYRWSKMDVHBN".iter().zip(b"TCGARYWSMKHBDVN".iter()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        comp
    };
}

fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

fn revcomp<C, T>(text: T) -> Vec<u8>
where
    C: Borrow<u8>,
    T: IntoIterator<Item = C>,
    T::IntoIter: DoubleEndedIterator,
{
    text.into_iter()
        .rev()
        .map(|a| complement(*a.borrow()))
        .collect()
}

#[derive(Debug, Serialize, Deserialize)]
struct FqGrepMetrics {
    total_alt_reads_found: usize,
    total_ref_reads_found: usize,
    total_both_reads_found: usize,
    total_reads_found: usize,
}

impl FqGrepMetrics {
    fn new(
        total_ref_reads_found: usize,
        total_alt_reads_found: usize,
        total_both_reads_found: usize,
    ) -> Self {
        Self {
            total_ref_reads_found,
            total_alt_reads_found,
            total_both_reads_found,
            total_reads_found: total_alt_reads_found
                + total_both_reads_found
                + total_ref_reads_found,
        }
    }
}

#[derive(Debug, Clone)]
struct Sample {
    name: String,
}

impl Sample {
    fn new(name: String) -> Self {
        Self { name }
    }

    fn from_r1_path<P: AsRef<Path>>(r1_path: P, prefix: &str) -> Result<Self> {
        let filename = r1_path
            .as_ref()
            .file_name()
            .expect("Unable to extract file name from R1")
            .to_string_lossy();
        let parts: Vec<&str> = R1_REGEX.splitn(&filename, 2).collect();
        if let Some(&part) = parts.first() {
            let name: &str = &part[0..part.len() - 1];
            Ok(Self::new(format!("{}_{}", prefix, name.to_owned())))
        } else {
            Err(Error::msg("Unable to extract a sample name from R1 path"))
        }
    }

    fn filenames<P: AsRef<Path>>(&self, dir: P, illumina: bool) -> (PathBuf, PathBuf) {
        if illumina {
            (
                dir.as_ref()
                    .join(format!("{}_S1_L001_R1_001.fastq.gz", self.name)),
                dir.as_ref()
                    .join(format!("{}_S1_L001_R2_001.fastq.gz", self.name)),
            )
        } else {
            (
                dir.as_ref().join(format!("{}.R1.fastq.gz", self.name)),
                dir.as_ref().join(format!("{}.R2.fastq.gz", self.name)),
            )
        }
    }
}

struct WriterStats {
    reads_written: usize,
}

impl WriterStats {
    fn new() -> Self {
        Self { reads_written: 0 }
    }
}

struct ThreadWriter {
    handle: JoinHandle<WriterStats>,
    tx: Option<Sender<Vec<OwnedRecord>>>,
}

impl ThreadWriter {
    fn new(mut writer: ZWriterWrapper) -> Self {
        // TODO: try making this unbounded
        let (tx, rx): (Sender<Vec<OwnedRecord>>, Receiver<Vec<OwnedRecord>>) =
            bounded(WRITER_CHANNEL_SIZE);
        let handle = std::thread::spawn(move || {
            let mut writer_stats = WriterStats::new();
            while let Ok(reads) = rx.recv() {
                for read in reads {
                    writer_stats.reads_written += 1;
                    fastq::write_to(&mut writer.0, &read.head, &read.seq, &read.qual)
                        .expect("failed writing read");
                }
            }
            writer.0.finish().expect("Error flushing writer");
            writer_stats
        });
        Self {
            handle,
            tx: Some(tx),
        }
    }
}
struct SampleWriter {
    r1: ThreadWriter, // + Send
    r2: ThreadWriter, // + Send
    lock: Mutex<()>,
}

impl SampleWriter {
    fn new<P: AsRef<Path>>(
        sample: &Sample,
        dir: P,
        compression_threads: usize,
        illumina: bool,
    ) -> Result<Self> {
        let (r1, r2) = sample.filenames(dir, illumina);

        let r1_writer = ZBuilder::<Bgzf, _>::new()
            .num_threads(compression_threads)
            .compression_level(Compression::new(2))
            .from_writer(BufWriter::with_capacity(BUFSIZE, File::create(r1)?));
        let r1 = ThreadWriter::new(ZWriterWrapper(r1_writer));

        let r2_writer = ZBuilder::<Bgzf, _>::new()
            .num_threads(compression_threads)
            .compression_level(Compression::new(2))
            .from_writer(BufWriter::with_capacity(BUFSIZE, File::create(r2)?));
        let r2 = ThreadWriter::new(ZWriterWrapper(r2_writer));
        let lock = Mutex::new(());

        Ok(Self { r1, r2, lock })
    }
}

struct ZWriterWrapper(Box<dyn ZWriter>);

// TODO: make ZWriter's return type be Send
unsafe impl Send for ZWriterWrapper {}
unsafe impl Sync for ZWriterWrapper {}

struct ThreadReader {
    handle: JoinHandle<()>,
    rx: Receiver<Vec<OwnedRecord>>,
}

impl ThreadReader {
    fn new(file: PathBuf) -> Self {
        let (tx, rx) = bounded(READER_CHANNEL_SIZE);
        let handle = std::thread::spawn(move || {
            let reader = fastq::Reader::with_capacity(
                MultiGzDecoder::new(BufReader::with_capacity(
                    BUFSIZE,
                    File::open(&file).expect("error in opening file"),
                )),
                BUFSIZE,
            )
            .into_records()
            .chunks(CHUNKSIZE * num_cpus::get());
            let chunks = reader.into_iter();

            for chunk in chunks {
                // let now = std::time::Instant::now();
                tx.send(chunk.map(|r| r.expect("Error reading")).collect())
                    .expect("Error sending");
            }
        });

        Self { handle, rx }
    }
}

#[derive(Debug, PartialEq, PartialOrd)]
enum Matches {
    Ref,
    Alt,
    Both,
    None,
}

/// A small utility to "grep" a pair of gzipped FASTQ files, outputting the read pair if the sequence
/// matches either the given ref sequence, the given alt sequence, or both.
#[derive(StructOpt, Debug)]
#[structopt(name = "fqgrep", global_setting(ColoredHelp), version = built_info::VERSION.as_str())]
struct Opts {
    /// Path to the input R1 gzipped FASTQ
    #[structopt(long, short = "1", display_order = 1)]
    r1_fastq: PathBuf,

    /// Path to the input R2 gzipped FASTQ
    #[structopt(long, short = "2", display_order = 2)]
    r2_fastq: PathBuf,

    /// The output directory to write to.
    ///
    /// `fqgrep` will write a pair of fastqs for the alt matched reads to: {output_dir}/ALT_{sample_name}R{1,2}.fastq.gz
    /// and the fastqs for the ref matched reads to {output_dir}/REF_{sample_name}R{1,2}.fastq.gz
    #[structopt(long, short = "o", display_order = 3)]
    output_dir: PathBuf,

    /// The "ref" search sequence to grep for, for R2 reads the search sequence is reverse complemented.
    #[structopt(long, short = "r", display_order = 4)]
    ref_search_sequence: String,

    /// The "alt" search sequence to grep for, for R2 reads the search sequence is reverse complemented.
    #[structopt(long, short = "a", display_order = 5)]
    alt_search_sequence: String,

    /// The number of threads to use for matching reads against pattern
    ///
    /// Keep in mind that there is are two reader threads and two writer threads created by default.
    #[structopt(long, short = "t", default_value = NUM_CPU.as_str(), display_order = 6)]
    threads_for_matching: usize,

    /// Number of threads to use for compression.
    ///
    /// Zero is probably the correct answer if the searched sequence is relatively rare.
    /// If it's common and many reads will be written, it may be worth decreasing the number of
    /// `treads_for_matching` by a small number and making this number > 2
    #[structopt(long, short = "c", default_value = "0", display_order = 7)]
    compression_threads: usize,

    /// True to name output FASTQs based on Illumin's FASTQ file name convetions.
    ///
    /// If true, will use the file name pattern "<Sample Name>_S1_L001_R<1 or 2>_001.fastq.gz".
    /// Otherwise, will use "<Sample Name>.R<1 or 2>.fastq.gz"
    #[structopt(long, short = "I", display_order = 8)]
    illumina: bool,
}

#[derive(Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
struct Barcode<'a>(&'a [u8], &'a [u8]);

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
    let opts = setup();
    let ref_sample = Sample::from_r1_path(&opts.r1_fastq, REF_PREFIX)?;
    let alt_sample = Sample::from_r1_path(&opts.r1_fastq, ALT_PREFIX)?;
    let both_sample = Sample::from_r1_path(&opts.r1_fastq, BOTH_PREFIX)?;

    // Create the output directory
    fs::create_dir_all(&opts.output_dir)?;

    let mut ref_writer = SampleWriter::new(
        &ref_sample,
        &opts.output_dir,
        opts.compression_threads,
        opts.illumina,
    )?;
    let mut alt_writer = SampleWriter::new(
        &alt_sample,
        &opts.output_dir,
        opts.compression_threads,
        opts.illumina,
    )?;
    let mut both_writer = SampleWriter::new(
        &both_sample,
        &opts.output_dir,
        opts.compression_threads,
        opts.illumina,
    )?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.threads_for_matching)
        .build()
        .unwrap();

    // TODO: make sure errors aren't being eaten
    info!("Creating reader threads");
    let r1_reader = ThreadReader::new(opts.r1_fastq.clone());
    let r2_reader = ThreadReader::new(opts.r2_fastq.clone());

    let ref_forward_seq = opts.ref_search_sequence.as_bytes();
    let ref_revcomp_seq = revcomp(ref_forward_seq);
    let alt_forward_seq = opts.alt_search_sequence.as_bytes();
    let alt_revcomp_seq = revcomp(alt_forward_seq);

    info!("Processing reads");
    pool.install(|| {
        for (r1_super, r2_super) in izip!(r1_reader.rx.iter(), r2_reader.rx.iter()) {
            (r1_super, r2_super)
                .into_par_iter()
                .map(|(r1, r2)| {
                    let ref_match = r1.seq.find(ref_forward_seq).is_some()
                        || r1.seq.find(&ref_revcomp_seq).is_some()
                        || r2.seq.find(ref_forward_seq).is_some()
                        || r2.seq.find(&ref_revcomp_seq).is_some();
                    let alt_match = r1.seq.find(alt_forward_seq).is_some()
                        || r1.seq.find(&alt_revcomp_seq).is_some()
                        || r2.seq.find(alt_forward_seq).is_some()
                        || r2.seq.find(&alt_revcomp_seq).is_some();
                    let m = match (ref_match, alt_match) {
                        (true, true) => Matches::Both,
                        (true, false) => Matches::Ref,
                        (false, true) => Matches::Alt,
                        (false, false) => Matches::None,
                    };
                    (m, (r1, r2))
                })
                .filter(|(m, _)| *m != Matches::None)
                .fold(
                    || ((vec![], vec![]), (vec![], vec![]), (vec![], vec![])),
                    |(mut ref_reads, mut alt_reads, mut both_reads): MatchedReads,
                     (match_type, (r1, r2)): (Matches, (OwnedRecord, OwnedRecord))| {
                        match match_type {
                            Matches::Ref => {
                                ref_reads.0.push(r1);
                                ref_reads.1.push(r2);
                            }
                            Matches::Alt => {
                                alt_reads.0.push(r1);
                                alt_reads.1.push(r2);
                            }
                            Matches::Both => {
                                both_reads.0.push(r1);
                                both_reads.1.push(r2);
                            }
                            _ => unreachable!()
                        }
                        (ref_reads, alt_reads, both_reads)
                    },
                )
                .for_each(
                    |(ref_reads, alt_reads, both_reads): MatchedReads| {
                        // Aquire lock to ensure R1 and R2 are enqueued at the same time
                        {
                            let _lock = ref_writer.lock.lock();
                            ref_writer.r1.tx.as_ref().unwrap().send(ref_reads.0).expect("Failed to send R1s");
                            ref_writer.r2.tx.as_ref().unwrap().send(ref_reads.1).expect("Failed to send R1s");
                        }
                        {
                            let _lock = alt_writer.lock.lock();
                            alt_writer.r1.tx.as_ref().unwrap().send(alt_reads.0).expect("Failed to send R1s");
                            alt_writer.r2.tx.as_ref().unwrap().send(alt_reads.1).expect("Failed to send R1s");
                        }
                        {
                            let _lock = both_writer.lock.lock();
                            both_writer.r1.tx.as_ref().unwrap().send(both_reads.0).expect("Failed to send R1s");
                            both_writer.r2.tx.as_ref().unwrap().send(both_reads.1).expect("Failed to send R1s");
                        }
                    },
                );
        }
    });

    info!("Joining reader threads");
    r1_reader.handle.join().expect("Failed to join r1");
    r2_reader.handle.join().expect("Failed to join r2");

    // Flush and close writers
    info!("Joining writer threads and collecting stats");
    while !ref_writer.r1.tx.as_ref().unwrap().is_empty() {}
    while !ref_writer.r2.tx.as_ref().unwrap().is_empty() {}
    ref_writer.r1.tx.take();
    ref_writer.r2.tx.take();
    let ref_r1_stats = ref_writer
        .r1
        .handle
        .join()
        .expect("Error joining r1 handle");
    let _ref_r2_stats = ref_writer
        .r2
        .handle
        .join()
        .expect("Error joining r2 handle");
    while !alt_writer.r1.tx.as_ref().unwrap().is_empty() {}
    while !alt_writer.r2.tx.as_ref().unwrap().is_empty() {}
    alt_writer.r1.tx.take();
    alt_writer.r2.tx.take();
    let alt_r1_stats = alt_writer
        .r1
        .handle
        .join()
        .expect("Error joining r1 handle");
    let _alt_r2_stats = alt_writer
        .r2
        .handle
        .join()
        .expect("Error joining r2 handle");
    while !both_writer.r1.tx.as_ref().unwrap().is_empty() {}
    while !both_writer.r2.tx.as_ref().unwrap().is_empty() {}
    both_writer.r1.tx.take();
    both_writer.r2.tx.take();
    let both_r1_stats = both_writer
        .r1
        .handle
        .join()
        .expect("Error joining r1 handle");
    let _both_r2_stats = both_writer
        .r2
        .handle
        .join()
        .expect("Error joining r2 handle");

    info!(
        "{} reads matched input ref sequence",
        ref_r1_stats.reads_written
    );
    info!(
        "{} reads matched input alt sequence",
        alt_r1_stats.reads_written
    );
    info!(
        "{} reads matched both input ref and alt sequences",
        both_r1_stats.reads_written
    );

    info!("Writing stats");
    let fq_grep_metrics = FqGrepMetrics::new(
        ref_r1_stats.reads_written,
        alt_r1_stats.reads_written,
        both_r1_stats.reads_written,
    );

    let writer = BufWriter::new(File::create(opts.output_dir.join("fqgrep_stats.tsv"))?);
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(writer);
    writer.serialize(fq_grep_metrics)?;
    writer.flush()?;

    Ok(())
}

/// Parse args and set up logging / tracing
fn setup() -> Opts {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    Opts::from_args()
}
