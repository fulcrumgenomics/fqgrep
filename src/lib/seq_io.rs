/// Additions to `seq_io::parallel` for interleaved paired end FASTQ files.
///
/// Relies on: <https://github.com/markschl/seq_io/pull/23>.
///
/// Adds support to use the `seq_io::parallel` module to read interleaved paired end FASTQ files as
/// well as read pairs across pairs of FASTQ files.
///
/// For interleaved paired end FASTQ files, a custom reader `InterleavedFastqReader` is required to
/// read an even number of records, since consecutive pairs are assumed to be read pairs (mates).
///
/// For pairs across pairs of FASTQ files, a custom reader `PairsAcrossPairsReader` is required to
/// synchronize records from the input files, ensuring the same number records are read from each
/// file.
use anyhow::Result;
use itertools::{self, Itertools};
use seq_io::fastq::{self, Record, RecordSet, RefRecord};
use seq_io::parallel::{ReusableReader, read_parallel_init};
use seq_io::parallel_record_impl;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io;

parallel_record_impl!(
    parallel_interleaved_fastq,
    parallel_interleaved_fastq_init,
    R,
    InterleavedFastqReader<R>,
    InterleavedRecordSet,
    (fastq::RefRecord, fastq::RefRecord),
    fastq::Error
);

/// Read that reads an interleaved paired end FASTQ file.
pub struct InterleavedFastqReader<R: std::io::Read, P = seq_io::policy::StdPolicy> {
    pub reader: fastq::Reader<R, P>,
    rset_len: Option<usize>,
}

impl<R, P> InterleavedFastqReader<R, P>
where
    R: std::io::Read,
    P: seq_io::policy::BufPolicy + Send,
{
    pub fn new(reader: fastq::Reader<R, P>) -> Self {
        Self {
            reader,
            rset_len: None,
        }
    }
}

impl<R, P> seq_io::parallel::Reader for InterleavedFastqReader<R, P>
where
    R: std::io::Read,
    P: seq_io::policy::BufPolicy + Send,
{
    type DataSet = InterleavedRecordSet;
    type Err = fastq::Error;

    fn fill_data(
        &mut self,
        record: &mut Self::DataSet,
    ) -> Option<std::result::Result<(), Self::Err>> {
        // self.rset_len will be None when fill_data is called for the first time, and from then
        // on it will have a value.
        let result = self
            .reader
            .read_record_set_exact(&mut record.set, self.rset_len);
        if let Some(Ok(())) = result {
            if self.rset_len.is_none() {
                self.rset_len = Some(record.set.len());
            }
            if record.set.len() % 2 != 0 {
                return Some(Err(fastq::Error::Io(std::io::Error::other(
                    "FASTQ file does not have an even number of records.",
                ))));
            }
        }
        result
    }
}

/// Thin wrapper around a ``RecordSet`` that contains an even number of records, with interleaved
/// read pairs.
#[derive(Default, Clone, Debug, Serialize, Deserialize)]
pub struct InterleavedRecordSet {
    set: RecordSet,
}

#[allow(clippy::into_iter_without_iter)]
impl<'a> std::iter::IntoIterator for &'a InterleavedRecordSet {
    type Item = (RefRecord<'a>, RefRecord<'a>);
    type IntoIter = PairedRecordSetIterator<'a>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        let iter = self.set.into_iter().tuples();
        PairedRecordSetIterator {
            iter: Box::new(iter),
        }
    }
}

impl InterleavedRecordSet {}

parallel_record_impl!(
    parallel_paired_fastq,
    parallel_paired_fastq_init,
    R,
    PairedFastqReader<R>,
    RecordSetTuple,
    (fastq::RefRecord, fastq::RefRecord),
    fastq::Error
);

/// Reader for paired ends that are spread across two FASTQ files.
/// The number of records read from each file is synchronized, but the number is not known until
/// after reading the first time, when we determine how many reads can fit into the fixed size
/// buffer.  While the buffer may grow later for each `RecordSet` later (if later reads are
/// longer), this is a good approximation given the buffer size.
pub struct PairedFastqReader<R: std::io::Read, P = seq_io::policy::StdPolicy> {
    reader1: fastq::Reader<R, P>,
    reader2: fastq::Reader<R, P>,
    rset_len: Option<usize>,
}

impl<R, P> PairedFastqReader<R, P>
where
    R: std::io::Read,
    P: seq_io::policy::BufPolicy + Send,
{
    pub fn new(reader1: fastq::Reader<R, P>, reader2: fastq::Reader<R, P>) -> Self {
        Self {
            reader1,
            reader2,
            rset_len: None,
        }
    }
}

impl<R, P> seq_io::parallel::Reader for PairedFastqReader<R, P>
where
    R: std::io::Read,
    P: seq_io::policy::BufPolicy + Send,
{
    type DataSet = RecordSetTuple;
    type Err = fastq::Error;

    fn fill_data(
        &mut self,
        record: &mut Self::DataSet,
    ) -> Option<std::result::Result<(), Self::Err>> {
        // self.rset_len will be None when fill_data is called for the first time, and from then
        // on it will have a value.
        let result1 = self
            .reader1
            .read_record_set_exact(&mut record.first, self.rset_len);
        if result1.is_some() {
            self.rset_len = Some(record.first.len());
        }
        let result2 = self
            .reader2
            .read_record_set_exact(&mut record.second, self.rset_len);

        match (result1, result2) {
            (None, None) => None,
            (Some(Ok(())), Some(Ok(()))) => {
                if record.first.len() == record.second.len() {
                    Some(Ok(()))
                } else {
                    let head1: String = record.first.into_iter().last().map_or_else(
                        || "No more records".to_string(),
                        |r| String::from_utf8_lossy(r.head()).to_string(),
                    );
                    let head2 = record.second.into_iter().last().map_or_else(
                        || "No more records".to_string(),
                        |r| String::from_utf8_lossy(r.head()).to_string(),
                    );
                    Some(Err(fastq::Error::Io(std::io::Error::other(format!(
                        "FASTQ files out of sync.  Last records:\n\t{head1}\n\t{head2}"
                    )))))
                }
            }
            (_, Some(Err(e))) | (Some(Err(e)), _) => Some(Err(e)),
            (None, _) => Some(Err(fastq::Error::UnexpectedEnd {
                pos: fastq::ErrorPosition {
                    line: self.reader2.position().line(),
                    id: None,
                },
            })),
            (_, None) => Some(Err(fastq::Error::UnexpectedEnd {
                pos: fastq::ErrorPosition {
                    line: self.reader1.position().line(),
                    id: None,
                },
            })),
        }
    }
}

/// Iterator over paired end reads from two FASTQ files.
pub struct PairedRecordSetIterator<'a> {
    iter: Box<dyn Iterator<Item = (RefRecord<'a>, RefRecord<'a>)> + 'a>,
}

impl<'a> Iterator for PairedRecordSetIterator<'a> {
    type Item = (RefRecord<'a>, RefRecord<'a>);
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

/// Stores two record sets with the same number of records for paired end reads.  Each record set
/// is used to read from the corresponding reader, one per end of a pair.
#[derive(Default, Clone, Debug, Serialize, Deserialize)]
pub struct RecordSetTuple {
    first: RecordSet,
    second: RecordSet,
}

#[allow(clippy::into_iter_without_iter)]
impl<'a> std::iter::IntoIterator for &'a RecordSetTuple {
    type Item = (RefRecord<'a>, RefRecord<'a>);
    type IntoIter = PairedRecordSetIterator<'a>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        #[allow(clippy::useless_conversion)]
        let iter = self.first.into_iter().zip(self.second.into_iter());
        PairedRecordSetIterator {
            iter: Box::new(iter),
        }
    }
}

/// A data set wrapper that carries a dispatch sequence number for ordered output.
pub(crate) struct OrderedDataSet<D> {
    inner: D,
    seq: u64,
}

impl<D: Default> Default for OrderedDataSet<D> {
    fn default() -> Self {
        Self {
            inner: D::default(),
            seq: 0,
        }
    }
}

/// Wraps a reader to assign monotonically increasing sequence numbers to each batch.
/// Sequence numbers are assigned in `fill_data`, which runs sequentially in the reader thread.
pub(crate) struct OrderedReader<R: seq_io::parallel::Reader> {
    inner: R,
    next_seq: u64,
}

impl<R: seq_io::parallel::Reader> OrderedReader<R> {
    pub fn new(inner: R) -> Self {
        Self { inner, next_seq: 0 }
    }
}

impl<R: seq_io::parallel::Reader> seq_io::parallel::Reader for OrderedReader<R>
where
    R::DataSet: Send,
{
    type DataSet = OrderedDataSet<R::DataSet>;
    type Err = R::Err;

    fn fill_data(
        &mut self,
        data: &mut Self::DataSet,
    ) -> Option<std::result::Result<(), Self::Err>> {
        let result = self.inner.fill_data(&mut data.inner);
        if let Some(Ok(())) = &result {
            data.seq = self.next_seq;
            self.next_seq += 1;
        }
        result
    }
}

/// Generates an ordered parallel processing function that preserves input order in output.
/// Uses `OrderedReader` to tag batches with sequence numbers and a `BTreeMap` reorder buffer
/// to yield results in dispatch order rather than completion order.
macro_rules! ordered_parallel_record_impl {
    ($name:ident, $reader:ty, $dataset:ty, $record:ty, $err:ty) => {
        /// Processes records in parallel while preserving input order in the output.
        /// When a batch arrives out of order, its record set and per-record data are moved
        /// into a reorder buffer via `mem::take` and yielded later in the correct sequence.
        pub fn $name<R, D, W, F, Out>(
            reader: $reader,
            n_threads: u32,
            queue_len: usize,
            work: W,
            mut func: F,
        ) -> std::result::Result<Option<Out>, $err>
        where
            R: io::Read + Send,
            D: Default + Send,
            W: Send + Sync + Fn($record, &mut D),
            F: FnMut($record, &mut D) -> Option<Out>,
        {
            read_parallel_init::<_, $err, _, _, _, _, $err, _, _, _>(
                n_threads,
                queue_len,
                move || {
                    Ok::<_, $err>(ReusableReader::<OrderedReader<$reader>, (Vec<D>, ())>::new(
                        OrderedReader::new(reader),
                    ))
                },
                || Ok::<_, $err>((OrderedDataSet::<$dataset>::default(), (Vec::<D>::new(), ()))),
                |data: &mut (OrderedDataSet<$dataset>, (Vec<D>, ()))| {
                    let recordset = &data.0.inner;
                    let out = &mut (data.1).0;
                    let mut record_iter = recordset.into_iter();
                    for (d, record) in out.iter_mut().zip(&mut record_iter) {
                        work(record, d);
                    }
                    for record in record_iter {
                        out.push(D::default());
                        work(record, out.last_mut().unwrap());
                    }
                    Ok::<_, $err>(())
                },
                |records| {
                    let mut next_seq: u64 = 0;
                    let mut buffer: BTreeMap<u64, ($dataset, Vec<D>)> = BTreeMap::new();

                    while let Some(result) = records.next() {
                        let (data, work_result) = result?;
                        work_result?;
                        let seq = data.0.seq;

                        if seq == next_seq {
                            for (record, d) in
                                (&data.0.inner).into_iter().zip((data.1).0.iter_mut())
                            {
                                if let Some(out) = func(record, d) {
                                    return Ok(Some(out));
                                }
                            }
                            next_seq += 1;
                            while let Some((rset, mut outs)) = buffer.remove(&next_seq) {
                                for (record, d) in (&rset).into_iter().zip(outs.iter_mut()) {
                                    if let Some(out) = func(record, d) {
                                        return Ok(Some(out));
                                    }
                                }
                                next_seq += 1;
                            }
                        } else {
                            buffer.insert(
                                seq,
                                (
                                    std::mem::take(&mut data.0.inner),
                                    std::mem::take(&mut (data.1).0),
                                ),
                            );
                        }
                    }
                    Ok(None)
                },
            )?
        }
    };
}

ordered_parallel_record_impl!(
    ordered_parallel_fastq,
    fastq::Reader<R>,
    RecordSet,
    fastq::RefRecord,
    fastq::Error
);

ordered_parallel_record_impl!(
    ordered_parallel_interleaved_fastq,
    InterleavedFastqReader<R>,
    InterleavedRecordSet,
    (fastq::RefRecord, fastq::RefRecord),
    fastq::Error
);

ordered_parallel_record_impl!(
    ordered_parallel_paired_fastq,
    PairedFastqReader<R>,
    RecordSetTuple,
    (fastq::RefRecord, fastq::RefRecord),
    fastq::Error
);
