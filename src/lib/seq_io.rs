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
use seq_io::parallel_record_impl;
use serde::{Deserialize, Serialize};
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
        let result = self.reader.read_record_set_limited(&mut record.set, 65536);
        if let Some(Ok(())) = result {
            if record.set.len() % 2 != 0 {
                return Some(Err(fastq::Error::Io(std::io::Error::other(
                    "FASTQ file does not have an even number of records.",
                ))));
            }
        }
        result
    }
}

/// Thin wrapper around a RecordSet that contains an even number of records, with interleaved read
/// pairs.
#[derive(Default, Clone, Debug, Serialize, Deserialize)]
pub struct InterleavedRecordSet {
    set: RecordSet,
}

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
pub struct PairedFastqReader<R: std::io::Read, P = seq_io::policy::StdPolicy> {
    pub reader1: fastq::Reader<R, P>,
    pub reader2: fastq::Reader<R, P>,
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
        let result1 = self
            .reader1
            .read_record_set_limited(&mut record.first, 65536);
        let result2 = self
            .reader2
            .read_record_set_limited(&mut record.second, 65536);

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
                    let head2 = record.first.into_iter().last().map_or_else(
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
