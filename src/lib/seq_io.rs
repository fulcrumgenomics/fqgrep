/// Additions to `seq_io::parallel` for interleaved paired end FASTQ files
use anyhow::Result;
use itertools::{self, Itertools};
use seq_io::fastq::{self, RecordSet, RefRecord};
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
        self.reader.read_record_set_limited(&mut record.set, 65536)
    }
}

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
        // TODO: check even #

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
            (Some(Ok(())), Some(Ok(()))) => Some(Ok(())),
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
        // TODO: check the same # of reads
    }
}

pub struct PairedRecordSetIterator<'a> {
    iter: Box<dyn Iterator<Item = (RefRecord<'a>, RefRecord<'a>)> + 'a>,
}

impl<'a> Iterator for PairedRecordSetIterator<'a> {
    type Item = (RefRecord<'a>, RefRecord<'a>);
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

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
