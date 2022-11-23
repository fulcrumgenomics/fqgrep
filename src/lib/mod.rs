#![deny(unsafe_code)]
#![allow(
    clippy::must_use_candidate,
    clippy::missing_panics_doc,
    clippy::missing_errors_doc,
    clippy::module_name_repetitions
)]
pub mod color;
pub mod matcher;
use lazy_static::lazy_static;
use std::{borrow::Borrow, path::Path};

pub const DNA_BASES: [u8; 5] = *b"ACGTN";
pub const IUPAC_BASES: [u8; 15] = *b"AGCTYRWSKMDVHBN";
pub const IUPAC_BASES_COMPLEMENT: [u8; 15] = *b"TCGARYWSMKHBDVN";

lazy_static! {
    pub static ref COMPLEMENT: [u8; 256] = {
        let mut comp = [0; 256];
        for (v, a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in IUPAC_BASES.iter().zip(IUPAC_BASES_COMPLEMENT.iter()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        comp
    };
}

fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

fn reverse_complement<C, T>(text: T) -> Vec<u8>
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

/// Returns true if the path ends with a recognized GZIP file extension
fn is_path_with_extension<P: AsRef<Path>>(p: &P, extensions: [&str; 2]) -> bool {
    if let Some(ext) = p.as_ref().extension() {
        match ext.to_str() {
            Some(x) => extensions.contains(&x),
            None => false,
        }
    } else {
        false
    }
}

/// The set of file extensions to treat as GZIPPED
const GZIP_EXTENSIONS: [&str; 2] = ["gz", "bgz"];

/// Returns true if the path ends with a recognized GZIP file extension
pub fn is_gzip_path<P: AsRef<Path>>(p: &P) -> bool {
    is_path_with_extension(p, GZIP_EXTENSIONS)
}

/// The set of file extensions to treat as FASTQ
const FASTQ_EXTENSIONS: [&str; 2] = ["fastq", "fq"];

/// Returns true if the path ends with a recognized FASTQ file extension
pub fn is_fastq_path<P: AsRef<Path>>(p: &P) -> bool {
    is_path_with_extension(p, FASTQ_EXTENSIONS)
}

// Tests
#[cfg(test)]
pub mod tests {
    use crate::*;
    use rstest::rstest;
    use seq_io::fastq::OwnedRecord;
    use std::str;
    use tempfile::TempDir;

    /// Helper function takes a sequence and returns a seq_io::fastq::OwnedRecord
    ///
    fn write_owned_record(seq: &str) -> OwnedRecord {
        let read = OwnedRecord {
            head: ("@Sample").as_bytes().to_vec(),
            seq: seq.as_bytes().to_vec(),
            qual: vec![b'X'; seq.len()],
        };
        read
    }

    // ############################################################################################
    // Tests reverse_complement()
    // ############################################################################################

    #[test]
    fn test_reverse_complement() {
        let read = write_owned_record("ACTG");
        let result = reverse_complement(read.seq);
        let string_result = str::from_utf8(&result).unwrap();

        // Correct
        assert_eq!(&string_result, &"CAGT");

        // Incorrect
        assert_ne!(&string_result, &"TGAC");
    }

    // ############################################################################################
    // Tests is_gzip_path()
    // ############################################################################################

    #[rstest]
    #[case("test_fastq.fq.gz", true)] // .fq.gz is valid gzip
    #[case("test_fastq.fq.bgz", true)] // .fq.bgz is valid gzip
    #[case("test_fastq.fq.tar", false)] // .fq.tar is invalid gzip
    fn test_is_gzip_path(#[case] file_name: &str, #[case] expected: bool) {
        let dir = TempDir::new().unwrap();
        let file_path = dir.path().join(file_name);
        let result = is_gzip_path(&file_path);
        assert_eq!(result, expected);
    }
    // ############################################################################################
    // Tests is_fastq_path()
    // ############################################################################################

    #[rstest]
    #[case("test_fastq.fq", true)] // .fq is valid fastq
    #[case("test_fastq.fastq", true)] // .fastq is valid fastq
    #[case("test_fastq.sam", false)] // .sam is invalid fastq
    fn test_is_fastq_path(#[case] file_name: &str, #[case] expected: bool) {
        let dir = TempDir::new().unwrap();
        let file_path = dir.path().join(file_name);
        let result = is_fastq_path(&file_path);
        assert_eq!(result, expected);
    }
}
