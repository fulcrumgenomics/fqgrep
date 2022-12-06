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
    use std::str;
    use tempfile::TempDir;

    // ############################################################################################
    // Tests reverse_complement()
    // ############################################################################################

    #[test]
    #[rstest]
    #[case("ACGT", "ACGT")] // Reverse complement with even length string
    #[case("ACG", "CGT")] // Reverse complement with odd length string (tests for off by one error)
    fn test_reverse_complement(#[case] seq: &str, #[case] expected: &str) {
        let result = reverse_complement(seq.as_bytes());
        let string_result = str::from_utf8(&result).unwrap();
        assert_eq!(&string_result, &expected);
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
