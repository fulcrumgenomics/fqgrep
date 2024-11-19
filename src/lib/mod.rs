#![deny(unsafe_code)]
#![allow(
    clippy::must_use_candidate,
    clippy::missing_panics_doc,
    clippy::missing_errors_doc,
    clippy::module_name_repetitions
)]
use bio::data_structures::bitenc::BitEnc;

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

    pub static ref BASE_A: usize = 1;
    pub static ref BASE_C: usize = 2;
    pub static ref BASE_G: usize = 4;
    pub static ref BASE_T: usize = 8;

    pub static ref IUPAC_MASKS: [u8; 256] = {
        let mut masks = [0; 256];
        let (a, c, g, t) = (1, 2, 4, 8);
        masks['A' as usize] = a;
        masks['C' as usize] = c;
        masks['G' as usize] = g;
        masks['T' as usize] = t;
        masks['U' as usize] = t;
        masks['M' as usize] = a | c;
        masks['R' as usize] = a | g;
        masks['W' as usize] = a | t;
        masks['S' as usize] = c | g;
        masks['Y' as usize] = c | t;
        masks['K' as usize] = g | t;
        masks['V' as usize] = a | c | g;
        masks['H' as usize] = a | c | t;
        masks['D' as usize] = a | g | t;
        masks['B' as usize] = c | g | t;
        masks['N' as usize] = a | c | g | t;
        masks
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

/// Expands the pattern containing IUPAC bases into one or more patterns.   For example,
/// GATK will be expanded to two fixed patterns GATG and GATT.
pub fn expand_iupac_fixed_pattern(
    pattern: &[u8],
    pattern_index: usize,
    prefix: &[u8],
    expanded: &mut Vec<Vec<u8>>,
) {
    if pattern_index == pattern.len() {
        // no more bases in the pattern, append the prefix, and recurse
        expanded.push(prefix.to_vec());
    } else {
        let mask = IUPAC_MASKS[pattern[pattern_index] as usize];
        for base in DNA_BASES.iter().take(4) {
            if (IUPAC_MASKS[*base as usize] & mask) != 0 {
                let new_prefix = [prefix, &[*base]].concat();
                expand_iupac_fixed_pattern(pattern, pattern_index + 1, &new_prefix, expanded);
            }
        }
    }
}

pub fn expand_iupac_regex(pattern: &[u8]) -> Vec<u8> {
    let mut new_pattern: Vec<u8> = Vec::new();
    for base in pattern {
        let mask = IUPAC_MASKS[*base as usize];
        let mut bases = Vec::new();
        for base in DNA_BASES.iter().take(4) {
            if (IUPAC_MASKS[*base as usize] & mask) != 0 {
                bases.push(*base);
            }
        }
        if bases.len() == 1 {
            new_pattern.push(*base);
        } else {
            new_pattern.push(b'[');
            new_pattern.extend_from_slice(&bases);
            new_pattern.push(b']');
        }
    }
    new_pattern
}

pub fn encode(bases: &[u8]) -> BitEnc {
    let mut vec = BitEnc::with_capacity(4, bases.len());
    for base in bases {
        vec.push(IUPAC_MASKS[*base as usize]);
    }
    vec
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

    #[rstest]
    #[case("GATTACA", vec!["GATTACA"])]
    #[case("GATMACA", vec!["GATAACA", "GATCACA"])]
    #[case("GRTBANA", vec!["GATCAAA", "GATCACA", "GATCAGA", "GATCATA", "GATGAAA", "GATGACA", "GATGAGA", "GATGATA", "GATTAAA", "GATTACA", "GATTAGA", "GATTATA", "GGTCAAA", "GGTCACA", "GGTCAGA", "GGTCATA", "GGTGAAA", "GGTGACA", "GGTGAGA", "GGTGATA", "GGTTAAA", "GGTTACA", "GGTTAGA", "GGTTATA"])]
    fn test_expand_iupac_fixed_pattern(#[case] pattern: &str, #[case] expected: Vec<&str>) {
        let mut actual = Vec::new();
        expand_iupac_fixed_pattern(pattern.as_bytes(), 0, &[], &mut actual);
        let actual_strings: Vec<String> = actual
            .iter()
            .map(|bytes| str::from_utf8(bytes).unwrap().to_string())
            .collect();
        assert_eq!(actual_strings, expected);
    }

    #[rstest]
    #[case("GATTACA", "GATTACA")]
    #[case("GATMACA", "GAT[AC]ACA")]
    #[case("GRTBANA", "G[AG]T[CGT]A[ACGT]A")]
    fn test_expand_iupac_regex(#[case] pattern: &str, #[case] expected: &str) {
        let expanded = expand_iupac_regex(pattern.as_bytes());
        let actual = str::from_utf8(&expanded).unwrap();
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_encode() {
        for base in IUPAC_BASES {
            let actual: u8 = encode(&[base]).get(0).unwrap();
            assert_eq!(actual, IUPAC_MASKS[base as usize]);
        }
    }
}
