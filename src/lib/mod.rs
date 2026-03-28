#![deny(unsafe_code)]
#![allow(
    clippy::must_use_candidate,
    clippy::missing_panics_doc,
    clippy::missing_errors_doc,
    clippy::module_name_repetitions
)]
pub mod color;
pub mod matcher;
pub mod seq_io;
use bio::data_structures::bitenc::BitEnc;
use std::{borrow::Borrow, path::Path, sync::LazyLock};

/// Lookup table mapping ASCII byte values to their 4-bit IUPAC bitmask encodings.
/// A=1, C=2, G=4, T=8, and ambiguity codes are the bitwise OR of their constituents.
pub const IUPAC_MASKS: [u8; 256] = {
    let mut masks = [0u8; 256];
    let (a, c, g, t) = (1, 2, 4, 8);
    masks[b'A' as usize] = a;
    masks[b'C' as usize] = c;
    masks[b'G' as usize] = g;
    masks[b'T' as usize] = t;
    masks[b'U' as usize] = t;
    masks[b'M' as usize] = a | c;
    masks[b'R' as usize] = a | g;
    masks[b'W' as usize] = a | t;
    masks[b'S' as usize] = c | g;
    masks[b'Y' as usize] = c | t;
    masks[b'K' as usize] = g | t;
    masks[b'V' as usize] = a | c | g;
    masks[b'H' as usize] = a | c | t;
    masks[b'D' as usize] = a | g | t;
    masks[b'B' as usize] = c | g | t;
    masks[b'N' as usize] = a | c | g | t;
    // Lowercase mappings
    masks[b'a' as usize] = a;
    masks[b'c' as usize] = c;
    masks[b'g' as usize] = g;
    masks[b't' as usize] = t;
    masks[b'u' as usize] = t;
    masks[b'm' as usize] = a | c;
    masks[b'r' as usize] = a | g;
    masks[b'w' as usize] = a | t;
    masks[b's' as usize] = c | g;
    masks[b'y' as usize] = c | t;
    masks[b'k' as usize] = g | t;
    masks[b'v' as usize] = a | c | g;
    masks[b'h' as usize] = a | c | t;
    masks[b'd' as usize] = a | g | t;
    masks[b'b' as usize] = c | g | t;
    masks[b'n' as usize] = a | c | g | t;
    masks
};

pub const DNA_BASES: [u8; 5] = *b"ACGTN";
pub const AMINO_ACIDS: [u8; 26] = *b"ARNDCEQGHILKMFPSTWYVBJOUXZ";
pub const IUPAC_BASES: [u8; 15] = *b"AGCTYRWSKMDVHBN";
pub const IUPAC_BASES_COMPLEMENT: [u8; 15] = *b"TCGARYWSMKHBDVN";

/// The encoded 4-bit mask values for the four standard DNA bases (A=1, C=2, G=4, T=8).
/// Used to distinguish DNA-only patterns from those containing IUPAC ambiguity codes.
pub const DNA_MASK_VALUES: [u8; 4] = [1, 2, 4, 8];

/// The maximum number of expanded patterns allowed before returning an error.
pub const MAX_IUPAC_EXPANSIONS: usize = 10_000;

/// Expands the pattern containing IUPAC bases into one or more patterns.
///
/// Returns an error if the expansion would exceed [`MAX_IUPAC_EXPANSIONS`] patterns.
pub fn expand_iupac_fixed_pattern(
    pattern: &[u8],
    pattern_index: usize,
    prefix: &mut Vec<u8>,
    expanded: &mut Vec<Vec<u8>>,
) -> Result<(), String> {
    if expanded.len() >= MAX_IUPAC_EXPANSIONS {
        return Err(format!(
            "IUPAC expansion exceeded {MAX_IUPAC_EXPANSIONS} patterns for input '{}'",
            String::from_utf8_lossy(pattern)
        ));
    }
    if pattern_index == pattern.len() {
        expanded.push(prefix.clone());
    } else {
        let mask = IUPAC_MASKS[pattern[pattern_index] as usize];
        for base in DNA_BASES.iter().take(4) {
            if (IUPAC_MASKS[*base as usize] & mask) != 0 {
                prefix.push(*base);
                expand_iupac_fixed_pattern(pattern, pattern_index + 1, prefix, expanded)?;
                prefix.pop();
            }
        }
    }
    Ok(())
}

/// Converts an IUPAC pattern into a regex pattern by replacing each ambiguity code with a
/// character class of its constituent DNA bases (e.g., `R` becomes `[AG]`).
pub fn expand_iupac_regex(pattern: &[u8]) -> Vec<u8> {
    let mut new_pattern: Vec<u8> = Vec::new();
    for base in pattern {
        let mask = IUPAC_MASKS[*base as usize];
        let mut bases = Vec::new();
        for dna_base in DNA_BASES.iter().take(4) {
            if (IUPAC_MASKS[*dna_base as usize] & mask) != 0 {
                bases.push(*dna_base);
            }
        }
        if bases.len() == 1 {
            new_pattern.push(bases[0]);
        } else {
            new_pattern.push(b'[');
            new_pattern.extend_from_slice(&bases);
            new_pattern.push(b']');
        }
    }
    new_pattern
}

/// Encodes a sequence of ASCII bases into a `BitEnc` using 4-bit IUPAC mask values.
pub fn encode(bases: &[u8]) -> BitEnc {
    let mut vec = BitEnc::with_capacity(4, bases.len());
    for base in bases {
        vec.push(IUPAC_MASKS[*base as usize]);
    }
    vec
}

static COMPLEMENT: LazyLock<[u8; 256]> = LazyLock::new(|| {
    let mut comp = [0; 256];
    for (v, a) in comp.iter_mut().enumerate() {
        *a = v as u8;
    }
    for (&a, &b) in IUPAC_BASES.iter().zip(IUPAC_BASES_COMPLEMENT.iter()) {
        comp[a as usize] = b;
        comp[a as usize + 32] = b + 32; // lowercase variants
    }
    comp
});

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

    // ############################################################################################
    // Tests expand_iupac_fixed_pattern()
    // ############################################################################################

    #[rstest]
    #[case("GATTACA", vec!["GATTACA"])]
    #[case("GATMACA", vec!["GATAACA", "GATCACA"])]
    #[case("GRTBANA", vec!["GATCAAA", "GATCACA", "GATCAGA", "GATCATA", "GATGAAA", "GATGACA", "GATGAGA", "GATGATA", "GATTAAA", "GATTACA", "GATTAGA", "GATTATA", "GGTCAAA", "GGTCACA", "GGTCAGA", "GGTCATA", "GGTGAAA", "GGTGACA", "GGTGAGA", "GGTGATA", "GGTTAAA", "GGTTACA", "GGTTAGA", "GGTTATA"])]
    fn test_expand_iupac_fixed_pattern(#[case] pattern: &str, #[case] expected: Vec<&str>) {
        let mut actual = Vec::new();
        expand_iupac_fixed_pattern(pattern.as_bytes(), 0, &mut Vec::new(), &mut actual).unwrap();
        let actual_strings: Vec<String> = actual
            .iter()
            .map(|bytes| str::from_utf8(bytes).unwrap().to_string())
            .collect();
        assert_eq!(actual_strings, expected);
    }

    // ############################################################################################
    // Tests expand_iupac_regex()
    // ############################################################################################

    #[rstest]
    #[case("GATTACA", "GATTACA")]
    #[case("GATMACA", "GAT[AC]ACA")]
    #[case("GRTBANA", "G[AG]T[CGT]A[ACGT]A")]
    #[case("GATUCA", "GATTCA")]
    fn test_expand_iupac_regex(#[case] pattern: &str, #[case] expected: &str) {
        let expanded = expand_iupac_regex(pattern.as_bytes());
        let actual = str::from_utf8(&expanded).unwrap();
        assert_eq!(actual, expected);
    }

    // ############################################################################################
    // Tests encode()
    // ############################################################################################

    #[rstest]
    fn test_encode() {
        for base in IUPAC_BASES {
            let actual: u8 = encode(&[base]).get(0).unwrap();
            assert_eq!(actual, IUPAC_MASKS[base as usize]);
        }
    }

    // ############################################################################################
    // Tests lowercase IUPAC masks
    // ############################################################################################

    #[rstest]
    fn test_iupac_masks_lowercase() {
        for &base in b"ACGTUMRWSYKVHDBNacgtumrwsykvhdbn" {
            let upper = (base as char).to_ascii_uppercase() as u8;
            assert_eq!(
                IUPAC_MASKS[base as usize], IUPAC_MASKS[upper as usize],
                "Lowercase '{}' should have the same mask as uppercase '{}'",
                base as char, upper as char
            );
        }
    }

    #[rstest]
    fn test_encode_lowercase() {
        for &base in b"acgt" {
            let upper = (base as char).to_ascii_uppercase() as u8;
            let lower_enc: u8 = encode(&[base]).get(0).unwrap();
            let upper_enc: u8 = encode(&[upper]).get(0).unwrap();
            assert_eq!(lower_enc, upper_enc);
        }
    }

    // ############################################################################################
    // Tests expand_iupac_fixed_pattern expansion limit
    // ############################################################################################

    #[rstest]
    fn test_expand_iupac_fixed_pattern_exceeds_limit() {
        // 'N' expands to 4 bases; N repeated 8 times = 4^8 = 65536 > MAX_IUPAC_EXPANSIONS
        let pattern = b"NNNNNNNN";
        let mut expanded = Vec::new();
        let result = expand_iupac_fixed_pattern(pattern, 0, &mut Vec::new(), &mut expanded);
        assert!(
            result.is_err(),
            "Should error when expansion exceeds the limit"
        );
        assert!(
            result
                .unwrap_err()
                .contains(&MAX_IUPAC_EXPANSIONS.to_string()),
            "Error message should mention the limit"
        );
    }
}
