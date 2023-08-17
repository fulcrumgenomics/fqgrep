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
use std::{borrow::Borrow};

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


// Tests
#[cfg(test)]
pub mod tests {
    use crate::*;
    use rstest::rstest;
    use std::str;

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
}
