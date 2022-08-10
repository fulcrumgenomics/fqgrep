#![deny(unsafe_code)]
#![allow(
    clippy::must_use_candidate,
    clippy::missing_panics_doc,
    clippy::missing_errors_doc,
    clippy::module_name_repetitions
)]
pub mod matcher;
use lazy_static::lazy_static;
use std::borrow::Borrow;

lazy_static! {
    pub static ref DNA_BASES: [u8; 5] = {
        let mut bases = [0; 5];
        for (i, base) in b"ACGTN".iter().enumerate() {
            bases[i] = *base;
        }
        bases
    };

    pub static ref IUPAC_BASES: [u8;15] = {
        let mut bases = [0; 15];
        for (i, base) in b"AGCTYRWSKMDVHBN".iter().enumerate() {
            bases[i] = *base;
        }
        bases
    };

    pub static ref IUPAC_BASES_COMPLEMENT: [u8;15] = {
        let mut bases = [0; 15];
        for (i, base) in b"TCGARYWSMKHBDVN".iter().enumerate() {
            bases[i] = *base;
        }
        bases
    };

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
