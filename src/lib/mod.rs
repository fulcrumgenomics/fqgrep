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
