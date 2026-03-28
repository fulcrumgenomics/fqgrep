use ahash::AHashSet;
use aho_corasick::AhoCorasick;
use bio::data_structures::bitenc::BitEnc;
use bitvec::prelude::*;
use seq_io::fastq::RefRecord;
use std::ops::Range;

use crate::AMINO_ACIDS;
use crate::DNA_BASES;
use crate::DNA_MASK_VALUES;
use crate::color::{COLOR_BACKGROUND, COLOR_BASES, COLOR_QUALS};
use crate::color::{color_background, color_head};
use crate::encode;
use crate::reverse_complement;
use anyhow::{Context, Result, bail};
use bstr::ByteSlice;
use regex::bytes::{Regex, RegexBuilder, RegexSet, RegexSetBuilder};
use seq_io::fastq::{OwnedRecord, Record};

/// Common options for pattern matchers
#[derive(Copy, Clone, Debug)]
pub struct MatcherOpts {
    /// Invert the matching. The bases are said to match if they do not contain the pattern
    pub invert_match: bool,
    /// Include the reverse complement of the bases.  The bases are said to match if they or their
    /// reverse complement contains the pattern.
    pub reverse_complement: bool,
}

/// Builds a bit vector from an iterator of ranges.  The ranges may overlap each other.
fn to_bitvec(ranges: impl Iterator<Item = Range<usize>>, len: usize) -> BitVec {
    let mut vec = bitvec![0; len];
    ranges.for_each(|range| {
        for index in range {
            vec.set(index, true);
        }
    });
    vec
}

/// Color the bases and qualities based on ranges specifying where a pattern matched bases.  The
/// ranges may overlap each other.
fn bases_colored(
    bases: &[u8],
    quals: &[u8],
    ranges: impl Iterator<Item = Range<usize>>,
) -> (Vec<u8>, Vec<u8>) {
    // The resulting colored bases
    let mut colored_bases = Vec::with_capacity(bases.len());
    let mut colored_quals = Vec::with_capacity(bases.len());

    // Merge the ranges into a bit mask, with 1 indicating that base is part of a pattern match
    let bits = to_bitvec(ranges, bases.len());

    // Iterate over the bit mask, finding stretches of matching and non-matching bases.  Color both
    // in both the bases and qualities
    let mut last_color_on = false;
    let mut last_bases_index = 0;
    let mut cur_bases_index = 0;
    for base_color_on in bits.iter() {
        if *base_color_on {
            // this base is to be colored
            if !last_color_on {
                // add up to but not including this base to the colored vector **as uncolored**
                if last_bases_index < cur_bases_index {
                    COLOR_BACKGROUND
                        .paint(&bases[last_bases_index..cur_bases_index])
                        .write_to(&mut colored_bases)
                        .unwrap();
                    COLOR_BACKGROUND
                        .paint(&quals[last_bases_index..cur_bases_index])
                        .write_to(&mut colored_quals)
                        .unwrap();
                }
                // first base in a run of bases to be colored
                last_bases_index = cur_bases_index;
            }

            last_color_on = true;
        } else {
            // this base is not to be colored
            if last_color_on {
                // add up to but not including this base to the colored vector **as colored**
                if last_bases_index < cur_bases_index {
                    COLOR_BASES
                        .paint(&bases[last_bases_index..cur_bases_index])
                        .write_to(&mut colored_bases)
                        .unwrap();
                    COLOR_QUALS
                        .paint(&quals[last_bases_index..cur_bases_index])
                        .write_to(&mut colored_quals)
                        .unwrap();
                }
                // first base in a run of bases to be colored
                last_bases_index = cur_bases_index;
            }
            last_color_on = false;
        }
        cur_bases_index += 1;
    }
    // Color to the end
    if last_bases_index < cur_bases_index {
        if last_color_on {
            COLOR_BASES
                .paint(&bases[last_bases_index..cur_bases_index])
                .write_to(&mut colored_bases)
                .unwrap();
            COLOR_QUALS
                .paint(&quals[last_bases_index..cur_bases_index])
                .write_to(&mut colored_quals)
                .unwrap();
        } else {
            COLOR_BACKGROUND
                .paint(&bases[last_bases_index..cur_bases_index])
                .write_to(&mut colored_bases)
                .unwrap();
            COLOR_BACKGROUND
                .paint(&quals[last_bases_index..cur_bases_index])
                .write_to(&mut colored_quals)
                .unwrap();
        }
    }

    (colored_bases, colored_quals)
}

/// Validates that a given FIXED pattern contains only valid bases (DNA or amino acid)
pub fn validate_fixed_pattern(pattern: &str, protein: bool) -> Result<()> {
    let valid_chars = if protein {
        &AMINO_ACIDS[..]
    } else {
        &DNA_BASES[..]
    };
    let kind = if protein { "amino acids" } else { "DNA bases" };
    for (index, base) in pattern.char_indices() {
        if !base.is_ascii() || !valid_chars.contains(&(base as u8)) {
            let next_index = index + base.len_utf8();
            bail!(
                "Fixed pattern must contain only {kind}: {} .. [{}] .. {}",
                &pattern[0..index],
                &pattern[index..next_index],
                &pattern[next_index..],
            )
        }
    }
    Ok(())
}

/// Base trait for all pattern matchers
pub trait Matcher {
    /// The options for the pattern matcher
    fn opts(&self) -> MatcherOpts;

    /// Returns true if the bases match the pattern, false otherwise
    fn bases_match(&self, bases: &[u8]) -> bool;

    /// Colors the bases and qualities based on where they match the pattern.  All bases
    /// that match the pattern are colored.  Colored in this case means adding ANSI color
    /// codes for printing to a terminal.
    fn color_matched_bases(&self, bases: &[u8], quals: &[u8]) -> (Vec<u8>, Vec<u8>);

    /// Returns true if the read's bases match the pattern, false otherwise
    #[inline]
    fn read_match(&self, read: &RefRecord) -> bool {
        let bases_match = self.bases_match(read.seq());
        if self.opts().invert_match {
            bases_match
                && (!self.opts().reverse_complement
                    || self.bases_match(&reverse_complement(read.seq())))
        } else {
            bases_match
                || (self.opts().reverse_complement
                    && self.bases_match(&reverse_complement(read.seq())))
        }
    }

    /// Adds ANSI color codes to the read's header, sequence, and quality based on where they
    /// match the pattern(s).
    #[inline]
    fn color(&self, read: &mut OwnedRecord, match_found: bool) {
        if match_found {
            let (seq, qual) = self.color_matched_bases(&read.seq, &read.qual);
            read.head = color_head(&read.head);
            read.seq = seq;
            read.qual = qual;
        } else {
            // always color, in case the read is paired
            read.head = color_background(&read.head);
            read.seq = color_background(&read.seq);
            read.qual = color_background(&read.qual);
        }
    }
}

/// Matcher for a fixed string pattern
pub struct FixedStringMatcher {
    pattern: Vec<u8>,
    opts: MatcherOpts,
}

impl Matcher for FixedStringMatcher {
    #[inline]
    fn bases_match(&self, bases: &[u8]) -> bool {
        bases.find(&self.pattern).is_some() != self.opts.invert_match
    }

    fn color_matched_bases(&self, bases: &[u8], quals: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let ranges = bases.find_iter(&self.pattern).map(|start| Range {
            start,
            end: start + self.pattern.len(),
        });
        if self.opts().reverse_complement {
            let bases_revcomp = &reverse_complement(bases);
            let ranges_revcomp = bases_revcomp
                .find_iter(&self.pattern)
                .map(|start| bases.len() - start - self.pattern.len())
                .map(|start| Range {
                    start,
                    end: start + self.pattern.len(),
                });
            bases_colored(bases, quals, ranges.chain(ranges_revcomp))
        } else {
            bases_colored(bases, quals, ranges)
        }
    }

    #[inline]
    fn opts(&self) -> MatcherOpts {
        self.opts
    }
}

impl FixedStringMatcher {
    pub fn new(pattern: &str, opts: MatcherOpts) -> Self {
        let pattern = pattern.as_bytes().to_vec();
        Self { pattern, opts }
    }
}

/// Matcher for a set of fixed string patterns
pub struct FixedStringSetMatcher {
    patterns: Vec<Vec<u8>>,
    aho_corasick: AhoCorasick,
    opts: MatcherOpts,
}

impl Matcher for FixedStringSetMatcher {
    #[inline]
    fn bases_match(&self, bases: &[u8]) -> bool {
        self.aho_corasick.is_match(bases) != self.opts.invert_match
    }

    fn color_matched_bases(&self, bases: &[u8], quals: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let ranges = self.patterns.iter().flat_map(|pattern| {
            bases
                .find_iter(&pattern)
                .map(|start| Range {
                    start,
                    end: start + pattern.len(),
                })
                .collect::<Vec<_>>()
        });
        if self.opts().reverse_complement {
            let bases_revcomp = &reverse_complement(bases);
            let ranges_revcomp = self.patterns.iter().flat_map(|pattern| {
                bases_revcomp
                    .find_iter(&pattern)
                    .map(|start| bases.len() - start - pattern.len())
                    .map(|start| Range {
                        start,
                        end: start + pattern.len(),
                    })
                    .collect::<Vec<_>>()
            });
            bases_colored(bases, quals, ranges.chain(ranges_revcomp))
        } else {
            bases_colored(bases, quals, ranges)
        }
    }

    #[inline]
    fn opts(&self) -> MatcherOpts {
        self.opts
    }
}

impl FixedStringSetMatcher {
    pub fn new<I, S>(patterns: I, opts: MatcherOpts) -> Result<Self>
    where
        S: AsRef<str>,
        I: IntoIterator<Item = S>,
    {
        let patterns: Vec<Vec<u8>> = patterns
            .into_iter()
            .map(|pattern| pattern.as_ref().to_owned().as_bytes().to_vec())
            .collect();
        let aho_corasick =
            AhoCorasick::new(&patterns).context("Failed to build Aho-Corasick automaton")?;
        Ok(Self {
            patterns,
            aho_corasick,
            opts,
        })
    }
}

/// Iterates over matching offsets where `needle` matches within `bases` using bitwise AND.
/// When `has_iupac` is false, uses a rolling additive hash to skip non-matching offsets.
/// Calls `on_match` for each matching offset; if it returns `false`, iteration stops early.
fn bitenc_match_offsets(
    bases: &BitEnc,
    needle: &BitEnc,
    has_iupac: bool,
    mut on_match: impl FnMut(usize) -> bool,
) {
    if bases.nr_symbols() < needle.nr_symbols() {
        return;
    }

    let mut bases_hash: usize = 0;
    let needle_hash = if has_iupac {
        0
    } else {
        (0..needle.nr_symbols())
            .map(|i| needle.get(i).unwrap() as usize)
            .sum()
    };

    'outer: for offset in 0..=bases.nr_symbols() - needle.nr_symbols() {
        if !has_iupac {
            if offset == 0 {
                bases_hash = (0..needle.nr_symbols())
                    .map(|i| bases.get(i).unwrap() as usize)
                    .sum();
            } else {
                bases_hash += bases.get(offset + needle.nr_symbols() - 1).unwrap() as usize;
                bases_hash -= bases.get(offset - 1).unwrap() as usize;
            }
            if bases_hash != needle_hash {
                continue 'outer;
            }
        }
        for i in 0..needle.nr_symbols() {
            let base = bases.get(offset + i).unwrap();
            if needle.get(i).unwrap() & base != base {
                continue 'outer;
            }
        }
        if !on_match(offset) {
            return;
        }
    }
}

/// Finds all positions where `needle` matches within `bases` using bitwise AND matching.
/// When `has_iupac` is false, uses a rolling additive hash to skip non-matching offsets.
fn bitenc_find_all(bases: &BitEnc, needle: &BitEnc, has_iupac: bool) -> Vec<Range<usize>> {
    let mut ranges: Vec<Range<usize>> = Vec::new();
    let needle_len = needle.nr_symbols();
    bitenc_match_offsets(bases, needle, has_iupac, |offset| {
        ranges.push(offset..offset + needle_len);
        true
    });
    ranges
}

/// Returns true if `needle` matches anywhere within `bases` using bitwise AND matching.
/// When `has_iupac` is false, uses a rolling additive hash to skip non-matching offsets.
fn bitenc_find(bases: &BitEnc, needle: &BitEnc, has_iupac: bool) -> bool {
    let mut found = false;
    bitenc_match_offsets(bases, needle, has_iupac, |_| {
        found = true;
        false // stop after first match
    });
    found
}

/// Returns true if the encoded pattern contains any IUPAC ambiguity codes
/// (i.e., any symbol that is not a single DNA base mask value).
fn bitenc_has_iupac(bitenc: &BitEnc) -> bool {
    !bitenc.iter().all(|value| DNA_MASK_VALUES.contains(&value))
}

/// Matcher for a single IUPAC bitmask pattern
pub struct BitMaskMatcher {
    bitenc: BitEnc,
    has_iupac: bool,
    opts: MatcherOpts,
}

impl Matcher for BitMaskMatcher {
    #[inline]
    fn bases_match(&self, bases: &[u8]) -> bool {
        let bases = encode(bases);
        bitenc_find(&bases, &self.bitenc, self.has_iupac) != self.opts.invert_match
    }

    fn color_matched_bases(&self, bases: &[u8], quals: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let encoded_bases = encode(bases);
        let ranges = bitenc_find_all(&encoded_bases, &self.bitenc, self.has_iupac);
        if self.opts().reverse_complement {
            let bases_revcomp = &reverse_complement(bases);
            let encoded_revcomp = encode(bases_revcomp);
            let ranges_revcomp = bitenc_find_all(&encoded_revcomp, &self.bitenc, self.has_iupac)
                .into_iter()
                .map(|range| Range {
                    start: bases.len() - range.start - range.len(),
                    end: bases.len() - range.start,
                });
            bases_colored(bases, quals, ranges.into_iter().chain(ranges_revcomp))
        } else {
            bases_colored(bases, quals, ranges.into_iter())
        }
    }

    #[inline]
    fn opts(&self) -> MatcherOpts {
        self.opts
    }
}

impl BitMaskMatcher {
    pub fn new(bitenc: BitEnc, opts: MatcherOpts) -> Self {
        let has_iupac = bitenc_has_iupac(&bitenc);
        Self {
            bitenc,
            has_iupac,
            opts,
        }
    }
}

/// Matcher for a set of IUPAC bitmask patterns
pub struct BitMaskSetMatcher {
    bitencs: Vec<BitEnc>,
    has_iupac: Vec<bool>,
    opts: MatcherOpts,
}

impl Matcher for BitMaskSetMatcher {
    fn bases_match(&self, bases: &[u8]) -> bool {
        let bases = encode(bases);
        self.bitencs
            .iter()
            .enumerate()
            .any(|(index, needle)| bitenc_find(&bases, needle, self.has_iupac[index]))
            != self.opts.invert_match
    }

    fn color_matched_bases(&self, bases: &[u8], quals: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let encoded_bases = encode(bases);
        let ranges = self.bitencs.iter().enumerate().flat_map(|(index, needle)| {
            bitenc_find_all(&encoded_bases, needle, self.has_iupac[index])
        });
        if self.opts().reverse_complement {
            let bases_revcomp = &reverse_complement(bases);
            let encoded_revcomp = encode(bases_revcomp);
            let ranges_revcomp = self.bitencs.iter().enumerate().flat_map(|(index, needle)| {
                bitenc_find_all(&encoded_revcomp, needle, self.has_iupac[index])
                    .into_iter()
                    .map(|range| Range {
                        start: bases.len() - range.start - range.len(),
                        end: bases.len() - range.start,
                    })
            });
            bases_colored(bases, quals, ranges.chain(ranges_revcomp))
        } else {
            bases_colored(bases, quals, ranges)
        }
    }

    #[inline]
    fn opts(&self) -> MatcherOpts {
        self.opts
    }
}

impl BitMaskSetMatcher {
    pub fn new(bitencs: Vec<BitEnc>, opts: MatcherOpts) -> Self {
        let has_iupac = bitencs.iter().map(bitenc_has_iupac).collect();
        Self {
            bitencs,
            has_iupac,
            opts,
        }
    }
}

/// Matcher for a regular expression pattern
pub struct RegexMatcher {
    regex: Regex,
    opts: MatcherOpts,
}

impl RegexMatcher {
    pub fn new(pattern: &str, opts: MatcherOpts) -> Result<Self> {
        let regex = RegexBuilder::new(pattern)
            .build()
            .with_context(|| format!("Invalid regular expression: '{}'", pattern))?;
        Ok(Self { regex, opts })
    }
}

impl Matcher for RegexMatcher {
    #[inline]
    fn bases_match(&self, bases: &[u8]) -> bool {
        self.regex.is_match(bases) != self.opts.invert_match
    }

    fn color_matched_bases(&self, bases: &[u8], quals: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let ranges = self.regex.find_iter(bases).map(|m| m.range());
        if self.opts().reverse_complement {
            let bases_revcomp = &reverse_complement(bases);
            let ranges_revcomp =
                self.regex
                    .find_iter(bases_revcomp)
                    .map(|m| m.range())
                    .map(|range| Range {
                        start: bases.len() - range.start - range.len(),
                        end: bases.len() - range.start,
                    });
            bases_colored(bases, quals, ranges.chain(ranges_revcomp))
        } else {
            bases_colored(bases, quals, ranges)
        }
    }

    #[inline]
    fn opts(&self) -> MatcherOpts {
        self.opts
    }
}

pub struct RegexSetMatcher {
    regex_set: RegexSet,
    regex_matchers: Vec<RegexMatcher>,
    opts: MatcherOpts,
}

/// Matcher for a set of regular expression patterns
impl RegexSetMatcher {
    pub fn new<I, S>(patterns: I, opts: MatcherOpts) -> Result<Self>
    where
        S: AsRef<str>,
        I: IntoIterator<Item = S>,
    {
        let string_patterns: Vec<String> = patterns
            .into_iter()
            .map(|p| p.as_ref().to_string())
            .collect();
        let regex_set = RegexSetBuilder::new(string_patterns.clone())
            .dfa_size_limit(usize::MAX)
            .build()
            .context("Failed to build regex set from patterns")?;
        let regex_matchers: Vec<RegexMatcher> = string_patterns
            .into_iter()
            .map(|pattern| RegexMatcher::new(pattern.as_ref(), opts))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            regex_set,
            regex_matchers,
            opts,
        })
    }
}

impl Matcher for RegexSetMatcher {
    #[inline]
    fn bases_match(&self, bases: &[u8]) -> bool {
        self.regex_set.is_match(bases) != self.opts.invert_match
    }

    fn color_matched_bases(&self, bases: &[u8], quals: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let ranges = self
            .regex_matchers
            .iter()
            .flat_map(|r| r.regex.find_iter(bases).map(|m| m.range()));
        if self.opts().reverse_complement {
            let bases_revcomp = &reverse_complement(bases);
            let ranges_revcomp = self.regex_matchers.iter().flat_map(|r| {
                r.regex
                    .find_iter(bases_revcomp)
                    .map(|m| m.range())
                    .map(|range| Range {
                        start: bases.len() - range.start - range.len(),
                        end: bases.len() - range.start,
                    })
            });
            bases_colored(bases, quals, ranges.chain(ranges_revcomp))
        } else {
            bases_colored(bases, quals, ranges)
        }
    }

    #[inline]
    fn opts(&self) -> MatcherOpts {
        self.opts
    }
}

/// Matcher for query names (read IDs).  Matches reads whose ID (the portion of the header before
/// the first whitespace) exactly matches one of the query names in the provided set.
pub struct QueryNameMatcher {
    query_names: AHashSet<Vec<u8>>,
    opts: MatcherOpts,
}

impl QueryNameMatcher {
    pub fn new(query_names: AHashSet<Vec<u8>>, opts: MatcherOpts) -> Self {
        Self { query_names, opts }
    }
}

impl Matcher for QueryNameMatcher {
    /// Not used — matching is by query name, not sequence. See `read_match()`.
    #[inline]
    fn bases_match(&self, _bases: &[u8]) -> bool {
        !self.opts.invert_match
    }

    fn color_matched_bases(&self, bases: &[u8], quals: &[u8]) -> (Vec<u8>, Vec<u8>) {
        (bases.to_vec(), quals.to_vec())
    }

    #[inline]
    fn opts(&self) -> MatcherOpts {
        self.opts
    }

    #[inline]
    fn read_match(&self, read: &RefRecord) -> bool {
        let id_match = self.query_names.contains(read.id_bytes());
        if self.opts.invert_match {
            !id_match
        } else {
            id_match
        }
    }
}

/// Factory for building a matcher
pub struct MatcherFactory;

impl MatcherFactory {
    pub fn new_matcher(
        pattern: &Option<String>,
        fixed_strings: bool,
        regexp: &Vec<String>,
        bitencs: Option<Vec<BitEnc>>,
        match_opts: MatcherOpts,
    ) -> Result<Box<dyn Matcher + Sync + Send>> {
        let use_bitmask = bitencs.is_some();
        let is_single = pattern.is_some();
        match (use_bitmask, is_single, fixed_strings) {
            (true, true, _) => {
                let bitencs = bitencs.unwrap();
                Ok(Box::new(BitMaskMatcher::new(
                    bitencs.into_iter().next().unwrap(),
                    match_opts,
                )))
            }
            (true, false, _) => Ok(Box::new(BitMaskSetMatcher::new(
                bitencs.unwrap(),
                match_opts,
            ))),
            (false, true, true) => Ok(Box::new(FixedStringMatcher::new(
                pattern.as_ref().unwrap(),
                match_opts,
            ))),
            (false, true, false) => Ok(Box::new(RegexMatcher::new(
                pattern.as_ref().unwrap(),
                match_opts,
            )?)),
            (false, false, true) => Ok(Box::new(FixedStringSetMatcher::new(regexp, match_opts)?)),
            (false, false, false) => Ok(Box::new(RegexSetMatcher::new(regexp, match_opts)?)),
        }
    }

    /// Create a matcher for query names (read IDs)
    pub fn new_query_name_matcher(
        query_names: AHashSet<Vec<u8>>,
        match_opts: MatcherOpts,
    ) -> Box<dyn Matcher + Sync + Send> {
        Box::new(QueryNameMatcher::new(query_names, match_opts))
    }
}

// Tests
#[cfg(test)]
pub mod tests {
    use crate::encode;
    use crate::matcher::*;
    use ahash::AHashSet;
    use rstest::rstest;

    // ############################################################################################
    // Tests to_bitvec()
    // ############################################################################################

    #[rstest]
    #[case(vec![(0, 1)], "AGG",  bitvec![1, 0, 0])] // single range that starts at the beginning of the read
    #[case(vec![(2, 3)], "AGG", bitvec![0, 0, 1])] // single range that ends at the ends of the read
    #[case(vec![(1, 2)], "AGG", bitvec![0, 1, 0])] // single range of a single base
    #[case(vec![(0 ,0)], "AGG", bitvec![0, 0, 0])] // empty set of ranges
    #[case(vec![(0, 3)], "AGG", bitvec![1, 1, 1])] // single range that encompasses the full read
    #[case(vec![(1, 4)], "AGGTC", bitvec![0, 1, 1, 1, 0])] // single range in the "middle" of the read
    #[case(vec![(0, 2), (3, 5)], "AGGTC", bitvec![1, 1, 0, 1, 1])] // two ranges non-overlapping
    #[case(vec![(0, 3), (3, 5)], "AGGTC", bitvec![1, 1, 1, 1, 1])] // two ranges abutting
    #[case(vec![(0, 4), (3, 5)], "AGGTC", bitvec![1, 1, 1, 1, 1])] // two ranges overlapping
    #[case(vec![(0, 3), (0, 5)], "AGGTC", bitvec![1, 1, 1, 1, 1])] // two ranges, one containing the other
    #[case(vec![(4, 5), (0, 2)], "AGGTC", bitvec![1, 1, 0, 0, 1])] // multiple ranges, but not in sorted order
    fn test_to_bitvec(
        #[case] ranges: Vec<(usize, usize)>,
        #[case] bases: &str,
        #[case] expected: BitVec,
    ) {
        let ranges = ranges
            .into_iter()
            .map(|(start, end)| std::ops::Range { start, end });
        let result_bitvec = to_bitvec(ranges, bases.len());
        assert_eq!(result_bitvec, expected);
    }

    // ############################################################################################
    // Tests bases_colored() single-base boundary
    // ############################################################################################

    /// Regression test: a single-base match at the start, middle, and end of a sequence must be
    /// colored.  Before the off-by-one fix (`last_bases_index + 1 < cur_bases_index` →
    /// `last_bases_index < cur_bases_index`), single-base matches produced empty output.
    #[test]
    fn test_single_base_match_colored() {
        let bases = b"ACGT";
        let quals = b"IIII";

        // Single-base match at position 0 ("A" in "ACGT")
        let ranges_start = vec![Range { start: 0, end: 1 }];
        let (colored_bases, _) = bases_colored(bases, quals, ranges_start.into_iter());
        assert!(
            !colored_bases.is_empty(),
            "single-base match at start must produce colored output"
        );

        // Single-base match in the middle (position 2, "G" in "ACGT")
        let ranges_mid = vec![Range { start: 2, end: 3 }];
        let (colored_bases, _) = bases_colored(bases, quals, ranges_mid.into_iter());
        assert!(
            !colored_bases.is_empty(),
            "single-base match in middle must produce colored output"
        );

        // Single-base match at the end (position 3, "T" in "ACGT")
        let ranges_end = vec![Range { start: 3, end: 4 }];
        let (colored_bases, _) = bases_colored(bases, quals, ranges_end.into_iter());
        assert!(
            !colored_bases.is_empty(),
            "single-base match at end must produce colored output"
        );

        // Verify the colored output contains the actual base bytes (not just ANSI escapes)
        // A single-base match ("G" at position 2 in "ACGT") should produce output containing
        // all four bases with coloring applied
        let ranges = vec![Range { start: 2, end: 3 }];
        let (colored_bases, colored_quals) = bases_colored(bases, quals, ranges.into_iter());
        let colored_str = String::from_utf8_lossy(&colored_bases);
        let colored_qual_str = String::from_utf8_lossy(&colored_quals);
        assert!(
            colored_str.contains("AC"),
            "uncolored prefix must be present"
        );
        assert!(colored_str.contains("G"), "colored base must be present");
        assert!(
            colored_str.contains("T"),
            "uncolored suffix must be present"
        );
        assert!(
            colored_qual_str.contains("II"),
            "uncolored qual prefix must be present"
        );
    }

    // ############################################################################################
    // Test FixedStringMatcher::read_match()
    // ############################################################################################

    #[rstest]
    #[case(false, "AG", "AGG", true)] // fixed string with match when reverse_complement is false
    #[case(false, "CC", "AGG", false)] // fixed string with no match when reverse_complement is false
    #[case(true, "CC", "AGG", true)] // fixed string with match when reverse_complement is true
    #[case(true, "TT", "AGG", false)] // fixed string with no match when reverse_complement is true
    #[case(false, "AT", "ATGAT", true)] // fixed string with multiple non-overlapping matches when reverse_complement is false
    #[case(true, "CG", "GCCG", true)] // fixed string with multiple non-overlapping matches when reverse_complement is true
    #[case(false, "AGAG", "AGAGAGAG", true)] // fixed string with overlapping matches when reverse_complemet is false
    #[case(true, "TCTC", "AGAGAGAG", true)] // fixed string with overlapping matches when reverse_complemet is true
    fn test_fixed_string_matcher_read_match(
        #[case] reverse_complement: bool,
        #[case] pattern: &str,
        #[case] seq: &str,
        #[case] expected: bool,
    ) {
        let invert_matches = [true, false];
        for invert_match in IntoIterator::into_iter(invert_matches) {
            let opts = MatcherOpts {
                invert_match,
                reverse_complement,
            };
            let matcher = FixedStringMatcher::new(pattern, opts);
            let qual = (0..seq.len()).map(|_| "X").collect::<String>();
            let record = format!("@id\n{seq}\n+\n{qual}\n");
            let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
            let read_record = reader.next().unwrap().unwrap();
            let result = matcher.read_match(&read_record);
            if invert_match {
                assert_ne!(result, expected);
            } else {
                assert_eq!(result, expected);
            }
        }
    }

    // ############################################################################################
    // Tests FixedStringSetMatcher::read_match()
    // ############################################################################################

    #[rstest]
    #[case(false, vec!["A", "AGG", "G"], "AGGG", true)] // match is true when reverse_complement is false
    #[case(true, vec!["A", "AGG", "G"], "TCCC", true)] // match is true when reverse_complement is true
    #[case(false, vec!["A", "AGG", "G"], "TTTT", false)] // match is false when reverse_complement is false
    #[case(true, vec!["T", "AAA"], "CCCCC", false)] // match is false when reverse_complement is true
    #[case(false, vec!["AGG", "C", "TT"], "AGGTT",  true)] // match is true but not all patterns match when reverse_complement is false
    #[case(true, vec!["AGG", "C", "TT"], "GGGGG",  true)] // match is true but not all patterns match when reverse_complement is true
    #[case(false, vec!["AC", "TT"], "TTACGTT",  true)] // match is true but one pattern matches multiple times when reverse_complement is false
    #[case(true, vec!["GT", "AA"], "TTACGTT",  true)] // match is true but one pattern matches multiple times when reverse_complement is true
    #[case(false, vec!["GAGA","AGTT"], "GAGAGTT",  true)] // match is true with overlapping matches when reverse_complment is false
    #[case(true, vec!["CTCT","AACT"], "GAGAGTT",  true)] // match is true with overlapping matches when reverse_complment is true
    fn test_fixed_string_set_metcher_read_match(
        #[case] reverse_complement: bool,
        #[case] patterns: Vec<&str>,
        #[case] seq: &str,
        #[case] expected: bool,
    ) {
        let invert_matches = [true, false];
        for invert_match in IntoIterator::into_iter(invert_matches) {
            let opts = MatcherOpts {
                invert_match,
                reverse_complement,
            };
            let matcher = FixedStringSetMatcher::new(patterns.iter(), opts).unwrap();
            let qual = (0..seq.len()).map(|_| "X").collect::<String>();
            let record = format!("@id\n{seq}\n+\n{qual}\n");
            let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
            let read_record = reader.next().unwrap().unwrap();
            let result = matcher.read_match(&read_record);
            if invert_match {
                assert_ne!(result, expected);
            } else {
                assert_eq!(result, expected);
            }
        }
    }

    // ############################################################################################
    // Test RegexMatcher::read_match()
    // ############################################################################################

    #[rstest]
    #[case(false, "^A", "AGG", true)] // regex with one match when reverse_complement is false
    #[case(false, "^T", "AGG", false)] // regex with no matches when reverse_complement is false
    #[case(true, "^C", "AGG", true)] // regex with one match when reverse_complement is true
    #[case(true, "^T", "AGG", false)] // regex with no matches when reverse_complement is true
    #[case(false, "A.A", "ATATA", true)] // regex with overlapping matches when reverse_complement is false
    #[case(true, "T.G", "CACACA", false)] // regex with overlapping matches when reverse_complement is true
    fn test_regex_matcher_read_match(
        #[case] reverse_complement: bool,
        #[case] pattern: &str,
        #[case] seq: &str,
        #[case] expected: bool,
    ) {
        let invert_matches = [true, false];
        for invert_match in IntoIterator::into_iter(invert_matches) {
            let opts = MatcherOpts {
                invert_match,
                reverse_complement,
            };

            let matcher = RegexMatcher::new(pattern, opts).unwrap();
            let qual = (0..seq.len()).map(|_| "X").collect::<String>();
            let record = format!("@id\n{seq}\n+\n{qual}\n");
            let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
            let read_record = reader.next().unwrap().unwrap();
            let result = matcher.read_match(&read_record);
            if invert_match {
                assert_ne!(result, expected);
            } else {
                assert_eq!(result, expected);
            }
        }
    }

    // ############################################################################################
    // Tests RegexSetMatcher::read_match()
    // ############################################################################################

    #[rstest]
    #[case(false, vec!["^A.G", "C..", "$T"], "AGGCTT", true)] // match is true when reverse_complement is false
    #[case(true, vec!["^T.C", "..G", "$A"], "AGGCTT", true)] // match is true when reverse_complement is true
    #[case(false, vec!["^A.G", "G..", "$T"], "CCTCA", false)] // match is false when reverse_complemet is false
    #[case(true, vec!["$A", "C.CC"], "CCTCA", false)] // match is false when reverse_complemet is true
    #[case(false, vec!["^T", ".GG", "A.+G"],  "ATCTACTACG",  true)] // match is true but not all patterns match when reverse_complement is false
    #[case(true, vec!["^C", ".CC", "C+.T"],  "ATCTACTACG",  true)] // match is true when reverse_complement is true
    #[case(false, vec!["^T", "T.A"], "TTAATAA", true)] // match is true but one pattern matches multiple times when reverse_complement is false
    #[case(true, vec!["^T", "T.A"], "AATA", true)] // match is true but one pattern matches multiple times when reverse_complemetn is true
    #[case(false, vec!["^T","T.+G"], "TAGAGTG",  true)] // match is true with overlapping matches when reverse_complement is false
    #[case(true, vec!["^A","A.+C"], "TAGAGTG",  true)] // match is true with overlapping matches when reverse_complement is true
    fn test_regex_set_metcher_read_match(
        #[case] reverse_complement: bool,
        #[case] patterns: Vec<&str>,
        #[case] seq: &str,
        #[case] expected: bool,
    ) {
        let invert_matches = [true, false];
        for invert_match in IntoIterator::into_iter(invert_matches) {
            let opts = MatcherOpts {
                invert_match,
                reverse_complement,
            };

            let matcher = RegexSetMatcher::new(patterns.iter(), opts).unwrap();
            let qual = (0..seq.len()).map(|_| "X").collect::<String>();
            let record = format!("@id\n{seq}\n+\n{qual}\n");
            let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
            let read_record = reader.next().unwrap().unwrap();
            let result = matcher.read_match(&read_record);
            if invert_match {
                assert_ne!(result, expected);
            } else {
                assert_eq!(result, expected);
            }
        }
    }

    // ############################################################################################
    // Tests QueryNameMatcher::read_match()
    // ############################################################################################

    #[rstest]
    #[case(false, vec!["read1", "read2"], "read1", true)] // exact match found
    #[case(false, vec!["read1", "read2"], "read3", false)] // no match
    #[case(false, vec!["read1"], "read1_extra", false)] // partial match should fail (query name is a prefix of read ID)
    #[case(false, vec!["read1_extra"], "read1", false)] // partial match should fail (read ID is a prefix of query name)
    #[case(true, vec!["read1", "read2"], "read1", false)] // invert_match: found becomes false
    #[case(true, vec!["read1", "read2"], "read3", true)] // invert_match: not found becomes true
    fn test_query_name_matcher_read_match(
        #[case] invert_match: bool,
        #[case] query_names: Vec<&str>,
        #[case] read_id: &str,
        #[case] expected: bool,
    ) {
        let opts = MatcherOpts {
            invert_match,
            reverse_complement: false,
        };
        let names: AHashSet<Vec<u8>> = query_names.iter().map(|s| s.as_bytes().to_vec()).collect();
        let matcher = QueryNameMatcher::new(names, opts);
        let seq = "ACGT";
        let qual = "XXXX";
        let record = format!("@{read_id} description\n{seq}\n+\n{qual}\n");
        let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
        let read_record = reader.next().unwrap().unwrap();
        assert_eq!(matcher.read_match(&read_record), expected);
    }

    #[test]
    fn test_query_name_matcher_with_header_description() {
        // Tests that only the read ID (before first whitespace) is matched, not the full header
        let opts = MatcherOpts {
            invert_match: false,
            reverse_complement: false,
        };
        let names: AHashSet<Vec<u8>> = vec!["SRR001666.1".as_bytes().to_vec()]
            .into_iter()
            .collect();
        let matcher = QueryNameMatcher::new(names, opts);

        // Record with description after read ID
        let record = "@SRR001666.1 071112_SLXA description\nACGT\n+\nXXXX\n";
        let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
        let read_record = reader.next().unwrap().unwrap();
        assert!(matcher.read_match(&read_record));
    }

    // ############################################################################################
    // Tests validate_fixed_pattern()
    // ############################################################################################

    #[rstest]
    #[case("AGTGTGATG", false)]
    #[case("QRNQRNQRN", true)]
    #[case("ARNDCEQGHILKMFPSTWYV", true)] // all standard 20
    #[case("BJOUXZ", true)] // extended IUPAC: B, J, O, U, X, Z
    fn test_validate_fixed_pattern_is_ok(#[case] pattern: &str, #[case] protein: bool) {
        let result = validate_fixed_pattern(pattern, protein);
        assert!(result.is_ok());
    }

    #[rstest]
    #[case(
        "AXGTGTGATG",
        false,
        "Fixed pattern must contain only DNA bases: A .. [X] .. GTGTGATG"
    )]
    #[case(
        "QRN1QRN",
        true,
        "Fixed pattern must contain only amino acids: QRN .. [1] .. QRN"
    )]
    fn test_validate_fixed_pattern_error(
        #[case] pattern: &str,
        #[case] protein: bool,
        #[case] msg: &str,
    ) {
        let result = validate_fixed_pattern(pattern, protein);
        let inner = result.unwrap_err().to_string();
        assert_eq!(inner, msg);
    }

    #[test]
    fn test_validate_fixed_pattern_rejects_non_ascii() {
        // A non-ASCII character whose low byte happens to match 'A' (0x41)
        // U+0141 (Latin capital letter L with stroke) has value 321, low byte = 0x41 = 'A'
        let pattern = "AC\u{0141}T";
        let result = validate_fixed_pattern(pattern, false);
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("Fixed pattern must contain only DNA bases")
        );
    }

    // ############################################################################################
    // Tests BitMaskMatcher: has_iupac detection
    // ############################################################################################

    #[test]
    fn test_bitmask_matcher_dna_only_has_no_iupac() {
        let opts = MatcherOpts {
            invert_match: false,
            reverse_complement: false,
        };
        let bitenc = encode(b"ACGT");
        let matcher = BitMaskMatcher::new(bitenc, opts);
        assert!(
            !matcher.has_iupac,
            "DNA-only pattern should have has_iupac=false"
        );
    }

    #[test]
    fn test_bitmask_matcher_iupac_pattern_has_iupac() {
        let opts = MatcherOpts {
            invert_match: false,
            reverse_complement: false,
        };
        let bitenc = encode(b"ACRT");
        let matcher = BitMaskMatcher::new(bitenc, opts);
        assert!(
            matcher.has_iupac,
            "IUPAC pattern should have has_iupac=true"
        );
    }

    // ############################################################################################
    // Tests BitMaskMatcher: edge cases
    // ############################################################################################

    #[test]
    fn test_bitmask_matcher_read_shorter_than_pattern() {
        let opts = MatcherOpts {
            invert_match: false,
            reverse_complement: false,
        };
        let bitenc = encode(b"ACGTACGT");
        let matcher = BitMaskMatcher::new(bitenc, opts);
        let seq = "ACG";
        let qual = (0..seq.len()).map(|_| "X").collect::<String>();
        let record = format!("@id\n{seq}\n+\n{qual}\n");
        let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
        let read_record = reader.next().unwrap().unwrap();
        assert!(!matcher.read_match(&read_record));
    }

    #[test]
    fn test_bitmask_set_matcher_read_shorter_than_pattern() {
        let opts = MatcherOpts {
            invert_match: false,
            reverse_complement: false,
        };
        let bitencs = vec![encode(b"ACGTACGT")];
        let matcher = BitMaskSetMatcher::new(bitencs, opts);
        let seq = "ACG";
        let qual = (0..seq.len()).map(|_| "X").collect::<String>();
        let record = format!("@id\n{seq}\n+\n{qual}\n");
        let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
        let read_record = reader.next().unwrap().unwrap();
        assert!(!matcher.read_match(&read_record));
    }

    // ############################################################################################
    // Test BitMaskMatcher::read_match() with DNA-only patterns
    // ############################################################################################

    #[rstest]
    #[case(false, "AG", "AGG", true)]
    #[case(false, "CC", "AGG", false)]
    #[case(true, "CC", "AGG", true)]
    #[case(true, "TT", "AGG", false)]
    #[case(false, "AT", "ATGAT", true)]
    #[case(true, "CG", "GCCG", true)]
    #[case(false, "AGAG", "AGAGAGAG", true)]
    #[case(true, "TCTC", "AGAGAGAG", true)]
    fn test_bitmask_matcher_dna_read_match(
        #[case] reverse_complement: bool,
        #[case] pattern: &str,
        #[case] seq: &str,
        #[case] expected: bool,
    ) {
        let invert_matches = [true, false];
        for invert_match in IntoIterator::into_iter(invert_matches) {
            let opts = MatcherOpts {
                invert_match,
                reverse_complement,
            };
            let matcher = BitMaskMatcher::new(encode(pattern.as_bytes()), opts);
            let qual = (0..seq.len()).map(|_| "X").collect::<String>();
            let record = format!("@id\n{seq}\n+\n{qual}\n");
            let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
            let read_record = reader.next().unwrap().unwrap();
            let result = matcher.read_match(&read_record);
            if invert_match {
                assert_ne!(result, expected);
            } else {
                assert_eq!(result, expected);
            }
        }
    }

    // ############################################################################################
    // Test BitMaskMatcher::read_match() with IUPAC patterns
    // ############################################################################################

    #[rstest]
    #[case(false, "GATK", "GATG", true)]
    #[case(false, "GATK", "GATT", true)]
    #[case(false, "GATK", "GATA", false)]
    #[case(false, "GATK", "GATC", false)]
    #[case(false, "GATR", "GATGA", true)]
    #[case(false, "GATR", "GATA", true)]
    #[case(false, "GATR", "GATGT", true)]
    #[case(false, "N", "A", true)]
    #[case(false, "N", "C", true)]
    #[case(false, "N", "G", true)]
    #[case(false, "N", "T", true)]
    #[case(false, "ANA", "ACA", true)]
    #[case(false, "ANA", "AAA", true)]
    #[case(true, "GATK", "CATC", true)]
    fn test_bitmask_matcher_iupac_read_match(
        #[case] reverse_complement: bool,
        #[case] pattern: &str,
        #[case] seq: &str,
        #[case] expected: bool,
    ) {
        let opts = MatcherOpts {
            invert_match: false,
            reverse_complement,
        };
        let matcher = BitMaskMatcher::new(encode(pattern.as_bytes()), opts);
        let qual = (0..seq.len()).map(|_| "X").collect::<String>();
        let record = format!("@id\n{seq}\n+\n{qual}\n");
        let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
        let read_record = reader.next().unwrap().unwrap();
        assert_eq!(
            matcher.read_match(&read_record),
            expected,
            "pattern={} seq={} revcomp={}",
            pattern,
            seq,
            reverse_complement
        );
    }

    // ############################################################################################
    // Tests BitMaskSetMatcher::read_match() with DNA-only patterns
    // ############################################################################################

    #[rstest]
    #[case(false, vec!["A", "AGG", "G"], "AGGG", true)]
    #[case(true, vec!["A", "AGG", "G"], "TCCC", true)]
    #[case(false, vec!["A", "AGG", "G"], "TTTT", false)]
    #[case(true, vec!["T", "AAA"], "CCCCC", false)]
    #[case(false, vec!["AGG", "C", "TT"], "AGGTT", true)]
    #[case(true, vec!["AGG", "C", "TT"], "GGGGG", true)]
    #[case(false, vec!["AC", "TT"], "TTACGTT", true)]
    #[case(true, vec!["GT", "AA"], "TTACGTT", true)]
    #[case(false, vec!["GAGA","AGTT"], "GAGAGTT", true)]
    #[case(true, vec!["CTCT","AACT"], "GAGAGTT", true)]
    fn test_bitmask_set_matcher_dna_read_match(
        #[case] reverse_complement: bool,
        #[case] patterns: Vec<&str>,
        #[case] seq: &str,
        #[case] expected: bool,
    ) {
        let invert_matches = [true, false];
        for invert_match in IntoIterator::into_iter(invert_matches) {
            let opts = MatcherOpts {
                invert_match,
                reverse_complement,
            };
            let bitencs: Vec<_> = patterns.iter().map(|p| encode(p.as_bytes())).collect();
            let matcher = BitMaskSetMatcher::new(bitencs, opts);
            let qual = (0..seq.len()).map(|_| "X").collect::<String>();
            let record = format!("@id\n{seq}\n+\n{qual}\n");
            let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
            let read_record = reader.next().unwrap().unwrap();
            let result = matcher.read_match(&read_record);
            if invert_match {
                assert_ne!(result, expected);
            } else {
                assert_eq!(result, expected);
            }
        }
    }

    // ############################################################################################
    // Tests BitMaskSetMatcher::read_match() with IUPAC patterns
    // ############################################################################################

    #[rstest]
    #[case(false, vec!["GATK", "ANA"], "GATG", true)]
    #[case(false, vec!["GATK", "ANA"], "ACA", true)]
    #[case(false, vec!["GATK", "ANA"], "CCCC", false)]
    #[case(false, vec!["R", "Y"], "A", true)]
    #[case(false, vec!["R", "Y"], "C", true)]
    #[case(false, vec!["R", "Y"], "AAAA", true)]
    fn test_bitmask_set_matcher_iupac_read_match(
        #[case] reverse_complement: bool,
        #[case] patterns: Vec<&str>,
        #[case] seq: &str,
        #[case] expected: bool,
    ) {
        let opts = MatcherOpts {
            invert_match: false,
            reverse_complement,
        };
        let bitencs: Vec<_> = patterns.iter().map(|p| encode(p.as_bytes())).collect();
        let matcher = BitMaskSetMatcher::new(bitencs, opts);
        let qual = (0..seq.len()).map(|_| "X").collect::<String>();
        let record = format!("@id\n{seq}\n+\n{qual}\n");
        let mut reader = seq_io::fastq::Reader::new(record.as_bytes());
        let read_record = reader.next().unwrap().unwrap();
        assert_eq!(
            matcher.read_match(&read_record),
            expected,
            "patterns={:?} seq={} revcomp={}",
            patterns,
            seq,
            reverse_complement
        );
    }
}
