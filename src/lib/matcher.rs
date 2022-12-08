use bitvec::prelude::*;
use std::ops::Range;

use crate::color::{color_background, color_head};
use crate::color::{COLOR_BACKGROUND, COLOR_BASES, COLOR_QUALS};
use crate::reverse_complement;
use crate::DNA_BASES;
use anyhow::{bail, Context, Result};
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
    /// Color the read based on the match
    pub color: bool,
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
                if last_bases_index + 1 < cur_bases_index {
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
                if last_bases_index + 1 < cur_bases_index {
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
    if last_bases_index + 1 < cur_bases_index {
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

/// Validates that a given FIXED pattern contains only valid DNA bases (ACGTN)
pub fn validate_fixed_pattern(pattern: &str) -> Result<()> {
    for (index, base) in pattern.chars().enumerate() {
        if !DNA_BASES.contains(&(base as u8)) {
            bail!(
                "Fixed pattern must contain only DNA bases: {} .. [{}] .. {}",
                &pattern[0..index],
                &pattern[index..=index],
                &pattern[index + 1..],
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
    fn read_match(&self, read: &mut OwnedRecord) -> bool {
        let match_found = if self.opts().invert_match {
            self.bases_match(read.seq())
                && (!self.opts().reverse_complement
                    || self.bases_match(&reverse_complement(read.seq())))
        } else {
            self.bases_match(read.seq())
                || (self.opts().reverse_complement
                    && self.bases_match(&reverse_complement(read.seq())))
        };

        if self.opts().color {
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

        match_found
    }
}

/// Matcher for a fixed string pattern
pub struct FixedStringMatcher {
    pattern: Vec<u8>,
    opts: MatcherOpts,
}

impl Matcher for FixedStringMatcher {
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
    opts: MatcherOpts,
}

impl Matcher for FixedStringSetMatcher {
    fn bases_match(&self, bases: &[u8]) -> bool {
        self.patterns
            .iter()
            .any(|pattern| bases.find(pattern).is_some())
            != self.opts.invert_match
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

    fn opts(&self) -> MatcherOpts {
        self.opts
    }
}

impl FixedStringSetMatcher {
    pub fn new<I, S>(patterns: I, opts: MatcherOpts) -> Self
    where
        S: AsRef<str>,
        I: IntoIterator<Item = S>,
    {
        let patterns: Vec<Vec<u8>> = patterns
            .into_iter()
            .map(|pattern| pattern.as_ref().to_owned().as_bytes().to_vec())
            .collect();
        Self { patterns, opts }
    }
}

/// Matcher for a regular expression pattern
pub struct RegexMatcher {
    regex: Regex,
    opts: MatcherOpts,
}

impl RegexMatcher {
    pub fn new(pattern: &str, opts: MatcherOpts) -> Self {
        let regex = RegexBuilder::new(pattern)
            .build()
            .context(format!("Invalid regular expression: {}", pattern))
            .unwrap();
        Self { regex, opts }
    }
}

impl Matcher for RegexMatcher {
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
    pub fn new<I, S>(patterns: I, opts: MatcherOpts) -> Self
    where
        S: AsRef<str>,
        I: IntoIterator<Item = S>,
    {
        let string_patterns: Vec<String> = patterns
            .into_iter()
            .map(|p| p.as_ref().to_string())
            .collect();
        let regex_set = RegexSetBuilder::new(string_patterns.clone())
            .build()
            .unwrap();
        let regex_matchers: Vec<RegexMatcher> = string_patterns
            .into_iter()
            .map(|pattern| RegexMatcher::new(pattern.as_ref(), opts))
            .collect();
        Self {
            regex_set,
            regex_matchers,
            opts,
        }
    }
}

impl Matcher for RegexSetMatcher {
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

    fn opts(&self) -> MatcherOpts {
        self.opts
    }
}

/// Factory for building a matcher
pub struct MatcherFactory;

impl MatcherFactory {
    pub fn new_matcher(
        pattern: &Option<String>,
        fixed_strings: bool,
        regexp: &Vec<String>,
        match_opts: MatcherOpts,
    ) -> Box<dyn Matcher + Sync + Send> {
        match (fixed_strings, &pattern) {
            (true, Some(pattern)) => Box::new(FixedStringMatcher::new(pattern, match_opts)),
            (false, Some(pattern)) => Box::new(RegexMatcher::new(pattern, match_opts)),
            (true, None) => Box::new(FixedStringSetMatcher::new(regexp, match_opts)),
            (false, None) => Box::new(RegexSetMatcher::new(regexp, match_opts)),
        }
    }
}

// Tests
#[cfg(test)]
pub mod tests {
    use crate::matcher::*;
    use rstest::rstest;

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
        for invert_match in IntoIterator::into_iter(invert_matches).to_owned() {
            let opts = MatcherOpts {
                invert_match,
                reverse_complement,
                color: false,
            };
            let matcher = FixedStringMatcher::new(pattern, opts);
            let mut read_record = write_owned_record(seq);
            let result = matcher.read_match(&mut read_record);
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
        for invert_match in IntoIterator::into_iter(invert_matches).to_owned() {
            let opts = MatcherOpts {
                invert_match,
                reverse_complement,
                color: false,
            };
            let matcher = FixedStringSetMatcher::new(patterns.iter(), opts);
            let mut read_record = write_owned_record(seq);
            let result = matcher.read_match(&mut read_record);
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
        for invert_match in IntoIterator::into_iter(invert_matches).to_owned() {
            let opts = MatcherOpts {
                invert_match,
                reverse_complement,
                color: false,
            };

            let matcher = RegexMatcher::new(&pattern, opts);
            let mut read_record = write_owned_record(seq);
            let result = matcher.read_match(&mut read_record);
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
        for invert_match in IntoIterator::into_iter(invert_matches).to_owned() {
            let opts = MatcherOpts {
                invert_match,
                reverse_complement,
                color: false,
            };

            let matcher = RegexSetMatcher::new(patterns.iter(), opts);
            let mut read_record = write_owned_record(seq);
            let result = matcher.read_match(&mut read_record);
            if invert_match {
                assert_ne!(result, expected);
            } else {
                assert_eq!(result, expected);
            }
        }
    }

    // ############################################################################################
    // Tests validate_fixed_pattern()
    // ############################################################################################

    #[test]
    fn test_validate_fixed_pattern_is_ok() {
        let pattern = "AGTGTGATG";
        let result = validate_fixed_pattern(&pattern);
        assert!(result.is_ok())
    }
    #[test]
    fn test_validate_fixed_pattern_error() {
        let pattern = "A.GTGTGATG";
        let msg = String::from("Fixed pattern must contain only DNA bases: A .. [.] .. GTGTGATG");
        let result = validate_fixed_pattern(&pattern);
        let inner = result.unwrap_err().to_string();
        assert_eq!(inner, msg);
    }
}
