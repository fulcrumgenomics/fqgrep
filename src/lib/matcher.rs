use aho_corasick::AhoCorasick;
use bitvec::prelude::*;
use seq_io::fastq::RefRecord;
use std::ops::Range;

use crate::DNA_BASES;
use crate::color::{COLOR_BACKGROUND, COLOR_BASES, COLOR_QUALS};
use crate::color::{color_background, color_head};
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

/// Validates that a given FIXED pattern contains only valid DNA bases (ACGTN)
pub fn validate_fixed_pattern(pattern: &str) -> Result<()> {
    for (index, base) in pattern.char_indices() {
        if !DNA_BASES.contains(&(base as u8)) {
            let end_index = index + base.len_utf8();
            bail!(
                "Fixed pattern must contain only DNA bases: {} .. [{}] .. {}",
                &pattern[0..index],
                &pattern[index..end_index],
                &pattern[end_index..],
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
    pub fn new<I, S>(patterns: I, opts: MatcherOpts) -> Self
    where
        S: AsRef<str>,
        I: IntoIterator<Item = S>,
    {
        let patterns: Vec<Vec<u8>> = patterns
            .into_iter()
            .map(|pattern| pattern.as_ref().to_owned().as_bytes().to_vec())
            .collect();
        let aho_corasick =
            AhoCorasick::new(&patterns).expect("Failed to build Aho-Corasick automaton");
        Self {
            patterns,
            aho_corasick,
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
    pub fn new(pattern: &str, opts: MatcherOpts) -> anyhow::Result<Self> {
        let regex = RegexBuilder::new(pattern)
            .build()
            .with_context(|| format!(
                "Invalid regular expression: '{}'. \
                Hint: Use --fixed-strings (-F) if you want to search for literal text containing regex special characters.",
                pattern
            ))?;
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
    pub fn new<I, S>(patterns: I, opts: MatcherOpts) -> anyhow::Result<Self>
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
            .collect::<anyhow::Result<Vec<_>>>()?;
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

/// Factory for building a matcher
pub struct MatcherFactory;

impl MatcherFactory {
    pub fn new_matcher(
        pattern: &Option<String>,
        fixed_strings: bool,
        regexp: &Vec<String>,
        match_opts: MatcherOpts,
    ) -> anyhow::Result<Box<dyn Matcher + Sync + Send>> {
        match (fixed_strings, &pattern) {
            (true, Some(pattern)) => Ok(Box::new(FixedStringMatcher::new(pattern, match_opts))),
            (false, Some(pattern)) => Ok(Box::new(RegexMatcher::new(pattern, match_opts)?)),
            (true, None) => Ok(Box::new(FixedStringSetMatcher::new(regexp, match_opts))),
            (false, None) => Ok(Box::new(RegexSetMatcher::new(regexp, match_opts)?)),
        }
    }
}

// Tests
#[cfg(test)]
pub mod tests {
    use crate::matcher::*;
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
            let matcher = FixedStringSetMatcher::new(patterns.iter(), opts);
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
    // Tests validate_fixed_pattern()
    // ############################################################################################

    #[test]
    fn test_validate_fixed_pattern_is_ok() {
        let pattern = "AGTGTGATG";
        let result = validate_fixed_pattern(pattern);
        assert!(result.is_ok());
    }
    #[test]
    fn test_validate_fixed_pattern_error() {
        let pattern = "AXGTGTGATG";
        let msg = String::from("Fixed pattern must contain only DNA bases: A .. [X] .. GTGTGATG");
        let result = validate_fixed_pattern(pattern);
        let inner = result.unwrap_err().to_string();
        assert_eq!(inner, msg);
    }

    // ############################################################################################
    // Tests for reverse complement position calculations
    // ############################################################################################

    #[test]
    fn test_rc_position_simple() {
        // Sequence: AAATTTCCC (9 bases, indices 0-8)
        // RC:       GGGAAATTT
        // Pattern: TTT (3 bases) - looking for this in forward
        // RC of TTT is AAA
        // In RC: AAA appears at position 3 (indices 3-5 in RC)
        // Forward position should be: 9 - 3 - 3 = 3 (indices 3-5 in forward)

        let seq = b"AAATTTCCC";
        let pattern = b"TTT"; // This is what we're looking for in forward

        // Verify the RC
        let rc = reverse_complement(seq);
        assert_eq!(&rc, b"GGGAAATTT");

        // Find RC of pattern in RC (AAA in RC)
        let pattern_rc = reverse_complement(pattern);
        let rc_start = rc.find(&pattern_rc).unwrap();
        assert_eq!(rc_start, 3);

        // Calculate forward position using the formula
        let forward_start = seq.len() - rc_start - pattern.len();
        assert_eq!(forward_start, 3); // 9 - 3 - 3 = 3

        // Verify the pattern actually exists at that position in forward
        assert_eq!(&seq[forward_start..forward_start + pattern.len()], pattern);
    }

    #[test]
    fn test_rc_position_at_start() {
        // Sequence: AAACCCGGG (9 bases)
        // RC:       CCCGGGΤΤΤ
        // Pattern: AAA (we want to find this in forward at position 0)
        // RC of AAA is TTT
        // In RC: TTT at position 6 (end of RC)
        // Forward: 9 - 6 - 3 = 0 (start of forward)

        let seq = b"AAACCCGGG";
        let pattern = b"AAA";
        let rc = reverse_complement(seq);
        assert_eq!(&rc, b"CCCGGGTTT");

        // Find RC of pattern in RC
        let pattern_rc = reverse_complement(pattern);
        let rc_start = rc.find(&pattern_rc).unwrap();
        assert_eq!(rc_start, 6);

        let forward_start = seq.len() - rc_start - pattern.len();
        assert_eq!(forward_start, 0);
        assert_eq!(&seq[forward_start..forward_start + pattern.len()], pattern);
    }

    #[test]
    fn test_rc_position_at_end() {
        // Sequence: CCCGGGTTT (9 bases)
        // RC:       AAACCCGGG
        // Pattern: AAA
        // In RC: AAA at position 0 (start of RC)
        // Forward: 9 - 0 - 3 = 6 (end of forward)

        let seq = b"CCCGGGTTT";
        let pattern = b"AAA";
        let rc = reverse_complement(seq);
        assert_eq!(&rc, b"AAACCCGGG");

        let rc_start = rc.find(pattern).unwrap();
        assert_eq!(rc_start, 0);

        let forward_start = seq.len() - rc_start - pattern.len();
        assert_eq!(forward_start, 6);

        // Check it matches in forward (TTT is RC of AAA)
        let expected_forward = reverse_complement(pattern);
        assert_eq!(
            &seq[forward_start..forward_start + pattern.len()],
            expected_forward.as_slice()
        );
    }

    #[test]
    fn test_rc_position_longer_pattern() {
        // Sequence: ACGTACGTACGT (12 bases)
        // Pattern: ACGT (4 bases)
        // RC:       ACGTACGTACGT (palindrome!)
        // Pattern appears at positions 0, 4, 8 in both forward and RC

        let seq = b"ACGTACGTACGT";
        let pattern = b"ACGT";
        let rc = reverse_complement(seq);

        // This sequence is a palindrome
        assert_eq!(seq, rc.as_slice());

        // Test one position
        let rc_start = 4;
        let forward_start = seq.len() - rc_start - pattern.len();
        assert_eq!(forward_start, 4); // 12 - 4 - 4 = 4 (palindrome property)
    }

    #[test]
    fn test_rc_position_single_base() {
        // Edge case: single base pattern
        let seq = b"ACGTACGT";
        let pattern = b"T";
        let rc = reverse_complement(seq);
        // RC: ACGTACGT (A<->T, C<->G, so ACGT reversed = ACGT)

        // Find first T in RC
        let rc_start = rc.find(pattern).unwrap();
        let forward_start = seq.len() - rc_start - pattern.len();

        // Verify the base at forward position is the RC of the RC base
        assert_eq!(seq[forward_start], b'A'); // T in RC corresponds to A in forward
    }

    // ############################################################################################
    // Property-based tests
    // ############################################################################################

    use proptest::prelude::*;

    /// Strategy for generating valid DNA sequences
    fn dna_sequence() -> impl Strategy<Value = Vec<u8>> {
        prop::collection::vec(prop::sample::select(b"ACGT".to_vec()), 1..100)
    }

    /// Strategy for generating DNA patterns
    fn dna_pattern() -> impl Strategy<Value = Vec<u8>> {
        prop::collection::vec(prop::sample::select(b"ACGT".to_vec()), 1..20)
    }

    proptest! {
        // Reduce test cases from default 256 to 20 for faster tests
        #![proptest_config(ProptestConfig::with_cases(20))]

        /// Property: Reverse complement is involutive - rc(rc(seq)) == seq
        #[test]
        fn prop_reverse_complement_involutive(seq in dna_sequence()) {
            let rc = reverse_complement(&seq);
            let rc_rc = reverse_complement(&rc);
            prop_assert_eq!(seq, rc_rc);
        }

        /// Property: A sequence always matches itself as a fixed string
        #[test]
        fn prop_sequence_matches_itself(seq in dna_sequence()) {
            let opts = MatcherOpts {
                invert_match: false,
                reverse_complement: false,
            };
            let matcher = FixedStringMatcher::new(
                std::str::from_utf8(&seq).unwrap(),
                opts
            );
            prop_assert!(matcher.bases_match(&seq));
        }

        /// Property: Pattern matching is consistent between forward and RC
        #[test]
        fn prop_rc_matching_consistent(seq in dna_sequence(), pattern in dna_pattern()) {
            let opts_forward = MatcherOpts {
                invert_match: false,
                reverse_complement: false,
            };
            let opts_rc = MatcherOpts {
                invert_match: false,
                reverse_complement: true,
            };

            let pattern_str = std::str::from_utf8(&pattern).unwrap();
            let matcher_forward = FixedStringMatcher::new(pattern_str, opts_forward);
            let matcher_rc = FixedStringMatcher::new(pattern_str, opts_rc);

            let matches_forward = matcher_forward.bases_match(&seq);
            let matches_with_rc = matcher_rc.bases_match(&seq);

            // If pattern matches forward, it should also match with RC enabled
            if matches_forward {
                prop_assert!(matches_with_rc);
            }
        }

        /// Property: Invert match negates the result
        #[test]
        fn prop_invert_match_negates(seq in dna_sequence(), pattern in dna_pattern()) {
            let opts_normal = MatcherOpts {
                invert_match: false,
                reverse_complement: false,
            };
            let opts_inverted = MatcherOpts {
                invert_match: true,
                reverse_complement: false,
            };

            let pattern_str = std::str::from_utf8(&pattern).unwrap();
            let matcher_normal = FixedStringMatcher::new(pattern_str, opts_normal);
            let matcher_inverted = FixedStringMatcher::new(pattern_str, opts_inverted);

            let matches_normal = matcher_normal.bases_match(&seq);
            let matches_inverted = matcher_inverted.bases_match(&seq);

            // Invert should always negate the result
            prop_assert_ne!(matches_normal, matches_inverted);
        }

        /// Property: Fixed string matcher and regex matcher agree on literal patterns
        #[test]
        fn prop_fixed_and_regex_agree(seq in dna_sequence(), pattern in dna_pattern()) {
            let opts = MatcherOpts {
                invert_match: false,
                reverse_complement: false,
            };

            let pattern_str = std::str::from_utf8(&pattern).unwrap();
            let fixed_matcher = FixedStringMatcher::new(pattern_str, opts);

            // Only test if the pattern is valid for regex (no special chars, but DNA is safe)
            if let Ok(regex_matcher) = RegexMatcher::new(pattern_str, opts) {
                let fixed_result = fixed_matcher.bases_match(&seq);
                let regex_result = regex_matcher.bases_match(&seq);

                // They should always agree on literal patterns
                prop_assert_eq!(fixed_result, regex_result);
            }
        }
    }
}
