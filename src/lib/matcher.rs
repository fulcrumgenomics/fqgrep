use crate::reverse_complement;
use crate::DNA_BASES;
use anyhow::{bail, Context, Result};
use bstr::ByteSlice;
use regex::bytes::{Regex, RegexBuilder, RegexSet, RegexSetBuilder};
use seq_io::fastq::{OwnedRecord, Record};

#[derive(Copy, Clone, Debug)]
pub struct MatcherOpts {
    pub invert_match: bool,
    pub reverse_complement: bool,
}

pub trait Matcher {
    fn opts(&self) -> MatcherOpts;

    fn bases_match(&self, bases: &[u8]) -> bool;

    fn read_match(&self, read: &OwnedRecord) -> bool {
        if self.opts().invert_match {
            self.bases_match(read.seq())
                && (self.opts().reverse_complement
                    && self.bases_match(&reverse_complement(read.seq())))
        } else {
            self.bases_match(read.seq())
                || (self.opts().reverse_complement
                    && self.bases_match(&reverse_complement(read.seq())))
        }
    }
}

pub struct FixedStringMatcher {
    pattern: Vec<u8>,
    opts: MatcherOpts,
}

impl Matcher for FixedStringMatcher {
    fn bases_match(&self, bases: &[u8]) -> bool {
        bases.find(&self.pattern).is_some() != self.opts.invert_match
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

    pub fn validate(pattern: &str) -> Result<()> {
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
}

pub struct FixedStringSetMatcher {
    patterns: Vec<Vec<u8>>,
    opts: MatcherOpts,
}

impl Matcher for FixedStringSetMatcher {
    fn bases_match(&self, bases: &[u8]) -> bool {
        self.patterns
            .iter()
            .any(|pattern| bases.find(&pattern).is_some())
            != self.opts.invert_match
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
    fn opts(&self) -> MatcherOpts {
        self.opts
    }
}

pub struct RegexSetMatcher {
    regexes: RegexSet,
    opts: MatcherOpts,
}

impl RegexSetMatcher {
    pub fn new<I, S>(patterns: I, opts: MatcherOpts) -> Self
    where
        S: AsRef<str>,
        I: IntoIterator<Item = S>,
    {
        let regexes = RegexSetBuilder::new(patterns).build().unwrap();
        Self { regexes, opts }
    }
}

impl Matcher for RegexSetMatcher {
    fn bases_match(&self, bases: &[u8]) -> bool {
        self.regexes.is_match(bases) != self.opts.invert_match
    }
    fn opts(&self) -> MatcherOpts {
        self.opts
    }
}
