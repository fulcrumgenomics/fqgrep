# fqgrep

<p align="center">
  <a href="https://github.com/fulcrumgenomics/fqgrep/actions?query=workflow%3ACheck"><img src="https://github.com/fulcrumgenomics/fqgrep/workflows/Check/badge.svg" alt="Build Status"></a>
  <img src="https://img.shields.io/crates/l/fqgrep.svg" alt="license">
  <a href="https://crates.io/crates/fqgrep"><img src="https://img.shields.io/crates/v/fqgrep.svg?colorB=319e8c" alt="Version info"></a>
  <a href="http://bioconda.github.io/recipes/fqgrep/README.html"><img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat" alt="Install with bioconda"></a>
  <br>
  Grep for FASTQ files.
</p>

Search fastq files for reads that match a given sequence.

## Install

### From bioconda

```console
conda install -c bioconda fqgrep
```

### From Source

```console 
git clone ... && cd fqgrep
cargo install --path .
```

## Usage

```console
fqgrep [FLAGS] [OPTIONS] [--] [args]...

FLAGS:
    -c, --count                 Only a count of selected lines is written to standard output
    -F, --fixed-strings         Interpret pattern as a set of fixed strings
    -v                          Selected lines are those not matching any of the specified patterns
    -Z, --decompress            Assume all unrecognized inputs are GZIP compressed
        --paired                Treat the input files as paired.  The number of input files must be a multiple of two,
                                with the first file being R1, second R2, third R1, fourth R2, and so on.  If the pattern
                                matches either R1 or R2, then both R1 and R2 will be output (interleaved).  If the input
                                is standard input, then treat the input as interlaved paired end reads
        --reverse-complement    Search the reverse complement for matches
        --progress              Write progress information
    -h, --help                  Prints help information
    -V, --version               Prints version information

OPTIONS:
    -t, --threads <threads>     The number of threads to use for matching reads against pattern.  See the full usage for
                                threads specific to reading and writing [default: 8]
        --color <color>         Mark up the matching text.  The possible values of when are “never”, “always” and “auto”
                                [default: never]
    -e, --regexp <regexp>...    Specify a pattern used during the search of the input: an input line is selected if it
                                matches any of the specified patterns.  This option is most useful when multiple `-e`
                                options are used to specify multiple patterns
    -f, --file <file>           Read one or more newline separated patterns from file.  Empty pattern lines match every
                                input line.  Newlines are not considered part of a pattern.  If file is empty, nothing
                                is matched

ARGS:
    <args>...    The first argument is the pattern to match, with the remaining arguments containing the files to
                 match.  If `-e` is given, then all the arguments are files to match. Use standard input if either
                 no files are given or `-` is given
```

## Example Invocation

```console
fqgrep --progress --threads 16 --color always <PATTERN> <TEST_1.fastq> <TEST_2.fq> <TEST_3>.fq.gz <TEST_3>.fastq.gz ...
```

## Help

See the following for usage:

```console
fqgrep -h
```
