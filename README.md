# fqgrep

<p align="center">
  <a href="https://github.com/fulcrumgenomics/fqgrep/actions?query=workflow%3ACheck"><img src="https://github.com/fulcrumgenomics/fqgrep/workflows/Check/badge.svg" alt="Build Status"></a>
  <img src="https://img.shields.io/crates/l/fqgrep.svg" alt="license">
  <a href="https://crates.io/crates/fqgrep"><img src="https://img.shields.io/crates/v/fqgrep.svg?colorB=319e8c" alt="Version info"></a>
  <a href="http://bioconda.github.io/recipes/fqgrep/README.html"><img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat" alt="Install with bioconda"></a>
  <a href="https://doi.org/10.5281/zenodo.14985002"><img src="https://zenodo.org/badge/416465549.svg" alt="zenodo"></a>
  <br>
  Grep for FASTQ files.
</p>

Search a pair of fastq files for reads that match a given ref or alt sequence.

<p>
<a href="https://fulcrumgenomics.com"><img src=".github/logos/fulcrumgenomics.svg" alt="Fulcrum Genomics" height="100"/></a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with fqgrep and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

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

To match a single pattern in a single FASTQ:

```console
fqgrep 'GACGAGATTA' /data/test.fastq.gz
```

To match two patterns across two FASTQs treating them as paired end:

```console
fqgrep  --paired --threads 10 --regexp 'GACGAGATTA' --regexp 'GACGTGATTA' /data/testR1.fastq.gz /data/testR2.fastq.gz
```

## Help

See the following for usage:

<!-- start usage -->
```console

OVERVIEW

The fqgrep utility searches any given input FASTQ files, selecting records whose bases match one or more patterns.  By default, a pattern matches the bases in a FASTQ record if the regular expression (RE) in the pattern matches the bases.  An empty expression matches every line.  Each FASTQ record that matches at least one of the patterns is written to the standard output.

INPUT COMPRESSION

By default, the input files are assumed to be uncompressed with the following exceptions: (1) If the input files are real files and end with .gz or .bgz, they are assumed to be GZIP compressed, or (2) if they end with .fastq or .fq, they are assumed to be uncompressed, or (3) if the -Z/--decompress option is specified then any unrecongized inputs (including standard input) are assumed to be GZIP compressed.

THREADS

The --threads option controls the number of threads used to search the reads. Independently, for single end reads or interleaved paired end reads, a single thread will be used to read each input FASTQ.  For paired end reads across pairs of FASTQs, two threads will be used to read the FASTQs for each end of a pair.  Finally, a single thread will be created for the writer.

EXIT STATUS

The fqgrep utility exits with one of the following values: 0 if one or more lines were selected, 1 if no lines were selected, and >1 if an error occurred.

Usage: fqgrep [OPTIONS] [ARGS]...

Arguments:
  [ARGS]...
          The first argument is the pattern to match, with the remaining arguments containing the files to match.  If -e is given, then all the arguments are files to match. Use standard input if either no files are given or - is given

Options:
  -t, --threads <THREADS>
          The number of threads to use for matching reads against pattern.  See the full usage for threads specific to reading and writing
          
          [default: 10]

      --color <COLOR>
          Mark up the matching text.  The possible values of when are “never”, “always” and “auto”
          
          [default: never]
          [possible values: never, always, auto]

  -c, --count
          Only a count of selected lines is written to standard output

  -e, --regexp <REGEXP>
          Specify a pattern used during the search of the input: an input line is selected if it matches any of the specified patterns.  This option is most useful when multiple -e options are used to specify multiple patterns

  -F, --fixed-strings
          Interpret pattern as a set of fixed strings

  -f, --file <FILE>
          Read one or more newline separated patterns from file.  Empty pattern lines match every input line.  Newlines are not considered part of a pattern.  If file is empty, nothing is matched

  -v
          Selected lines are those not matching any of the specified patterns

  -Z, --decompress
          Assume all unrecognized inputs are GZIP compressed

      --paired
          Treat the input files as paired.  The number of input files must be a multiple of two, with the first file being R1, second R2, third R1, fourth R2, and so on.  If the pattern matches either R1 or R2, then both R1 and R2 will be output (interleaved).  If the input is standard input, then treat the input as interlaved paired end reads

      --reverse-complement
          Search the reverse complement for matches

      --progress
          Write progress information

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```
<!-- end usage -->
