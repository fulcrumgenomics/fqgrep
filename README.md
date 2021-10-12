# fqgrep

<p align="center">
  <a href="https://github.com/fulcrumgenomics/fqgrep/actions?query=workflow%3ACheck"><img src="https://github.com/fulcrumgenomics/fqgrep/workflows/Check/badge.svg" alt="Build Status"></a>
  <img src="https://img.shields.io/crates/l/fqgrep.svg" alt="license">
  <a href="https://crates.io/crates/fqgrep"><img src="https://img.shields.io/crates/v/fqgrep.svg?colorB=319e8c" alt="Version info"></a><br>
  Grep for FASTQ files.
</p>

Search a pair of fastq files for reads that match a given ref or alt sequence.

## Install

```bash 
git clone ... && cd fqgrep
cargo install --path .
```

## Usage

```
fqgrep 0.1.0
A small utility to "grep" a pair of gzipped FASTQ files, outputting the read pair if the sequence given matches either

USAGE:
    fqgrep [OPTIONS] --alt-search-sequence <alt-search-sequence> --output-dir <output-dir> --r1-fastq <r1-fastq> --r2-fastq <r2-fastq> --ref-search-sequence <ref-search-sequence>

FLAGS:
    -h, --help       
            Prints help information

    -V, --version    
            Prints version information


OPTIONS:
    -a, --alt-search-sequence <alt-search-sequence>
            The "alt" search sequence to grep for, for R2 reads the search sequence is reverse complemented

    -c, --compression-threads <compression-threads>      
            Number of threads to use for compression.
            
            Zero is probably the correct answer if the searched sequence is relatively rare. If it's common and many
            reads will be written, it may be worth decreasing the number of `treads_for_matching` by a small number and
            making this number > 2 [default: 0]
    -o, --output-dir <output-dir>                        
            The output directory to write to (must exist)
            
            `fqgrep` will write a pair of fastqs for the alt matched reads to:
            {output_dir}/ALT_{sample_name}R{1,2}.fastq.gz and the fastqs for the ref matched reads to
            {output_dir}/REF_{sample_name}R{1,2}.fastq.gz
    -1, --r1-fastq <r1-fastq>                            
            Path to the input R1 gzipped FASTQ

    -2, --r2-fastq <r2-fastq>                            
            Path to the input R2 gzipped FASTQ

    -r, --ref-search-sequence <ref-search-sequence>
            The "ref" search sequence to grep for, for R2 reads the search sequence is reverse complemented

    -t, --threads-for-matching <threads-for-matching>
            The number of threads to use for matching reads against pattern
            
            Keep in mind that there is are two reader threads and two writer threads created by default. [default: 32]
```