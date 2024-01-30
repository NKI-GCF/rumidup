# rumidup
UMI aware mark (optical) duplicates in NGS BAM

## Installation
### From source
Create a working rust environment using e.g. [rustup](https://rustup.rs/).
Clone the repo
Run `cargo build --release`

The resulting (statically linked) binary can be found in ./target/release/rumidup

### Using bioconda
It is our intention to publish to bioconda soon.


### Input
Coordinate sorted BAM (file or stdin). When paired end reads are used the BAM
file should hve the `MC` and `ms` tags added. This can usually be done by
piping the usually name-sorted bam/sam output from the aligner through
`samtools fixmate`.

The UMI should be present in the BAM `RX` tag. Because the standard Illumina
BCLconvert program attaches the UMI sequence to the read name it is also
possible to extract the UMI from the read name and (optionally) add a `RX` AUX
tags.

### Output
All records are written in the same order to a file or stdout with the BAM flag
adjusted to reflect duplication status. Additional AUX tags can be added
depending on configuration.

### Options
```
    -b, --bam <FILE>              The input bam file. rumidup reads from stdin when omitted
    -o, --output <FILE>           The output bam file. rumidup writes to stdout when omitted
    -m, --metrics <FILE>          The duplication metrics file, if missing metrics will be written
                                  to stderr
    -p, --pixel-distance <INT>    Optical duplicate pixel distance. Maximum distance between
                                  clusters to consider them optical duplicates use 100 for
                                  HiSeq/NextSeq, 2500 for NovaSeq. use 0 to disable optical
                                  duplicate counting. Only affects metrics [default: 100]
    -d, --umi-distance <INT>      UMI distance. The maximin hamming distance between the UMI
                                  seqeunces used to consider read(pairs) to be duplicates [default:
                                  1]
    -x, --no-original-tag         When a UMI is corrected the original tag is writtento the OX
                                  field. use --no-original-tag to suppress this
    -k, --keep-readname           Don't remove UMI from read name
    -f, --force                   Ignore previous duplicate marking applied to BAM file. This
                                  information is extracted from the header. Use --force to redoing
                                  duplicate marking
        --no-pg                   Do not add PG tag to header
        --force-compression       Force compression when writing to stdout. Normally when writing to
                                  a pipe the BAM stream is written without compression. Toggle this
                                  flag when redirecting to a file and force compression of the
                                  stream
    -h, --help                    Print help information
    -V, --version                 Print version information
```

### Example
Because marking duplicates only adds information it's convenient to pipe
`rumidup` after the sorting step. Default output for `rumidup` is currently BAM,
but can of course be piped into `samtools view` for conversion to CRAM.


```
bwa mem ref r1.fastq.gz r2.fastq.gz | \
   samtools fixmate -m - - | \
   samtools sort -@4 -m4g - | \
   rumidup --pixel-distance 2500 --metrics metrics.txt -o mdup.bam
```


## Method
`rumidup` collects all reads that map to the same genomic location. It then determines
the reads fragment coordinates with the help of the mapped position and the
CIGAR line information:
 - Unpaired reads take the 5' sequencing start and the mapped fragment end base
 - Paired-end reads take the 5' sequencing start of both reads

The reads are then grouped on these calculated coordinate pairs and for each
group continues to a markduplicate phase. In this phase the UMI sequences are
collected. Starting with the most common UMI sequence all other UMIs that are
within `umi-distance` hamming distance are considered identical. The read with
the highest score (sum of base qualities, including mate score if available)
will be the candidate. All others will be marked as duplicates. In the default
case when a UMI is within `umi-distane` from the candidate the `RX` tag will be
fixed and an `OX` tag will be added to the record that contains the original
UMI sequence.


## Limitations
Because the program works by parsing CIGAR lines and only considers the 5'
sequencing start coordinates as the fragment coordinates in paired end mode
`rumidup` is unaffected by quality trimming and soft-clipping. For unpaired reads
the read end is also used for grouping. Soft clipping would likely be
identical across duplicates, but quality trimming would not. This could result
in shifting of coordinates and failing to group possible PCR duplicates so
quality trimming is not recommended.

## Why
Currently the only tool capable of doing UMI aware markduplicates is in the [Picard suite](https://broadinstitute.github.io/picard/command-line-overview.html#UmiAwareMarkDuplicatesWithMateCigar). This tools has been listed as experimental for a long time and suffers from multiple issues:
  - A bug that leads to inconsistent flagging of duplicates of both end of a paired-end sequence [1449](https://github.com/broadinstitute/picard/issues/1449).
  - Poor performance
  - Since version 1.16 `samtools markdup` also supports UMI tags, but only allows exact matches in stead of a configurable hamming distance

**Why rumidup?**
Naming is hard. If you squint you can see letters from. Rust/UMI/duplicate/Illumina.

