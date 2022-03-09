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
file should hve the `MC` and `ms` tags added. This can be done by piping the
unsorted bam through `samtools fixmate`.

The UMI should either be present in the RX tag. But because the standard Illumina BCLconvert program attaches the UMI sequence to the readname it is also possible to 

### Output
All records are written in the same order to a file or stdout with the flag
adjusted to reflect duplication status. Additional aux tags can be added
depending on configuration.

### Options
```
    -b, --bam <BAM>
            The input bam file. rumidup reads from stdin when omitted

    -d, --umi-distance <UMI_DISTANCE>
            UMI distance. The maximin hamming distance between the UMI seqeunces used to consider
            read(pairs) to be duplicates [default: 1]

    -f, --force
            Ignore previous duplicate marking applied to BAM file. This information is extracted
            from the header. Use --force to redoing duplicate marking

    -h, --help
            Print help information

    -k, --keep-readname
            Don't remove UMI from read name

    -m, --metrics <METRICS>
            The duplication metrics file, if missing metrics will be written to stderr

    -o, --output <OUTPUT>
            The output bam file. rumidup writes to stdout when omitted

    -p, --pixel-distance <PIXEL_DISTANCE>
            Optical duplicate pixel distance. Maximum distance between clusters to consider them
            optical duplicates use 100 for HiSeq/NextSeq, 2500 for NovaSeq. use 0 to disable optical
            duplicate counting Only affects metrics [default: 100]

    -V, --version
            Print version information

    -x, --no-original-tag
            When a UMI is corrected the original tag is writtento the OX field. use --no-original-
            tag to suppress this

```

### Example
Because marking duplicates only adds information it's convenient to pipe `rumidup` after the sorting step. Default output for `rumidup` is current BAM, but can of course be piped into `samtools view` for conversion to CRAM.


```
bwa mem ref r1.fastq.gz r2.fastq.gz | \
   samtools fixmate -m - - | \
   samtools sort -@4 -m4g - | \
   rumidup --pixel-distance 2500 --metrics metrics.txt -o mdup.bam
```


## Method
rumidup reads all reads that map to the same genomic location. It then determines
the reads fragment coordinates with the help of the mapped position and the
CIGAR line information:
 - Unpaired reads take the 5' sequencing start and the mapped fragment end based
 - Paired-end reads take the 5' sequecning start of both reads

The reads are then grouped on these coordinate pairs for each group continues
to the markduplicate phase. In this phase the UMI sequences are collected.
Starting with the most common UMI sequence all others UMIs that are within
`umi-distance` hamming distance are considered identical. The read with the
highest score (base qualities) will be the candiate. All others will be marked
as duplicates. In the default case when a UMI is different from the candidate
the `RX` tag will be fixed and an `OX` tag will be added to the record.


## Limitations
Because the program works by parsing CIGAR lines and only considering the 5'
sequencing start coordinate as the read position in paired end mode rumidup is
unaffacted by quality trimming and soft-clipping. For unpaired reads it can be
assumed that soft clipping is identical across duplicates, but quality trimming
would not. This could result in shifting of read start position so quality
trimming is not recommended.

## Why
Currently the only tool capable of doing UMI aware markduplicates is in the [Picard suite](https://broadinstitute.github.io/picard/command-line-overview.html#UmiAwareMarkDuplicatesWithMateCigar). This tools has been listed as experimental for a long time and suffers from multiple issues:
  - A bug that leads to inconsistent flagging of duplicates of both end of a paired-end sequence.
  - Poor performance

**Why rumidup?**
Naming is hard. If you squint you can see letters from. Rust/UMI/duplicate/illumina.

