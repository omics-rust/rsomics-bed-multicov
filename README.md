# rsomics-bed-multicov

Count reads from multiple BAM files overlapping BED intervals — a `bedtools multicov` equivalent written in Rust.

## Usage

```
rsomics-bed-multicov --bams s1.bam s2.bam s3.bam --bed peaks.bed
```

Output is the BED file with one read-count column appended per BAM.

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--bams` | required | Input BAM files (indexed, one or more) |
| `--bed` | required | BED file defining the regions |
| `--min-mapq` | 0 | Minimum mapping quality |
| `-D` | off | Include duplicate reads (FLAG 0x400) |
| `-F` | off | Include failed-QC reads (FLAG 0x200) |
| `-p` | off | Only count properly-paired reads |
| `-f` | 1e-9 | Minimum overlap as fraction of BED region length |
| `-r` | off | Require reciprocal overlap |
| `-s` | off | Require same strandedness (BED col 6 vs read) |
| `-S` | off | Require opposite strandedness |

## Requirements

Each input BAM must be coordinate-sorted and indexed (`.bam.bai`).

## Origin

This crate is a Rust reimplementation of `bedtools multicov`, informed by:
- The `BEDTools` documentation and man page
- Quinlan & Hall (2010). *BEDTools: a flexible suite of utilities for comparing genomic features.* Bioinformatics 26(6): 841–842. DOI: 10.1093/bioinformatics/btq033
- Black-box behavior testing against bedtools 2.31.1

License: MIT OR Apache-2.0
Upstream credit: bedtools <https://github.com/arq5x/bedtools2> (MIT)
