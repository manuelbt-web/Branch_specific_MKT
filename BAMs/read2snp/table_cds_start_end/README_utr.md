# UTR Length Extractor

A script to extract 5' and 3' UTR lengths from GFF3 files, with transcript IDs and strand information.

## Features

- Extracts UTR lengths for all mRNA features
- Handles both 5' and 3' UTRs
- Preserves strand information
- Identifies transcripts missing UTR annotations
- Provides validation statistics

## Usage

```bash
./get_utr_lengths.sh <input.gff3> <output.tab>
