# CDS Position Calculator

## Description
This Bash script calculates Coding Sequence (CDS) positions for transcripts based on:
- Transcript length statistics
- Annotated UTR (Untranslated Region) lengths
- Strand orientation information

## Key Features
- **Strand-aware calculations**: Correctly handles + and - strand transcripts
- **Position validation**: Checks for invalid CDS ranges (start ≥ end)
- **Error reporting**: Identifies transcripts missing from UTR annotations
- **Flexible input**: Handles transcript IDs with pipe (`|`) delimiters
- **Header preservation**: Maintains original full transcript IDs

## Usage
```bash
./CDS_positions.sh <transcript_stats> <utr_lengths> <output_file>


