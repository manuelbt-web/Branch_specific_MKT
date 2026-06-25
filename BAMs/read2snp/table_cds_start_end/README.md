# CDS Position Calculator

Computes the start and end positions of coding sequences (CDS) within
consensus transcript sequences.

Given a **transcript FASTA** and a **GFF3 annotation**, the script
automatically:
1. measures each transcript's length from the FASTA,
2. extracts UTR lengths from the GFF3, and
3. combines both to produce a CDS position table.

No pre-processing step is needed — everything runs in one command.

---

## Background — the flanking region

In this project, consensus transcripts are assembled with a fixed **flanking
region** added on both the 5' and 3' ends of the annotated gene boundaries.
This flanking provides sequence context for read mapping and SNP calling.

Because of this, the CDS does **not** start at position 1 of the consensus
sequence. Its position is calculated as:

```
cds_start = FLANK + five_prime_UTR_length + 1   (1-based, inclusive)
cds_end   = transcript_length − FLANK − three_prime_UTR_length
```

The default flanking is **200 bp** (change with `--flank` if your assembly
used a different value).

> **Note on strand:** Consensus transcriptome assemblies are oriented 5'→3'
> regardless of the genomic strand, so the formula above applies uniformly.
> The strand column in the output is extracted from the GFF3 for reference.

---

## Script: `compute_cds_positions.sh`

### Usage

```bash
# Standard usage — all inputs provided:
bash compute_cds_positions.sh \
    --fasta   sequences.fasta \
    --gff3    annotation.gff3 \
    --outdir  results/cds/

# Re-use an existing UTR table (e.g., to re-run with a different --flank):
bash compute_cds_positions.sh \
    --fasta   sequences.fasta \
    --utr     results/cds/utr_lengths.tab \
    --outdir  results/cds/

# Show all options:
bash compute_cds_positions.sh --help
```

### Arguments

| Argument | Required | Description |
|---|---|---|
| `--fasta FILE` | Always | Transcript FASTA (single- or multi-line). Lengths computed from this file. |
| `--outdir DIR` | Always | Output directory (created if absent). |
| `--gff3 FILE` | Yes (unless `--utr`) | GFF3 annotation — used to extract UTR lengths. |
| `--utr FILE` | Alternative to `--gff3` | Use an existing UTR table (skips step 2). |
| `--flank N` | No | Flanking region in bp (default: **200**). |
| `--help` | No | Show usage and exit. |

### GFF3 requirements

The GFF3 must contain `mRNA`, `five_prime_UTR`, and `three_prime_UTR` feature
types with standard attributes (`ID=`, `Parent=`).
The `transcript:` prefix in IDs is stripped automatically.
Pipe-delimited IDs (e.g. `GENE001|HOG12345`) are supported.

---

## Outputs (in `--outdir`)

| File | Description |
|---|---|
| `transcript_lengths.tab` | Step 1: sequence ID → length in bp |
| `utr_lengths.tab` | Step 2: UTR lengths per transcript |
| `cds_positions.tab` | **Final output**: CDS coordinates per transcript |

### `utr_lengths.tab` columns

| mRNA_ID | five_prime_UTR_len | three_prime_UTR_len | strand |
|---|---|---|---|
| GENE001 | 45 | 120 | + |

### `cds_positions.tab` columns

| full_transcript_id | transcript_id | cds_start | cds_end | strand | five_prime_utr | three_prime_utr | original_length |
|---|---|---|---|---|---|---|---|
| GENE001\|HOG12345 | GENE001 | 246 | 1180 | + | 45 | 120 | 1500 |

> `cds_start` and `cds_end` are **1-based, inclusive** coordinates within
> the consensus sequence.

---

## Example

```bash
bash compute_cds_positions.sh \
    --fasta   data/speltoides_consensus.fasta \
    --gff3    data/speltoides_annotation.gff3 \
    --outdir  results/cds/
```

Expected output:

```
======================================================
  compute_cds_positions.sh
======================================================
  FASTA file   : data/speltoides_consensus.fasta
  GFF3 file    : data/speltoides_annotation.gff3
  Flanking     : 200 bp
  Output dir   : results/cds/
======================================================

[10:32:01] Step 1 — Computing transcript lengths from FASTA...
[10:32:03]   Lengths table: results/cds/transcript_lengths.tab  (24581 sequences)

[10:32:03] Step 2 — Extracting UTR lengths from GFF3...
[10:32:06]   UTR table    : results/cds/utr_lengths.tab  (24581 transcripts)

[10:32:06] Step 3 — Computing CDS positions...
[10:32:07]   CDS table    : results/cds/cds_positions.tab  (24312 transcripts)

======================================================
  Done  (6s)
  Outputs: transcript_lengths.tab / utr_lengths.tab / cds_positions.tab
======================================================
```

---

## Dependencies

| Tool | Notes |
|---|---|
| GNU awk (gawk) | Required for 3-argument `match()` used in GFF3 parsing. Standard on Linux. |

Verify: `awk --version | head -1` — should print `GNU Awk …`

---

## Next step

Feed `cds_positions.tab` to
[`../extract_cds/extraction_of_covered_contigs.py`](../extract_cds/) to
extract CDS sequences and filter by completeness threshold.

```bash
python ../extract_cds/extraction_of_covered_contigs.py \
    --fasta            data/speltoides_consensus.fasta \
    --cds-table        results/cds/cds_positions.tab \
    --min-completeness 0.7 \
    --min-fraction     0.5 \
    --gene-list        results/retained_genes.txt \
    --cds-dir          results/cds_sequences/
```
