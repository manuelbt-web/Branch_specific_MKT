# CDS Extraction and Coverage Filtering

Extracts CDS sequences from a multi-FASTA consensus file and retains only
genes that are sufficiently complete across individuals.

**Script:** `extraction_of_covered_contigs.py`

---

## What it does

1. Reads the CDS position table from
   [`../table_cds_start_end/compute_cds_positions.sh`](../table_cds_start_end/).
2. For each gene in the multi-FASTA, extracts the CDS region from **every
   individual's sequence** (one FASTA entry per individual).
3. Computes **sequence completeness** = fraction of non-N bases in the CDS.
4. Keeps a gene if at least `--min-fraction` of its individual sequences meet
   `--min-completeness`.
5. Outputs a gene list and one CDS FASTA file per retained gene.

> **"Coverage" here means sequence completeness (% non-N bases), not
> sequencing depth.**  N bases are gaps or low-confidence positions in the
> consensus sequence, typically arising from insufficient read support at
> that position during consensus calling.

---

## Installation

Requires Python ≥ 3.7 and [Biopython](https://biopython.org/).

```bash
# With pip:
pip install biopython

# With conda / mamba:
conda install -c conda-forge biopython
```

---

## Usage

```bash
python extraction_of_covered_contigs.py \
    --fasta            speltoides_consensus.fasta \
    --cds-table        results/cds/cds_positions.tab \
    --min-completeness 0.7 \
    --min-fraction     0.5 \
    --gene-list        results/retained_genes.txt \
    --cds-dir          results/cds_sequences/
```

Show help:

```bash
python extraction_of_covered_contigs.py --help
```

---

## Arguments

| Argument | Required | Description |
|---|---|---|
| `--fasta FILE` | Yes | Multi-FASTA of consensus sequences (all individuals) |
| `--cds-table FILE` | Yes | CDS position table from `compute_cds_positions.sh` |
| `--min-completeness FLOAT` | Yes | Min fraction of non-N bases per sequence (0–1) |
| `--min-fraction FLOAT` | Yes | Min fraction of complete sequences per gene (0–1) |
| `--gene-list FILE` | Yes | Output file: retained gene IDs (one per line) |
| `--cds-dir DIR` | Yes | Output directory: one `<GENE_ID>_CDS.fasta` per gene |

---

## Parameters explained

### `--min-completeness`

Threshold for a single individual's sequence to be considered "complete".

- `0.7` → at least 70% of the CDS positions must be non-N
- `1.0` → no missing bases at all (very strict)
- `0.5` → at least half the positions are called

Choose based on your tolerance for missing data.  Values between 0.5 and 0.8
are typical for RNA-seq consensus data.

### `--min-fraction`

Fraction of individuals that must have a complete sequence for the gene to be
retained.

- `0.5` → at least half the individuals must meet `--min-completeness`
- `1.0` → all individuals must have a complete CDS (very strict)
- `0.0` → any gene with at least one complete sequence is kept

Adjust based on how many individuals you require for downstream analysis
(e.g., MKT requires genotypes across a sufficient number of individuals).

---

## FASTA format expected

The multi-FASTA must group sequences by gene using the pipe `|` separator in
the header:

```
>GENE001|sample_Ae_spe_1
ATGCATGC...
>GENE001|sample_Ae_spe_2
ATGNATGC...
>GENE002|sample_Ae_spe_1
ATGCATGC...
```

The first `|`-delimited token (`GENE001`) is used to group sequences per gene.

---

## Outputs

| Output | Description |
|---|---|
| `<gene-list>` | Plain text file listing retained gene IDs, one per line |
| `<cds-dir>/<GENE_ID>_CDS.fasta` | CDS sequences for each retained gene |

**CDS FASTA header format:**

```
>GENE001|sample_Ae_spe_1|CDS
```

The `|CDS` suffix is appended to mark that the sequence is the extracted
CDS region (not the full transcript).

---

## Example output

```
Loading CDS table : results/cds/cds_positions.tab
  24312 genes with CDS annotations
Parsing FASTA     : speltoides_consensus.fasta
  218808 sequences across 24312 genes

Results:
  Genes with CDS annotations : 24312
  Genes found in FASTA       : 24312
  Genes evaluated            : 23891
  Genes retained             : 17245

Filters applied:
  --min-completeness 70%  (each sequence must have >= 70% non-N bases in CDS)
  --min-fraction     50%  (at least 50% of sequences per gene must be complete)

Gene list written : results/retained_genes.txt
CDS FASTA files   : 17245 files written to results/cds_sequences/
```

---

## Connection to the MKT pipeline

The retained genes list and CDS FASTA files are the input to the **read2snp**
SNP-calling step ([`../read2snp.sbatch`](../read2snp.sbatch)), which generates
per-gene SNP tables used in the McDonald-Kreitman test.

---

## Dependencies

| Package | Version | Install |
|---|---|---|
| Python | ≥ 3.7 | `python --version` |
| Biopython | ≥ 1.79 | `pip install biopython` |
