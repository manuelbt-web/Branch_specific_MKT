# read2snp — SNP Calling and CDS Extraction Pipeline

This folder contains scripts to call SNPs from BAM files and extract CDS
(coding sequence) regions from consensus transcripts, as preparation for
McDonald-Kreitman test (MKT) analysis.

---

## Pipeline Overview

```
BAM files + reference transcriptome
         │
         ▼  read2snp.sbatch                               [SLURM]
         SNP tables per transcript (reads2snp format)

Transcript FASTA + GFF3 annotation
         │
         ▼  Step 1 — table_cds_start_end/compute_cds_positions.sh
         │   transcript_lengths.tab   (sequence ID → length in bp)
         │   utr_lengths.tab          (5'/3' UTR lengths per transcript)
         │   cds_positions.tab        (CDS start/end per transcript)
         │
         ▼  Step 2 — extract_cds/extraction_of_covered_contigs.py
             retained_genes.txt       (genes passing completeness filter)
             cds_sequences/           (one CDS FASTA per retained gene)
```

---

## Step 0 — SNP calling with read2snp

**Script:** `read2snp.sbatch`
**Cluster:** SLURM — edit the cluster-specific lines before submission.

```bash
sbatch read2snp.sbatch \
    /path/to/bam_list.txt \
    /path/to/reference.fasta \
    /path/to/output/
```

**Inputs:**
- `bam_list.txt` — plain-text list of BAM file paths (one per line)
- `reference.fasta` — reference transcriptome used for mapping
- `output/` — directory for SNP tables

**Key parameter inside the script:**
- `FIS_VALUE` — inbreeding coefficient; edit this variable at the top of
  `read2snp.sbatch` before submission.

**Note:** `reads2snp` must be installed and accessible in `$PATH`.
See [https://lbbe-software.univ-lyon1.fr/reads2snp](https://lbbe-software.univ-lyon1.fr/reads2snp).

---

## Step 1 — Compute CDS positions

**Script:** `table_cds_start_end/compute_cds_positions.sh`

Computes CDS start/end positions from the transcript FASTA and GFF3
annotation in a single command.  No pre-processing needed.

```bash
bash table_cds_start_end/compute_cds_positions.sh \
    --fasta   reference.fasta \
    --gff3    annotation.gff3 \
    --outdir  results/cds/
```

**What it does internally:**
1. Measures every transcript length from the FASTA.
2. Extracts 5' and 3' UTR lengths from the GFF3.
3. Applies `cds_start = FLANK + 5'UTR + 1` and
   `cds_end = length − FLANK − 3'UTR` (default FLANK = 200 bp).

**Outputs:**
- `results/cds/transcript_lengths.tab` — sequence lengths (intermediate)
- `results/cds/utr_lengths.tab` — UTR lengths (intermediate)
- `results/cds/cds_positions.tab` — **CDS coordinates per transcript**

See [table_cds_start_end/README.md](table_cds_start_end/README.md) for full
documentation, including how to change the flanking value.

---

## Step 2 — Extract CDS sequences and filter by completeness

**Script:** `extract_cds/extraction_of_covered_contigs.py`

Filters genes based on the completeness of their CDS sequences (fraction of
non-N bases) across individuals.

```bash
python extract_cds/extraction_of_covered_contigs.py \
    --fasta            reference.fasta \
    --cds-table        results/cds/cds_positions.tab \
    --min-completeness 0.7 \
    --min-fraction     0.5 \
    --gene-list        results/retained_genes.txt \
    --cds-dir          results/cds_sequences/
```

**Parameters to adjust:**

| Parameter | Meaning | Example |
|---|---|---|
| `--min-completeness` | Min fraction of non-N bases in CDS per sequence | `0.7` (70%) |
| `--min-fraction` | Min fraction of individuals with a complete CDS | `0.5` (50%) |

**Outputs:**
- `results/retained_genes.txt` — list of genes passing the filter
- `results/cds_sequences/<GENE_ID>_CDS.fasta` — one file per retained gene

See [extract_cds/README.md](extract_cds/README.md) for full documentation.

---

## Full example

```bash
# 0. Run SNP calling (on SLURM cluster):
sbatch read2snp.sbatch bam_list.txt reference.fasta snp_output/

# 1. Compute CDS positions:
bash table_cds_start_end/compute_cds_positions.sh \
    --fasta   reference.fasta \
    --gff3    annotation.gff3 \
    --outdir  results/cds/

# 2. Extract CDS and filter by completeness:
python extract_cds/extraction_of_covered_contigs.py \
    --fasta            reference.fasta \
    --cds-table        results/cds/cds_positions.tab \
    --min-completeness 0.7 \
    --min-fraction     0.5 \
    --gene-list        results/retained_genes.txt \
    --cds-dir          results/cds_sequences/
```

---

## Dependencies

| Tool | Version | Used in |
|---|---|---|
| reads2snp | 2.0 | Step 0 — SNP calling |
| GNU awk (gawk) | any | Step 1 — GFF3 parsing |
| Python | ≥ 3.7 | Step 2 |
| Biopython | ≥ 1.79 | Step 2 |

Install Biopython:
```bash
pip install biopython
# or
conda install -c conda-forge biopython
```

---

## Folder structure

```
read2snp/
├── read2snp.sbatch                              Step 0: SNP calling (SLURM)
├── table_cds_start_end/
│   ├── compute_cds_positions.sh                 Step 1: CDS coordinates
│   └── README.md
└── extract_cds/
    ├── extraction_of_covered_contigs.py         Step 2: CDS extraction
    └── README.md
```
