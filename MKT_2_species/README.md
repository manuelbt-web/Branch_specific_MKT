# MKT 2-Species Pipeline

This pipeline implements the two-species McDonald-Kreitman Test (MKT) using
the dNdSpNpS framework.

Starting from per-gene CDS FASTA files containing multiple individuals of a
focal species, the pipeline adds an outgroup ortholog, removes non-homologous
sequence fragments, produces codon-aware multiple sequence alignments,
standardises FASTA headers, computes divergence statistics with dNdSpiNpiS,
and computes polymorphism statistics with EggLib.

---

## Pipeline overview

```
Input: per-gene CDS FASTA files
       (one file per gene, N individuals of the focal species)

  step1_add_outgroup/
  └── add_outgroup.sh
        ↓  appends outgroup ortholog via RBH table
        ↓  output: one FASTA per gene with outgroup at end

  step2_MACSE_trim/
  └── MACSE_trim_non_homologous.sbatch     [SLURM array]
        ↓  removes non-homologous sequence fragments (trimNonHomologousFragments)
        ↓  output: trimmed NT and AA sequences

  step3_MACSE_align/
  └── MACSE_aligning_sequences.sbatch     [SLURM array]
        ↓  codon-aware multiple sequence alignment (alignSequences)
        ↓  output: aligned NT (with ! for frameshifts) and AA sequences

  step4_prep_dNdS/
  └── input_preparation_for_dNdSpNpS.py
        ↓  reorders sequences (focal first, outgroup last)
        ↓  standardises FASTA headers for dNdSpNpS
        ↓  output: analysis-ready FASTA files

  ┌─────────────────── step4_out/ ────────────────────┐
  │  (Steps 5 and 6 are independent — run in parallel) │
  └──────────────┬────────────────────────┬────────────┘
                 ↓                        ↓

  step5_run_dNdSpNpS/               step6_polymorphism/
  └── dNdSpiNpiS.sbatch  [SLURM]    ├── remove_outgroup_sequence.py
        ↓  divergence per gene       │     ↓  strips outgroup sequence
        ↓  dN, dS, pN, pS,          └── polymorphisms_stats.py
        ↓  dN/dS, pN/pS, DoS, NI         ↓  S, Pi, D, thetaW, pi_S, pi_NS
        ↓  output: per-gene .out          ↓  folded SFS (EggLib CodingDiversity)
                                          ↓  SDM estimation (integrated)
                                          ↓  PnNeutral, deleterious
```

---

## Directory structure

```
MKT_2_species/
├── README.md                                   ← this file
├── step1_add_outgroup/
│   ├── README.md
│   └── add_outgroup.sh
├── step2_MACSE_trim/
│   ├── README.md
│   └── MACSE_trim_non_homologous.sbatch
├── step3_MACSE_align/
│   ├── README.md
│   └── MACSE_aligning_sequences.sbatch
├── step4_prep_dNdS/
│   ├── README.md
│   └── input_preparation_for_dNdSpNpS.py
├── step5_run_dNdSpNpS/
│   ├── README.md
│   └── dNdSpiNpiS.sbatch
├── step6_polymorphism/
│   ├── README.md
│   ├── remove_outgroup_sequence.py
│   ├── polymorphisms_stats.py
│   └── estimation_of_SDM.py
└── step7_MKT/
    ├── README.md
    ├── merge_tables.py
    └── MKT_analysis.Rmd
```

---

## Prerequisites

| Tool | Version | Used in | Installation |
|---|---|---|---|
| **seqkit** | ≥ 2.0 | Step 1 | `conda install -c bioconda seqkit` |
| **Java** | ≥ 8 | Steps 2–3 | `conda install -c conda-forge openjdk` or cluster module |
| **MACSE** | v2.07 (inside MACSE_V2_PIPELINES v12.01) | Steps 2–3 | see below |
| **Python** | ≥ 3.7 | Steps 4, 6 | system / conda |
| **BioPython** | ≥ 1.79 | Step 4 | `pip install biopython` |
| **dNdSpiNpiS** | v1.0 | Step 5 | Download from PopPhyl (see Step 5 README) |
| **EggLib** | 3.5.2 | Steps 6b | Pixi — see Step 6 README (optional) |
| **pandas** | any | Step 7a | `pip install pandas` |

### Install MACSE (Steps 2 and 3)

MACSE is distributed as a JAR file bundled inside the MACSE_V2_PIPELINES
repository.  Download the v12.01 archive and extract the JAR once:

```bash
wget https://github.com/ranwez/MACSE_V2_PIPELINES/archive/refs/tags/v12.01.tar.gz
tar xzf v12.01.tar.gz
MACSE_JAR=$(find MACSE_V2_PIPELINES-12.01/ -name "*.jar" | head -1)
echo "MACSE_JAR=$MACSE_JAR"   # note this path for Steps 2 and 3
```

---

## Data flow between steps

| Step | Input directory | Output directory | Key files |
|---|---|---|---|
| 1 | `cds_per_gene/` | `step1_out/` | `<gene>.fasta` (+ outgroup) |
| 2 | `step1_out/` | `step2_out/` | `NT/<gene>_NT.fasta`, `AA/<gene>_AA.fasta` |
| 3 | `step2_out/NT/` | `step3_out/` | `<gene>_aligned_NT.fasta`, `<gene>_aligned_AA.fasta` |
| 4 | `step3_out/` | `step4_out/` | `<gene>.fasta` (standardised headers) |
| 5 | `step4_out/` | `step5_out/` | `<gene>.out` (dN, dS, pN, pS, DoS, NI) |
| 6a | `step4_out/` | `step6a_no_outgroup/` | `<gene>.fasta` (focal only) |
| 6b | `step6a_no_outgroup/` | `step6_stats/` | `polymorphism_stats.tsv` (stats + SDM), `folded_sfs.tsv` (optional) |

---

## Quick-start (complete example)

```bash
# ── 0. Prerequisites ──────────────────────────────────────────────────────────
# Download MACSE once and note the JAR path:
wget https://github.com/ranwez/MACSE_V2_PIPELINES/archive/refs/tags/v12.01.tar.gz
tar xzf v12.01.tar.gz
MACSE_JAR=$(find MACSE_V2_PIPELINES-12.01/ -name "*.jar" | head -1)

# ── Step 1 ────────────────────────────────────────────────────────────────────
bash step1_add_outgroup/add_outgroup.sh \
    -i cds_per_gene/ \
    -o step1_out/ \
    -g H_vulgare_CDS.fasta \
    -r RBH_Ae_speltoides_H_vulgare.tab

# ── Step 2 ────────────────────────────────────────────────────────────────────
# Edit [CLUSTER] lines in step2_MACSE_trim/MACSE_trim_non_homologous.sbatch first
N=$(ls step1_out/*.fasta | wc -l)
sbatch --array=0-$((N-1))%15 step2_MACSE_trim/MACSE_trim_non_homologous.sbatch \
    -i step1_out/ \
    -o step2_out/ \
    -j "$MACSE_JAR"

# Wait for all Step 2 jobs to complete (squeue -u $USER), then:

# ── Step 3 ────────────────────────────────────────────────────────────────────
# Edit [CLUSTER] lines in step3_MACSE_align/MACSE_aligning_sequences.sbatch first
N=$(ls step2_out/NT/*.fasta | wc -l)
sbatch --array=0-$((N-1))%15 step3_MACSE_align/MACSE_aligning_sequences.sbatch \
    -i step2_out/NT/ \
    -o step3_out/ \
    -j "$MACSE_JAR"

# Wait for all Step 3 jobs to complete, then:

# ── Step 4 ────────────────────────────────────────────────────────────────────
python step4_prep_dNdS/input_preparation_for_dNdSpNpS.py \
    --input-dir  step3_out/ \
    --species    Ae_speltoides \
    --outgroup   H_vulgare \
    --output-dir step4_out/ \
    --threads    4

# ── Step 5 — Divergence (dNdSpiNpiS) ────────────────────────────────────────
# Edit [CLUSTER] lines in step5_run_dNdSpNpS/dNdSpiNpiS.sbatch first
N=$(ls step4_out/*.fasta | wc -l)
sbatch --array=0-$((N-1))%15 step5_run_dNdSpNpS/dNdSpiNpiS.sbatch \
    -i step4_out/ \
    -o step5_out/ \
    -b /path/to/dNdSpiNpiS_1.0 \
    -g Ae_speltoides \
    -k H_vulgare

# Merge per-gene results into one TSV:
head -1 "$(ls step5_out/*.out | head -1)" > all_genes_dNdSpNpS.tsv
for f in step5_out/*.out; do tail -n +2 "$f"; done >> all_genes_dNdSpNpS.tsv

# ── Step 6 — Polymorphism (EggLib) ───────────────────────────────────────────
# Step 6a: remove outgroup
python step6_polymorphism/remove_outgroup_sequence.py \
    --input-dir  step4_out/ \
    --output-dir step6a_no_outgroup/

# Step 6b: compute statistics + SDM estimation (all-in-one)
# SDM estimation (T=0.15) runs automatically at the end.
# --sfs-output keeps the intermediate folded SFS for traceability.
python step6_polymorphism/polymorphisms_stats.py \
    --input-dir  step6a_no_outgroup/ \
    --output     step6_stats/polymorphism_stats.tsv \
    --sfs-output step6_stats/folded_sfs.tsv \
    --freq-cutoff 0.15 \
    --threads    4
```

Each step's subdirectory contains a dedicated `README.md` with full details
on arguments, input/output formats, and troubleshooting.
