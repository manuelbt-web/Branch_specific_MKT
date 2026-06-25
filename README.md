# Branch-Specific McDonald-Kreitman Test Analysis

**Population genomics of positive selection in wild *Aegilops* and *Triticum* species**

---

## Overview

This repository contains all scripts and analysis pipelines for a **branch-specific McDonald-Kreitman test** (MKT) to detect genes under positive selection in four wild  species:
- *Aegilops speltoides*
- *Aegilops mutica*
- *Aegilops tauschii*
- *Triticum urartu*

The analysis identifies genes evolving under positive selection at the branch leading
to each focal species, while using *T. aestivum* orthologs as the outgroup for
divergence estimation.

---

## Context

The McDonald-Kreitman test compares rates of non-synonymous (N) and synonymous (S)
substitutions at the divergence (D) and polymorphism (P) levels. Under strict
neutrality, the ratio D_N/D_S = P_N/P_S. An excess of non-synonymous divergence
(D_N/D_S > P_N/P_S) signals positive selection.

The **branch-specific** extension isolates divergence on a specific phylogenetic
branch using a three-taxon scheme (focal species, outgroup 1, outgroup 2), allowing
gene-by-gene inference of positive selection on branches leading to each wild species.

The **imputed MKT** variant (Murga-Moreno *et al.* 2022) corrects for slightly deleterious mutations that inflate P_N, improving sensitivity.

---

## Repository Structure

```
Branch_specific_MKT/
│
├── BAMs/                              # Step 1 — BAM processing
│   ├── coverage/                      # Sequencing depth and coverage filters
│   └── read2snp/                      # Polymorphism calling (read2snp)
│
├── Orthology/                         # Step 2 — Orthology
│   ├── reciprocal_best_hits.py        # BLAST-based pairwise orthologs
│   ├── orthofinder.sbatch             # OrthoFinder (HOG-based orthogroups)
│   └── reciprocal_blast.sbatch
│
├── MKT_2_species/                     # Step 3a — Standard 2-species MKT (reference)
│   ├── step1_add_outgroup/
│   ├── step2_MACSE_trim/
│   ├── step3_MACSE_align/
│   ├── step4_prep_dNdS/
│   ├── step5_run_dNdSpNpS/
│   ├── step6_polymorphism/
│   └── step7_MKT/
│
├── MKT_3_species/                     # Step 3b — Branch-specific 3-species MKT (MAIN)
│   ├── Alignment_divergence/          # MAFFT alignment → HMMcleaner → PAL2NAL
│   ├── Alignment_polymorphism/        # MACSE alignment for focal species
│   ├── EnrichAlignment/               # Codon alignment enrichment + PAML setup
│   │   ├── orthogroup_table.py        # Build HOG → species gene table
│   │   ├── MACSE_enrichment.py
│   │   ├── filter_codon_alignment.py
│   │   └── filter_codon_alignment/
│   │       ├── divergence/            # PAML codeml branch model (dN/dS)
│   │       │   ├── codeml_inputs/
│   │       │   └── run_codeml/
│   │       └── polymorphism/          # Polymorphism statistics (egglib)
│   └── branch_specific_MKT/          # R analysis package (MAIN ENTRY POINT)
│       ├── render_analysis.R          # ← RUN THIS to generate all MKT results
│       └── branch_specific_MKT_analysis.Rmd
│
├── Recombinaison_analysis/            # Step 4 — Recombination landscape
│
├── Expression_analysis/               # Step 5 — Differential expression (DESeq2)
│   ├── render_expression.R            # ← RUN THIS per experiment
│   └── expression_analysis.Rmd
│
├── protein_protein_interactions/      # Step 6 — STRING network analysis
│   ├── render_ppi.R                   # ← RUN THIS
│   └── protein_protein_interactions.Rmd
│
├── 2DNS/                              # Step 7 — Gene set MK test (KEGG pathways)
│   ├── render_2dns.R                  # ← RUN THIS
│   └── 2dns_analysis.Rmd
│
└── Accelerated_Protein_Evolution_analysis/
```

---

## Prerequisites

### System requirements

| Tool | Version | Purpose |
|------|---------|---------|
| R | ≥ 4.2 | Statistical analyses |
| Python | ≥ 3.8 | Pre-processing scripts |
| MACSE | 2.x | Codon-aware alignment |
| MAFFT | ≥ 7 | Multiple sequence alignment |
| PAML | 4.x | Branch-specific dN/dS (codeml) |
| HMMcleaner | 1.x | Alignment cleaning |
| read2snp | 2.x | Polymorphism calling from BAMs |
| BLAST+ | ≥ 2.12 | Reciprocal best-hit orthologs |
| OrthoFinder | ≥ 2.5 | Hierarchical orthogroup inference |

### R packages

```r
# CRAN
install.packages(c(
  "tidyverse",    # data manipulation and visualization
  "kableExtra",   # HTML table formatting
  "patchwork",    # figure composition
  "rmarkdown",    # report rendering
  "UpSetR",       # UpSet plots (branch_specific_MKT)
  "DT"            # interactive tables (branch_specific_MKT)
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")   # differential expression (Expression_analysis)
```

### Python packages

```bash
pip install biopython pandas numpy scipy
```

---

## Workflow Overview

The complete analysis proceeds in the following order:

```
Raw BAMs
    │
    ▼
1. BAMs/           — depth/coverage QC, polymorphism calling (read2snp)
    │
    ▼
2. Orthology/      — compute reciprocal best hits OR OrthoFinder HOGs
    │
    ▼
3. MKT_3_species/  — align sequences (divergence + polymorphism),
                     run codeml branch model, merge tables
    │
    ▼
4. branch_specific_MKT/render_analysis.R   ← PRODUCES branch_specific_MKT_results.tsv
    │
    ├──▶ 5. Recombinaison_analysis/
    ├──▶ 6. Expression_analysis/render_expression.R
    ├──▶ 7. protein_protein_interactions/render_ppi.R
    └──▶ 8. 2DNS/render_2dns.R
```

All downstream analyses (steps 5–8) read **`branch_specific_MKT_results.tsv`**
produced by `render_analysis.R`. Run `render_analysis.R` **before** any downstream
analysis.

---

## Quick Start

### Step 1 — Run the branch-specific MKT

```r
# In RStudio or R console, from the MKT_3_species/branch_specific_MKT/ directory:
source("render_analysis.R")

# Or from the terminal:
Rscript MKT_3_species/branch_specific_MKT/render_analysis.R
```

Fill in `SPECIES_LIST` in `render_analysis.R` with species-specific paths before running.
See `MKT_3_species/branch_specific_MKT/README.md` for full instructions.

### Step 2 — Run downstream analyses

Each downstream analysis has its own launcher (`render_*.R`):

```r
# Expression enrichment (requires DESeq2)
source("Expression_analysis/render_expression.R")

# Protein–protein interaction network
source("protein_protein_interactions/render_ppi.R")

# Gene set MK test (KEGG pathways)
source("2DNS/render_2dns.R")
```

### Step 3 — View reports

Each launcher generates an HTML report in its `results/` subdirectory.
Open `results/<analysis>/<analysis>.html` in a browser.

---

## Input Data

The pre-processed input files required by each analysis are included in this repository:

| Analysis | Input files |
|----------|-------------|
| `Expression_analysis/` | `*_count.tsv` — RNA-seq count matrices per experiment |
| `2DNS/` | `input_2DNS_*.csv` — Pathway-level PN/PS/DN/DS per species |
| `protein_protein_interactions/` | STRING edge/node CSVs, candidate ID lists |

Raw sequencing data (FASTQ/BAM files) are deposited at:
> **[NCBI/ENA accession — to be added upon publication]**

---

## Output Structure

Each analysis produces self-contained HTML reports and TSV result tables:

```
results/
├── speltoides/
│   ├── branch_specific_MKT_analysis.html
│   ├── branch_specific_MKT_results.tsv   ← used by all downstream analyses
│   ├── candidates_all.tsv
│   └── candidates_branch_specific.tsv
├── mutica/          (same structure)
├── tauschii/        (same structure)
├── urartu/          (same structure)
├── combined/        (cross-species candidate tables)
├── expression/      (from render_expression.R)
├── ppi/             (from render_ppi.R)
└── 2dns/            (from render_2dns.R)
```

---

## Reproducibility

All analyses are fully automated through the launcher scripts. The only manual
input is specifying file paths in the `SPECIES_LIST` configuration block at the top
of each launcher.

All R analyses set a fixed random seed (`seed = 123` by default) and report the
full session information (`sessionInfo()`) at the end of each HTML report.

The `data/` directory for large raw files is excluded from version control
(see `.gitignore`). Input summary statistics and processed pathway tables
are versioned directly in the repository.

---


## License

This code is released under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## Contact

Manuel Barrientos Tecun
manuel.barrientos@inrae.fr

