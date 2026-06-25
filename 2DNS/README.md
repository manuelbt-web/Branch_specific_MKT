# 2DNS — Gene Set MK Analysis

## Objective

Test whether KEGG pathways show an excess of non-synonymous divergence relative
to polymorphism (positive selection signal) using a four-layer evidence framework
applied to four wild *Aegilops*/*Triticum* species.

**Statistical approach:**

| Layer | Method | Input |
|-------|--------|-------|
| 1 | Fisher exact test (two-sided + one-sided) | Pathway-level PN/PS/DN/DS |
| 2 | Benjamini-Hochberg FDR correction | Layer 1 p-values |
| 3 | Stratified permutation (length-matched null) | MKT gene pool |
| 4 | Jackknife leave-one-gene-out | Gene-level data (optional) |

---

## Contents

| File | Description |
|------|-------------|
| `render_2dns.R` | **Main launcher** — fill in species configs and run |
| `2dns_analysis.Rmd` | Analysis template (called by the launcher) |
| `input_2DNS_speltoides.csv` | Pathway-level MK counts — *Ae. speltoides* |
| `input_2DNS_mutica.csv` | Pathway-level MK counts — *Ae. mutica* |
| `input_2DNS_tauschii.csv` | Pathway-level MK counts — *Ae. tauschii* |
| `input_2DNS_urartu.csv` | Pathway-level MK counts — *T. urartu* |

---

## Prerequisites

Run `branch_specific_MKT/render_analysis.R` first for each species. This
analysis reads `branch_specific_MKT_results.tsv` from each `output_dir` to build
the permutation gene pool and identify candidate genes.

---

## Dependencies

```r
install.packages(c("tidyverse", "kableExtra", "patchwork", "rmarkdown"))
```

---

## How to run

### Step 1 — Set paths in `render_2dns.R`

```r
SPECIES_LIST <- list(

  speltoides = list(
    label      = "Ae. speltoides",
    input_file = "input_2DNS_speltoides.csv",  # in this folder
    output_dir = "results/speltoides"           # from render_analysis.R
  ),
  mutica = list(
    label      = "Ae. mutica",
    input_file = "input_2DNS_mutica.csv",
    output_dir = "results/mutica"
  ),
  tauschii = list(
    label      = "Ae. tauschii",
    input_file = "input_2DNS_tauschii.csv",
    output_dir = "results/tauschii"
  ),
  urartu = list(
    label      = "T. urartu",
    input_file = "input_2DNS_urartu.csv",
    output_dir = "results/urartu"
  )

)
```

### Step 2 — Adjust parameters (optional)

```r
MIN_GENES      <- 8L      # minimum genes per pathway
N_PERMUTATIONS <- 1000L   # 200 for quick check, 1000+ for publication
CANDIDATE_TYPE <- "bs_any"
```

### Step 3 — Run

```r
source("render_2dns.R")   # RStudio / R console
Rscript render_2dns.R     # terminal
```

---

## Input file format

### `input_2DNS_*.csv` (required, in this folder)

One row per KEGG pathway. Required columns:

| Column | Description |
|--------|-------------|
| `pathway_name` | KEGG pathway name |
| `PN` | Non-synonymous polymorphism count |
| `PS` | Synonymous polymorphism count |
| `DN` | Non-synonymous divergence count |
| `DS` | Synonymous divergence count |
| `n_genes` | Number of genes in the pathway |
| `set_length` | Total CDS length (optional; used for display) |

### `branch_specific_MKT_results.tsv` (from `render_analysis.R`)

Used to:
- Build the permutation gene pool (all MKT-analyzable genes with their PN/PS/DN/DS)
- Identify BS candidate genes (`pos_branch_specific_*` columns)

### `gene_level_file` (optional, for Layer 4 jackknife)

Not included by default. If provided, must contain one row per gene per pathway:

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier |
| `pathway_name` | KEGG pathway name |
| `PN`, `PS`, `DN`, `DS` | Per-gene counts |

Add to `SPECIES_LIST` as `gene_level_file = "genes_stats_speltoides_pathways.csv"`.

---

## Outputs

All written to `results/2dns/`:

| File | Content |
|------|---------|
| `2dns_analysis.html` | Interactive HTML report |
| `fisher_results.tsv` | Layer 1 Fisher results for all pathways |
| `permutation_results.tsv` | Layer 3 permutation results (Fisher-significant pathways) |
| `comprehensive_summary.tsv` | All layers combined per pathway |
| `candidates.tsv` | BS candidate genes identified from MKT results |

---

## Candidate type options

| Value | Meaning |
|-------|---------|
| `"bs_any"` | Union of branch-specific standard + imputed MKT (recommended) |
| `"bs_impMKT"` | Branch-specific imputed MKT only |
| `"bs_mkt"` | Branch-specific standard MKT only |
| `"any"` | Union of all four MKT tests |

---

## Permutation null model

The permutation null tests whether the observed pathway-level MK signal is
stronger than expected when sampling the same number of genes from the background
pool of all MKT-analyzable genes. Sampling is stratified by gene length (proxy:
total coding sites from `branch_specific_MKT_results.tsv`) using
N quantile-based strata (default: 5; set via `LENGTH_STRATA_N` in the launcher).

Because gene-to-pathway membership is not encoded in
`branch_specific_MKT_results.tsv`, the stratification uses the overall pool
length distribution (not pathway-specific). This is a valid conservative null
model. A gene-level pathway CSV (`gene_level_file`) would enable exact
pathway-specific stratification and the jackknife layer.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `File not found: branch_specific_MKT_results.tsv` | Pipeline not run | Run `render_analysis.R` first |
| Missing columns in input CSV | Wrong file format | Check columns: `pathway_name`, `PN`, `PS`, `DN`, `DS`, `n_genes` |
| Layer 4 skipped | No `gene_level_file` | Set `gene_level_file` in `SPECIES_LIST` |
| 0 candidates | Filter too strict | Check that MKT results contain `pos_branch_specific_*` columns |

---

← [Scripts — overview](../README.md)
← [branch_specific_MKT — run this first](../MKT_3_species/branch_specific_MKT/README.md)
