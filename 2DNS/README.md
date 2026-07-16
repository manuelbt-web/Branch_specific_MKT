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
| 3 | Stratified permutation (length- and pathway-overlap-matched null) | MKT gene pool, or `gene_level_file` if provided (recommended) |
| 4 | Jackknife leave-one-gene-out | `gene_level_file` (recommended — see below) |
| — | Candidate gene × significant-pathway enrichment test | `gene_level_file` + branch-specific candidates |

**Providing a `gene_level_file` per species is strongly recommended** — see
[Input file format](#input-file-format) below. Without it, Layer 3 falls back
to a conservative length-only null (no pathway-overlap control) and Layers 4
and the enrichment test are skipped entirely.

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
    label           = "Ae. speltoides",
    input_file      = "input_2DNS_speltoides.csv",  # in this folder
    output_dir      = "results/speltoides",          # from render_analysis.R
    gene_level_file = "genes_stats_speltoides_pathways.csv"  # recommended, see below
  ),
  mutica = list(
    label           = "Ae. mutica",
    input_file      = "input_2DNS_mutica.csv",
    output_dir      = "results/mutica",
    gene_level_file = "genes_stats_mutica_pathways.csv"
  ),
  tauschii = list(
    label           = "Ae. tauschii",
    input_file      = "input_2DNS_tauschii.csv",
    output_dir      = "results/tauschii",
    gene_level_file = "genes_stats_tauschii_pathways.csv"
  ),
  urartu = list(
    label           = "T. urartu",
    input_file      = "input_2DNS_urartu.csv",
    output_dir      = "results/urartu",
    gene_level_file = "genes_stats_urartu_pathways.csv"
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

### `gene_level_file` (optional, but recommended for every species)

Not included by default — this repo ships one per species anyway
(`genes_stats_speltoides_pathways.csv`, etc.) that already has the right
format; point `gene_level_file` at it. One row per gene per pathway a gene
belongs to (genes in >1 pathway repeat over multiple rows):

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier |
| `pathway_name` | KEGG pathway name |
| `PN`, `PS`, `DN`, `DS` | Per-gene counts |
| `length` | Coding length (used for pathway-specific length stratification in Layer 3) |

Add to `SPECIES_LIST` as `gene_level_file = "genes_stats_speltoides_pathways.csv"`.

Supplying this file changes three things, not just the jackknife:
1. **Layer 3 permutation** draws its length-stratum targets from the
   *observed pathway's own* gene-length composition (not the overall pool),
   and resamples genes weighted by how many pathways each one belongs to
   ("pathway degree") — so a pathway sharing many genes with other pathways
   isn't treated as if all its genes were equally exchangeable with any
   random gene. This is what "controlling for pathway size, gene-length
   composition, and pathway overlap structure" means in practice, and it is
   required to reproduce the published pathway counts (without it, Layer 3
   over-estimates how many pathways survive permutation).
2. **Layer 4 jackknife** becomes available (as before).
3. **The candidate × pathway enrichment test** becomes available: whether
   branch-specific candidate genes (`pos_branch_specific_*` in
   `branch_specific_MKT_results.tsv`) are over-represented among the genes
   that make up pathways significant at each layer, restricted to genes with
   a KEGG annotation (a candidate gene with no pathway annotation can never
   appear in a pathway, and is correctly excluded from the test rather than
   silently counted as "not enriched").

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
| `candidate_pathway_enrichment.tsv` | Enrichment test: candidate genes vs. significant-pathway membership, per layer (written only if `gene_level_file` was supplied for >= 1 species) |
| `candidate_pathway_presence.tsv` | Per-candidate-gene detail: which pathway(s) it belongs to and whether those pass each layer |

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
stronger than expected when sampling the same number of genes from a
background pool, stratified by gene length using N quantile-based strata
(default: 5; set via `LENGTH_STRATA_N` in the launcher).

**With a `gene_level_file` (recommended, used whenever supplied):** the
background pool and length strata are built from that file. The permutation
target for each stratum comes from the *observed pathway's own* gene-length
composition, and genes are resampled with probability proportional to how
many pathways they belong to (their "degree" in the gene↔pathway bipartite
graph) — controlling for pathway size, gene-length composition, **and**
pathway overlap structure simultaneously. This is the method that reproduces
the published pathway counts.

**Without a `gene_level_file`:** gene-to-pathway membership isn't available
from `branch_specific_MKT_results.tsv` alone, so the pool is built from that
file instead, using the *overall* pool length distribution (not
pathway-specific) with no degree-weighting. This is a valid but more
conservative null model — expect it to retain more pathways as significant
than the method above, since it doesn't discount hub genes shared across
many pathways. Layer 4 (jackknife) and the candidate/pathway enrichment test
also require `gene_level_file` and are skipped without it.

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
