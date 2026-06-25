# Expression_analysis

## Objective

Test whether genes under positive selection in wild *Aegilops*/*Triticum* species
(identified by branch-specific MKT) are enriched among differentially expressed
genes in *T. aestivum* under biotic and abiotic stress conditions.

**Statistical approach:** Fisher's exact test (two-sided) comparing the proportion
of DE genes (padj < 0.05, |LFC| > 0.5) between positively selected candidates and
the neutral background, for each stress condition independently.

**Orthologue mapping:** focal-species genes are mapped to *T. aestivum* via the
`aestivum_ortholog` column in `branch_specific_MKT_results.tsv`.

---

## Contents

| File | Description |
|------|-------------|
| `render_expression.R` | **Main launcher** — fill in species configs and run |
| `expression_analysis.Rmd` | Analysis template (called by the launcher) |
| `SRP064598_count.tsv` | Cold stress RNA-seq counts |
| `SRP045409_count.tsv` | Heat / drought stress RNA-seq counts |
| `ERP013829_count.tsv` | *Fusarium graminearum* inoculation counts |
| `SRP041017_count.tsv` | Powdery mildew + stripe rust counts |
| `ERP009837_count.tsv` | *Zymoseptoria tritici* (Septoria) counts |
| `DRP000768_count.tsv` | Phosphorus starvation counts |
| `PAMP_Triggered_Imune_Response_count.tsv` | PAMP-triggered immune response counts |

---

## Prerequisites

Run `branch_specific_MKT/render_analysis.R` first. This analysis reads
`branch_specific_MKT_results.tsv` from each species output directory.

---

## Dependencies

```r
install.packages(c("tidyverse", "patchwork", "kableExtra", "rmarkdown"))

# DESeq2 (Bioconductor)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
```

---

## How to run

### Step 1 — Fill in `SPECIES_LIST` in `render_expression.R`

```r
SPECIES_LIST <- list(

  Aegilopsspeltoides = list(
    output_dir = "results/speltoides"   # from branch_specific_MKT/render_analysis.R
  ),
  Aegilopsmutica = list(
    output_dir = "results/mutica"
  ),
  Aegilopstauschii = list(
    output_dir = "results/tauschii"
  ),
  Triticumurartu = list(
    output_dir = "results/urartu"
  )

)
```

### Step 2 — Choose the candidate type (optional)

```r
CANDIDATE_TYPE <- "bs_any"   # union of branch-specific standard + imputed MKT
```

| Value | Meaning |
|-------|---------|
| `"bs_any"` | Union of `pos_branch_specific_MKT` and `pos_branch_specific_impMKT` (recommended) |
| `"bs_impMKT"` | Branch-specific imputed MKT only |
| `"bs_mkt"` | Branch-specific standard MKT only |
| `"any"` | Union of all 4 MKT tests |

### Step 3 — Run

```r
source("render_expression.R")   # RStudio / R console
# or
Rscript render_expression.R     # terminal
```

To run only specific experiments, set `enabled = FALSE` in the other entries
of `EXPERIMENTS` in `render_expression.R`.

---

## Outputs

All written to `results/expression/`:

| File | Content |
|------|---------|
| `<experiment>.html` | HTML report per experiment |
| `<experiment>/<experiment>_<species>_fisher_tests.tsv` | Per-species Fisher test table |
| `<experiment>/<experiment>_cross_species_summary.tsv` | Combined species comparison |

---

## Experiments

| ID | Name | Count file | Control | Conditions |
|----|------|-----------|---------|------------|
| cold | Cold stress | `SRP064598_count.tsv` | none | 10-day cold, tissue culture + cold |
| heat | Heat/drought | `SRP045409_count.tsv` | control | 1h/6h heat, drought, combined |
| fusarium | *Fusarium* | `ERP013829_count.tsv` | mock_* | 6 time points (3h–48h) vs matched mock |
| powdery | Powdery + rust | `SRP041017_count.tsv` | none | Powdery mildew E09, stripe rust CYR31 |
| septoria | *Septoria* | `ERP009837_count.tsv` | mock_* | 5 time points (1d–21d) vs matched mock |
| phosphorus | P starvation | `DRP000768_count.tsv` | none | 10-day P starvation |
| pamp | PAMP response | `PAMP_Triggered_Imune_Response_count.tsv` | water | chitin (1 g/L), flg22 (500 nM) |

---

## Experimental design notes

### Fusarium & Septoria (time-matched design)

These experiments compare pathogen-inoculated vs mock-inoculated samples at the same
time point. The `contrasts` list in `render_expression.R` maps each treatment condition
to its corresponding mock reference:

```r
contrasts = list(
  "Fusarium_3 hours" = "mock_3 hours",
  "Fusarium_6 hours" = "mock_6 hours",
  ...
)
```

Conditions are named as `Treatment_time` (e.g. `"Fusarium_6 hours"`, `"mock_6 hours"`).

### Other experiments

Simple one-vs-control design. All non-control conditions are automatically tested
against `control_condition` using DESeq2's `contrast = c("condition", trt, ref)`.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `File not found: branch_specific_MKT_results.tsv` | Pipeline not run | Run `branch_specific_MKT/render_analysis.R` first |
| `sample_conditions length != count matrix columns` | Wrong vector length | Count columns in count TSV with `ncol(read_tsv(count_file)) - 1` |
| 0 candidates match expression data | Wrong ortholog mapping | Check `aestivum_ortholog` column in MKT results is populated |
| DESeq2 errors on contrasts | Typo in condition name | Check that all names in `contrasts` appear in `sample_conditions` |
| `No significant DE genes` | Normal result | Fewer candidates than neutral background → no enrichment |

---

← [Scripts — overview](../README.md)
← [branch_specific_MKT — run this first](../MKT_3_species/branch_specific_MKT/README.md)
