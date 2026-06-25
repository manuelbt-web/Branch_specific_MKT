# Step 7 — McDonald-Kreitman Test (Standard and Imputed MKT)

## Objective

Compute the Standard MKT and the Imputed MKT (impMKT) per gene by merging
the divergence output (Step 5), the polymorphism statistics (Step 6b), and the
SDM estimates (Step 6c), then running Fisher's exact tests and BH p-value
adjustment.

**Two scripts:**

| Script | Language | Role |
|---|---|---|
| `merge_tables.py` | Python | Step 7a — merge tables, compute impMKT intermediate values |
| `MKT_analysis.Rmd` | R Markdown | Step 7b — apply filters, compute MKT statistics, produce report |

---

## Dependencies

### Python (merge_tables.py)

| Package | Installation |
|---|---|
| **pandas** | `pip install pandas` |

### R (MKT_analysis.Rmd)

| Package | Installation |
|---|---|
| **tidyverse** | `install.packages("tidyverse")` |
| **kableExtra** | `install.packages("kableExtra")` |
| **patchwork** | `install.packages("patchwork")` |
| **rmarkdown** | `install.packages("rmarkdown")` |

---

## Step 7a — merge_tables.py

### Input

| Input | Description | Source |
|---|---|---|
| `--dndspnps FILE` | Merged dNdSpiNpiS TSV | `all_genes_dNdSpNpS.tsv` (Step 5) |
| `--polymorphism FILE` | Polymorphism + SDM statistics TSV | `step6_stats/polymorphism_stats.tsv` (Step 6b) |

The polymorphism table must come from `polymorphisms_stats.py` (which now includes
SDM estimation). It already contains `PnMinus`, `PnGreater`, `PnNeutral`, etc.

### Output

| File | Description |
|---|---|
| `step7_merged.tsv` | Unified table with all MKT columns |

### What the merge adds

`merge_tables.py` only adds two columns from dNdSpiNpiS:

| Column | Source | Description |
|---|---|---|
| `Dn_counts` | dNdSpiNpiS `fixN` | NS fixed differences (divergence) |
| `Ds_counts` | dNdSpiNpiS `fixS` | S fixed differences (divergence) |

All other columns (polymorphism, SDM, impMKT intermediate values) already come
from `polymorphisms_stats.py`.

### Usage

```bash
python step7_MKT/merge_tables.py \
    --dndspnps     all_genes_dNdSpNpS.tsv \
    --polymorphism step6_stats/polymorphism_stats.tsv \
    --output       step7_merged.tsv \
    --species      Ae_speltoides
```

### All options

| Option | Required | Default | Description |
|---|---|---|---|
| `--dndspnps FILE` | Yes | — | Merged dNdSpiNpiS TSV |
| `--polymorphism FILE` | Yes | — | Polymorphism stats TSV (from Step 6b) |
| `--output FILE` | Yes | — | Output merged TSV |
| `--species NAME` | No | — | Focal species name (informational only) |

---

## Step 7b — MKT_analysis.Rmd

### Input

| Input | Description |
|---|---|
| `params$merged_table` | `step7_merged.tsv` from Step 7a |

### Output

| File | Description |
|---|---|
| `step7_MKT_results.tsv` | Final table with all MKT statistics and p-values |
| `MKT_analysis.html` | Interactive report with plots and summary tables |

### Statistics computed

**Standard MKT** (filter: Pn + Ps ≥ 5; Pn, Ps, Dn, Ds ≥ 1):

| Column | Formula |
|---|---|
| `Standard_MKT` | α = 1 − (Ds × Pn) / (Dn × Ps) |
| `DoS_MKT_standard` | DoS = Dn/(Dn+Ds) − Pn/(Pn+Ps)  (Stoletzki & Eyre-Walker 2011) |
| `standard_MKT_pvalue` | Fisher's exact test on [[Pn, Dn], [Ps, Ds]] |
| `standard_MKT_p_adj` | BH-adjusted p-value (computed on filtered genes only) |

**Imputed MKT** (filter: PnNeutral + Ps ≥ 5; PnNeutral, Ps, Dn, Ds ≥ 1; same filter as standard MKT but applied after subtracting SDMs from Pn):

| Column | Formula |
|---|---|
| `standard_impMKT` | α_imp = 1 − (Ds × PnNeutral) / (Dn × Ps) |
| `DoS_impMKT_standard` | DoS_imp = Dn/(Dn+Ds) − PnNeutral/(PnNeutral+Ps) |
| `standard_impMKT_pvalue` | Fisher's exact test on [[PnNeutral, Dn], [Ps, Ds]] |
| `standard_impMKT_p_adj` | BH-adjusted p-value (computed on filtered genes only) |

**Positively selected genes** are defined as DoS > 0 AND raw p-value < 0.05.

### Figures produced (Sections 7.1 – 7.11)

| Section | Figure | Description |
|---|---|---|
| 7.1 | Alignment stats | number_sites and sample_size distributions |
| 7.2 | Polymorphism stats | π, θ_W, Tajima's D, π_S vs π_NS |
| 7.3 | Count distributions | Pn, Ps, Dn, Ds histograms |
| 7.4 | MK scatter | Pn vs Ps; Dn vs Ds |
| 7.5 | SDM quality | ratioPs, deleterious, PnNeutral vs Pn |
| 7.6 | Filter summary | Gene counts passing each filter |
| 7.7 | DoS distributions | Standard MKT vs impMKT |
| 7.8 | α distributions | Standard MKT vs impMKT |
| 7.9 | p-value distributions | Standard MKT vs impMKT |
| 7.10 | Volcano plots | DoS vs −log₁₀(p), both tests |
| 7.11 | MK plot | Pn/Ps vs Dn/Ds scatter (log scale) |

### Usage

```r
rmarkdown::render(
  "step7_MKT/MKT_analysis.Rmd",
  params = list(
    merged_table = "step7_merged.tsv",
    output_table = "step7_MKT_results.tsv",
    species      = "Ae_speltoides"          # used in plot titles
  ),
  output_file = "MKT_analysis.html"
)
```

Or open `MKT_analysis.Rmd` in RStudio and click **Knit**.

---

## Full two-step example

```bash
# ── Step 7a — merge tables ────────────────────────────────────────────────────
python step7_MKT/merge_tables.py \
    --dndspnps     all_genes_dNdSpNpS.tsv \
    --polymorphism step6_stats/polymorphism_stats.tsv \
    --output       step7_merged.tsv \
    --species      Ae_speltoides
```

```r
# ── Step 7b — MKT analysis (in R) ────────────────────────────────────────────
rmarkdown::render(
  "step7_MKT/MKT_analysis.Rmd",
  params = list(
    merged_table = "step7_merged.tsv",
    output_table = "step7_MKT_results.tsv",
    species      = "Ae_speltoides"
  )
)
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `Missing required columns: ['PnMinus', ...]` | polymorphism table from old script version | Re-run `polymorphisms_stats.py` (new version integrates SDM) |
| Many genes excluded by impMKT filter | Very few high-freq S polymorphisms (small samples) | Inspect `ratioPs` distribution; re-run `polymorphisms_stats.py` with a lower `--freq-cutoff` |
| All `standard_MKT_p_adj` = NA | 0 genes passed filter | Relax filters or check counts in merged table |
| `PnNeutral` > `num_pol_NS` | SDM estimate is negative (more neutral than expected) | Capped automatically; `deleterious` will be negative |

---

## Previous steps

← [Step 5: dNdSpiNpiS divergence](../step5_run_dNdSpNpS/README.md)
← [Step 6: Polymorphism statistics](../step6_polymorphism/README.md)
