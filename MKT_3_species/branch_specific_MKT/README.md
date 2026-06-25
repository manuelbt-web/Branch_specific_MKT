# branch_specific_MKT — Branch-Specific McDonald-Kreitman Analysis

## Objective

Run four McDonald-Kreitman tests per focal species and compare results across
species. The branch-specific tests use dN and dS estimated **on the focal
species branch** by codeml (PAML), converting rates back to divergence counts:

$$D_N^{branch} = \text{round}(dN \times N_{codeml}), \quad D_S^{branch} = \text{round}(dS \times S_{codeml})$$

| Test | Pn used | Dn / Ds source |
|------|---------|----------------|
| Standard MKT | raw Pn | fixN / fixS from dNdSpiNpiS |
| Standard impMKT | PnNeutral (SDM-corrected) | fixN / fixS from dNdSpiNpiS |
| Branch-specific MKT | raw Pn | round(dN × N\_codeml) / round(dS × S\_codeml) |
| Branch-specific impMKT | PnNeutral (SDM-corrected) | round(dN × N\_codeml) / round(dS × S\_codeml) |

---

## Contents

| File | Description |
|------|-------------|
| `render_analysis.R` | **Main launcher** — edit species paths, run this file |
| `branch_specific_MKT_analysis.Rmd` | Per-species analysis (called by the launcher) |
| `merge_mkt_results.py` | Merges polymorphism + dNdSpiNpiS + codeml into one table |

---

## Dependencies

### Python (`merge_mkt_results.py`)

```bash
pip install pandas
```

### R (`branch_specific_MKT_analysis.Rmd`)

```r
install.packages(c("tidyverse", "kableExtra", "patchwork", "UpSetR", "DT", "rmarkdown"))
```

---

## Inputs

### Step 1 — `merge_mkt_results.py`

Three data sources are merged by `gene_id` into one table per focal species:

| Argument | Source | Description |
|----------|--------|-------------|
| `--polymorphism-dir DIR` | `polymorphism_stats.py` output | Per-gene TSV files with Pi, Pn, Ps, PnNeutral, … |
| `--divergence-dir DIR` | `dNdSpiNpiS.sbatch` output | Per-gene `.out` files; `fixN` → `dn_counts`, `fixS` → `ds_counts` |
| `--codeml-file FILE` | `parse_results.py` output | `focal_species_results.tsv` with dN, dS, N\_codeml, S\_codeml, omega, … |
| `-o / --output FILE` | — | Output merged TSV (input to the Rmd) |

```bash
python merge_mkt_results.py \
    --polymorphism-dir  poly_stats/ \
    --divergence-dir    dNdSpiNpiS_output/ \
    --codeml-file       codeml_results/focal_species_results.tsv \
    --focal-species     Aegilopsspeltoides \
    --output            results/speltoides/merged_mkt_results.tsv
```

### Step 2 — `render_analysis.R` + `branch_specific_MKT_analysis.Rmd`

| Input | Source | Description |
|-------|--------|-------------|
| `merged_table` | `merge_mkt_results.py` | One per focal species |
| `orthogroup_table` | `orthogroup_table.py` | HOG → gene ID mapping for all species |

---

## How to run

### Single or multi-species (recommended)

**1.** Open `render_analysis.R`.

**2.** Edit **Section 1** — one entry per focal species:

```r
SPECIES_LIST <- list(

  Aegilopsspeltoides = list(
    merged_table                 = "results/speltoides/merged_mkt_results.tsv",
    output_dir                   = "results/speltoides",
    focal_species_in_ortho_table = "Aegilops_speltoides"
  ),
  Aegilopsmutica = list(
    merged_table                 = "results/mutica/merged_mkt_results.tsv",
    output_dir                   = "results/mutica",
    focal_species_in_ortho_table = "Aegilops_mutica"
  )

)
```

**3.** Edit **Section 2** — shared orthogroup table and combined output directory:

```r
ORTHOGROUP_TABLE <- "data/orthogroup_table.tsv"
COMBINED_DIR     <- "results/combined"
```

**4.** Run:

```r
# RStudio: click Source
source("path/to/render_analysis.R")   # R console
Rscript path/to/render_analysis.R     # terminal
```

The launcher:
- Creates output directories automatically
- Runs the Rmd once per species → per-species HTML report + per-species TSVs
- When ≥2 species: builds combined cross-species tables in `COMBINED_DIR/`

---

### Direct `rmarkdown::render()` (single species, advanced use)

```r
rmarkdown::render(
  input       = "branch_specific_MKT_analysis.Rmd",
  output_file = "results/speltoides/branch_specific_MKT_analysis.html",
  params = list(
    merged_table                 = "results/speltoides/merged_mkt_results.tsv",
    orthogroup_table             = "data/orthogroup_table.tsv",
    output_table                 = "results/speltoides/branch_specific_MKT_results.tsv",
    candidates_table_all         = "results/speltoides/candidates_all.tsv",
    candidates_table_bs          = "results/speltoides/candidates_branch_specific.tsv",
    candidates_shared            = "results/speltoides/candidates_shared.tsv",
    species                      = "Aegilopsspeltoides",
    focal_species_in_ortho_table = "Aegilops_speltoides"
  )
)
```

---

## Outputs

### Per-species (one set per focal species in `output_dir/`)

| File | Description |
|------|-------------|
| `branch_specific_MKT_analysis.html` | Interactive HTML report: all figures and tables |
| `branch_specific_MKT_results.tsv` | All genes: 4 MKT tests + ortholog columns |
| `candidates_all.tsv` | Genes with DoS > 0 and p < 0.05 in ≥1 test |
| `candidates_branch_specific.tsv` | Candidates detected by ≥1 branch-specific test |
| `candidates_shared.tsv` | *(used internally by the Rmd for single-species cross-species option)* |

### Combined (only when ≥2 species — written to `COMBINED_DIR/`)

| File | Description |
|------|-------------|
| `candidates_all_combined.tsv` | All species' positive candidates in one file |
| `candidates_shared_combined.tsv` | HOGs detected in ≥2 focal species (any test) |
| `candidates_branch_specific_combined.tsv` | All species' branch-specific candidates |
| `candidates_shared_branch_specific_combined.tsv` | HOGs positive in ≥2 species (branch-specific tests) |

All combined tables have a `focal_species` column identifying which species each
row belongs to, and an `ortholog` (HOG) column as the shared cross-species key.

---

## Output layout (4 species example)

```
results/
├── speltoides/
│   ├── branch_specific_MKT_analysis.html
│   ├── branch_specific_MKT_results.tsv
│   ├── candidates_all.tsv
│   └── candidates_branch_specific.tsv
├── mutica/
│   └── (same structure)
├── tauschii/
│   └── (same structure)
├── urartu/
│   └── (same structure)
└── combined/
    ├── candidates_all_combined.tsv
    ├── candidates_shared_combined.tsv
    ├── candidates_branch_specific_combined.tsv
    └── candidates_shared_branch_specific_combined.tsv
```

---

## Output table columns (`branch_specific_MKT_results.tsv`)

**Identifiers and orthologs**

| Column | Description |
|--------|-------------|
| `gene_id` | Focal species gene identifier |
| `ortholog` | HOG identifier (e.g. HOG0005558) |
| `{species}_ortholog` | Gene ID in each species (from orthogroup table; one column per species) |

**Polymorphism statistics** (from `polymorphism_stats.py`)

| Column | Description |
|--------|-------------|
| `num_pol_NS` / `num_pol_S` | Non-synonymous / synonymous polymorphic sites (Pn / Ps) |
| `PnNeutral_estimation` | SDM-corrected Pn (float) |
| `deleterious` | Estimated SDM count |
| `Pi_per_site`, `pi_S`, `pi_NS` | Nucleotide diversity |
| `Tajimas_D` | Tajima's D |
| `ratioPs` | Neutral ratio |

**Divergence** (from `dNdSpiNpiS`)

| Column | Description |
|--------|-------------|
| `dn_counts` | Fixed non-synonymous differences (fixN) |
| `ds_counts` | Fixed synonymous differences (fixS) |

**Branch-specific divergence** (computed from codeml rates)

| Column | Formula |
|--------|---------|
| `Dn_branch_specific` | round(dN × N\_codeml) |
| `Ds_branch_specific` | round(dS × S\_codeml) |

**codeml parameters** (from `parse_results.py`)

| Column | Description |
|--------|-------------|
| `dN`, `dS` | Branch dN / dS (focal species) |
| `N_codeml`, `S_codeml` | Expected non-synonymous / synonymous sites (PAML) |
| `omega` | Foreground dN/dS (branch model, focal species) |

**MKT statistics** (one set per test × 4 tests)

| Column | Description |
|--------|-------------|
| `standard_MKT` | α from standard MKT |
| `DoS_MKT_standard` | DoS from standard MKT |
| `standard_MKT_pvalue` | Fisher's exact p-value |
| `standard_MKT_p_adj` | BH FDR-adjusted p-value |
| `pos_standard_MKT` | TRUE if DoS > 0 and p < 0.05 |
| *(same pattern for `standard_impMKT`, `branch_specific_MKT`, `branch_specific_impMKT`)* | |

---

## Figures in the HTML report

| Section | Figure |
|---------|--------|
| 8.1 | **UpSet plot** — intersections of positive candidates across 4 tests |
| 8.2 | MKT filter summary — gene counts at each filter stage |
| 8.3 | Alignment length, sample size, π, θ_W, Tajima's D distributions |
| 8.4 | Standard vs branch-specific divergence count distributions and correlation |
| 8.5 | SDM estimation diagnostics |
| 8.6–8.9 | DoS, α, p-value, and volcano plots (all 4 tests) |
| 8.10 | Classic McDonald-Kreitman plots (all 4 tests) |
| 8.11 | DoS: standard vs branch-specific scatter |
| 8.12 | omega (codeml) vs DoS from branch-specific tests |

---

## Full example — 4 focal species

```bash
# Step 1: merge inputs for each species
for SP in speltoides mutica tauschii urartu; do
  python merge_mkt_results.py \
      --polymorphism-dir  poly_stats_${SP}/ \
      --divergence-dir    dNdSpiNpiS_${SP}/ \
      --codeml-file       codeml_results/focal_species_results.tsv \
      --focal-species     Aegilops${SP} \
      --output            results/${SP}/merged_mkt_results.tsv
done
```

```r
# Step 2: edit SPECIES_LIST in render_analysis.R, then:
source("branch_specific_MKT/render_analysis.R")
# → 4 per-species HTML reports
# → results/combined/  with 4 cross-species tables
```

---

## Positive candidate definition

A gene is a **positive candidate** if, for at least one of the four tests:

> **DoS > 0** AND **raw Fisher p-value < 0.05**

FDR correction (Benjamini-Hochberg) is applied **per test** on genes passing
that test's filter (Pn + Ps ≥ 5 and all cells ≥ 1; PnNeutral + Ps ≥ 5 for
imputed). Adjusted p-values are in the output table for stricter filtering.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `Could not auto-detect focal species` | Species name in merged table doesn't match orthogroup header | Set `focal_species_in_ortho_table = "Aegilops_speltoides"` in `SPECIES_LIST` |
| `merged_table not found` | Wrong path | Check `merged_table` in `SPECIES_LIST`; use absolute path if needed |
| `Dn_branch_specific` all NA | `dN` or `N_codeml` missing in merged table | Verify `--codeml-file` was passed to `merge_mkt_results.py` |
| 0 genes in candidate tables | No genes pass DoS > 0 and p < 0.05 | Inspect volcano plot in HTML report; check divergence count distributions |
| UpSet plot empty | No genes positive in any test | See above |
| Combined tables not created | Fewer than 2 species in `SPECIES_LIST` | Add more species; combined step requires ≥2 |
| `ortholog` column missing in combined | Orthogroup table not joined | Check `ORTHOGROUP_TABLE` path; verify `orthogroup_table.py` was run |

---

← [MKT_3_species — overview](../README.md)
