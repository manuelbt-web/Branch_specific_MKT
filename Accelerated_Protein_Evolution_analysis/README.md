# Accelerated_Protein_Evolution_analysis

## Objective

Test whether genes under positive selection (branch-specific MKT candidates)
show convergent signals of protein sequence acceleration: elevated branch-
specific dN/dS and/or significant branch model likelihood ratio test (LRT).

The key question: do genes detected by branch-specific McDonald-Kreitman tests
also show elevated rates of non-synonymous substitution relative to background,
suggesting genuine adaptive evolution rather than relaxed constraint?

---

## Contents

| File | Description |
|------|-------------|
| `render_high_dn_ds.R` | **Main launcher** — edit species dirs here and run |
| `high_dn_ds.Rmd` | Analysis script (called by the launcher) |

---

## Prerequisites

This analysis reads outputs from `branch_specific_MKT/render_analysis.R`.
Run that step first. **No other inputs are needed.**

---

## Dependencies

```r
install.packages(c("tidyverse", "ComplexUpset", "kableExtra", "patchwork", "rmarkdown"))
```

---

## How to run

**1.** Open `render_high_dn_ds.R`.

**2.** Copy the `output_dir` values from the `render_analysis.R` you already
configured — they are the only paths needed here:

```r
SPECIES_LIST <- list(

  Aegilopsspeltoides = list(output_dir = "results/speltoides"),
  Aegilopsmutica     = list(output_dir = "results/mutica"),
  Aegilopstauschii   = list(output_dir = "results/tauschii"),
  Triticumurartu     = list(output_dir = "results/urartu")

)
```

**3.** Adjust analysis parameters if needed:

```r
Z_THRESHOLD     <- 1.5    # top ~7% of the dN/dS distribution
LRT_P_THRESHOLD <- 0.05   # branch model LRT significance threshold
OUTPUT_DIR      <- "results/high_dn_ds"
```

**4.** Run:

```r
source("render_high_dn_ds.R")   # R console / RStudio Source button
Rscript render_high_dn_ds.R     # terminal
```

---

## Files read (automatically from each `output_dir`)

| File | Content |
|------|---------|
| `branch_specific_MKT_results.tsv` | All genes with MKT statistics, codeml results, dN/dS |
| `candidates_all.tsv` | Positive candidates (DoS > 0, p < 0.05 in ≥1 test) |
| `candidates_branch_specific.tsv` | Candidates detected by ≥1 branch-specific test |

---

## Outputs

All written to `OUTPUT_DIR/` (default: `results/high_dn_ds/`):

| File | Content |
|------|---------|
| `high_dn_ds.html` | Interactive HTML report |
| `enrichment_high_dNdS.tsv` | Per-species Fisher/Wilcoxon enrichment in high dN/dS |
| `enrichment_lrt_acceleration.tsv` | Per-species Fisher enrichment in accelerated genes |
| `lrt_classifications.tsv` | Per-gene LRT results (all species combined) |
| `ortholog_pairs_{A}_vs_{B}.tsv` | Cross-species ortholog pair tables (≥2 species only) |

---

## Analyses

### 1 — High dN/dS identification

For each species:
- **Universe**: genes analyzable in ≥1 branch-specific MKT test (Pn+Ps ≥ 5, all counts ≥ 1)
  with a valid `DnDs_branch_specific_counts` value
- **z-score**: computed on log₁₊ₓ(dN/dS\_branch\_specific) across the universe
- **High dN/dS**: genes with z > `Z_THRESHOLD`
- **Fisher enrichment** (one-sided greater): branch-specific MKT candidates vs background
- **Wilcoxon test**: distribution of dN/dS in candidates vs background

### 2 — Branch model LRT

For each species, using codeml values already in the stats table:

$$LRT = 2 \times (lnL_{branch} - lnL_{M0}) \sim \chi^2_{df=1}$$

Genes classified as **Accelerated** (LRT p < threshold AND omega\_branch > omega\_M0),
**Decelerated** (significant AND omega\_branch < omega\_M0), or **Neutral**.

Fisher test: are branch-specific MKT candidates enriched in the accelerated category?

### 3 — Triple overlap (per species)

For each species: genes in high dN/dS ∩ accelerated (LRT) ∩ branch-specific MKT candidates.

### 4 — Cross-species patterns (≥2 species)

Ortholog pairs are identified from the `{species}_ortholog` columns in the
stats table (populated by `orthogroup_table.py`).

For each pair of focal species:
- **Lineage-specific high dN/dS**: z\_A > threshold AND z\_B < −threshold (or vice versa)
- **Shared high dN/dS**: both species high
- **UpSet plots** (using `ComplexUpset`):
  - Methods × species: Standard impMKT and Branch-specific impMKT for each species
  - Methods + high dN/dS × species: adds high dN/dS categories
- **Fisher enrichment**: BS candidates enriched in lineage-specific patterns?

---

## Column mapping (old → new pipeline)

| Old column / concept | New column / derivation |
|----------------------|------------------------|
| `dNdS_branch_specific` | `DnDs_branch_specific_counts` |
| `branch_specific_MKT == "ANALYZABLE"` | Derived from MKT filter on `Dn_branch_specific`, `Ds_branch_specific`, `num_pol_NS`, `num_pol_S` |
| `divergence_table_*.txt` | `lnL_branch_model`, `lnL_m0_model`, `omega`, `omega_m0` columns in stats table |
| Method `"impMKT 3sp"` | `pos_standard_impMKT == TRUE` |
| Method `"BS impMKT"` | `pos_branch_specific_impMKT == TRUE` |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `Files not found: branch_specific_MKT_results.tsv` | render_analysis.R not yet run | Run `render_analysis.R` first; check `output_dir` paths |
| `No ortholog column for '...' found` | Orthogroup table wasn't joined | Check that `ORTHOGROUP_TABLE` was set in render_analysis.R |
| `DnDs_branch_specific_counts` all NA | codeml results missing from merged table | Verify `--codeml-file` in merge_mkt_results.py |
| Empty UpSet plots | No candidates detected | Check results from branch_specific_MKT_analysis.Rmd |
| Cross-species section skipped | Only 1 species in SPECIES_LIST | Add more species — cross-species requires ≥2 |

---

← [Scripts — overview](../README.md)
← [branch_specific_MKT — run this first](../MKT_3_species/branch_specific_MKT/README.md)
