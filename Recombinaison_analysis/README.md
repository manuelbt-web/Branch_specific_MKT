# Recombinaison_analysis

## Objective

Test whether branch-specific positive selection candidates are enriched in
genomic regions with high recombination rates. Under Hill-Robertson
interference, linked selection reduces the efficacy of natural selection in
low-recombination regions, so we expect an excess of adaptive evolution in
high-recombination regions.

**Statistical approach:** binomial logistic GLM with binary selection status as
the response variable. Assignment uncertainty (gene-to-recombination-interval
distance) enters the model as a standardised covariate rather than being
handled with separate overlap-only sub-analyses.

---

## Contents

| File | Description |
|------|-------------|
| `render_recombinaison.R` | **Main launcher** — fill in species configs and run |
| `recombinaison.Rmd` | Analysis script (called by the launcher) |
| `prepare_gene_coords.py` | Convert GFF3 annotation → gene coordinates BED |
| `prepare_recmap.py` | Convert Glemin et al. recombination maps → BED |

---

## Prerequisites

1. **Run `branch_specific_MKT/render_analysis.R` first.** This analysis reads
   `branch_specific_MKT_results.tsv` and `candidates_branch_specific.tsv`
   produced there.
2. **Prepare gene coordinates BED** with `prepare_gene_coords.py` (one per species).
3. **Convert recombination maps** with `prepare_recmap.py` (one per species).

---

## Dependencies

### Python (data preparation scripts)

Standard library only — no `pip install` needed.

```bash
python --version   # requires Python ≥ 3.8
```

### R (analysis)

```r
install.packages(c("tidyverse", "kableExtra", "patchwork", "rmarkdown"))

# Bioconductor packages (genome interval operations)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "IRanges", "S4Vectors"))
```

---

## Step-by-step workflow

### Step 1 — Prepare gene coordinates BED

Run `prepare_gene_coords.py` once per focal species using the GFF3 annotation
file for that species.

**1a. Preview the GFF3 file** to check feature types and column names:

```bash
python prepare_gene_coords.py \
    --gff  annotation/Ae_speltoides.gff3 \
    --preview
```

This prints the available feature types (gene, mRNA, CDS, …) and the first
few lines, so you can identify the correct `--feature` and `--id-attr`.

**1b. Extract gene coordinates:**

```bash
python prepare_gene_coords.py \
    --gff      annotation/Ae_speltoides.gff3 \
    --feature  gene \
    --id-attr  ID \
    --out      data/speltoides_gene_coords.bed
```

**1c. Verify that gene IDs match the MKT pipeline:**

```bash
# Check a few gene IDs in the BED file
head -5 data/speltoides_gene_coords.bed

# Check gene IDs in the MKT results
cut -f1 results/speltoides/branch_specific_MKT_results.tsv | head -5
```

They must match exactly (e.g. both `EVM0000001.1` or both `TRAES_001`).

#### Common GFF3 ID issues

| Annotation | ID attribute value | Fix |
|------------|-------------------|-----|
| Ensembl plants | `gene:TRAES001` | Add `--strip-prefix` |
| MAKER/EVM with versions | `EVM.model.Chr1A.1.10` | Add `--strip-version` |
| AUGUSTUS / BRAKER | `gene1` | Use `--id-attr gene_id` or `--id-attr Name` |
| Chromosome named `chr1` vs `Chr1` | `chr1` | Use `--chr-prefix Chr --chr-strip chr` |

---

### Step 2 — Convert recombination maps to BED

**Source:** [Glemin et al. ms-rec-triticeae](https://github.com/sylvainglemin/ms-rec-triticeae/tree/main/outputs/recombination)

Download the relevant `.txt` file for each species (or the remapped version if
your genome differs from the reference). Then convert with `prepare_recmap.py`.

**2a. Preview the file** to check column names:

```bash
python prepare_recmap.py \
    --input  Ae_speltoides_genomeB_1cM_remapped.txt \
    --preview
```

Expected output for Glemin et al. files:

```
Columns (8):
  [0] poscM
  [1] Chromosome   ← default --chr-col
  [2] Start        ← default --start-col
  [3] End          ← default --end-col
  [4] recRate      ← default --rec-col
  [5] nb_complete_site
  [6] piSyn
  [7] f0
```

**2b. Convert to BED** — defaults match the Glemin et al. format, so no flags are needed:

```bash
python prepare_recmap.py \
    --input  Ae_speltoides_genomeB_1cM_remapped.txt \
    --out    data/Ae_speltoides_recmap.bed
```

NA values in `recRate` are automatically skipped. Float coordinates (`Start`/`End`)
are rounded to integers. The resulting chromosome names will match whatever
is in the `Chromosome` column (e.g. `1B`, `2A`, etc.).

**2c. Verify chromosome name consistency** (most common issue):

```bash
cut -f1 data/Ae_speltoides_recmap.bed | sort -u
cut -f1 data/speltoides_gene_coords.bed | sort -u
```

Both outputs must show the same chromosome names (e.g. both `Chr1A` or both
`1A`). Fix mismatches with `--chr-prefix` or `--chr-strip`:

```bash
# If recmap uses "1A" but gene BED uses "Chr1A":
python prepare_recmap.py \
    --input   Ae_speltoides_genomeB_1cM_remapped.txt \
    --rec-col mean_rho \
    --chr-prefix Chr \
    --out     data/Ae_speltoides_recmap.bed
```

#### Unit note

The `recRate` column in the Glemin et al. files is in **cM/bp** (very small
numbers, e.g. `7.69e-07`). The GLM uses `rec_rate × rec_scale` where
`rec_scale = 1e6` (default in `render_recombinaison.R`), which converts
cM/bp → cM/Mb at analysis time. **Do not scale the BED file** — leave values
as-is and keep `REC_SCALE = 1e6` in the launcher. The logistic regression
coefficient is then interpretable as the change in log-odds per 1 cM/Mb
increase in recombination rate.

---

### Step 3 — Run the analysis

**3a. Open `render_recombinaison.R`** and fill in `SPECIES_LIST`:

```r
SPECIES_LIST <- list(

  Aegilopsspeltoides = list(
    output_dir      = "results/speltoides",          # from render_analysis.R
    gene_coords_bed = "data/speltoides_gene_coords.bed",
    rec_map_bed     = "data/Ae_speltoides_recmap.bed"
  ),
  Aegilopsmutica = list(
    output_dir      = "results/mutica",
    gene_coords_bed = "data/mutica_gene_coords.bed",
    rec_map_bed     = "data/Ae_mutica_recmap.bed"
  )

)
```

**3b. Run:**

```r
source("render_recombinaison.R")   # R console / RStudio Source button
Rscript render_recombinaison.R     # terminal
```

---

## Full example — 2 species

```bash
# ── Step 0: get recombination maps from Glemin et al. ──────────────────────
# Download manually from:
# https://github.com/sylvainglemin/ms-rec-triticeae/tree/main/outputs/recombination

# ── Step 1: gene coordinates ───────────────────────────────────────────────
# Preview first
python Recombinaison_analysis/prepare_gene_coords.py \
    --gff annotation/Ae_speltoides.gff3 --preview

# Extract
python Recombinaison_analysis/prepare_gene_coords.py \
    --gff annotation/Ae_speltoides.gff3 \
    --out data/speltoides_gene_coords.bed

python Recombinaison_analysis/prepare_gene_coords.py \
    --gff annotation/Ae_mutica.gff3 \
    --out data/mutica_gene_coords.bed

# ── Step 2: recombination maps ─────────────────────────────────────────────
# Preview first (check column names match Glemin et al. format)
python Recombinaison_analysis/prepare_recmap.py \
    --input Ae_speltoides_genomeB_1cM_remapped.txt --preview

# Convert (no extra flags needed — defaults match Glemin et al. format)
python Recombinaison_analysis/prepare_recmap.py \
    --input  Ae_speltoides_genomeB_1cM_remapped.txt \
    --out    data/Ae_speltoides_recmap.bed

python Recombinaison_analysis/prepare_recmap.py \
    --input  Ae_mutica_1cM_remapped.txt \
    --out    data/Ae_mutica_recmap.bed

# ── Step 3: verify chromosome name match ───────────────────────────────────
diff <(cut -f1 data/speltoides_gene_coords.bed | sort -u) \
     <(cut -f1 data/Ae_speltoides_recmap.bed | sort -u)
# empty output = names match

# ── Step 4: run GLM analysis ───────────────────────────────────────────────
# Edit SPECIES_LIST in render_recombinaison.R, then:
Rscript Recombinaison_analysis/render_recombinaison.R
```

---

## What is read automatically from `branch_specific_MKT`

| File | Used for |
|------|---------|
| `{output_dir}/branch_specific_MKT_results.tsv` | Defines analyzable gene universe (MKT filter) |
| `{output_dir}/candidates_branch_specific.tsv` | Binary response variable: selected = 1 |

The `selection` column and the analyzable filter are both derived
**automatically** — you do not need to add them to the gene BED file.

---

## Outputs

All written to `OUTPUT_DIR/` (default: `results/recombination/`):

| File | Content |
|------|---------|
| `recombination.html` | Interactive HTML report |
| `gene_recombination_assignments.tsv` | Per-gene: rec rate, assignment method, distance_z, selection status |
| `glm_coefficients.tsv` | GLM coefficient table (all species combined) |
| `sensitivity_analysis.tsv` | Full vs distance-filtered model comparison |

---

## Analyses in the HTML report

| Section | Content |
|---------|---------|
| 1 — Input loading | Summary table: analyzable genes, candidates, rec intervals per species |
| 2 — Recombination assignment | Per-gene assignment method (overlap / nearest / population mean) and quality |
| 3 — Rec rate vs selection | Density plots: candidates vs background; median comparison table |
| 4 — Binomial logistic GLM | Main model + interaction model + LRT; OR plot with 95% CI |
| 5 — Sensitivity analysis | Exclude |distance_z| > threshold; compare full vs filtered estimates |
| 6 — Genomic distribution | Per-chromosome scatter: rec rate vs position, candidates highlighted |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `File not found: branch_specific_MKT_results.tsv` | render_analysis.R not run | Run `render_analysis.R` first |
| `No genes remain after restricting to chromosomes with rec data` | Chr name mismatch | Run the `diff` check above; use `--chr-prefix`/`--chr-strip` |
| All genes assigned by `population_mean` | No overlap/nearest match | Same as above |
| 0 selected genes | Wrong `CANDIDATE_TYPE` or empty candidates table | Check `candidates_branch_specific.tsv` exists and is non-empty |
| `mean_rho column not found` | Different column name in your recmap file | Run `prepare_recmap.py --preview`; use correct `--rec-col` |
| GFF extracts 0 genes | Wrong `--feature` or `--id-attr` | Run `prepare_gene_coords.py --preview`; adjust flags |
| Gene IDs don't match MKT results | GFF uses different ID format | Use `--strip-prefix`, `--strip-version`, or adjust `--id-attr` |
| Bioconductor install fails | BiocManager not available | `install.packages("BiocManager")` first |

---

← [Scripts — overview](../README.md)
← [branch_specific_MKT — run this first](../MKT_3_species/branch_specific_MKT/README.md)
