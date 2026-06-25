# Coverage Analysis Pipeline

This folder contains a set of scripts for computing and summarizing sequencing
coverage from BAM files.

## Running the Pipeline

### Option A — Master wrapper (recommended)

Use `run_coverage_pipeline.sh` to run all steps or any subset in one command:

```bash
# Full pipeline (steps 1–7):
bash run_coverage_pipeline.sh --bam-dir /path/to/bams --output /path/to/results

# Single step (e.g. redo step 2 with custom thresholds):
bash run_coverage_pipeline.sh --steps 2 --output /path/to/results \
     --min-depth 15 --min-ind 8

# See all options:
bash run_coverage_pipeline.sh --help
```

### Option B — Run individual scripts

Each script can also be run independently (see the step-by-step instructions
below). This is useful for debugging or re-running a specific step.

---

The scripts are meant to be run in the order described below.
Each step produces outputs that are required by the next step.

---

## Pipeline Overview

```
BAM files (one per individual)
    │
    ▼  Step 1 — depth_per_individual/depth_per_individual.sbatch      [SLURM]
depth_per_position_all_individuals_BAMs.depth
    │
    ▼  Step 2 — covered_positions/covered_positions.sh                [bash]
    ├── depth_per_position_covered_all_individuals.depth   (filtered positions)
    ├── contig_percentage_covered_positions.txt            (per-contig stats)
    └── number_covered_position.txt                        (genome-wide summary)
         │                     │
         │                     ▼  Step 3a — coverage_summary_contigs/percentage_coverage.sh
         │                     contig_coverage_summary.txt
         │
         ▼  Step 3b — mean_depth/calculate_mean_depth.sh             [bash]
         ├── mean_depth_per_position_covered.txt
         └── mean_depth_per_position.txt
              │
              ├── Step 4 — mean_depth_per_contig/mean_depth_per_contig.sh
              │   mean_depth_per_contig.txt
              │
              └── Step 5 — number_position_per_depth/number_positions_per_depth_per_contig.sh
                  number_position_per_depth_per_contig_sorted_with_total_positions.txt
                      │
                      ▼  Step 6 — table_contigs_per_depth/table_contigs_per_depth_coverage.sh
                      ├── coverage_depth_summary_counts.txt
                      ├── coverage_depth_summary_percentage.txt
                      └── coverage_depth_min_positions_counts.txt

BAM files (optional, independent)
    │
    ▼  Step 7 — whole_coverage_estimation/coverage_estimation_samtool.sbatch  [SLURM]
    coverage_results.tsv  (Sample | Reads_mapped | Mean_Read_Length | Coverage)

Coverage files from two species + RBH pairs file (optional, downstream MKT input)
    │
    ▼  Step 8 — count_paired_covered_contig_from_RBH/count_covered_pairs_contigs.sh
    (count of best-hit pairs where both contigs meet the coverage threshold)

Visualization (Step 10)
    └── Chromosome_coverage.Rmd  → step10_report/coverage_report.html
```

---

## Step-by-Step Instructions

### Step 1 — Generate depth matrix across individuals

**Script:** `depth_per_individual/depth_per_individual.sbatch`
**Cluster:** SLURM — edit the `[CLUSTER]` lines before submission

```bash
sbatch depth_per_individual/depth_per_individual.sbatch \
    /path/to/bam_files/ \
    /path/to/output/depth/
```

**Input:** directory of `.bam` files (one per individual, indexed)
**Output:** `depth_per_position_all_individuals_BAMs.depth`
— tab-delimited matrix `contig | pos | sample1 | sample2 | ...`

---

### Step 2 — Filter positions by coverage threshold

**Script:** `covered_positions/covered_positions.sh`

```bash
bash covered_positions/covered_positions.sh \
    /path/to/depth_per_position_all_individuals_BAMs.depth \
    /path/to/output/covered/ \
    [min_depth]        # default: 10 reads
    [min_individuals]  # default: 6 individuals
```

**Input:** depth matrix from Step 1
**Output:**
- `depth_per_position_covered_all_individuals.depth` — positions where ≥ N individuals have ≥ D reads
- `contig_percentage_covered_positions.txt` — per-contig coverage statistics
- `number_covered_position.txt` — genome-wide summary

---

### Step 3a — Summarize contigs by coverage bins

**Script:** `coverage_summary_contigs/percentage_coverage.sh`

```bash
bash coverage_summary_contigs/percentage_coverage.sh \
    /path/to/contig_percentage_covered_positions.txt \
    /path/to/output/summary/
```

**Input:** `contig_percentage_covered_positions.txt` from Step 2
**Output:** `contig_coverage_summary.txt` — count of contigs in each 10% coverage bin (0%, 10%, ..., 100%)

---

### Step 3b — Calculate mean depth per position

**Script:** `mean_depth/calculate_mean_depth.sh`

```bash
bash mean_depth/calculate_mean_depth.sh \
    /path/to/depth_per_position_covered_all_individuals.depth \
    /path/to/depth_per_position_all_individuals_BAMs.depth \
    /path/to/output/mean_depth/
```

**Input:** covered depth file and full depth file from Step 2 / Step 1
**Output:**
- `mean_depth_per_position_covered.txt` — mean depth across individuals for covered positions
- `mean_depth_per_position.txt` — mean depth across individuals for all positions

---

### Step 4 — Calculate mean depth per contig

**Script:** `mean_depth_per_contig/mean_depth_per_contig.sh`

```bash
bash mean_depth_per_contig/mean_depth_per_contig.sh \
    /path/to/mean_depth_per_position_covered.txt \
    /path/to/output/mean_depth_per_contig.txt
```

**Input:** mean depth per position file from Step 3b
**Output:** `mean_depth_per_contig.txt` — columns: `contig | mean_depth`

---

### Step 5 — Count positions per depth threshold per contig

**Script:** `number_position_per_depth/number_positions_per_depth_per_contig.sh`

```bash
bash number_position_per_depth/number_positions_per_depth_per_contig.sh \
    /path/to/mean_depth_per_position_covered.txt \
    /path/to/contig_percentage_covered_positions.txt \
    /path/to/output/per_depth/
```

**Input:** mean depth per position (Step 3b) and contig stats (Step 2)
**Output:** `number_position_per_depth_per_contig_sorted_with_total_positions.txt`
— columns: `Contig | Total_Positions | Pos_10X | Pos_15X | Pos_20X | ... | Pos_80X`

---

### Step 6 — Build summary tables by depth threshold

**Script:** `table_contigs_per_depth/table_contigs_per_depth_coverage.sh`

```bash
bash table_contigs_per_depth/table_contigs_per_depth_coverage.sh \
    /path/to/number_position_per_depth_per_contig_sorted_with_total_positions.txt \
    /path/to/output/tables/
```

**Input:** per-depth matrix from Step 5
**Output:**
- `coverage_depth_summary_counts.txt` — number of contigs meeting each % coverage threshold at each depth
- `coverage_depth_summary_percentage.txt` — same as above, expressed as % of total contigs
- `coverage_depth_min_positions_counts.txt` — number of contigs with at least N positions at each depth

---

### Step 7 — Estimate whole-genome sequencing coverage (independent)

**Script:** `whole_coverage_estimation/coverage_estimation_samtool.sbatch`
**Cluster:** SLURM — edit the `[CLUSTER]` lines before submission

```bash
sbatch whole_coverage_estimation/coverage_estimation_samtool.sbatch \
    /path/to/bam_files/ \
    /path/to/output/coverage_results.tsv \
    /path/to/reference.fasta
```

**Input:** BAM directory + reference FASTA
**Output:** tab-delimited file: `Sample | Reads_mapped | Mean_Read_Length | Coverage`

This step is independent of Steps 1–6 and can be run in parallel.

---

### Step 8 — Count covered reciprocal best-hit pairs (downstream MKT input)

**Script:** `count_paired_covered_contig_from_RBH/count_covered_pairs_contigs.sh`

```bash
bash count_paired_covered_contig_from_RBH/count_covered_pairs_contigs.sh \
    /path/to/covered_contigs_species1.txt \
    /path/to/covered_contigs_species2.txt \
    /path/to/best_hits_pairs.tsv
```

**Input:**
- A list of covered contigs for species 1 (one contig name per line)
- A list of covered contigs for species 2
- A tab-delimited file of reciprocal best-hit pairs (`contig_sp1 contig_sp2`)

**Output:** counts of pairs where both contigs meet the coverage threshold (printed to stdout)

---

## Dependencies

| Tool      | Minimum version | Used in       |
|-----------|----------------|---------------|
| samtools  | 1.14           | Steps 1, 7    |
| seqkit    | 2.8            | Step 7        |
| awk       | any            | Steps 2–6, 8  |
| R         | 4.0            | Visualization |

---

## Key Parameters

| Parameter        | Default | Where to change                  |
|------------------|---------|----------------------------------|
| min_depth        | 10      | Step 2 CLI argument              |
| min_individuals  | 6       | Step 2 CLI argument              |
| SLURM account    | —       | `[CLUSTER]` lines in .sbatch     |
| SLURM partition  | —       | `[CLUSTER]` lines in .sbatch     |
| module names     | —       | `[CLUSTER]` lines in .sbatch     |

---

## Notes

- All `.sbatch` scripts require SLURM. For non-SLURM systems, run the commands
  inside the script directly in bash (remove the `#SBATCH` header and
  `module load` lines, ensure tools are in `PATH`).
- The deprecated files `depth_per_individual_io.sbatch` and
  `coverage_estimation.sbatch` are kept for reference only. Use the
  cluster-agnostic versions described above.
- Step 10 requires R with the `rmarkdown`, `ggplot2`, `dplyr`, `tidyr`,
  `scales`, `patchwork`, `knitr`, and `kableExtra` packages. Missing packages
  are installed automatically when the Rmd is rendered.
