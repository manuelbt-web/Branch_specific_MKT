# protein_protein_interactions

## Objective

Test whether genes under branch-specific positive selection (from
`branch_specific_MKT/render_analysis.R`) form a more densely connected
subnetwork in STRING than expected by chance.

**Statistical approach:** four null models compared against the observed
candidate subgraph density:

| Null model | Description |
|-----------|-------------|
| Degree-matched | Random gene sets matched for network degree (±bin) |
| Degree + length | Degree AND gene length jointly matched |
| Network rewiring | Degree-sequence-preserving edge rewiring (Maslov-Sneppen) |
| Naive | Uniformly random gene sets of the same size |

---

## Contents

| File | Description |
|------|-------------|
| `render_ppi.R` | **Main launcher** — fill in species configs and run |
| `protein_protein_interactions.Rmd` | Analysis template (called by the launcher) |
| `string_id_*_bs_candidates.txt` | STRING IDs of branch-specific positive selection candidates |
| `STRING_network_*_default_edge.csv` | STRING interaction network (Cytoscape edge export) |
| `STRING_network_*_default_node.csv` | STRING node annotations (Cytoscape node export) |

---

## Prerequisites

1. **Run `branch_specific_MKT/render_analysis.R` first.**
   Reads `branch_specific_MKT_results.tsv` from each `output_dir`
   to obtain gene length data for the degree+length matched null.
   Analysis still runs without this file (falls back to degree-only null).

2. **STRING network files must be in this folder** (already present).

---

## Dependencies

```r
install.packages(c("tidyverse", "igraph", "kableExtra", "patchwork", "rmarkdown"))
```

---

## How to run

### Step 1 — Fill in `SPECIES_LIST` in `render_ppi.R`

```r
SPECIES_LIST <- list(

  Aegilopsspeltoides = list(
    label              = "Ae. speltoides",
    output_dir         = "results/speltoides",   # from render_analysis.R
    candidate_ids_file = "string_id_speltoides_bs_candidates.txt",
    network_edge_file  = "STRING_network_speltoides_default_edge.csv",
    network_node_file  = "STRING_network_speltoides_default_node.csv"
  )

)
```

### Step 2 — Adjust parameters (optional)

```r
SCORE_THRESHOLD <- 0.7    # STRING confidence threshold
N_PERMUTATIONS  <- 1000   # increase for publication (200 = quick check)
```

### Step 3 — Run

```r
source("render_ppi.R")   # RStudio / R console
Rscript render_ppi.R     # terminal
```

---

## Outputs

All written to `results/ppi/`:

| File | Content |
|------|---------|
| `ppi_analysis.html` | Interactive HTML report |
| `<Species>/network_summary.tsv` | Observed density + p-values per null model |
| `<Species>/candidate_candidate_interactions.tsv` | Candidate–candidate interaction pairs |
| `<Species>/candidate_connectivity.tsv` | Per-candidate degree breakdown |
| `<Species>/jackknife_influence.tsv` | Leave-one-out subgraph density contributions |
| `all_species_network_summary.tsv` | Cross-species comparison table |

---

## STRING network files

The network files were downloaded from [STRING](https://string-db.org) using
T. *aestivum* protein IDs (taxon 4565) and exported from Cytoscape.

- **Edge file format:** Cytoscape STRING edge CSV with `shared name` column
  (`"4565.XXXXX (pp) 4565.YYYYY"`) and `stringdb::score` (combined confidence).
- **Node file format:** Cytoscape STRING node CSV with `name` (full STRING ID),
  `display name` (short UniProt ID), and `stringdb::description`.

To add a new species: download the STRING network for that species' T. *aestivum*
orthologs, export edge + node CSVs from Cytoscape, and add an entry to
`SPECIES_LIST` in `render_ppi.R`.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `File not found: branch_specific_MKT_results.tsv` | Pipeline not run | Run `render_analysis.R`; analysis still works without it |
| 0 candidates in network | Wrong STRING IDs | Check that `.txt` IDs match the taxon prefix in the edge file (e.g. `4565.`) |
| `Cannot detect node columns` | Non-standard edge file format | Ensure `shared name` column exists, or use a file with `node1`/`node2` columns |
| Rewiring null very slow | Large network + high rewire_multiplier | Reduce `REWIRE_MULTIPLIER` to 5 or lower `N_PERMUTATIONS` |
| Score column not found | Non-standard column name | Rename to `combined_score` in the edge CSV |

---

← [Scripts — overview](../README.md)
← [branch_specific_MKT — run this first](../MKT_3_species/branch_specific_MKT/README.md)
