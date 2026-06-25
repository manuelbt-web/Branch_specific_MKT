# =============================================================================
# render_ppi.R — Launcher for protein_protein_interactions.Rmd
# =============================================================================
#
# PURPOSE
#   Test whether branch-specific positive selection candidates form a more
#   densely connected subnetwork in STRING than expected by chance.
#   Four null models are used:
#     1. Degree-matched permutation
#     2. Degree + gene-length matched permutation (gene length from MKT results)
#     3. Degree-sequence-preserving network rewiring
#     4. Naive random permutation
#
# PREREQUISITES
#   Run branch_specific_MKT/render_analysis.R first (provides gene length data).
#   STRING network files and candidate ID lists must be in this folder
#   (see SECTION 1 below for expected filenames).
#
# HOW TO USE
#   1. Fill in SECTION 1: SPECIES_LIST — one entry per species.
#   2. Adjust SECTION 2: analysis parameters if needed.
#   3. Run:
#        RStudio : open this file → click Source
#        Console : source("render_ppi.R")
#        Terminal: Rscript render_ppi.R
#
# INPUT FILES (must be in the protein_protein_interactions/ folder)
#   candidate_ids_file : plain-text file, one STRING ID per line
#                        (e.g. string_id_speltoides_bs_candidates.txt)
#   network_edge_file  : STRING / Cytoscape edge CSV export
#                        (e.g. STRING_network_speltoides_default_edge.csv)
#   network_node_file  : STRING / Cytoscape node CSV export — used for gene
#                        descriptions and aestivum_ortholog→STRING ID mapping
#                        (e.g. STRING_network_speltoides_default_node.csv)
#   output_dir         : output directory from branch_specific_MKT/render_analysis.R
#                        Used to read branch_specific_MKT_results.tsv for gene
#                        length data (improves degree+length matched null).
#
# OUTPUTS (in results/ppi/<Species>/)
#   network_summary.tsv                  — observed metrics + p-values per species
#   candidate_candidate_interactions.tsv — edges between candidate pairs
#   candidate_connectivity.tsv           — per-candidate degree breakdown
#   jackknife_influence.tsv              — leave-one-out density contribution
#   results/ppi/all_species_network_summary.tsv  — cross-species comparison
#   results/ppi/ppi_analysis.html        — full HTML report
# =============================================================================


# ── 0. Dependency check ───────────────────────────────────────────────────────

required_cran <- c("rmarkdown", "tidyverse", "kableExtra", "patchwork")
required_cran_other <- c("igraph")  # igraph is CRAN but large; listed separately

missing_cran <- setdiff(
  c(required_cran, required_cran_other),
  rownames(installed.packages())
)

if (length(missing_cran) > 0)
  stop(
    "Missing CRAN packages: ", paste(missing_cran, collapse = ", "),
    "\nInstall with:\n  install.packages(c(",
    paste0('"', missing_cran, '"', collapse = ", "), "))"
  )


# =============================================================================
# SECTION 1 — Species configurations
# =============================================================================
#
# Each entry specifies:
#   label             : display name used in the report
#   output_dir        : from branch_specific_MKT/render_analysis.R
#                       (the script reads branch_specific_MKT_results.tsv here
#                        for gene length data used in the degree+length null)
#   candidate_ids_file: plain-text file listing candidate STRING IDs (one per line)
#   network_edge_file : STRING Cytoscape edge CSV export
#   network_node_file : STRING Cytoscape node CSV export (for gene descriptions)
#
# All file paths are relative to the protein_protein_interactions/ folder.
# The engine resolves them as absolute paths automatically.

SPECIES_LIST <- list(

  Aegilopsspeltoides = list(
    label              = "Ae. speltoides",
    output_dir         = "results/speltoides",
    candidate_ids_file = "string_id_speltoides_bs_candidates.txt",
    network_edge_file  = "STRING_network_speltoides_default_edge.csv",
    network_node_file  = "STRING_network_speltoides_default_node.csv"
  )

  # Uncomment to add more species:
  # ,
  # Aegilopsmutica = list(
  #   label              = "Ae. mutica",
  #   output_dir         = "results/mutica",
  #   candidate_ids_file = "string_id_mutica_bs_candidates.txt",
  #   network_edge_file  = "STRING_network_mutica_default_edge.csv",
  #   network_node_file  = "STRING_network_mutica_default node.csv"
  # ),
  # Aegilopstauschii = list(
  #   label              = "Ae. tauschii",
  #   output_dir         = "results/tauschii",
  #   candidate_ids_file = "string_id_tauschii_bs_candidates.txt",
  #   network_edge_file  = "STRING_network_tauschii_default_edge.csv",
  #   network_node_file  = "STRING_network_tauschii_default_node.csv"
  # ),
  # Triticumurartu = list(
  #   label              = "T. urartu",
  #   output_dir         = "results/urartu",
  #   candidate_ids_file = "string_id_urartu_bs_candidates.txt",
  #   network_edge_file  = "STRING_network_urartu_default_edge.csv",
  #   network_node_file  = "STRING_network_urartu_default_node.csv"
  # )

)


# =============================================================================
# SECTION 2 — Analysis parameters
# =============================================================================

# STRING combined score threshold (0–1). Interactions below this are discarded.
# 0.4 = medium confidence, 0.7 = high confidence, 0.9 = very high confidence.
SCORE_THRESHOLD <- 0.7

# Number of permutations for each null model.
# 200 for quick checks; 1000+ recommended for publication.
N_PERMUTATIONS <- 200

# niter = ecount(graph) * REWIRE_MULTIPLIER for the degree-preserving rewiring.
REWIRE_MULTIPLIER <- 10

# Random seed for reproducibility.
SEED <- 123

# Output directory for HTML report and result TSVs.
OUTPUT_DIR <- "results/ppi"


# =============================================================================
# ENGINE — do not edit below this line
# =============================================================================

# ── Locate script and Rmd directories ─────────────────────────────────────────

.script_dir <- tryCatch(
  dirname(normalizePath(
    c(sys.frame(1)$ofile,
      sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
    )[1]
  )),
  error = function(e) {
    warning("Could not auto-detect script directory; using getwd().")
    getwd()
  }
)

.rmd_file <- file.path(.script_dir, "protein_protein_interactions.Rmd")
if (!file.exists(.rmd_file))
  stop("Rmd not found: ", .rmd_file)


# ── Validate species configs ───────────────────────────────────────────────────

if (length(SPECIES_LIST) == 0)
  stop("SPECIES_LIST is empty. Add at least one species in Section 1.")

.species_configs <- setNames(
  lapply(names(SPECIES_LIST), function(.sp) {
    .cfg <- SPECIES_LIST[[.sp]]

    # Resolve file paths relative to the script directory
    .cfg$candidate_ids_file <- normalizePath(
      file.path(.script_dir, .cfg$candidate_ids_file), mustWork = FALSE
    )
    .cfg$network_edge_file <- normalizePath(
      file.path(.script_dir, .cfg$network_edge_file), mustWork = FALSE
    )
    .cfg$network_node_file <- normalizePath(
      file.path(.script_dir, .cfg$network_node_file), mustWork = FALSE
    )

    # MKT results (gene length data)
    .mkt_file <- file.path(.cfg$output_dir, "branch_specific_MKT_results.tsv")
    if (!file.exists(.mkt_file))
      warning(sprintf(
        "[%s] branch_specific_MKT_results.tsv not found in '%s'.\n",
        .sp, .cfg$output_dir
      ), "  → Gene length matching will fall back to degree-only. Run render_analysis.R first.")

    # Required input files
    for (.key in c("candidate_ids_file", "network_edge_file")) {
      if (!file.exists(.cfg[[.key]]))
        stop(sprintf("[%s] File not found for '%s': %s", .sp, .key, .cfg[[.key]]))
    }

    .cfg
  }),
  names(SPECIES_LIST)
)


# ── Render ────────────────────────────────────────────────────────────────────

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
.out_html <- normalizePath(
  file.path(OUTPUT_DIR, "ppi_analysis.html"), mustWork = FALSE
)

message("=== PPI Network Analysis ===")
message(sprintf("Species        : %s", paste(names(SPECIES_LIST), collapse = ", ")))
message(sprintf("Score threshold: %.2f", SCORE_THRESHOLD))
message(sprintf("Permutations   : %d",   N_PERMUTATIONS))
message(sprintf("Output         : %s",   .out_html))

rmarkdown::render(
  input       = .rmd_file,
  output_file = .out_html,
  params = list(
    experiment_name = "Positive Selection Candidates",
    species_configs = .species_configs,
    score_threshold = SCORE_THRESHOLD,
    n_permutations  = N_PERMUTATIONS,
    rewire_multiplier = REWIRE_MULTIPLIER,
    seed            = SEED,
    output_dir      = normalizePath(OUTPUT_DIR, mustWork = FALSE)
  ),
  envir = new.env(parent = globalenv())
)

message("Done. Report: ", .out_html)
