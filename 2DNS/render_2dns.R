# =============================================================================
# render_2dns.R — Launcher for 2dns_analysis.Rmd
# =============================================================================
#
# PURPOSE
#   Run a four-layer gene set MK analysis (2DNS) for one or more species:
#     Layer 1 — Fisher exact test (two-sided + one-sided)
#     Layer 2 — Benjamini-Hochberg FDR correction
#     Layer 3 — Stratified permutation (length-matched null from MKT gene pool)
#     Layer 4 — Jackknife leave-one-gene-out (optional; requires gene_level_file)
#
# PREREQUISITES
#   Run branch_specific_MKT/render_analysis.R first for each species.
#   Pathway-level input CSVs (input_2DNS_*.csv) must be in this folder.
#
# HOW TO USE
#   1. Fill in SECTION 1 — SPECIES_LIST (one entry per focal species).
#   2. Adjust SECTION 2 — parameters if needed.
#   3. Run:
#        RStudio : open this file → click Source
#        Console : source("render_2dns.R")
#        Terminal: Rscript render_2dns.R
#
# INPUTS (relative paths are resolved from this script's directory)
#   input_file     : pathway-level CSV, e.g. input_2DNS_speltoides.csv
#                    Columns: pathway_name, PN, PS, DN, DS, n_genes [, set_length]
#   output_dir     : directory from branch_specific_MKT/render_analysis.R
#                    Used to read branch_specific_MKT_results.tsv for the gene
#                    pool (permutation null) and candidate gene identification.
#   gene_level_file: (OPTIONAL) gene-level pathway CSV for the jackknife layer.
#                    Columns: gene_id, pathway_name, PN, PS, DN, DS.
#                    Set to NULL to skip Layer 4.
#
# OUTPUTS (in OUTPUT_DIR)
#   2dns_analysis.html          — full HTML report
#   fisher_results.tsv          — Layer 1 results for all pathways and species
#   permutation_results.tsv     — Layer 3 results (Fisher-significant pathways)
#   comprehensive_summary.tsv   — all layers combined per pathway
#   candidates.tsv              — BS candidate genes from MKT results
# =============================================================================


# ── 0. Dependency check ───────────────────────────────────────────────────────

required_pkgs <- c("rmarkdown", "tidyverse", "kableExtra", "patchwork")
missing_pkgs  <- setdiff(required_pkgs, rownames(installed.packages()))

if (length(missing_pkgs) > 0)
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall with:\n  install.packages(c(",
    paste0('"', missing_pkgs, '"', collapse = ", "), "))"
  )


# =============================================================================
# SECTION 1 — Species configurations
# =============================================================================
#
# Each species entry requires:
#   label          : display name used in the report
#   input_file     : pathway-level CSV in this folder (input_2DNS_*.csv)
#   output_dir     : from branch_specific_MKT/render_analysis.R
#                    (reads branch_specific_MKT_results.tsv for gene pool +
#                     candidate identification)
#   gene_level_file: gene-level pathway CSV (gene_id, pathway_name, PN, PS,
#                    DN, DS, length). Strongly recommended for every species:
#                    without it, Layer 3's permutation null can't control for
#                    pathway overlap structure (over-estimates significant
#                    pathways), and Layer 4 (jackknife) plus the candidate x
#                    pathway enrichment test are skipped. Set to NULL to
#                    fall back to the conservative length-only null instead.
#
# File paths are relative to the 2DNS/ folder (resolved automatically).

SPECIES_LIST <- list(

  speltoides = list(
    label           = "Ae. speltoides",
    input_file      = "input_2DNS_speltoides.csv",
    output_dir      = "results/speltoides",        # from render_analysis.R
    gene_level_file = "genes_stats_speltoides_pathways.csv"
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


# =============================================================================
# SECTION 2 — Analysis parameters
# =============================================================================

# Minimum number of genes per pathway (pathways with fewer genes are excluded).
MIN_GENES <- 8L

# Number of permutations per pathway (Layer 3).
# 200 for quick checks; 1000+ recommended for publication.
N_PERMUTATIONS <- 1000L

# Number of gene-length strata for the permutation null (Layer 3).
LENGTH_STRATA_N <- 5L

# Significance thresholds.
FISHER_P_THRESHOLD <- 0.05    # Layer 1 raw p-value threshold
FDR_Q_THRESHOLD    <- 0.05    # Layer 2 BH q-value threshold
PERM_P_THRESHOLD   <- 0.05    # Layer 3 permutation p-value threshold

# Random seed for reproducibility.
SEED <- 123L

# Candidate type — which MKT test(s) define positive selection candidates.
#   "bs_any"    : branch-specific impMKT OR standard (recommended)
#   "bs_impMKT" : branch-specific imputed MKT only
#   "bs_mkt"    : branch-specific standard MKT only
#   "any"       : union of all four MKT tests
CANDIDATE_TYPE <- "bs_any"

# Output directory for the HTML report and result TSVs.
OUTPUT_DIR <- "results/2dns"


# =============================================================================
# ENGINE — do not edit below this line
# =============================================================================

# ── Locate script and Rmd ────────────────────────────────────────────────────

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

.rmd_file <- file.path(.script_dir, "2dns_analysis.Rmd")
if (!file.exists(.rmd_file))
  stop("Rmd not found: ", .rmd_file)


# ── Validate and resolve paths ────────────────────────────────────────────────

if (length(SPECIES_LIST) == 0)
  stop("SPECIES_LIST is empty. Add at least one species in Section 1.")

.species_configs <- setNames(
  lapply(names(SPECIES_LIST), function(.sp) {
    .cfg <- SPECIES_LIST[[.sp]]

    # Required: input_file
    if (is.null(.cfg$input_file))
      stop("'input_file' missing for species '", .sp, "'")
    .cfg$input_file <- normalizePath(
      file.path(.script_dir, .cfg$input_file), mustWork = FALSE
    )
    if (!file.exists(.cfg$input_file))
      stop(sprintf("[%s] input_file not found: %s", .sp, .cfg$input_file))

    # Required: output_dir (for branch_specific_MKT_results.tsv)
    if (is.null(.cfg$output_dir))
      stop("'output_dir' missing for species '", .sp, "'")
    .mkt_file <- file.path(.cfg$output_dir, "branch_specific_MKT_results.tsv")
    if (!file.exists(.mkt_file))
      warning(sprintf(
        "[%s] branch_specific_MKT_results.tsv not found in '%s'.\n",
        .sp, .cfg$output_dir
      ), "  → Gene pool will be unavailable; permutation and candidate steps skipped.",
        call. = FALSE)

    # Optional: gene_level_file (for jackknife Layer 4)
    if (!is.null(.cfg$gene_level_file)) {
      .cfg$gene_level_file <- normalizePath(
        file.path(.script_dir, .cfg$gene_level_file), mustWork = FALSE
      )
      if (!file.exists(.cfg$gene_level_file))
        warning(sprintf("[%s] gene_level_file not found: %s — jackknife skipped.",
                        .sp, .cfg$gene_level_file), call. = FALSE)
    }

    .cfg
  }),
  names(SPECIES_LIST)
)


# ── Render ────────────────────────────────────────────────────────────────────

.out_html <- normalizePath(
  file.path(.script_dir, OUTPUT_DIR, "2dns_analysis.html"), mustWork = FALSE
)
dir.create(dirname(.out_html), recursive = TRUE, showWarnings = FALSE)

message("=== 2DNS Gene Set MK Analysis ===")
message(sprintf("Species        : %s", paste(names(SPECIES_LIST), collapse = ", ")))
message(sprintf("Min genes      : %d",   MIN_GENES))
message(sprintf("Permutations   : %d",   N_PERMUTATIONS))
message(sprintf("Candidate type : %s",   CANDIDATE_TYPE))
message(sprintf("Output         : %s",   .out_html))

rmarkdown::render(
  input       = .rmd_file,
  output_file = .out_html,
  params = list(
    experiment_name   = "Gene Set MK Analysis (2DNS)",
    species_configs   = .species_configs,
    min_genes         = MIN_GENES,
    n_permutations    = N_PERMUTATIONS,
    length_strata_n   = LENGTH_STRATA_N,
    fisher_p_threshold = FISHER_P_THRESHOLD,
    fdr_q_threshold    = FDR_Q_THRESHOLD,
    perm_p_threshold   = PERM_P_THRESHOLD,
    seed              = SEED,
    candidate_type    = CANDIDATE_TYPE,
    output_dir        = normalizePath(
      file.path(.script_dir, OUTPUT_DIR), mustWork = FALSE
    )
  ),
  envir = new.env(parent = globalenv())
)

message("Done. Report: ", .out_html)
