# =============================================================================
# render_expression.R — Launcher for expression_analysis.Rmd
# =============================================================================
#
# PURPOSE
#   Test whether branch-specific positive selection candidates (from
#   branch_specific_MKT/render_analysis.R) are enriched among differentially
#   expressed genes under various stress conditions (cold, heat, Fusarium,
#   powdery mildew, Septoria, phosphorus starvation, PAMP response).
#
# PREREQUISITES
#   Run branch_specific_MKT/render_analysis.R first.
#   This script reads branch_specific_MKT_results.tsv from each output_dir.
#
# HOW TO USE
#   1. Fill in SECTION 1: SPECIES_LIST with your species and output directories.
#   2. Optionally adjust SECTION 2: ANALYSIS PARAMETERS.
#   3. Optionally enable/disable individual experiments in SECTION 3.
#   4. Run:
#        RStudio : open this file → click Source
#        Console : source("render_expression.R")
#        Terminal: Rscript render_expression.R
#
# OUTPUTS
#   One HTML report per experiment in:
#     results/expression/<experiment_name>.html
#   TSV result tables per experiment in:
#     results/expression/<experiment_name>_<species>_fisher_tests.tsv
#     results/expression/<experiment_name>_cross_species_summary.tsv
# =============================================================================


# ── 0. Dependency check ───────────────────────────────────────────────────────

required_pkgs <- c("rmarkdown", "tidyverse", "DESeq2", "patchwork", "kableExtra")
missing_cran  <- setdiff(
  setdiff(required_pkgs, "DESeq2"),
  rownames(installed.packages())
)
missing_bioc  <- if (!requireNamespace("DESeq2", quietly = TRUE)) "DESeq2" else character(0)

if (length(missing_cran) > 0)
  stop("Missing CRAN packages: ", paste(missing_cran, collapse = ", "),
       "\nInstall with:\n  install.packages(c(",
       paste0('"', missing_cran, '"', collapse = ", "), "))")

if (length(missing_bioc) > 0)
  stop("Missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "),
       "\nInstall with:\n",
       "  if (!requireNamespace('BiocManager', quietly=TRUE)) ",
       "install.packages('BiocManager')\n",
       "  BiocManager::install(c(",
       paste0('"', missing_bioc, '"', collapse = ", "), "))")


# =============================================================================
# SECTION 1 — Species configurations
# =============================================================================
#
# Each entry must have:
#   output_dir : the same output_dir used in branch_specific_MKT/render_analysis.R
#                The script reads branch_specific_MKT_results.tsv from here.
#
# T. aestivum orthologs are extracted from the aestivum_ortholog column in
# branch_specific_MKT_results.tsv to map focal-species genes to wheat expression data.

SPECIES_LIST <- list(

  Aegilopsspeltoides = list(
    output_dir = "results/speltoides"
  )

  # Uncomment to add more species:
  # ,
  # Aegilopsmutica = list(
  #   output_dir = "results/mutica"
  # ),
  # Aegilopstauschii = list(
  #   output_dir = "results/tauschii"
  # ),
  # Triticumurartu = list(
  #   output_dir = "results/urartu"
  # )

)


# =============================================================================
# SECTION 2 — Analysis parameters
# =============================================================================

# Which MKT candidate set defines "positively selected" genes:
#   "bs_any"    = union of branch-specific standard + imputed MKT (recommended)
#   "bs_impMKT" = branch-specific imputed MKT only
#   "bs_mkt"    = branch-specific standard MKT only
#   "any"       = union of all 4 MKT tests (standard, imputed, branch-specific)
CANDIDATE_TYPE <- "bs_any"

# |log2FoldChange| threshold for classifying genes as up- or down-regulated
LFC_THRESHOLD <- 0.5

# BH-adjusted p-value threshold for differential expression
PADJ_THRESHOLD <- 0.05

# Output directory for HTML reports and result tables
OUTPUT_DIR <- "results/expression"


# =============================================================================
# SECTION 3 — Experiment definitions
# =============================================================================
#
# Each experiment entry specifies:
#   name              : label used in filenames and report titles
#   count_file        : basename of the count TSV (in Expression_analysis/)
#   control_condition : reference level for DESeq2 (relevel)
#   sample_conditions : character vector — one element per sample column in
#                       the count TSV, in the SAME ORDER as the columns
#   contrasts         : (optional) named list treatment → reference for
#                       non-standard contrast specification (see Fusarium,
#                       Septoria examples below). If empty, all non-control
#                       conditions are tested vs control automatically.
#
# Set enabled = FALSE to skip an experiment without deleting its config.

EXPERIMENTS <- list(

  # ── Cold stress (SRP064598) ─────────────────────────────────────────────────
  cold = list(
    enabled           = TRUE,
    name              = "cold_stress",
    count_file        = "SRP064598_count.tsv",
    control_condition = "none",
    sample_conditions = c(
      "tissue culture following 10 days cold treatment",
      "tissue culture following 10 days cold treatment",
      "tissue culture following 10 days cold treatment",
      "10 days cold treatment",
      "10 days cold treatment",
      "10 days cold treatment",
      "none", "none", "none"
    ),
    contrasts = list()
  ),

  # ── Heat / drought stress (SRP045409) ──────────────────────────────────────
  heat = list(
    enabled           = TRUE,
    name              = "heat_stress",
    count_file        = "SRP045409_count.tsv",
    control_condition = "control",
    sample_conditions = c(
      "6_hour_drought_heat", "6_hour_drought_heat",
      "1_hour_drought_heat", "1_hour_drought_heat",
      "6_hour_heat",         "6_hour_heat",
      "1_hour_heat",         "1_hour_heat",
      "6_hour_drought",      "6_hour_drought",
      "1_hour_drought",      "1_hour_drought",
      "control",             "control"
    ),
    contrasts = list()
  ),

  # ── Fusarium graminearum (ERP013829) ────────────────────────────────────────
  # Time-matched contrasts: Fusarium_X hours vs mock_X hours
  fusarium = list(
    enabled           = TRUE,
    name              = "fusarium",
    count_file        = "ERP013829_count.tsv",
    control_condition = "mock_6 hours",   # reference for relevel (not used directly)
    sample_conditions = c(
      rep("mock_6 hours",      3), rep("mock_48 hours",     3),
      rep("mock_3 hours",      3), rep("mock_36 hours",     3),
      rep("mock_24 hours",     3), rep("mock_12 hours",     3),
      rep("Fusarium_6 hours",  3), rep("Fusarium_48 hours", 3),
      rep("Fusarium_3 hours",  3), rep("Fusarium_36 hours", 3),
      rep("Fusarium_24 hours", 3), rep("Fusarium_12 hours", 3),
      rep("mock_6 hours",      3), rep("mock_48 hours",     3),
      rep("mock_3 hours",      3), rep("mock_36 hours",     3),
      rep("mock_24 hours",     3), rep("mock_12 hours",     3),
      rep("Fusarium_6 hours",  3), rep("Fusarium_48 hours", 3),
      rep("Fusarium_3 hours",  3), rep("Fusarium_36 hours", 3),
      rep("Fusarium_24 hours", 3), rep("Fusarium_12 hours", 3)
    ),
    contrasts = list(
      "Fusarium_3 hours"  = "mock_3 hours",
      "Fusarium_6 hours"  = "mock_6 hours",
      "Fusarium_12 hours" = "mock_12 hours",
      "Fusarium_24 hours" = "mock_24 hours",
      "Fusarium_36 hours" = "mock_36 hours",
      "Fusarium_48 hours" = "mock_48 hours"
    )
  ),

  # ── Powdery mildew + stripe rust (SRP041017) ────────────────────────────────
  powdery = list(
    enabled           = TRUE,
    name              = "powdery_mildew",
    count_file        = "SRP041017_count.tsv",
    control_condition = "none",
    sample_conditions = c(
      rep("Powdery_mildew_E09_72h", 3),
      rep("Powdery_mildew_E09_48h", 3),
      rep("Powdery_mildew_E09_24h", 3),
      rep("none",                   3),
      rep("Stripe_rust_CYR31_72h",  3),
      rep("Stripe_rust_CYR31_48h",  3),
      rep("Stripe_rust_CYR31_24h",  3)
    ),
    contrasts = list()
  ),

  # ── Septoria (Zymoseptoria tritici, ERP009837) ──────────────────────────────
  # Time-matched contrasts: Zymoseptoria_X vs mock_X
  septoria = list(
    enabled           = TRUE,
    name              = "septoria",
    count_file        = "ERP009837_count.tsv",
    control_condition = "mock_1 day",
    sample_conditions = c(
      "Zymoseptoria_1 day",   "Zymoseptoria_1 day",
      "mock_21 days",         "mock_21 days",
      "mock_14 days",         "mock_9 days",
      "mock_9 days",          "mock_9 days",
      "mock_4 days",          "mock_4 days",         "mock_4 days",
      "mock_1 day",           "mock_1 day",          "mock_1 day",
      "Zymoseptoria_21 days", "Zymoseptoria_21 days",
      "Zymoseptoria_14 days", "Zymoseptoria_9 days",
      "Zymoseptoria_21 days", "Zymoseptoria_14 days",
      "Zymoseptoria_14 days", "Zymoseptoria_9 days",
      "Zymoseptoria_9 days",  "Zymoseptoria_4 days",
      "Zymoseptoria_4 days",  "Zymoseptoria_4 days",
      "Zymoseptoria_1 day",   "mock_21 days",
      "mock_14 days",         "mock_14 days"
    ),
    contrasts = list(
      "Zymoseptoria_1 day"   = "mock_1 day",
      "Zymoseptoria_4 days"  = "mock_4 days",
      "Zymoseptoria_9 days"  = "mock_9 days",
      "Zymoseptoria_14 days" = "mock_14 days",
      "Zymoseptoria_21 days" = "mock_21 days"
    )
  ),

  # ── Phosphorus starvation (DRP000768) ───────────────────────────────────────
  phosphorus = list(
    enabled           = TRUE,
    name              = "phosphorus",
    count_file        = "DRP000768_count.tsv",
    control_condition = "none",
    sample_conditions = c(
      "phosphorus_starvation_10d", "none",
      "phosphorus_starvation_10d", "phosphorus_starvation_10d",
      "none",                      "none",
      "phosphorus_starvation_10d", "none",
      "none",                      "phosphorus_starvation_10d",
      "phosphorus_starvation_10d", "none"
    ),
    contrasts = list()
  ),

  # ── PAMP-triggered immune response (PAMP_Triggered_Imune_Response) ───────────
  pamp = list(
    enabled           = TRUE,
    name              = "PAMP",
    count_file        = "PAMP_Triggered_Imune_Response_count.tsv",
    control_condition = "water",
    sample_conditions = c(
      rep("water",                  9),
      rep("chitin_1g_per_L",        6),
      rep("flg22_500nM",            6)
    ),
    contrasts = list()
  )

)


# =============================================================================
# ENGINE — do not edit below this line
# =============================================================================

# ── Locate Rmd ────────────────────────────────────────────────────────────────

.rmd_dir <- tryCatch(
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
.rmd_file <- file.path(.rmd_dir, "expression_analysis.Rmd")
if (!file.exists(.rmd_file))
  stop("Rmd not found: ", .rmd_file)


# ── Validate species ──────────────────────────────────────────────────────────

if (length(SPECIES_LIST) == 0)
  stop("SPECIES_LIST is empty. Add at least one species in Section 1.")

.species_output_dirs <- setNames(
  lapply(names(SPECIES_LIST), function(.sp) {
    .od <- SPECIES_LIST[[.sp]]$output_dir
    .f  <- file.path(.od, "branch_specific_MKT_results.tsv")
    if (!file.exists(.f))
      warning(sprintf("[%s] File not found: %s\n  → Run render_analysis.R first.", .sp, .f))
    normalizePath(.od, mustWork = FALSE)
  }),
  names(SPECIES_LIST)
)


# ── Run each experiment ───────────────────────────────────────────────────────

dir.create(normalizePath(OUTPUT_DIR, mustWork = FALSE), recursive = TRUE, showWarnings = FALSE)

.enabled <- Filter(function(e) isTRUE(e$enabled), EXPERIMENTS)
message(sprintf("=== Expression Analysis ==="))
message(sprintf("Species   : %s", paste(names(SPECIES_LIST), collapse = ", ")))
message(sprintf("Candidate : %s", CANDIDATE_TYPE))
message(sprintf("Experiments to run: %s", paste(names(.enabled), collapse = ", ")))

for (.exp_name in names(.enabled)) {
  .exp <- .enabled[[.exp_name]]

  # Resolve count file path (relative to this script's directory)
  .count_path <- normalizePath(
    file.path(.rmd_dir, .exp$count_file),
    mustWork = FALSE
  )
  if (!file.exists(.count_path))
    stop(sprintf("[%s] Count file not found: %s", .exp_name, .count_path))

  .out_html <- normalizePath(
    file.path(OUTPUT_DIR, paste0(.exp$name, ".html")),
    mustWork = FALSE
  )
  .out_subdir <- normalizePath(
    file.path(OUTPUT_DIR, .exp$name),
    mustWork = FALSE
  )

  message(sprintf("\n--- Running: %s ---", .exp$name))
  message(sprintf("    Count file : %s", basename(.count_path)))
  message(sprintf("    Output     : %s", .out_html))

  rmarkdown::render(
    input       = .rmd_file,
    output_file = .out_html,
    params = list(
      experiment_name    = .exp$name,
      count_file         = .count_path,
      control_condition  = .exp$control_condition,
      sample_conditions  = .exp$sample_conditions,
      contrasts          = .exp$contrasts,
      species_output_dirs = .species_output_dirs,
      candidate_type     = CANDIDATE_TYPE,
      lfc_threshold      = LFC_THRESHOLD,
      padj_threshold     = PADJ_THRESHOLD,
      output_dir         = .out_subdir
    ),
    envir = new.env(parent = globalenv())
  )

  message(sprintf("    Done. Report: %s", .out_html))
}

message("\n=== All experiments completed ===")
message(sprintf("Results in: %s", normalizePath(OUTPUT_DIR, mustWork = FALSE)))
