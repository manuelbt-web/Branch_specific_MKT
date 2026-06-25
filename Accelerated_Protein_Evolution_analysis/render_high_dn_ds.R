# =============================================================================
# render_high_dn_ds.R — Launcher for high_dn_ds.Rmd
# =============================================================================
#
# PURPOSE
#   Detect accelerated protein evolution and test whether branch-specific MKT
#   positive-selection candidates are enriched among genes with high branch-
#   specific dN/dS or significant branch model likelihood ratio (LRT).
#
# PREREQUISITES
#   Run render_analysis.R (branch_specific_MKT) first.
#   This script reads the output files produced by that step.
#
# HOW TO USE
#   1. Fill in SECTION 1 with the same species output directories you used
#      in render_analysis.R.  You can copy-paste them directly.
#   2. Adjust analysis parameters in SECTION 2 if needed.
#   3. Run:
#        a) RStudio  : open this file → click Source
#        b) R console: source("path/to/render_high_dn_ds.R")
#        c) Terminal : Rscript path/to/render_high_dn_ds.R
#
# FILES READ AUTOMATICALLY (from each species output_dir)
#   branch_specific_MKT_results.tsv   ← all genes, all MKT stats
#   candidates_all.tsv                ← positive candidates (any test)
#   candidates_branch_specific.tsv    ← branch-specific candidates
#
# OUTPUT
#   results/high_dn_ds/high_dn_ds.html               ← HTML report
#   results/high_dn_ds/enrichment_high_dNdS.tsv       ← per-species enrichment stats
#   results/high_dn_ds/enrichment_lrt_acceleration.tsv ← LRT enrichment stats
#   results/high_dn_ds/lrt_classifications.tsv         ← per-gene LRT results
#   results/high_dn_ds/ortholog_pairs_*.tsv            ← cross-species pair tables
# =============================================================================


# ── 0. Dependency check ───────────────────────────────────────────────────────

required_pkgs <- c("rmarkdown", "tidyverse", "ComplexUpset", "kableExtra", "patchwork")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")


# =============================================================================
# SECTION 1 — Focal species
# =============================================================================
#
# Copy the output_dir values from render_analysis.R.
# Only output_dir is needed here — all input files are derived from it.
#
# Format:
#   species_label = list(output_dir = "path/used/in/render_analysis.R")
#
# For a single species: provide one entry (cross-species section is skipped).
# For multiple species: provide all entries (cross-species analysis runs automatically).

SPECIES_LIST <- list(

  Aegilopsspeltoides = list(output_dir = "results/speltoides")

  # Uncomment and fill in to add more species (same dirs as in render_analysis.R):
  # ,
  # Aegilopsmutica = list(output_dir = "results/mutica"),
  # Aegilopstauschii = list(output_dir = "results/tauschii"),
  # Triticumurartu = list(output_dir = "results/urartu")

)


# =============================================================================
# SECTION 2 — Analysis parameters
# =============================================================================

# z-score threshold on log(dN/dS): genes above this are "high dN/dS"
# (z > 1.5 corresponds approximately to the top ~7% of a normal distribution)
Z_THRESHOLD <- 1.5

# p-value threshold for the branch model LRT (chi-squared, df=1)
LRT_P_THRESHOLD <- 0.05

# Output directory for the HTML report and result TSV tables
OUTPUT_DIR <- "results/high_dn_ds"


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
.rmd_file <- file.path(.rmd_dir, "high_dn_ds.Rmd")
if (!file.exists(.rmd_file))
  stop("Rmd file not found: ", .rmd_file)


# ── Validate inputs ───────────────────────────────────────────────────────────

if (length(SPECIES_LIST) == 0)
  stop("SPECIES_LIST is empty. Add at least one species in Section 1.")

for (.sp in names(SPECIES_LIST)) {
  .out <- SPECIES_LIST[[.sp]]$output_dir
  .req <- c("branch_specific_MKT_results.tsv",
            "candidates_all.tsv",
            "candidates_branch_specific.tsv")
  .missing <- .req[!file.exists(file.path(.out, .req))]
  if (length(.missing) > 0)
    warning(sprintf(
      "[%s] Required files not found in '%s': %s\n",
      .sp, .out, paste(.missing, collapse = ", "),
      "\nRun render_analysis.R first."
    ))
}


# ── Build param list for the Rmd ─────────────────────────────────────────────

.species_output_dirs <- setNames(
  lapply(SPECIES_LIST, function(cfg)
    normalizePath(cfg$output_dir, mustWork = FALSE)),
  names(SPECIES_LIST)
)


# ── Create output directory and render ───────────────────────────────────────

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
.out_html <- normalizePath(file.path(OUTPUT_DIR, "high_dn_ds.html"), mustWork = FALSE)

message(sprintf("=== Accelerated Protein Evolution Analysis ==="))
message(sprintf("Species : %s", paste(names(SPECIES_LIST), collapse = ", ")))
message(sprintf("Output  : %s", .out_html))

rmarkdown::render(
  input       = .rmd_file,
  output_file = .out_html,
  params = list(
    species_output_dirs = .species_output_dirs,
    z_threshold         = Z_THRESHOLD,
    lrt_p_threshold     = LRT_P_THRESHOLD,
    output_dir          = normalizePath(OUTPUT_DIR, mustWork = FALSE)
  ),
  envir = new.env(parent = globalenv())
)

message("Done. Report: ", .out_html)
