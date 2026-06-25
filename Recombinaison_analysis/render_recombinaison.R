# =============================================================================
# render_recombinaison.R — Launcher for recombinaison.Rmd
# =============================================================================
#
# PURPOSE
#   Test whether branch-specific positive selection candidates are enriched in
#   high-recombination regions (binomial logistic GLM, uncertainty-adjusted).
#
# PREREQUISITES
#   Run render_analysis.R (branch_specific_MKT) first.
#   This script reads the candidate and stats files it produces.
#
# HOW TO USE
#   1. Fill in SECTION 1 with your species configurations.
#      Each species needs the same output_dir used in render_analysis.R,
#      plus two BED files (gene coordinates and recombination map).
#   2. Adjust analysis parameters in SECTION 2 if needed.
#   3. Run:
#        a) RStudio  : open this file → click Source
#        b) R console: source("path/to/render_recombinaison.R")
#        c) Terminal : Rscript path/to/render_recombinaison.R
#
# OUTPUT
#   results/recombination/recombination.html               ← HTML report
#   results/recombination/gene_recombination_assignments.tsv
#   results/recombination/glm_coefficients.tsv
#   results/recombination/sensitivity_analysis.tsv
#
# CHROMOSOME NAMING NOTE
#   The chromosome names in gene_coords_bed and rec_map_bed must match exactly
#   (e.g., both use "Chr1A" or both use "1A"). Mismatches produce a warning
#   and genes on unmatched chromosomes are dropped from the analysis.
# =============================================================================


# ── 0. Dependency check ───────────────────────────────────────────────────────

required_pkgs <- c("rmarkdown", "tidyverse", "GenomicRanges", "IRanges",
                   "S4Vectors", "kableExtra", "patchwork")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

# GenomicRanges, IRanges, and S4Vectors are Bioconductor packages
bioc_pkgs   <- intersect(c("GenomicRanges", "IRanges", "S4Vectors"), missing_pkgs)
cran_pkgs   <- setdiff(missing_pkgs, bioc_pkgs)

if (length(cran_pkgs) > 0)
  stop("Missing CRAN packages: ", paste(cran_pkgs, collapse = ", "),
       "\nInstall with:\n  install.packages(c(",
       paste0('"', cran_pkgs, '"', collapse = ", "), "))")

if (length(bioc_pkgs) > 0)
  stop("Missing Bioconductor packages: ", paste(bioc_pkgs, collapse = ", "),
       "\nInstall with:\n",
       "  if (!requireNamespace('BiocManager', quietly=TRUE)) ",
       "install.packages('BiocManager')\n",
       "  BiocManager::install(c(",
       paste0('"', bioc_pkgs, '"', collapse = ", "), "))")


# =============================================================================
# SECTION 1 — Species configurations
# =============================================================================
#
# Add one entry per focal species. Each entry must have:
#
#   output_dir      : the same output_dir used in render_analysis.R — the script
#                     reads branch_specific_MKT_results.tsv and
#                     candidates_branch_specific.tsv from here automatically.
#
#   gene_coords_bed : BED file with gene genomic coordinates.
#                     NO header, 4 tab-separated columns:
#                       chr  start  end  gene_id
#                     Generate with: python prepare_gene_coords.py --gff annotation.gff3 --out data/species_gene_coords.bed
#                     Run with --preview first to check feature type and ID attribute.
#
#   rec_map_bed     : BED file with recombination rate per interval.
#                     NO header, 4 tab-separated columns:
#                       chr  start  end  rec_rate
#                     Generate from Glemin et al. ms-rec-triticeae TXT files with:
#                       python prepare_recmap.py --input recmap.txt --out data/species_recmap.bed
#                     Run with --preview first to check column names.

SPECIES_LIST <- list(

  Aegilopsspeltoides = list(
    output_dir      = "results/speltoides",
    gene_coords_bed = "data/speltoides_gene_coords.bed",
    rec_map_bed     = "data/Ae_speltoides_genomeB_1cM.bed"
  )

  # Uncomment to add more species:
  # ,
  # Aegilopsmutica = list(
  #   output_dir      = "results/mutica",
  #   gene_coords_bed = "data/mutica_gene_coords.bed",
  #   rec_map_bed     = "data/Ae_mutica_recombination.bed"
  # ),
  # Aegilopstauschii = list(
  #   output_dir      = "results/tauschii",
  #   gene_coords_bed = "data/tauschii_gene_coords.bed",
  #   rec_map_bed     = "data/Ae_tauschii_recombination.bed"
  # ),
  # Triticumurartu = list(
  #   output_dir      = "results/urartu",
  #   gene_coords_bed = "data/urartu_gene_coords.bed",
  #   rec_map_bed     = "data/T_urartu_recombination.bed"
  # )

)


# =============================================================================
# SECTION 2 — Analysis parameters
# =============================================================================

# Which candidate set defines "selected" (binary response in the GLM):
#   "bs"  = branch-specific candidates (pos_branch_specific_MKT or
#            pos_branch_specific_impMKT) — recommended for this analysis
#   "any" = positive in any of the 4 MKT tests
CANDIDATE_TYPE <- "bs"

# MKT filter: genes must pass this filter to enter the recombination analysis
# (same filter used in branch_specific_MKT_analysis.Rmd for branch-specific tests)
MIN_POLY  <- 5   # minimum Pn + Ps
MIN_COUNT <- 1   # minimum for each of Pn, Ps, Dn, Ds (branch-specific)

# Sensitivity analysis: exclude genes with |distance_z| > this threshold
# z = 3 removes the top ~0.3% of genes with very large assignment distances
DISTANCE_Z_THRESHOLD <- 3

# Scaling factor for rec_rate in the GLM
# Default 1e6 converts from per-bp to per-Mb, giving coefficients of reasonable scale
REC_SCALE <- 1e6

# Output directory for the HTML report and result tables
OUTPUT_DIR <- "results/recombination"


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
.rmd_file <- file.path(.rmd_dir, "recombinaison.Rmd")
if (!file.exists(.rmd_file))
  stop("Rmd file not found: ", .rmd_file)


# ── Validate inputs ───────────────────────────────────────────────────────────

if (length(SPECIES_LIST) == 0)
  stop("SPECIES_LIST is empty. Add at least one species in Section 1.")

for (.sp in names(SPECIES_LIST)) {
  .cfg <- SPECIES_LIST[[.sp]]

  # MKT result files (from render_analysis.R)
  .mkt_files <- c(
    file.path(.cfg$output_dir, "branch_specific_MKT_results.tsv"),
    file.path(.cfg$output_dir,
              if (CANDIDATE_TYPE == "bs") "candidates_branch_specific.tsv"
              else                         "candidates_all.tsv")
  )
  for (.f in .mkt_files)
    if (!file.exists(.f))
      warning(sprintf("[%s] MKT result file not found: %s\n",
                      .sp, .f),
              "  → Run render_analysis.R first.")

  # User-provided BED files
  for (.key in c("gene_coords_bed", "rec_map_bed")) {
    .f <- .cfg[[.key]]
    if (is.null(.f) || !nzchar(.f))
      stop(sprintf("[%s] '%s' is missing from SPECIES_LIST.", .sp, .key))
    if (!file.exists(.f))
      stop(sprintf("[%s] File not found for '%s': %s", .sp, .key, .f))
  }
}


# ── Build species_configs param ───────────────────────────────────────────────

.species_configs <- setNames(
  lapply(names(SPECIES_LIST), function(.sp) {
    .cfg <- SPECIES_LIST[[.sp]]
    list(
      output_dir      = normalizePath(.cfg$output_dir,      mustWork = FALSE),
      gene_coords_bed = normalizePath(.cfg$gene_coords_bed, mustWork = FALSE),
      rec_map_bed     = normalizePath(.cfg$rec_map_bed,     mustWork = FALSE)
    )
  }),
  names(SPECIES_LIST)
)


# ── Create output directory and render ───────────────────────────────────────

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
.out_html <- normalizePath(file.path(OUTPUT_DIR, "recombination.html"),
                            mustWork = FALSE)

message("=== Recombination Rate Analysis ===")
message(sprintf("Species         : %s", paste(names(SPECIES_LIST), collapse = ", ")))
message(sprintf("Candidate type  : %s", CANDIDATE_TYPE))
message(sprintf("Output          : %s", .out_html))

rmarkdown::render(
  input       = .rmd_file,
  output_file = .out_html,
  params = list(
    species_configs      = .species_configs,
    candidate_type       = CANDIDATE_TYPE,
    min_poly             = MIN_POLY,
    min_count            = MIN_COUNT,
    distance_z_threshold = DISTANCE_Z_THRESHOLD,
    rec_scale            = REC_SCALE,
    output_dir           = normalizePath(OUTPUT_DIR, mustWork = FALSE)
  ),
  envir = new.env(parent = globalenv())
)

message("Done. Report: ", .out_html)
