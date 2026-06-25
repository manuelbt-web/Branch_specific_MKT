# =============================================================================
# render_analysis.R — Multi-species Branch-Specific MKT Launcher
# =============================================================================
#
# PURPOSE
#   Run branch_specific_MKT_analysis.Rmd for one or several focal species.
#   When multiple species are provided, combined cross-species candidate tables
#   are built automatically after all per-species analyses are done.
#
# HOW TO USE
#   1. Edit SECTION 1 — provide one entry in SPECIES_LIST per focal species.
#   2. Edit SECTION 2 — set shared paths (orthogroup table, output root).
#   3. Run:
#        a) RStudio  : open this file → click Source
#        b) R console: source("path/to/render_analysis.R")
#        c) Terminal : Rscript path/to/render_analysis.R
#
# OUTPUT LAYOUT (example with 2 species)
#   results/
#   ├── speltoides/
#   │   ├── branch_specific_MKT_analysis.html        ← per-species report
#   │   ├── branch_specific_MKT_results.tsv          ← all genes + 4 MKT tests
#   │   ├── candidates_all.tsv                       ← positive in ≥1 test
#   │   └── candidates_branch_specific.tsv           ← positive in branch-sp. test
#   ├── mutica/
#   │   └── (same structure)
#   └── combined/                          ← created only when ≥2 species
#       ├── candidates_all_combined.tsv              ← all species, all tests
#       ├── candidates_branch_specific_combined.tsv  ← branch-specific only
#       ├── candidates_shared_combined.tsv           ← HOGs in ≥2 focal species
#       └── candidates_shared_branch_specific_combined.tsv
#
# NOTES
#   • File paths can be absolute or relative to your working directory when
#     you run this script — NOT relative to the .Rmd location.
#   • Output directories are created automatically.
#   • The .Rmd file is found automatically next to this render_analysis.R.
# =============================================================================


# ── 0.  Dependency check ─────────────────────────────────────────────────────

required_pkgs <- c("rmarkdown", "tidyverse", "kableExtra",
                   "patchwork", "UpSetR", "DT")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace,
                                        quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")


# =============================================================================
# SECTION 1 — Define focal species
# =============================================================================
#
# Add one entry per focal species.
#
# Each entry is a named list:
#   name                         : species label (as it appears in the
#                                  'species' column of the merged table,
#                                  e.g. "Aegilopsspeltoides")
#   merged_table                 : path to merge_mkt_results.py output  [REQUIRED]
#   output_dir                   : directory for per-species outputs     [REQUIRED]
#   focal_species_in_ortho_table : species name EXACTLY as in the
#                                  orthogroup table header column
#                                  (e.g. "Aegilops_speltoides").
#                                  Leave "" for automatic detection.
#
# For a single species: provide one entry and the combined step is skipped.
# For multiple species: add all entries; combined tables are built automatically.

SPECIES_LIST <- list(

  Aegilopsspeltoides = list(
    merged_table                 = "results/speltoides/merged_mkt_results.tsv",
    output_dir                   = "results/speltoides",
    focal_species_in_ortho_table = "Aegilops_speltoides"
  )

  # Uncomment and fill in to add more species:
  # ,
  # Aegilopsmutica = list(
  #   merged_table                 = "results/mutica/merged_mkt_results.tsv",
  #   output_dir                   = "results/mutica",
  #   focal_species_in_ortho_table = "Aegilops_mutica"
  # ),
  # Aegilopstauschii = list(
  #   merged_table                 = "results/tauschii/merged_mkt_results.tsv",
  #   output_dir                   = "results/tauschii",
  #   focal_species_in_ortho_table = "Aegilops_tauschii"
  # ),
  # Triticumurartu = list(
  #   merged_table                 = "results/urartu/merged_mkt_results.tsv",
  #   output_dir                   = "results/urartu",
  #   focal_species_in_ortho_table = "Triticum_urartu"
  # )

)


# =============================================================================
# SECTION 2 — Shared settings
# =============================================================================

# Orthogroup table from orthogroup_table.py
# Columns: orthogroup | Aegilops_speltoides | Aegilops_mutica | ...
ORTHOGROUP_TABLE <- "data/orthogroup_table.tsv"

# Combined cross-species tables go here (created only when ≥2 species)
COMBINED_DIR <- "results/combined"


# =============================================================================
# ENGINE — do not edit below this line
# =============================================================================

# ── Locate Rmd file ───────────────────────────────────────────────────────────

.rmd_dir <- tryCatch(
  dirname(normalizePath(
    c(sys.frame(1)$ofile,
      sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
    )[1]
  )),
  error = function(e) {
    warning("Could not auto-detect script directory; using getwd(). ",
            "Ensure render_analysis.R is in the same directory as the .Rmd.")
    getwd()
  }
)
.rmd_file <- file.path(.rmd_dir, "branch_specific_MKT_analysis.Rmd")
if (!file.exists(.rmd_file))
  stop("Rmd file not found: ", .rmd_file,
       "\nEnsure render_analysis.R and branch_specific_MKT_analysis.Rmd ",
       "are in the same directory.")


# ── Validate SPECIES_LIST ─────────────────────────────────────────────────────

if (length(SPECIES_LIST) == 0)
  stop("SPECIES_LIST is empty. Add at least one species in Section 1.")

for (.sp in names(SPECIES_LIST)) {
  .cfg <- SPECIES_LIST[[.sp]]
  if (is.null(.cfg$merged_table))
    stop("'merged_table' missing for species '", .sp, "'")
  if (!file.exists(.cfg$merged_table))
    stop("merged_table not found for '", .sp, "': ", .cfg$merged_table)
  if (is.null(.cfg$output_dir))
    stop("'output_dir' missing for species '", .sp, "'")
}

.n_species <- length(SPECIES_LIST)
message(sprintf("\n=== Branch-specific MKT: %d focal species ===", .n_species))
for (.sp in names(SPECIES_LIST))
  message(sprintf("  • %s → %s", .sp, SPECIES_LIST[[.sp]]$output_dir))


# ── Pass 1: per-species analysis ──────────────────────────────────────────────

.per_species <- list()  # collect output paths for the combined step

for (.sp in names(SPECIES_LIST)) {
  .cfg     <- SPECIES_LIST[[.sp]]
  .out_dir <- .cfg$output_dir

  dir.create(.out_dir, recursive = TRUE, showWarnings = FALSE)

  .out_tsv  <- file.path(.out_dir, "branch_specific_MKT_results.tsv")
  .out_html <- file.path(.out_dir, "branch_specific_MKT_analysis.html")
  .cand_all <- file.path(.out_dir, "candidates_all.tsv")
  .cand_bs  <- file.path(.out_dir, "candidates_branch_specific.tsv")
  .cand_sh  <- file.path(.out_dir, "candidates_shared.tsv")

  message(sprintf("\n[%s] Rendering...", .sp))

  rmarkdown::render(
    input       = .rmd_file,
    output_file = normalizePath(.out_html, mustWork = FALSE),
    params = list(
      merged_table                 = normalizePath(.cfg$merged_table, mustWork = FALSE),
      orthogroup_table             = normalizePath(ORTHOGROUP_TABLE,  mustWork = FALSE),
      output_table                 = normalizePath(.out_tsv,  mustWork = FALSE),
      candidates_table_all         = normalizePath(.cand_all, mustWork = FALSE),
      candidates_table_bs          = normalizePath(.cand_bs,  mustWork = FALSE),
      candidates_shared            = normalizePath(.cand_sh,  mustWork = FALSE),
      species                      = .sp,
      focal_species_in_ortho_table = .cfg$focal_species_in_ortho_table %||% "",
      other_results                = list()   # cross-species done in combined step
    ),
    envir = new.env(parent = globalenv())
  )

  .per_species[[.sp]] <- list(
    out_tsv   = .out_tsv,
    out_html  = .out_html,
    cand_all  = .cand_all,
    cand_bs   = .cand_bs,
    cand_sh   = .cand_sh,
    out_dir   = .out_dir
  )

  message(sprintf("[%s] Done — report: %s", .sp, .out_html))
}


# ── Combined cross-species tables (≥2 species only) ───────────────────────────

if (.n_species >= 2) {

  message("\n=== Building combined cross-species tables ===")
  dir.create(COMBINED_DIR, recursive = TRUE, showWarnings = FALSE)

  .load_tsv <- function(path, sp_label) {
    if (!file.exists(path)) {
      warning("File not found, skipping: ", path)
      return(NULL)
    }
    readr::read_tsv(
      path,
      col_types  = readr::cols(.default = readr::col_character()),
      na         = c("NA", "NaN", ""),
      show_col_types = FALSE
    )
  }

  # ── All candidates combined ─────────────────────────────────────────────────
  .all_rows <- lapply(names(.per_species), function(.sp)
    .load_tsv(.per_species[[.sp]]$cand_all, .sp)
  )
  .all_cands <- dplyr::bind_rows(Filter(Negate(is.null), .all_rows))

  .f_all <- file.path(COMBINED_DIR, "candidates_all_combined.tsv")
  readr::write_tsv(.all_cands, .f_all, na = "NA")
  message(sprintf("  candidates_all_combined.tsv      — %d rows across %d species",
                  nrow(.all_cands), .n_species))

  # ── Shared candidates (HOG detected in ≥2 focal species) ───────────────────
  if ("ortholog" %in% colnames(.all_cands)) {
    .shared <- .all_cands %>%
      dplyr::filter(!is.na(ortholog) & ortholog != "") %>%
      dplyr::group_by(ortholog) %>%
      dplyr::filter(dplyr::n_distinct(focal_species) >= 2) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(ortholog, focal_species)

    .f_shared <- file.path(COMBINED_DIR, "candidates_shared_combined.tsv")
    readr::write_tsv(.shared, .f_shared, na = "NA")
    message(sprintf("  candidates_shared_combined.tsv   — %d HOGs in ≥2 species",
                    dplyr::n_distinct(.shared$ortholog)))
  } else {
    warning("'ortholog' column not found in candidate tables — shared table not created.")
  }

  # ── Branch-specific candidates combined ────────────────────────────────────
  .bs_rows <- lapply(names(.per_species), function(.sp)
    .load_tsv(.per_species[[.sp]]$cand_bs, .sp)
  )
  .bs_cands <- dplyr::bind_rows(Filter(Negate(is.null), .bs_rows))

  .f_bs <- file.path(COMBINED_DIR, "candidates_branch_specific_combined.tsv")
  readr::write_tsv(.bs_cands, .f_bs, na = "NA")
  message(sprintf("  candidates_branch_specific_combined.tsv — %d rows", nrow(.bs_cands)))

  # ── Shared branch-specific candidates ──────────────────────────────────────
  if ("ortholog" %in% colnames(.bs_cands)) {
    .shared_bs <- .bs_cands %>%
      dplyr::filter(!is.na(ortholog) & ortholog != "") %>%
      dplyr::group_by(ortholog) %>%
      dplyr::filter(dplyr::n_distinct(focal_species) >= 2) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(ortholog, focal_species)

    .f_shared_bs <- file.path(COMBINED_DIR,
                               "candidates_shared_branch_specific_combined.tsv")
    readr::write_tsv(.shared_bs, .f_shared_bs, na = "NA")
    message(sprintf(
      "  candidates_shared_branch_specific_combined.tsv — %d HOGs in ≥2 species",
      dplyr::n_distinct(.shared_bs$ortholog)
    ))
  }

  message(sprintf("\nCombined tables written to: %s/", COMBINED_DIR))
}


# ── Summary ───────────────────────────────────────────────────────────────────

message("\n=== All done ===")
message("Per-species reports:")
for (.sp in names(.per_species))
  message(sprintf("  %s → %s", .sp, .per_species[[.sp]]$out_html))

if (.n_species >= 2)
  message(sprintf("Combined tables → %s/", COMBINED_DIR))
