#!/usr/bin/env Rscript
# Full run — Ae. speltoides branch-specific MKT (3560 genes)
#
# HOW TO RUN (from the Scripts/ directory):
#   Rscript test/speltoides/run_report.R
#
# The orthogroup table is auto-generated from the split-step divergence FASTAs
# (02_split/divergence/) via orthogroup_table.py — no manual pre-computation needed.

# ── Resolve base directory ─────────────────────────────────────────────────────
.this_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
BASE <- if (length(.this_file) == 1L && nchar(.this_file) > 0) {
  normalizePath(file.path(dirname(.this_file), "../.."), mustWork = FALSE)
} else {
  "C:/Users/barrientos/Documents/data/MKT_2024/article/Scripts"
}

# ── Settings ──────────────────────────────────────────────────────────────────
SPECIES_LIST <- list(
  Aegilopsspeltoides = list(
    merged_table                 = file.path(BASE, "test/speltoides/05_merged/merged_mkt_results.tsv"),
    output_dir                   = file.path(BASE, "test/speltoides/06_results"),
    focal_species_in_ortho_table = "Aegilops_speltoides"
  )
)

SPLIT_DIVERGENCE_DIR <- file.path(BASE, "test/speltoides/02_split/divergence")
ORTHOGROUP_TABLE     <- file.path(BASE, "test/speltoides/02_split/orthogroup_table.tsv")
COMBINED_DIR         <- file.path(BASE, "test/speltoides/06_results/combined")

.rmd_dir  <- file.path(BASE, "MKT_3_species/branch_specific_MKT")
.rmd_file <- file.path(.rmd_dir, "branch_specific_MKT_analysis.Rmd")
if (!file.exists(.rmd_file))
  stop("Rmd not found: ", .rmd_file)

# ── Dependency check ──────────────────────────────────────────────────────────
required_pkgs <- c("rmarkdown", "tidyverse", "kableExtra", "patchwork", "UpSetR", "DT")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with:  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")

# ── Auto-generate orthogroup table from split divergence FASTAs ───────────────
if (nchar(SPLIT_DIVERGENCE_DIR) > 0) {
  .ortho_py <- normalizePath(
    file.path(.rmd_dir, "..", "EnrichAlignment", "orthogroup_table.py"),
    mustWork = FALSE
  )
  if (!file.exists(.ortho_py))
    stop("orthogroup_table.py not found: ", .ortho_py)
  message(sprintf("\nGenerating orthogroup table from: %s", SPLIT_DIVERGENCE_DIR))
  .ret <- system2(
    "python3",
    args   = c(shQuote(.ortho_py),
               "--input-dir", shQuote(SPLIT_DIVERGENCE_DIR),
               "--output",    shQuote(normalizePath(ORTHOGROUP_TABLE, mustWork = FALSE))),
    stdout = TRUE, stderr = TRUE
  )
  message(paste(.ret, collapse = "\n"))
  if (!file.exists(ORTHOGROUP_TABLE))
    stop("orthogroup_table.py failed — output not found: ", ORTHOGROUP_TABLE)
}
if (!file.exists(ORTHOGROUP_TABLE))
  stop("Orthogroup table not found: ", ORTHOGROUP_TABLE)

# ── Render per species ────────────────────────────────────────────────────────
.n_species <- length(SPECIES_LIST)
message(sprintf("\n=== Branch-specific MKT: %d focal species ===", .n_species))
for (.sp in names(SPECIES_LIST))
  message(sprintf("  • %s → %s", .sp, SPECIES_LIST[[.sp]]$output_dir))

for (.sp in names(SPECIES_LIST)) {
  .cfg     <- SPECIES_LIST[[.sp]]
  .out_dir <- .cfg$output_dir
  dir.create(.out_dir, recursive = TRUE, showWarnings = FALSE)

  .out_html <- file.path(.out_dir, "branch_specific_MKT_analysis.html")
  .out_tsv  <- file.path(.out_dir, "branch_specific_MKT_results.tsv")
  .cand_all <- file.path(.out_dir, "candidates_all.tsv")
  .cand_bs  <- file.path(.out_dir, "candidates_branch_specific.tsv")
  .cand_sh  <- file.path(.out_dir, "candidates_shared.tsv")

  message(sprintf("\n[%s] Rendering...", .sp))
  rmarkdown::render(
    input       = .rmd_file,
    output_file = normalizePath(.out_html, mustWork = FALSE),
    params = list(
      merged_table                 = normalizePath(.cfg$merged_table,  mustWork = FALSE),
      orthogroup_table             = normalizePath(ORTHOGROUP_TABLE,   mustWork = FALSE),
      output_table                 = normalizePath(.out_tsv,  mustWork = FALSE),
      candidates_table_all         = normalizePath(.cand_all, mustWork = FALSE),
      candidates_table_bs          = normalizePath(.cand_bs,  mustWork = FALSE),
      candidates_shared            = normalizePath(.cand_sh,  mustWork = FALSE),
      species                      = .sp,
      focal_species_in_ortho_table = .cfg$focal_species_in_ortho_table,
      other_results                = list()
    ),
    envir = new.env(parent = globalenv())
  )
  message(sprintf("[%s] Done — %s", .sp, .out_html))
}

message("\n=== All done ===")
