#!/bin/bash
# ==============================================================================
# CLUSTER CONFIGURATION — read this section before running on a new cluster
# ==============================================================================
#
# macse_enrichment.sbatch — SLURM job for MACSE enrichAlignment
#
# Adds per-species CDS polymorphism sequences to multi-species orthogroup
# nucleotide alignments using MACSE enrichAlignment.
#
# This script submits a single SLURM job that runs MACSE_enrichment.py, which
# internally parallelises across genes using Python threads (--threads).
# There is no SLURM array because all genes are managed by one Python process.
#
# PRE-REQUISITE — build the orthogroup table (run once on the login node):
#
#   python orthogroup_table.py \
#       --input-dir  <nt_alignment_dir>   \
#       --output     orthogroup_table.tsv
#
#   python MACSE_enrichment.py \
#       --ortho-table orthogroup_table.tsv \
#       --list-species                       # shows available species column names
#
# THEN — edit the USER SETTINGS section below and submit:
#
#   sbatch macse_enrichment.sbatch
#
# Lines marked [CLUSTER] MUST be adapted before the first submission.
#
# Resource guidance:
#   --cpus-per-task : set to the same value as THREADS below
#   --mem           : MACSE uses ~1–2 GB per running instance;
#                     set to THREADS × 2G  (e.g. 8 threads → 16G)
#   --time          : allow generous wall time; large datasets can run many hours
#
# ==============================================================================

#SBATCH --job-name=macse_enrich
#SBATCH --output=macse_enrich_%j.out
#SBATCH --error=macse_enrich_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8              # [CLUSTER] match THREADS below
#SBATCH --mem=16G                      # [CLUSTER] THREADS × 2G recommended
#SBATCH --time=24:00:00               # [CLUSTER] adjust per dataset size
#SBATCH --partition=cpu-dedicated     # [CLUSTER] set to your cluster partition
##SBATCH --account=myaccount          # [CLUSTER] uncomment and set if required

# ==============================================================================
# MODULE LOADING — adapt to your cluster
# ==============================================================================
module purge
# module load bioinfo-ifb              # [CLUSTER] uncomment your meta-module
# module load python/3.9              # [CLUSTER] uncomment if Python needs loading
# module load java-jdk/17             # [CLUSTER] uncomment if Java needs loading
# conda activate macse_env            # [CLUSTER] uncomment if using a conda env
# ==============================================================================

set -euo pipefail

# ==============================================================================
# USER SETTINGS — edit these paths before submitting
# ==============================================================================

# Orthogroup-to-gene-ID table (output of orthogroup_table.py)
ORTHO_TABLE="/path/to/orthogroup_table.tsv"

# Directory containing per-HOG nucleotide alignment FASTAs
# (e.g. output of pal2nal.sbatch: HOG0000001_aln_hmm_nt.fasta)
ORTHO_DIR="/path/to/nt_alignments/"

# Directory containing per-gene CDS polymorphism FASTAs for the target species
# (e.g. *_CDS_aligned_NT.fasta from an earlier MACSE alignment run)
POLY_DIR="/path/to/cds_polymorphism/"

# Path to the MACSE JAR file
MACSE_JAR="/path/to/macse_v2.07.jar"

# Column name in ORTHO_TABLE for the target species
# Run:  python MACSE_enrichment.py --ortho-table "$ORTHO_TABLE" --list-species
# to see all available column names.
SPECIES_COL="Aegilops_speltoides"

# Output directory (sub-directories enriched/ and logs/ are created automatically)
OUTPUT_DIR="/path/to/macse_enriched/"

# Number of parallel MACSE workers (must equal --cpus-per-task above)
THREADS=8

# ==============================================================================
# VALIDATION — do not edit below this line
# ==============================================================================

# Locate MACSE_enrichment.py in the same directory as this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACSE_PY="${SCRIPT_DIR}/MACSE_enrichment.py"

show_help() {
cat <<'EOF'
macse_enrichment.sbatch — SLURM job for MACSE enrichAlignment

WORKFLOW
  Step 1 — Build the orthogroup table (run once on the login node):
    python orthogroup_table.py \
        --input-dir  <nt_alignment_dir> \
        --output     orthogroup_table.tsv

  Step 2 — List available species columns:
    python MACSE_enrichment.py \
        --ortho-table orthogroup_table.tsv \
        --list-species

  Step 3 — Edit USER SETTINGS in this file, then submit:
    sbatch macse_enrichment.sbatch

  Step 4 — Dry run (logs commands without executing):
    bash macse_enrichment.sbatch --dry-run

  Step 5 — Resume after partial run:
    sbatch macse_enrichment.sbatch       # already-done genes are skipped

OUTPUTS
  $OUTPUT_DIR/enriched/              — enriched alignments (*_NT.fasta, *_AA.fasta)
  $OUTPUT_DIR/logs/enrichment_summary.tsv   — per-gene result table
  $OUTPUT_DIR/logs/failed_alignments.txt    — genes that failed
  $OUTPUT_DIR/logs/rerun_failed.sh          — re-run script for failed genes
  $OUTPUT_DIR/logs/macse_commands.log       — all MACSE commands issued
EOF
}

# Parse optional --dry-run flag from command line
DRY_RUN_FLAG=""
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN_FLAG="--dry-run" ;;
        -h|--help) show_help; exit 0 ;;
    esac
done

# ── Banner ─────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  MACSE enrichAlignment — polymorphism integration"
echo "======================================================"
echo "  Job ID        : ${SLURM_JOB_ID:-N/A (local run)}"
echo "  Start         : $(date)"
echo "  Ortho table   : $ORTHO_TABLE"
echo "  Ortho dir     : $ORTHO_DIR"
echo "  Poly dir      : $POLY_DIR"
echo "  MACSE jar     : $MACSE_JAR"
echo "  Species col   : $SPECIES_COL"
echo "  Output dir    : $OUTPUT_DIR"
echo "  Threads       : $THREADS"
[ -n "$DRY_RUN_FLAG" ] && echo "  [DRY RUN — commands logged, not executed]"
echo "======================================================"
echo ""

# ── Pre-flight checks ──────────────────────────────────────────────────────────
ERRORS=0

[ -f "$ORTHO_TABLE" ] || { echo "ERROR: ORTHO_TABLE not found: $ORTHO_TABLE"; ERRORS=1; }
[ -d "$ORTHO_DIR"   ] || { echo "ERROR: ORTHO_DIR not found: $ORTHO_DIR";     ERRORS=1; }
[ -d "$POLY_DIR"    ] || { echo "ERROR: POLY_DIR not found: $POLY_DIR";       ERRORS=1; }
[ -f "$MACSE_JAR"   ] || { echo "ERROR: MACSE_JAR not found: $MACSE_JAR";    ERRORS=1; }
[ -f "$MACSE_PY"    ] || { echo "ERROR: MACSE_enrichment.py not found: $MACSE_PY"; ERRORS=1; }

command -v java    &>/dev/null || { echo "ERROR: java not found in PATH";    ERRORS=1; }
command -v python3 &>/dev/null || { echo "ERROR: python3 not found in PATH"; ERRORS=1; }

if [ "$ERRORS" -gt 0 ]; then
    echo ""
    echo "Fix the errors above and resubmit."
    exit 1
fi

# Report input counts
N_ORTHO=$(ls -1 "$ORTHO_DIR"/*HOG*.fasta 2>/dev/null | wc -l || echo 0)
N_POLY=$(ls  -1 "$POLY_DIR"/*.fasta      2>/dev/null | wc -l || echo 0)
echo "  Ortho alignment files : $N_ORTHO"
echo "  Poly CDS files        : $N_POLY"
echo ""

# ── Run MACSE_enrichment.py ────────────────────────────────────────────────────
echo "Launching MACSE_enrichment.py..."
python3 "$MACSE_PY" \
    --ortho-table  "$ORTHO_TABLE" \
    --ortho-dir    "$ORTHO_DIR"   \
    --poly-dir     "$POLY_DIR"    \
    --output-dir   "$OUTPUT_DIR"  \
    --macse-jar    "$MACSE_JAR"   \
    --species-col  "$SPECIES_COL" \
    --threads      "$THREADS"     \
    $DRY_RUN_FLAG

EXIT_CODE=$?

echo ""
echo "======================================================"
echo "  Done  ($(date))"
echo "  Python exit code : $EXIT_CODE"
echo "  Results in       : $OUTPUT_DIR"
echo "======================================================"
exit "$EXIT_CODE"
