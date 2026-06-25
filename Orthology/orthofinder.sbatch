#!/bin/bash
# ==============================================================================
# CLUSTER CONFIGURATION — read this section before running on a new cluster
# ==============================================================================
#
# OrthoFinder is a comprehensive platform for comparative genomics that finds
# orthogroups and orthologs, infers rooted gene trees and a rooted species
# tree, and identifies all gene duplication events.
#
# IMPORTANT — INPUT FORMAT
# -------------------------
# OrthoFinder requires ONE PROTEIN FASTA FILE PER SPECIES, all placed in a
# single directory.  The filename (without extension) becomes the species name
# in the output.
#
#   proteomes/
#   ├── Ae_speltoides.fasta     ← species 1
#   ├── Ae_tauschii.fasta       ← species 2
#   ├── Ae_urartu.fasta         ← species 3
#   └── T_aestivum.fasta        ← species 4
#
# Do NOT combine all species into a single file — OrthoFinder needs to know
# which sequences belong to which species.
#
# Lines marked [CLUSTER] MUST be adapted before submission.
#
# Resource guidance:
#   --cpus-per-task : passed to OrthoFinder as -t (DIAMOND search threads)
#                     and -a (analysis threads). 32 CPUs is a good default.
#   --mem           : ~2 GB per species is a safe estimate for typical
#                     transcriptome-scale proteomes (5,000–30,000 proteins).
#                     Increase for whole-genome proteomes.
#   --time          : highly variable (minutes to hours depending on
#                     number of species and proteome size).
#
# Cluster examples:
#   IFB:   --partition=cpu-dedicated  --account=<project>@cpu
#          module load bioinfo-ifb
#          module load orthofinder/3.0.1b1   (if available as a module)
#
#   CIRAD: --partition=cpu-dedicated  --account=dedicated-cpu@cirad-normal
#          module load bioinfo-cirad
#          module load orthofinder/3.0.1b1
#
#   Conda installation (recommended — see INSTALLATION section in README.md):
#          conda activate orthofinder_env
#          (remove or comment out the module load lines below)
#
# INSTALLATION (if not available as a module)
# --------------------------------------------
# Option A — conda (recommended):
#   conda create -n orthofinder_env -c bioconda orthofinder diamond
#   conda activate orthofinder_env
#
# Option B — manual (OrthoFinder 3.x):
#   wget https://github.com/davidemms/OrthoFinder/releases/download/3.0.1b1/OrthoFinder.tar.gz
#   tar xzf OrthoFinder.tar.gz
#   pip install scipy biopython
#   # Add to PATH or call with full path
#
# ==============================================================================

#SBATCH --job-name=orthofinder
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --account=dedicated-cpu@cirad-normal   # [CLUSTER] change to your account
#SBATCH --partition=cpu-dedicated              # [CLUSTER] change to your partition

# ==============================================================================
# MODULE LOADING — adapt to your cluster
# ==============================================================================
module load bioinfo-ifb          # [CLUSTER] change or remove if using conda
# module load orthofinder/3.0.1b1  # [CLUSTER] uncomment if available as module
# conda activate orthofinder_env   # [CLUSTER] uncomment if using conda
# ==============================================================================

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-32}"

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
Usage: sbatch orthofinder.sbatch [OPTIONS]

Required:
  --proteomes DIR   Directory containing one protein FASTA per species.
                    File names (without extension) are used as species names.
                    Example: Ae_speltoides.fasta → species "Ae_speltoides"

  --outdir    DIR   Output directory for OrthoFinder results.

Optional:
  --search    TOOL  Search tool: diamond (default), diamond_ultra_sens,
                    blast, mmseqs2_default
  --threads   N     Number of threads (default: SLURM_CPUS_PER_TASK or 32)
  --help            Show this message

IMPORTANT — directory structure:
  proteomes/
  ├── species1.fasta    ← one file per species
  ├── species2.fasta
  └── species3.fasta

Do NOT put all sequences in a single file.

Outputs (in --outdir):
  Orthogroups/
    Orthogroups.tsv                  — genes grouped into orthogroups
    Orthogroups_SingleCopyOrthologues.txt  — strict single-copy OGs
  Orthologues/
    Orthologues_<species>/           — pairwise ortholog tables
  Species_Tree/
    SpeciesTree_rooted.txt           — rooted species tree (Newick)
  Gene_Trees/                        — one rooted gene tree per orthogroup

Examples:
  sbatch orthofinder.sbatch \
      --proteomes proteomes/ \
      --outdir    results/orthofinder/

  # Re-run from existing DIAMOND results (skip search step):
  # Use -b to point to a previous OrthoFinder results directory
EOF
}

# ── Argument parsing ──────────────────────────────────────────────────────────
PROTEOMES_DIR=""
OUTPUT_DIR=""
SEARCH_TOOL="diamond"

[ $# -eq 0 ] && { show_help; exit 0; }

while [ $# -gt 0 ]; do
    case "$1" in
        --proteomes) PROTEOMES_DIR="$2"; shift 2 ;;
        --outdir)    OUTPUT_DIR="$2";    shift 2 ;;
        --search)    SEARCH_TOOL="$2";   shift 2 ;;
        --threads)   THREADS="$2";       shift 2 ;;
        --help|-h)   show_help; exit 0 ;;
        *) echo "ERROR: Unknown argument: '$1'" >&2; show_help >&2; exit 1 ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$PROTEOMES_DIR" ] || die "--proteomes is required"
[ -n "$OUTPUT_DIR" ]    || die "--outdir is required"
[ -d "$PROTEOMES_DIR" ] || die "Proteomes directory not found: $PROTEOMES_DIR"

N_SPECIES=$(find "$PROTEOMES_DIR" -maxdepth 1 \
    \( -name "*.fa" -o -name "*.fasta" -o -name "*.faa" \) | wc -l)
[ "$N_SPECIES" -ge 2 ] || \
    die "Need at least 2 protein FASTA files in $PROTEOMES_DIR, found $N_SPECIES"

command -v orthofinder &>/dev/null || \
    die "OrthoFinder not found. Check module loading, conda env, or PATH."

mkdir -p "$OUTPUT_DIR"

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  OrthoFinder"
echo "======================================================"
echo "  Job ID      : ${SLURM_JOB_ID:-N/A}"
echo "  Start       : $(date)"
echo "  Proteomes   : $PROTEOMES_DIR  ($N_SPECIES species)"
echo "  Output      : $OUTPUT_DIR"
echo "  Search      : $SEARCH_TOOL"
echo "  Threads     : $THREADS"
echo "======================================================"
echo ""

# List input files for reproducibility
echo "Species detected:"
find "$PROTEOMES_DIR" -maxdepth 1 \
    \( -name "*.fa" -o -name "*.fasta" -o -name "*.faa" \) | sort | \
    while read -r f; do
        n_seqs=$(grep -c "^>" "$f" || true)
        printf "  %-40s  %d sequences\n" "$(basename "$f")" "$n_seqs"
    done
echo ""

# ── Run OrthoFinder ───────────────────────────────────────────────────────────
orthofinder \
    -f  "$PROTEOMES_DIR" \
    -o  "$OUTPUT_DIR" \
    -S  "$SEARCH_TOOL" \
    -t  "$THREADS" \
    -a  "$THREADS"

echo ""
echo "======================================================"
echo "  Done  ($(date))"
echo ""
echo "  Key output files:"
echo "    $OUTPUT_DIR/Orthogroups/Orthogroups.tsv"
echo "    $OUTPUT_DIR/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"
echo "    $OUTPUT_DIR/Species_Tree/SpeciesTree_rooted.txt"
echo "    $OUTPUT_DIR/Orthologues/"
echo "======================================================"
