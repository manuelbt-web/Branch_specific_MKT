#!/bin/bash
# ==============================================================================
# CLUSTER CONFIGURATION — read this section before running on a new cluster
# ==============================================================================
#
# mafft_alignments.sbatch — SLURM array job for multiple-sequence alignment
# of OrthoFinder single-copy orthologs using MAFFT.
#
# INPUT  : Single_Copy_Orthologue_Sequences/ directory from OrthoFinder
#          (one .fa file per single-copy orthogroup; protein sequences)
# OUTPUT : one aligned FASTA per orthogroup  (*_aln.fa)
#
# Lines marked [CLUSTER] MUST be adapted before the first submission.
#
# Resource guidance:
#   --cpus-per-task : passed to MAFFT as --thread; 4-8 CPUs recommended
#                     (MAFFT parallelises over columns, benefit saturates ~8)
#   --mem           : 2–4 GB per job is sufficient for typical orthogroups
#   --time          : 2–4 h is enough for most transcriptome-scale datasets
#
# Cluster examples:
#   IFB:   --partition=cpu-dedicated  --account=<project>@cpu
#          module load bioinfo-ifb; module load mafft/7.515
#
#   CIRAD: --partition=cpu-dedicated  --account=dedicated-cpu@cirad-normal
#          module load bioinfo-cirad; module load mafft/7.515
#
#   Conda (recommended if module not available):
#          conda create -n mafft_env -c bioconda mafft
#          conda activate mafft_env
#          (remove or comment out the module load lines below)
#
# SUBMISSION EXAMPLE:
#   N=$(ls -1 <input_dir>/*.fa | wc -l)
#   sbatch --array=0-$((N-1))%15 mafft_alignments.sbatch \
#          --input-dir  <input_dir>  \
#          --output-dir <output_dir>
#
# ==============================================================================

#SBATCH --job-name=mafft_align
#SBATCH --output=mafft_%A_%a.out
#SBATCH --error=mafft_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4           # [CLUSTER] adjust; passed to MAFFT --thread
#SBATCH --mem=4G                    # [CLUSTER] adjust per dataset size
#SBATCH --time=04:00:00             # [CLUSTER] adjust per dataset size
#SBATCH --partition=cpu-dedicated   # [CLUSTER] set to your cluster partition
##SBATCH --account=myaccount        # [CLUSTER] uncomment and set if required

# ==============================================================================
# MODULE LOADING — adapt to your cluster
# ==============================================================================
module load mafft/7.515             # [CLUSTER] change or remove if using conda
# conda activate mafft_env          # [CLUSTER] uncomment if using conda
# ==============================================================================

set -euo pipefail

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
mafft_alignments.sbatch — SLURM array job for MAFFT protein alignment

USAGE
  # 1. Count input files:
  N=$(ls -1 <input_dir>/*.fa | wc -l)

  # 2. Submit as array:
  sbatch --array=0-$((N-1))%15 mafft_alignments.sbatch \
         --input-dir  <input_dir>  \
         --output-dir <output_dir>

  # 3. Or run sequentially without SLURM:
  bash mafft_alignments.sbatch --input-dir <input_dir> --output-dir <output_dir>

REQUIRED
  -i | --input-dir  DIR    OrthoFinder Single_Copy_Orthologue_Sequences/ path
  -o | --output-dir DIR    Output directory for aligned FASTA files

OPTIONAL
  -e | --ext        EXT    File extension to look for (default: fa)
  -t | --threads    N      MAFFT --thread value (default: SLURM_CPUS_PER_TASK or 1)
  --mafft-opts      OPTS   Extra MAFFT options (default: "--auto")
                           Example: --mafft-opts "--linsi"
  -h | --help              Show this message

INPUT
  OrthoFinder Single_Copy_Orthologue_Sequences/*.fa
  Each file contains one protein sequence per species for one orthogroup.
  Typical path:
    results/orthofinder/OrthoFinder/Results_*/Single_Copy_Orthologue_Sequences/

OUTPUT
  <output_dir>/OG0000001_aln.fa   — aligned FASTA (MAFFT --auto)
  <output_dir>/<gene>_mafft.err   — per-gene MAFFT stderr (only if non-empty)

NOTE
  Already-aligned files are skipped automatically (idempotent).
  Submit the same command again to resume after failures.
EOF
}

# ── Argument parsing ──────────────────────────────────────────────────────────
INPUT_DIR=""
OUTPUT_DIR=""
EXT="fa"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
MAFFT_OPTS="--auto"

[ $# -eq 0 ] && { show_help; exit 0; }

while [ $# -gt 0 ]; do
    case "$1" in
        -i|--input-dir)   INPUT_DIR="$2";   shift 2 ;;
        -o|--output-dir)  OUTPUT_DIR="$2";  shift 2 ;;
        -e|--ext)         EXT="$2";         shift 2 ;;
        -t|--threads)     THREADS="$2";     shift 2 ;;
        --mafft-opts)     MAFFT_OPTS="$2";  shift 2 ;;
        -h|--help)        show_help; exit 0 ;;
        *) echo "ERROR: Unknown argument: '$1'" >&2; show_help >&2; exit 1 ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$INPUT_DIR"  ] || die "--input-dir is required"
[ -n "$OUTPUT_DIR" ] || die "--output-dir is required"
[ -d "$INPUT_DIR"  ] || die "Input directory not found: $INPUT_DIR"

command -v mafft &>/dev/null || \
    die "MAFFT not found. Check module loading, conda environment, or PATH."

mkdir -p "$OUTPUT_DIR"

# ── Collect input files ───────────────────────────────────────────────────────
shopt -s nullglob
FILES=("$INPUT_DIR"/*."$EXT")
shopt -u nullglob
N=${#FILES[@]}

[ "$N" -gt 0 ] || die "No .$EXT files found in $INPUT_DIR"

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  MAFFT alignment — single-copy ortholog sequences"
echo "======================================================"
echo "  Job ID      : ${SLURM_JOB_ID:-N/A}"
echo "  Array task  : ${SLURM_ARRAY_TASK_ID:-sequential}"
echo "  Start       : $(date)"
echo "  Input dir   : $INPUT_DIR"
echo "  Output dir  : $OUTPUT_DIR"
echo "  Extension   : .$EXT"
echo "  Files found : $N"
echo "  MAFFT opts  : $MAFFT_OPTS --thread $THREADS"
echo "======================================================"

# ── Per-file function ─────────────────────────────────────────────────────────
run_mafft() {
    local FILE="$1"
    local BASE
    BASE=$(basename "$FILE" ."$EXT")
    local OUT="${OUTPUT_DIR}/${BASE}_aln.fa"
    local ERRLOG="${OUTPUT_DIR}/${BASE}_mafft.err"

    if [ -s "$OUT" ]; then
        echo "[SKIP] ${BASE}.${EXT}  →  $OUT already exists"
        return 0
    fi

    echo "[RUN ] ${BASE}.${EXT}"
    # shellcheck disable=SC2086
    mafft $MAFFT_OPTS --thread "$THREADS" "$FILE" > "$OUT" 2> "$ERRLOG"

    # Remove empty error log; keep non-empty ones for inspection
    if [ ! -s "$ERRLOG" ]; then
        rm -f "$ERRLOG"
    else
        echo "[WARN] MAFFT produced stderr for $BASE — see $ERRLOG"
    fi
}

# ── Run: array mode or sequential fallback ────────────────────────────────────
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    IDX="$SLURM_ARRAY_TASK_ID"
    if (( IDX >= N )); then
        echo "Array task $IDX is out of range (0..$((N-1))). Nothing to do."
        exit 0
    fi
    run_mafft "${FILES[$IDX]}"
else
    echo "No SLURM_ARRAY_TASK_ID set — running all $N files sequentially."
    echo "To parallelise, resubmit as:"
    echo "  sbatch --array=0-$((N-1))%15 mafft_alignments.sbatch \\"
    echo "         --input-dir '$INPUT_DIR' --output-dir '$OUTPUT_DIR'"
    echo ""
    for FILE in "${FILES[@]}"; do
        run_mafft "$FILE"
    done
fi

echo ""
echo "======================================================"
echo "  Done  ($(date))"
echo "  Aligned files written to $OUTPUT_DIR/"
echo "======================================================"
