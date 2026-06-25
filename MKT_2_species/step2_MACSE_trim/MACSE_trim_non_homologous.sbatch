#!/bin/bash
# ==============================================================================
# CLUSTER CONFIGURATION — read this section before running on a new cluster
# ==============================================================================
#
# SLURM array job: runs MACSE trimNonHomologousFragments on one FASTA file
# per array task.  Each task trims non-homologous regions from a multi-sequence
# FASTA (one gene per file, all individuals + outgroup).
#
# MACSE writes output files in the same directory as the input; this script
# uses a per-task temporary directory to avoid touching the input directory.
#
# Lines marked [CLUSTER] MUST be adapted before submission.
#
# Resource guidance:
#   --mem     : 30G covers most transcriptome-scale genes. Reduce to 10G for
#               small genes; increase to 60G+ for large genes/alignments.
#   --time    : trimNonHomologousFragments typically takes 1-60 min per gene.
#               24h is a safe upper limit for a large batch.
#   --array%N : maximum N simultaneous tasks. 15 avoids overloading the I/O.
#
# MACSE installation (see README.md):
#   wget https://github.com/ranwez/MACSE_V2_PIPELINES/archive/refs/tags/v12.01.tar.gz
#   tar xzf v12.01.tar.gz
#   MACSE_JAR=$(find MACSE_V2_PIPELINES-12.01/ -name "*.jar" | head -1)
#
# Submission command:
#   N=$(ls INPUT_DIR/*.fasta | wc -l)
#   sbatch --array=0-$((N-1))%15 MACSE_trim_non_homologous.sbatch \
#       -i INPUT_DIR -o OUTPUT_DIR -j /path/to/macse.jar
#
# Cluster examples:
#   IFB:   --partition=cpu-dedicated  --account=<project>@cpu
#          module load bioinfo-ifb; module load java-jdk
#
#   CIRAD: --partition=cpu-dedicated  --account=dedicated-cpu@cirad-long
#          module load bioinfo-cirad; module load java-jdk
#
# ==============================================================================

#SBATCH --job-name=MACSE_trim
#SBATCH --output=log_%A_%a_out.txt
#SBATCH --error=log_%A_%a_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --account=dedicated-cpu@cirad-long    # [CLUSTER] change to your account
#SBATCH --partition=cpu-dedicated             # [CLUSTER] change to your partition
#SBATCH --array=0-0%15                        # set at submission (see header)

# ==============================================================================
# MODULE LOADING — adapt to your cluster
# ==============================================================================
module load bioinfo-ifb   # [CLUSTER] change or remove if using conda / system Java
module load java-jdk      # [CLUSTER] change to your Java module (e.g. java/17)
# ==============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT_DIR=""
OUTPUT_DIR=""
MACSE_JAR=""
MIN_HOMOLOGY="0.6"
JAVA_MEM="24g"

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
Usage: sbatch --array=0-$((N-1))%15 MACSE_trim_non_homologous.sbatch \
         -i INPUT_DIR -o OUTPUT_DIR -j MACSE_JAR [-m MIN_HOMOLOGY]

Required:
  -i DIR    Input directory with per-gene FASTA files (.fasta)
            Each file must contain unaligned nucleotide sequences for
            one gene: focal individuals + outgroup (added by add_outgroup.sh)
  -o DIR    Output directory (created if absent). Subdirectories:
              NT/            — trimmed nucleotide alignments
              AA/            — trimmed amino acid alignments
              mask_detail/   — NT with masking details
              trim_stats/    — per-gene trim statistics CSV
  -j FILE   Path to MACSE JAR file (e.g. /path/to/macse_v2.07.jar)

Optional:
  -m FLOAT  Minimum fraction of homologous sites to keep a sequence (default: 0.6)
  -h        Show this message

Submission:
  N=$(ls INPUT_DIR/*.fasta | wc -l)
  sbatch --array=0-$((N-1))%15 MACSE_trim_non_homologous.sbatch \
      -i INPUT_DIR -o OUTPUT_DIR -j /path/to/macse.jar

MACSE installation (MACSE_V2_PIPELINES v12.01):
  wget https://github.com/ranwez/MACSE_V2_PIPELINES/archive/refs/tags/v12.01.tar.gz
  tar xzf v12.01.tar.gz
  MACSE_JAR=$(find MACSE_V2_PIPELINES-12.01/ -name "*.jar" | head -1)
EOF
}

# ── Parse arguments ───────────────────────────────────────────────────────────
[ $# -eq 0 ] && { show_help; exit 0; }

while getopts ":i:o:j:m:h" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG"   ;;
        o) OUTPUT_DIR="$OPTARG"  ;;
        j) MACSE_JAR="$OPTARG"   ;;
        m) MIN_HOMOLOGY="$OPTARG" ;;
        h) show_help; exit 0     ;;
        :) echo "ERROR: Option -$OPTARG requires an argument." >&2; exit 1 ;;
       \?) echo "ERROR: Unknown option: -$OPTARG." >&2; show_help >&2; exit 1 ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$INPUT_DIR" ]  || die "-i INPUT_DIR is required"
[ -n "$OUTPUT_DIR" ] || die "-o OUTPUT_DIR is required"
[ -n "$MACSE_JAR" ]  || die "-j MACSE_JAR is required"
[ -d "$INPUT_DIR" ]  || die "Input directory not found: $INPUT_DIR"
[ -f "$MACSE_JAR" ]  || die "MACSE JAR not found: $MACSE_JAR"

command -v java &>/dev/null || die "java not found — load the Java module."

# ── Output directory structure ────────────────────────────────────────────────
mkdir -p "$OUTPUT_DIR/NT"
mkdir -p "$OUTPUT_DIR/AA"
mkdir -p "$OUTPUT_DIR/mask_detail"
mkdir -p "$OUTPUT_DIR/trim_stats"

# ── Select FASTA for this array task ─────────────────────────────────────────
mapfile -t FASTA_FILES < <(find "$INPUT_DIR" -maxdepth 1 -name "*.fasta" | sort)
NUM_FASTA=${#FASTA_FILES[@]}

[ "$NUM_FASTA" -gt 0 ] || die "No .fasta files found in $INPUT_DIR"

if [ "${SLURM_ARRAY_TASK_ID:-0}" -ge "$NUM_FASTA" ]; then
    echo "[INFO] Task ID ${SLURM_ARRAY_TASK_ID} >= number of files ($NUM_FASTA) — nothing to do."
    exit 0
fi

FASTA_FILE="${FASTA_FILES[${SLURM_ARRAY_TASK_ID:-0}]}"
BASENAME=$(basename "$FASTA_FILE" .fasta)

# Expected outputs
OUT_NT="${OUTPUT_DIR}/NT/${BASENAME}_NT.fasta"
OUT_AA="${OUTPUT_DIR}/AA/${BASENAME}_AA.fasta"
OUT_MASK="${OUTPUT_DIR}/mask_detail/${BASENAME}_mask_detail_NT.fasta"
TRIM_STATS="${OUTPUT_DIR}/trim_stats/${BASENAME}_trim_stats.csv"

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  MACSE trimNonHomologousFragments"
echo "======================================================"
echo "  Job ID        : ${SLURM_JOB_ID:-N/A}"
echo "  Array task    : ${SLURM_ARRAY_TASK_ID:-0} / $((NUM_FASTA - 1))"
echo "  Start         : $(date)"
echo "  Gene          : $BASENAME"
echo "  Input FASTA   : $FASTA_FILE"
echo "  MACSE JAR     : $MACSE_JAR"
echo "  Min homology  : $MIN_HOMOLOGY"
echo "  Output dir    : $OUTPUT_DIR"
echo "======================================================"

# Skip if outputs already exist (allows re-submission after partial failure)
if [ -f "$OUT_NT" ] && [ -f "$OUT_AA" ]; then
    echo "[INFO] Outputs already exist for $BASENAME — skipping."
    exit 0
fi

# ── Run MACSE in a temporary directory ────────────────────────────────────────
# MACSE writes outputs next to its input file; we use a temp dir to avoid
# modifying the input directory.
WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT

WORK_FASTA="$WORK_DIR/${BASENAME}.fasta"
cp "$FASTA_FILE" "$WORK_FASTA"

echo "[INFO] Running trimNonHomologousFragments..."
java -Xmx${JAVA_MEM} -jar "$MACSE_JAR" \
    -prog   trimNonHomologousFragments \
    -seq    "$WORK_FASTA" \
    -out_trim_info         "$TRIM_STATS" \
    -min_homology_to_keep_seq "$MIN_HOMOLOGY"

# ── Move MACSE outputs to structured output directory ─────────────────────────
for ext in _NT.fasta _AA.fasta _mask_detail_NT.fasta; do
    src="$WORK_DIR/${BASENAME}${ext}"
    [ -f "$src" ] || continue
    case "$ext" in
        _NT.fasta)              dst="$OUT_NT"   ;;
        _AA.fasta)              dst="$OUT_AA"   ;;
        _mask_detail_NT.fasta)  dst="$OUT_MASK" ;;
    esac
    mv "$src" "$dst"
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
if [ -s "$OUT_NT" ]; then
    N_SEQS=$(grep -c '^>' "$OUT_NT" || true)
    echo "[INFO] NT output: $N_SEQS sequences in $OUT_NT"
else
    echo "[WARN] NT output is empty — all sequences may have been filtered out."
fi
[ -f "$TRIM_STATS" ] && echo "[INFO] Trim stats: $TRIM_STATS"
echo "[INFO] Finished on $(date)"
