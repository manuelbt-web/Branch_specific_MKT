#!/bin/bash
# ==============================================================================
# CLUSTER CONFIGURATION — read this section before running on a new cluster
# ==============================================================================
#
# SLURM array job: runs MACSE alignSequences on one FASTA file per array task.
# Input is the NT output directory from MACSE_trim_non_homologous.sbatch
# (trimmed, unaligned nucleotide sequences).
#
# Unlike trimNonHomologousFragments, MACSE alignSequences accepts explicit
# output paths (-out_NT, -out_AA) so no temporary directory is needed.
#
# Lines marked [CLUSTER] MUST be adapted before submission.
#
# Resource guidance:
#   --mem     : 30G covers most genes. MACSE alignment is memory-intensive
#               for alignments with >100 sequences; increase to 60G if needed.
#   --time    : alignSequences can take minutes to hours per gene depending
#               on gene length and number of sequences.
#   --array%N : maximum N simultaneous tasks (15 avoids I/O saturation).
#
# MACSE installation (see README.md):
#   wget https://github.com/ranwez/MACSE_V2_PIPELINES/archive/refs/tags/v12.01.tar.gz
#   tar xzf v12.01.tar.gz
#   MACSE_JAR=$(find MACSE_V2_PIPELINES-12.01/ -name "*.jar" | head -1)
#
# Submission command:
#   N=$(ls INPUT_DIR/NT/*.fasta | wc -l)
#   sbatch --array=0-$((N-1))%15 MACSE_aligning_sequences.sbatch \
#       -i INPUT_DIR/NT -o OUTPUT_DIR -j /path/to/macse.jar
#
# Cluster examples:
#   IFB:   --partition=cpu-dedicated  --account=<project>@cpu
#          module load bioinfo-ifb; module load java-jdk
#
#   CIRAD: --partition=cpu-dedicated  --account=dedicated-cpu@cirad-long
#          module load bioinfo-cirad; module load java-jdk
#
# ==============================================================================

#SBATCH --job-name=MACSE_align
#SBATCH --output=log_%A_%a_out.txt
#SBATCH --error=log_%A_%a_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=84:00:00
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
JAVA_MEM="24g"

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
Usage: sbatch --array=0-$((N-1))%15 MACSE_aligning_sequences.sbatch \
         -i INPUT_DIR -o OUTPUT_DIR -j MACSE_JAR

Required:
  -i DIR    Input directory with trimmed nucleotide FASTA files
            (the NT/ subdirectory from MACSE_trim_non_homologous.sbatch)
  -o DIR    Output directory for aligned sequences (created if absent):
              <gene>_aligned_AA.fasta   — amino acid alignment
              <gene>_aligned_NT.fasta   — nucleotide alignment (with ! for
                                          frameshifts)
  -j FILE   Path to MACSE JAR file (e.g. /path/to/macse_v2.07.jar)

Optional:
  -h        Show this message

Submission:
  N=$(ls INPUT_DIR/NT/*.fasta | wc -l)
  sbatch --array=0-$((N-1))%15 MACSE_aligning_sequences.sbatch \
      -i INPUT_DIR/NT -o aligned/ -j /path/to/macse.jar

MACSE installation (MACSE_V2_PIPELINES v12.01):
  wget https://github.com/ranwez/MACSE_V2_PIPELINES/archive/refs/tags/v12.01.tar.gz
  tar xzf v12.01.tar.gz
  MACSE_JAR=$(find MACSE_V2_PIPELINES-12.01/ -name "*.jar" | head -1)

Next step:
  python input_preparation_for_dNdSpNpS.py \
      --input-dir aligned/ \
      --species   <species_name> \
      --outgroup  <outgroup_name>
EOF
}

# ── Parse arguments ───────────────────────────────────────────────────────────
[ $# -eq 0 ] && { show_help; exit 0; }

while getopts ":i:o:j:h" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG"  ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        j) MACSE_JAR="$OPTARG"  ;;
        h) show_help; exit 0    ;;
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

mkdir -p "$OUTPUT_DIR"

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

OUT_AA="${OUTPUT_DIR}/${BASENAME}_aligned_AA.fasta"
OUT_NT="${OUTPUT_DIR}/${BASENAME}_aligned_NT.fasta"

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  MACSE alignSequences"
echo "======================================================"
echo "  Job ID      : ${SLURM_JOB_ID:-N/A}"
echo "  Array task  : ${SLURM_ARRAY_TASK_ID:-0} / $((NUM_FASTA - 1))"
echo "  Start       : $(date)"
echo "  Gene        : $BASENAME"
echo "  Input FASTA : $FASTA_FILE"
echo "  MACSE JAR   : $MACSE_JAR"
echo "  Output dir  : $OUTPUT_DIR"
echo "======================================================"

# Skip if outputs already exist
if [ -f "$OUT_NT" ] && [ -f "$OUT_AA" ]; then
    echo "[INFO] Outputs already exist for $BASENAME — skipping."
    exit 0
fi

# ── Run MACSE alignment ───────────────────────────────────────────────────────
echo "[INFO] Running alignSequences..."
java -Xmx${JAVA_MEM} -jar "$MACSE_JAR" \
    -prog   alignSequences \
    -seq    "$FASTA_FILE" \
    -out_AA "$OUT_AA" \
    -out_NT "$OUT_NT"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
if [ -s "$OUT_NT" ]; then
    N_SEQS=$(grep -c '^>' "$OUT_NT" || true)
    echo "[INFO] NT alignment: $N_SEQS sequences in $OUT_NT"
else
    echo "[WARN] NT alignment output is empty for $BASENAME."
fi
echo "[INFO] Finished on $(date)"
