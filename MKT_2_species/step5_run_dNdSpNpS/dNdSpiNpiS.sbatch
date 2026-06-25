#!/bin/bash
# ==============================================================================
# CLUSTER CONFIGURATION — read this section before running on a new cluster
# ==============================================================================
#
# SLURM array job: runs dNdSpiNpiS v1.0 on one codon-aligned FASTA file per
# array task.  Input is the directory produced by step 4 (headers standardised,
# focal individuals first, outgroup last).
#
# Lines marked [CLUSTER] MUST be adapted before submission.
#
# dNdSpiNpiS v1.0 binary:
#   Download from https://kimura.univ-montp2.fr/PopPhyl/index.php?section=tools
#   and pass the path with -b.
#
# Resource guidance:
#   --mem   : 10G is sufficient for most genes.
#   --time  : a few minutes per gene in most cases.
#   %N      : maximum simultaneous tasks; 15 is a safe default.
#
# Submission command:
#   N=$(ls INPUT_DIR/*.fasta | wc -l)
#   sbatch --array=0-$((N-1))%15 dNdSpiNpiS.sbatch \
#       -i INPUT_DIR -o OUTPUT_DIR -b /path/to/dNdSpiNpiS_1.0 \
#       -g Ae_speltoides -k H_vulgare
#
# Cluster examples:
#   IFB:   --partition=cpu-dedicated  --account=<project>@cpu
#   CIRAD: --partition=cpu-dedicated  --account=dedicated-cpu@cirad-normal
#
# ==============================================================================

#SBATCH --job-name=dNdSpiNpiS
#SBATCH --output=log_%A_%a_out.txt
#SBATCH --error=log_%A_%a_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=06:00:00
#SBATCH --account=dedicated-cpu@cirad-normal   # [CLUSTER] change to your account
#SBATCH --partition=cpu-dedicated              # [CLUSTER] change to your partition
#SBATCH --array=0-0%15                         # set at submission (see header)

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT_DIR=""
OUTPUT_DIR=""
BINARY=""
INGROUP=""
OUTGROUP=""
GC_BIN=8
GAP_N=10
KAPPA=1

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
Usage: sbatch --array=0-$((N-1))%15 dNdSpiNpiS.sbatch \
         -i INPUT_DIR -o OUTPUT_DIR -b BINARY -g INGROUP -k OUTGROUP

Required:
  -i DIR    Input directory with per-gene codon-aligned FASTA files
            (output of step4_prep_dNdS, headers standardised)
  -o DIR    Output directory — one result file per gene (created if absent)
  -b FILE   Path to the dNdSpiNpiS_1.0 binary
  -g NAME   Ingroup species name (must match the |species| field in headers)
  -k NAME   Outgroup species name (must match the |outgroup| field in headers)

Optional:
  -c INT    GC bins for analysis (default: 8)
  -n INT    Max % of N/gap sites per position (default: 10)
  -t FLOAT  Transition/transversion ratio kappa (default: 1)
  -h        Show this message

Submission:
  N=$(ls INPUT_DIR/*.fasta | wc -l)
  sbatch --array=0-$((N-1))%15 dNdSpiNpiS.sbatch \
      -i step4_out/ -o step5_out/ \
      -b /path/to/dNdSpiNpiS_1.0 \
      -g Ae_speltoides -k H_vulgare

Binary download:
  https://kimura.univ-montp2.fr/PopPhyl/index.php?section=tools

Outputs (in OUTPUT_DIR):
  <gene>.out       Raw dNdSpiNpiS result (dN, dS, pN, pS, dN/dS, pN/pS, DoS, NI)
  log_<job>.txt    SLURM log per task

Merge all results:
  head -1 $(ls step5_out/*.out | head -1) > all_genes_dNdSpNpS.tsv
  for f in step5_out/*.out; do tail -n +2 "$f"; done >> all_genes_dNdSpNpS.tsv
EOF
}

# ── Parse arguments ───────────────────────────────────────────────────────────
[ $# -eq 0 ] && { show_help; exit 0; }

while getopts ":i:o:b:g:k:c:n:t:h" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG"  ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        b) BINARY="$OPTARG"     ;;
        g) INGROUP="$OPTARG"    ;;
        k) OUTGROUP="$OPTARG"   ;;
        c) GC_BIN="$OPTARG"     ;;
        n) GAP_N="$OPTARG"      ;;
        t) KAPPA="$OPTARG"      ;;
        h) show_help; exit 0    ;;
        :) echo "ERROR: Option -$OPTARG requires an argument." >&2; exit 1 ;;
       \?) echo "ERROR: Unknown option: -$OPTARG." >&2; show_help >&2; exit 1 ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$INPUT_DIR" ]  || die "-i INPUT_DIR is required"
[ -n "$OUTPUT_DIR" ] || die "-o OUTPUT_DIR is required"
[ -n "$BINARY" ]     || die "-b BINARY (path to dNdSpiNpiS_1.0) is required"
[ -n "$INGROUP" ]    || die "-g INGROUP species name is required"
[ -n "$OUTGROUP" ]   || die "-k OUTGROUP species name is required"
[ -d "$INPUT_DIR" ]  || die "Input directory not found: $INPUT_DIR"
[ -f "$BINARY" ]     || die "dNdSpiNpiS binary not found: $BINARY"
[ -x "$BINARY" ]     || die "dNdSpiNpiS binary is not executable: $BINARY (run: chmod +x $BINARY)"

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
OUT_FILE="${OUTPUT_DIR}/${BASENAME}.out"

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  dNdSpiNpiS v1.0"
echo "======================================================"
echo "  Job ID       : ${SLURM_JOB_ID:-N/A}"
echo "  Array task   : ${SLURM_ARRAY_TASK_ID:-0} / $((NUM_FASTA - 1))"
echo "  Start        : $(date)"
echo "  Gene         : $BASENAME"
echo "  Input FASTA  : $FASTA_FILE"
echo "  Ingroup      : $INGROUP"
echo "  Outgroup     : $OUTGROUP"
echo "  GC bins      : $GC_BIN"
echo "  Gap/N cutoff : $GAP_N"
echo "  Kappa        : $KAPPA"
echo "  Output       : $OUT_FILE"
echo "======================================================"

# Skip if output already exists
if [ -f "$OUT_FILE" ] && [ -s "$OUT_FILE" ]; then
    echo "[INFO] Output already exists for $BASENAME — skipping."
    exit 0
fi

# ── Run dNdSpiNpiS ────────────────────────────────────────────────────────────
"$BINARY" \
    -alignment_file="$FASTA_FILE" \
    -ingroup="$INGROUP" \
    -outgroup="$OUTGROUP" \
    -gc_bin="$GC_BIN" \
    -gapN_site="$GAP_N" \
    -kappa="$KAPPA" \
    -compute_distances=yes \
    -out="${OUTPUT_DIR}/${BASENAME}" > "$OUT_FILE"

# ── Summary ───────────────────────────────────────────────────────────────────
if [ -s "$OUT_FILE" ]; then
    echo "[INFO] Output written: $OUT_FILE"
else
    echo "[WARN] Output is empty for $BASENAME"
fi
echo "[INFO] Finished on $(date)"
