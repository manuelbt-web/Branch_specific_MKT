#!/bin/bash
# ==============================================================================
# CLUSTER CONFIGURATION — read this section before running on a new cluster
# ==============================================================================
#
# Complete pairwise ortholog detection pipeline in one command:
#
#   Step 1  BLAST: query A → subject B      → <A>_vs_<B>.tab
#   Step 2  BLAST: query B → subject A      → <B>_vs_<A>.tab
#   Step 3  Reciprocal Best Hits (RBH)      → RBH_<A>_<B>.tab
#
# Two modes are supported:
#   nucl  — blastn   (nucleotide sequences, e.g. CDS)
#   prot  — blastp   (protein sequences / proteomes)
#
# Lines marked [CLUSTER] MUST be adapted before submission.
#
#   1. --account   : your project/billing account
#   2. --partition : the queue/partition name
#   3. module load : replace with your cluster's BLAST module
#
# Resource guidance:
#   --cpus-per-task : passed as --num_threads to BLAST. 8 is suitable for
#                     transcriptome-scale data. Use 16-32 for whole proteomes.
#   --mem           : 10G is sufficient for most transcriptome comparisons.
#                     Use 20G+ for large proteome blastp jobs.
#   --time          : blastn on CDS: ~minutes. blastp on proteomes: hours.
#
# Python 3 must be accessible for Step 3 (RBH).  On most clusters it is
# available by default after loading a bioinfo module.  If not, add:
#   module load python/3.x   (uncomment and adapt below)
#
# Cluster examples:
#   IFB:   --partition=cpu-dedicated  --account=<project>@cpu
#          module load bioinfo-ifb; module load blast/2.16.0
#
#   CIRAD: --partition=cpu-dedicated  --account=dedicated-cpu@cirad-normal
#          module load bioinfo-cirad; module load blast/2.16.0
#
#   Without modules (BLAST and Python already in PATH):
#          Comment out the module load lines below.
#
# ==============================================================================

#SBATCH --job-name=reciprocal_blast
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --account=dedicated-cpu@cirad-normal   # [CLUSTER] change to your account
#SBATCH --partition=cpu-dedicated              # [CLUSTER] change to your partition

# ==============================================================================
# MODULE LOADING — adapt to your cluster
# ==============================================================================
module load bioinfo-ifb          # [CLUSTER] change or remove if not needed
module load blast/2.16.0         # [CLUSTER] change to your available version
# module load python/3.9.13      # [CLUSTER] uncomment if Python 3 is not in PATH
# ==============================================================================

set -euo pipefail

# Locate reciprocal_best_hits.py next to this script (works with SLURM sbatch)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RBH_SCRIPT="$SCRIPT_DIR/reciprocal_best_hits.py"

# ── Defaults ──────────────────────────────────────────────────────────────────
MODE=""
QUERY=""
SUBJECT=""
OUTPUT_DIR=""
EVALUE=""
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
Usage: sbatch reciprocal_blast.sbatch [OPTIONS]

Runs a complete pairwise ortholog detection pipeline:
  Step 1: BLAST A→B
  Step 2: BLAST B→A
  Step 3: Reciprocal Best Hits (RBH)

Required:
  --mode    nucl|prot   BLAST type: nucleotide (blastn) or protein (blastp)
  --query   FILE        FASTA file for species A
  --subject FILE        FASTA file for species B
  --outdir  DIR         Output directory

Optional:
  --evalue  FLOAT       E-value cutoff (default: 1e-20 for nucl, 1e-10 for prot)
  --threads N           BLAST threads (default: SLURM_CPUS_PER_TASK or 8)
  --help                Show this message

Outputs (in --outdir):
  <A>_vs_<B>.tab           BLAST hits: A as query
  <B>_vs_<A>.tab           BLAST hits: B as query
  RBH_<A>_<B>.tab          Reciprocal Best Hit pairs
                           Columns: A_id | B_id | score_A→B | score_B→A

Examples:
  # Nucleotide (CDS sequences):
  sbatch reciprocal_blast.sbatch \
      --mode    nucl \
      --query   speltoides_cds.fasta \
      --subject tauschii_cds.fasta \
      --outdir  results/blast/

  # Protein (proteomes):
  sbatch reciprocal_blast.sbatch \
      --mode    prot \
      --query   speltoides_proteome.fasta \
      --subject tauschii_proteome.fasta \
      --outdir  results/blast/

Note: reciprocal_best_hits.py must be in the same directory as this script.
      To re-run only the RBH step on existing BLAST tables:
        python reciprocal_best_hits.py --help
EOF
}

# ── Argument parsing ──────────────────────────────────────────────────────────
[ $# -eq 0 ] && { show_help; exit 0; }

while [ $# -gt 0 ]; do
    case "$1" in
        --mode)    MODE="$2";       shift 2 ;;
        --query)   QUERY="$2";      shift 2 ;;
        --subject) SUBJECT="$2";    shift 2 ;;
        --outdir)  OUTPUT_DIR="$2"; shift 2 ;;
        --evalue)  EVALUE="$2";     shift 2 ;;
        --threads) THREADS="$2";    shift 2 ;;
        --help|-h) show_help; exit 0 ;;
        *) echo "ERROR: Unknown argument: '$1'" >&2; show_help >&2; exit 1 ;;
    esac
done

# ── Validate arguments ────────────────────────────────────────────────────────
die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$MODE" ]       || die "--mode is required (nucl or prot)"
[ -n "$QUERY" ]      || die "--query is required"
[ -n "$SUBJECT" ]    || die "--subject is required"
[ -n "$OUTPUT_DIR" ] || die "--outdir is required"
[ -f "$QUERY" ]      || die "Query file not found: $QUERY"
[ -f "$SUBJECT" ]    || die "Subject file not found: $SUBJECT"
[ -f "$RBH_SCRIPT" ] || die "reciprocal_best_hits.py not found in: $SCRIPT_DIR"

case "$MODE" in
    nucl)
        BLAST_CMD="blastn"
        [ -z "$EVALUE" ] && EVALUE="1e-20"
        ;;
    prot)
        BLAST_CMD="blastp"
        [ -z "$EVALUE" ] && EVALUE="1e-10"
        ;;
    *) die "--mode must be 'nucl' (blastn) or 'prot' (blastp), got: '$MODE'" ;;
esac

command -v "$BLAST_CMD" &>/dev/null || \
    die "$BLAST_CMD not found. Check module loading or PATH."
command -v python3 &>/dev/null || \
    die "python3 not found. Check module loading or PATH."

mkdir -p "$OUTPUT_DIR"

QUERY_ID=$(basename "$QUERY"   | sed 's/\.[^.]*$//')
SUBJECT_ID=$(basename "$SUBJECT" | sed 's/\.[^.]*$//')
OUT_AvsB="${OUTPUT_DIR}/${QUERY_ID}_vs_${SUBJECT_ID}.tab"
OUT_BvsA="${OUTPUT_DIR}/${SUBJECT_ID}_vs_${QUERY_ID}.tab"
OUT_RBH="${OUTPUT_DIR}/RBH_${QUERY_ID}_${SUBJECT_ID}.tab"

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  Reciprocal BLAST + RBH Pipeline"
echo "======================================================"
echo "  Job ID    : ${SLURM_JOB_ID:-N/A}"
echo "  Start     : $(date)"
echo "  Mode      : $MODE ($BLAST_CMD)"
echo "  Query     : $QUERY  [$QUERY_ID]"
echo "  Subject   : $SUBJECT  [$SUBJECT_ID]"
echo "  E-value   : $EVALUE"
echo "  Threads   : $THREADS"
echo "  Output    : $OUTPUT_DIR"
echo "======================================================"

# ==============================================================================
# STEP 1 — BLAST: A → B
# ==============================================================================
echo ""
echo "[1/3] BLAST: $QUERY_ID → $SUBJECT_ID ..."
"$BLAST_CMD" \
    -query       "$QUERY" \
    -subject     "$SUBJECT" \
    -evalue      "$EVALUE" \
    -outfmt      6 \
    -num_threads "$THREADS" \
    -out         "$OUT_AvsB"

echo "      $(wc -l < "$OUT_AvsB") hits  →  $(basename "$OUT_AvsB")"

# ==============================================================================
# STEP 2 — BLAST: B → A
# ==============================================================================
echo ""
echo "[2/3] BLAST: $SUBJECT_ID → $QUERY_ID ..."
"$BLAST_CMD" \
    -query       "$SUBJECT" \
    -subject     "$QUERY" \
    -evalue      "$EVALUE" \
    -outfmt      6 \
    -num_threads "$THREADS" \
    -out         "$OUT_BvsA"

echo "      $(wc -l < "$OUT_BvsA") hits  →  $(basename "$OUT_BvsA")"

[ -s "$OUT_AvsB" ] || die "Step 1 output empty: $OUT_AvsB"
[ -s "$OUT_BvsA" ] || die "Step 2 output empty: $OUT_BvsA"

# ==============================================================================
# STEP 3 — Reciprocal Best Hits
# ==============================================================================
echo ""
echo "[3/3] Finding reciprocal best hits ..."
python3 "$RBH_SCRIPT" \
    --a-vs-b "$OUT_AvsB" \
    --b-vs-a "$OUT_BvsA" \
    --output "$OUT_RBH"

RBH_COUNT=$(( $(wc -l < "$OUT_RBH") - 1 ))
echo "      $RBH_COUNT reciprocal best hits  →  $(basename "$OUT_RBH")"

# ==============================================================================
# SUMMARY
# ==============================================================================
echo ""
echo "======================================================"
echo "  Done  ($(date))"
echo ""
echo "  Output files:"
printf "    %-45s  %s hits\n"  "$(basename "$OUT_AvsB")"   "$(wc -l < "$OUT_AvsB")"
printf "    %-45s  %s hits\n"  "$(basename "$OUT_BvsA")"   "$(wc -l < "$OUT_BvsA")"
printf "    %-45s  %s pairs\n" "$(basename "$OUT_RBH")"    "$RBH_COUNT"
echo "======================================================"
