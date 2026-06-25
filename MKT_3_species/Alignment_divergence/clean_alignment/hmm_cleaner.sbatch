#!/bin/bash
# ==============================================================================
# CLUSTER CONFIGURATION — read this section before running on a new cluster
# ==============================================================================
#
# hmm_cleaner.sbatch — SLURM array job for HmmCleaner post-processing of
# MAFFT protein alignments.
#
# HmmCleaner detects and masks unreliable alignment columns by building a
# profile HMM and scoring each position; columns below the threshold are
# replaced with "?" characters (or removed, depending on downstream tool).
#
# INPUT  : MAFFT-aligned FASTA files (*_aln.fa) from mafft_alignments.sbatch
# OUTPUT : one cleaned FASTA + log + score file per orthogroup
#
# DEPENDENCIES
#   - HmmCleaner (Perl script; see INSTALLATION below)
#   - HMMER ≥ 3.3  (hmmbuild, hmmalign)
#   - Perl ≥ 5.20
#
# Lines marked [CLUSTER] MUST be adapted before the first submission.
#
# Resource guidance:
#   --cpus-per-task : HmmCleaner is single-threaded; set to 1
#   --mem           : 2–4 GB per job is sufficient
#   --time          : 1–4 h depending on alignment length and number of genes
#
# Cluster examples:
#   IFB:   --partition=cpu-dedicated  --account=<project>@cpu
#          module load bioinfo-ifb; module load hmmer/3.3.2
#
#   CIRAD: --partition=cpu-dedicated  --account=dedicated-cpu@cirad-normal
#          module load bioinfo-cirad; module load hmmer/3.3.2-singularity
#
# INSTALLATION — HmmCleaner
#   Option A — conda/bioconda (recommended, adds HmmCleaner.pl to PATH):
#     conda create -n hmmcleaner_env -c bioconda hmmcleaner hmmer
#     conda activate hmmcleaner_env
#     (then use --hmmcleaner-bin "$(which HmmCleaner.pl)")
#
#   Option B — source install:
#     git clone https://github.com/ranwez/HmmCleaner.git
#     cd HmmCleaner
#     # Note the full path to bin/HmmCleaner.pl and lib/
#     # Pass them as --hmmcleaner-bin and --hmmcleaner-lib below
#
# SUBMISSION EXAMPLE:
#   N=$(ls -1 <input_dir>/*_aln.fa | wc -l)
#   sbatch --array=0-$((N-1))%15 hmm_cleaner.sbatch \
#          --input-dir      <mafft_output_dir>  \
#          --output-dir     <hmmcleaner_output_dir> \
#          --hmmcleaner-bin /path/to/HmmCleaner.pl \
#          --hmmcleaner-lib /path/to/HmmCleaner/lib   # only for source install
#
# ==============================================================================

#SBATCH --job-name=hmmcleaner
#SBATCH --output=hmm_%A_%a.out
#SBATCH --error=hmm_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G                    # [CLUSTER] adjust per alignment size
#SBATCH --time=04:00:00             # [CLUSTER] adjust per dataset size
#SBATCH --partition=cpu-dedicated   # [CLUSTER] set to your cluster partition
##SBATCH --account=myaccount        # [CLUSTER] uncomment and set if required

# ==============================================================================
# MODULE LOADING — adapt to your cluster
# ==============================================================================
# [CLUSTER] Load HMMER (required by HmmCleaner); change module name if needed
module load hmmer/3.3.2             # [CLUSTER] change or remove if using conda
# [CLUSTER] Load HmmCleaner if available as a module:
# module load hmmcleaner            # [CLUSTER] uncomment if available
# conda activate hmmcleaner_env     # [CLUSTER] uncomment if using conda
# ==============================================================================

set -euo pipefail

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
hmm_cleaner.sbatch — SLURM array job for HmmCleaner alignment cleaning

USAGE
  # 1. Count input files:
  N=$(ls -1 <mafft_output_dir>/*_aln.fa | wc -l)

  # 2. Submit as array:
  sbatch --array=0-$((N-1))%15 hmm_cleaner.sbatch \
         --input-dir      <mafft_output_dir> \
         --output-dir     <output_dir> \
         --hmmcleaner-bin /path/to/HmmCleaner.pl

  # 3. Or run sequentially (no SLURM):
  bash hmm_cleaner.sbatch \
       --input-dir  <mafft_output_dir> \
       --output-dir <output_dir> \
       --hmmcleaner-bin /path/to/HmmCleaner.pl

REQUIRED
  -i | --input-dir      DIR    MAFFT output directory (contains *_aln.fa files)
  -o | --output-dir     DIR    Output directory for cleaned files
  -b | --hmmcleaner-bin PATH   Full path to HmmCleaner.pl

OPTIONAL
  -l | --hmmcleaner-lib DIR    HmmCleaner lib/ directory (source install only)
                               Not needed for conda/bioconda installations
  -e | --ext            EXT    Input file extension (default: fa)
       --pattern         PAT   Input filename pattern (default: *_aln)
                               Example: --pattern "*_aligned_NT"
  -p | --perl           PATH   Perl interpreter to use (default: perl)
  -h | --help                  Show this message

INPUT
  <input_dir>/*_aln.fa
  Protein alignments produced by mafft_alignments.sbatch (one file per
  orthogroup, all sequences already aligned).

OUTPUT (per orthogroup, in <output_dir>/)
  <gene>_aln_hmm.fasta   — cleaned alignment (unreliable columns masked)
  <gene>_aln_hmm.log     — positions removed per sequence
  <gene>_aln_hmm.score   — per-position HMM alignment scores

NOTE
  Already-cleaned files are skipped automatically (idempotent).
  Submit the same command again to resume after failures.
EOF
}

# ── Argument parsing ──────────────────────────────────────────────────────────
INPUT_DIR=""
OUTPUT_DIR=""
HMMCLEANER_BIN=""
HMMCLEANER_LIB=""
EXT="fa"
PATTERN="*_aln"
PERL_BIN="perl"

[ $# -eq 0 ] && { show_help; exit 0; }

while [ $# -gt 0 ]; do
    case "$1" in
        -i|--input-dir)       INPUT_DIR="$2";       shift 2 ;;
        -o|--output-dir)      OUTPUT_DIR="$2";      shift 2 ;;
        -b|--hmmcleaner-bin)  HMMCLEANER_BIN="$2";  shift 2 ;;
        -l|--hmmcleaner-lib)  HMMCLEANER_LIB="$2";  shift 2 ;;
        -e|--ext)             EXT="$2";              shift 2 ;;
        --pattern)            PATTERN="$2";          shift 2 ;;
        -p|--perl)            PERL_BIN="$2";         shift 2 ;;
        -h|--help)            show_help; exit 0 ;;
        *) echo "ERROR: Unknown argument: '$1'" >&2; show_help >&2; exit 1 ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$INPUT_DIR"      ] || die "--input-dir is required"
[ -n "$OUTPUT_DIR"     ] || die "--output-dir is required"
[ -n "$HMMCLEANER_BIN" ] || die "--hmmcleaner-bin is required (full path to HmmCleaner.pl)"

[ -d "$INPUT_DIR"      ] || die "Input directory not found: $INPUT_DIR"
[ -f "$HMMCLEANER_BIN" ] || die "HmmCleaner.pl not found: $HMMCLEANER_BIN"

# Verify Perl can find HmmCleaner (basic syntax check)
command -v "$PERL_BIN" &>/dev/null || die "Perl not found: $PERL_BIN"

# Verify HMMER is available (HmmCleaner calls hmmbuild and hmmalign internally)
command -v hmmbuild  &>/dev/null || die "hmmbuild not found. Check HMMER module loading."
command -v hmmalign  &>/dev/null || die "hmmalign not found. Check HMMER module loading."

# Build Perl include-path flag (-Ilib) for source installations
PERL_ILIB=""
if [ -n "$HMMCLEANER_LIB" ]; then
    [ -d "$HMMCLEANER_LIB" ] || die "HmmCleaner lib directory not found: $HMMCLEANER_LIB"
    PERL_ILIB="-I${HMMCLEANER_LIB}"
fi

mkdir -p "$OUTPUT_DIR"

# ── Collect input files ───────────────────────────────────────────────────────
shopt -s nullglob
FILES=("$INPUT_DIR"/${PATTERN}."$EXT")
shopt -u nullglob
N=${#FILES[@]}

[ "$N" -gt 0 ] || \
    die "No files matching '${PATTERN}.${EXT}' found in $INPUT_DIR"

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  HmmCleaner — alignment cleaning"
echo "======================================================"
echo "  Job ID       : ${SLURM_JOB_ID:-N/A}"
echo "  Array task   : ${SLURM_ARRAY_TASK_ID:-sequential}"
echo "  Start        : $(date)"
echo "  Input dir    : $INPUT_DIR"
echo "  Output dir   : $OUTPUT_DIR"
echo "  Pattern      : ${PATTERN}.${EXT}"
echo "  Files found  : $N"
echo "  HmmCleaner   : $HMMCLEANER_BIN"
[ -n "$HMMCLEANER_LIB" ] && echo "  Perl lib     : $HMMCLEANER_LIB"
echo "  Perl         : $PERL_BIN"
echo "======================================================"

# ── Per-file function ─────────────────────────────────────────────────────────
run_hmmcleaner() {
    local FILE="$1"
    local BASENAME
    BASENAME=$(basename "$FILE")                # e.g. OG0000001_aln.fa
    local STEM="${BASENAME%.*}"                 # e.g. OG0000001_aln
    local OUT_FASTA="${OUTPUT_DIR}/${STEM}_hmm.fasta"

    if [ -s "$OUT_FASTA" ]; then
        echo "[SKIP] $BASENAME  →  $OUT_FASTA already exists"
        return 0
    fi

    echo "[RUN ] $BASENAME"

    # HmmCleaner writes output to the directory of its input file.
    # To control output placement, we copy the input to OUTPUT_DIR,
    # run HmmCleaner from there (outputs land in OUTPUT_DIR), then
    # remove the working copy.
    local TMP_INPUT="${OUTPUT_DIR}/${BASENAME}"
    cp "$FILE" "$TMP_INPUT"

    # Run HmmCleaner in a subshell so the cd stays local
    (
        cd "$OUTPUT_DIR"
        # shellcheck disable=SC2086
        "$PERL_BIN" $PERL_ILIB "$HMMCLEANER_BIN" "$BASENAME"
    )

    # Remove the working copy (keep only HmmCleaner outputs)
    rm -f "$TMP_INPUT"
}

# ── Run: array mode or sequential fallback ────────────────────────────────────
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    IDX="$SLURM_ARRAY_TASK_ID"
    if (( IDX >= N )); then
        echo "Array task $IDX is out of range (0..$((N-1))). Nothing to do."
        exit 0
    fi
    run_hmmcleaner "${FILES[$IDX]}"
else
    echo "No SLURM_ARRAY_TASK_ID set — running all $N files sequentially."
    echo "To parallelise, resubmit as:"
    echo "  sbatch --array=0-$((N-1))%15 hmm_cleaner.sbatch \\"
    echo "         --input-dir '$INPUT_DIR' --output-dir '$OUTPUT_DIR' \\"
    echo "         --hmmcleaner-bin '$HMMCLEANER_BIN'"
    echo ""
    for FILE in "${FILES[@]}"; do
        run_hmmcleaner "$FILE"
    done
fi

echo ""
echo "======================================================"
echo "  Done  ($(date))"
echo "  Cleaned files written to $OUTPUT_DIR/"
echo "======================================================"
