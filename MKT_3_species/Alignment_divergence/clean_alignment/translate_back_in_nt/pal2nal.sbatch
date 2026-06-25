#!/bin/bash
# ==============================================================================
# CLUSTER CONFIGURATION — read this section before running on a new cluster
# ==============================================================================
#
# pal2nal.sbatch — Convert protein alignments to codon-level nucleotide
# alignments using PAL2NAL (Suyama et al. 2006; fork by liaochenlanruo).
#
# PAL2NAL back-translates a protein multiple-sequence alignment (from MAFFT /
# HmmCleaner) into a nucleotide alignment using unaligned CDS sequences.
# The protein alignment structure is preserved exactly.
#
# pal2nal.pl is downloaded automatically from GitHub on first run.
# No manual installation is required.
#
# INPUT  A: protein alignment directory   (e.g., HmmCleaner output)
# INPUT  B: nucleotide CDS directory      (one .fa per orthogroup)
# OUTPUT  : nucleotide alignment directory (one .fasta per orthogroup)
#
# FILE MATCHING
#   For each protein alignment (e.g., OG0000001_aln_hmm.fasta), the script
#   extracts the orthogroup ID (OG0000001) and looks for the corresponding
#   nucleotide file in --nt-dir:
#     OG0000001.fa, OG0000001.fasta, OG0000001_nt.fa, OG0000001_nt.fasta
#   Genes with no matching nucleotide file are skipped (with a warning).
#
# NUCLEOTIDE CDS DIRECTORY
#   Must contain one FASTA file per orthogroup, with sequences in the same
#   order and with matching IDs as the protein alignment.
#   Typical way to prepare:
#     • Run OrthoFinder in nucleotide mode (-d flag) — not recommended
#     • Or use a CDS-extraction script to pull CDS per orthogroup from
#       the original transcriptomes, using OrthoFinder's Orthogroups.tsv
#
# Lines marked [CLUSTER] MUST be adapted before the first submission.
#
# Resource guidance:
#   --cpus-per-task : PAL2NAL is single-threaded; set to 1
#   --mem           : 1–2 GB per job is sufficient
#   --time          : <1 h for typical datasets
#
# SUBMISSION EXAMPLE:
#   N=$(ls -1 <protein_aln_dir>/*.fasta | wc -l)
#   sbatch --array=0-$((N-1))%20 pal2nal.sbatch \
#          --protein-aln-dir  <protein_aln_dir> \
#          --nt-dir           <nt_cds_dir> \
#          --output-dir       <output_dir>
#
# ==============================================================================

#SBATCH --job-name=pal2nal
#SBATCH --output=pal2nal_%A_%a.out
#SBATCH --error=pal2nal_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G                    # [CLUSTER] adjust per alignment size
#SBATCH --time=01:00:00             # [CLUSTER] adjust per dataset size
#SBATCH --partition=cpu-dedicated   # [CLUSTER] set to your cluster partition
##SBATCH --account=myaccount        # [CLUSTER] uncomment and set if required

# ==============================================================================
# MODULE LOADING — adapt to your cluster
# ==============================================================================
# PAL2NAL is a Perl script — Perl must be in PATH (usually available by default)
# module load perl/5.36               # [CLUSTER] uncomment if Perl needs loading
# conda activate pal2nal_env          # [CLUSTER] uncomment if using a conda env
# ==============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
# Default cache directory for the auto-downloaded pal2nal.pl
PAL2NAL_CACHE_DEFAULT="${HOME}/.local/share/pal2nal"
PAL2NAL_URL="https://raw.githubusercontent.com/liaochenlanruo/PAL2NAL/master/pal2nal.pl"

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
pal2nal.sbatch — Convert protein alignments to nucleotide alignments (PAL2NAL)

PAL2NAL is downloaded automatically from GitHub on first run and cached
locally (~/.local/share/pal2nal/pal2nal.pl).

USAGE
  # 1. Count input files:
  N=$(ls -1 <protein_aln_dir>/*.fasta | wc -l)

  # 2. Submit as SLURM array:
  sbatch --array=0-$((N-1))%20 pal2nal.sbatch \
         --protein-aln-dir <hmmcleaner_out/> \
         --nt-dir           <nt_cds_per_orthogroup/> \
         --output-dir       <nt_aln_out/>

  # 3. Or run sequentially without SLURM:
  bash pal2nal.sbatch \
       --protein-aln-dir <hmmcleaner_out/> \
       --nt-dir           <nt_cds_per_orthogroup/> \
       --output-dir       <nt_aln_out/>

REQUIRED
  -a | --protein-aln-dir DIR   Protein alignment directory (HmmCleaner output)
  -n | --nt-dir           DIR   Nucleotide CDS directory (per-orthogroup FASTAs)
  -o | --output-dir       DIR   Output directory for nucleotide alignments

OPTIONAL
  --pal2nal-bin       PATH   Path to an existing pal2nal.pl (skips download)
  --pal2nal-cache     DIR    Directory to cache pal2nal.pl (default: ~/.local/share/pal2nal)
  --prot-ext          EXT    Protein file extension (default: fasta)
  --nt-ext            EXT    Nucleotide file extension (default: fa,fasta)
  --format            FMT    PAL2NAL output format: fasta|paml|clustal|codon
                             (default: fasta)
  --nogap                    Remove alignment columns that contain gaps
  --codontable        N      Codon table number (default: 1 = standard)
                             1=universal, 2=vertebrate mitochondrial, 5=invertebrate mt...
  -h | --help                Show this message

INPUT FILE MATCHING
  For each protein alignment file (e.g. OG0000001_aln_hmm.fasta), the script
  extracts the prefix before the first underscore-delimited non-alphanumeric
  part. Matching is tried against these nucleotide file patterns (in order):
    OG0000001.fa, OG0000001.fasta, OG0000001_nt.fa, OG0000001_nt.fasta
  Unmatched genes are skipped with a warning in the SLURM log.

OUTPUT
  <output_dir>/OG0000001_aln_hmm_nt.fasta   — codon-level nucleotide alignment
  Unaligned columns (from PAL2NAL) are removed when --nogap is used.

NOTE
  Already-converted files are skipped automatically (idempotent).
  Submit the same command again to resume after partial failures.
EOF
}

# ── Argument parsing ──────────────────────────────────────────────────────────
PROTEIN_ALN_DIR=""
NT_DIR=""
OUTPUT_DIR=""
PAL2NAL_BIN=""
PAL2NAL_CACHE="$PAL2NAL_CACHE_DEFAULT"
PROT_EXT="fasta"
NT_EXT_LIST="fa fasta"
FORMAT="fasta"
NOGAP=""
CODONTABLE="1"

[ $# -eq 0 ] && { show_help; exit 0; }

while [ $# -gt 0 ]; do
    case "$1" in
        -a|--protein-aln-dir) PROTEIN_ALN_DIR="$2"; shift 2 ;;
        -n|--nt-dir)          NT_DIR="$2";           shift 2 ;;
        -o|--output-dir)      OUTPUT_DIR="$2";       shift 2 ;;
        --pal2nal-bin)        PAL2NAL_BIN="$2";      shift 2 ;;
        --pal2nal-cache)      PAL2NAL_CACHE="$2";    shift 2 ;;
        --prot-ext)           PROT_EXT="$2";         shift 2 ;;
        --nt-ext)             NT_EXT_LIST="$2";      shift 2 ;;
        --format)             FORMAT="$2";           shift 2 ;;
        --nogap)              NOGAP="-nogap";        shift   ;;
        --codontable)         CODONTABLE="$2";       shift 2 ;;
        -h|--help)            show_help; exit 0 ;;
        *) echo "ERROR: Unknown argument: '$1'" >&2; show_help >&2; exit 1 ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$PROTEIN_ALN_DIR" ] || die "--protein-aln-dir is required"
[ -n "$NT_DIR"          ] || die "--nt-dir is required"
[ -n "$OUTPUT_DIR"      ] || die "--output-dir is required"

[ -d "$PROTEIN_ALN_DIR" ] || die "Protein alignment directory not found: $PROTEIN_ALN_DIR"
[ -d "$NT_DIR"          ] || die "Nucleotide CDS directory not found: $NT_DIR"

command -v perl &>/dev/null || die "Perl not found. Please load the Perl module or add it to PATH."

# ── Resolve PAL2NAL ───────────────────────────────────────────────────────────
find_or_download_pal2nal() {
    # 1. User-specified path
    if [ -n "$PAL2NAL_BIN" ]; then
        [ -f "$PAL2NAL_BIN" ] || die "--pal2nal-bin: file not found: $PAL2NAL_BIN"
        echo "$PAL2NAL_BIN"
        return
    fi

    # 2. Already in PATH
    if command -v pal2nal.pl &>/dev/null; then
        echo "$(command -v pal2nal.pl)"
        return
    fi

    # 3. Cached copy
    local CACHED="${PAL2NAL_CACHE}/pal2nal.pl"
    if [ -f "$CACHED" ]; then
        echo "$CACHED"
        return
    fi

    # 4. Download from GitHub (liaochenlanruo fork — newest version)
    echo ""
    echo "[INFO] pal2nal.pl not found. Downloading from GitHub..."
    echo "       Source: $PAL2NAL_URL"
    mkdir -p "$PAL2NAL_CACHE"

    local DOWNLOAD_OK=0
    if command -v wget &>/dev/null; then
        wget -q -O "$CACHED" "$PAL2NAL_URL" && DOWNLOAD_OK=1
    elif command -v curl &>/dev/null; then
        curl -fsSL -o "$CACHED" "$PAL2NAL_URL" && DOWNLOAD_OK=1
    else
        die "Cannot download pal2nal.pl: neither wget nor curl is available."
    fi

    if [ "$DOWNLOAD_OK" -eq 1 ] && [ -s "$CACHED" ]; then
        chmod +x "$CACHED"
        echo "[INFO] pal2nal.pl downloaded and cached at: $CACHED"
        echo "       Subsequent runs will use the cached copy."
        echo "$CACHED"
    else
        rm -f "$CACHED"
        die "Download failed. Check internet connectivity or provide --pal2nal-bin."
    fi
}

PAL2NAL_BIN=$(find_or_download_pal2nal)

# ── Collect protein alignment files ──────────────────────────────────────────
shopt -s nullglob
FILES=("$PROTEIN_ALN_DIR"/*."$PROT_EXT")
shopt -u nullglob
N=${#FILES[@]}

[ "$N" -gt 0 ] || die "No .$PROT_EXT files found in $PROTEIN_ALN_DIR"

mkdir -p "$OUTPUT_DIR"

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  PAL2NAL — protein-to-nucleotide alignment"
echo "======================================================"
echo "  Job ID         : ${SLURM_JOB_ID:-N/A}"
echo "  Array task     : ${SLURM_ARRAY_TASK_ID:-sequential}"
echo "  Start          : $(date)"
echo "  Protein dir    : $PROTEIN_ALN_DIR  ($N files)"
echo "  NT CDS dir     : $NT_DIR"
echo "  Output dir     : $OUTPUT_DIR"
echo "  PAL2NAL        : $PAL2NAL_BIN"
echo "  Format         : $FORMAT"
echo "  Codon table    : $CODONTABLE"
echo "  Remove gaps    : ${NOGAP:-no}"
echo "======================================================"

# ── File matching ─────────────────────────────────────────────────────────────
# Extract the orthogroup ID from a protein alignment filename.
# Strategy: strip known suffixes (_aln, _hmm, _aligned, etc.) to get the base ID.
get_orthogroup_id() {
    local fname
    fname=$(basename "$1" ."$PROT_EXT")
    # Remove common pipeline suffixes (longest match wins because we strip iteratively)
    fname="${fname%_hmm}"
    fname="${fname%_aln}"
    fname="${fname%_aligned}"
    fname="${fname%_filtered}"
    fname="${fname%_nt}"
    echo "$fname"
}

find_nt_file() {
    local OG_ID="$1"
    for ext in $NT_EXT_LIST; do
        for suffix in "" "_nt"; do
            local candidate="${NT_DIR}/${OG_ID}${suffix}.${ext}"
            if [ -f "$candidate" ]; then
                echo "$candidate"
                return 0
            fi
        done
    done
    return 1
}

# ── Per-file function ─────────────────────────────────────────────────────────
run_pal2nal() {
    local PROT_FILE="$1"
    local BASENAME
    BASENAME=$(basename "$PROT_FILE" ."$PROT_EXT")
    local OUT="${OUTPUT_DIR}/${BASENAME}_nt.fasta"

    if [ -s "$OUT" ]; then
        echo "[SKIP] $BASENAME (output already exists)"
        return 0
    fi

    local OG_ID
    OG_ID=$(get_orthogroup_id "$PROT_FILE")

    local NT_FILE
    if ! NT_FILE=$(find_nt_file "$OG_ID"); then
        echo "[SKIP] $BASENAME — no nucleotide file found for orthogroup '$OG_ID' in $NT_DIR"
        echo "       Tried patterns: ${OG_ID}.{$NT_EXT_LIST}, ${OG_ID}_nt.{$NT_EXT_LIST}"
        return 0
    fi

    echo "[RUN ] $BASENAME  +  $(basename "$NT_FILE")  →  $(basename "$OUT")"

    local TMP_OUT="${OUTPUT_DIR}/${BASENAME}_nt_tmp.fasta"

    # shellcheck disable=SC2086
    perl "$PAL2NAL_BIN" \
        "$PROT_FILE" \
        "$NT_FILE" \
        -output "$FORMAT" \
        -codontable "$CODONTABLE" \
        $NOGAP \
        > "$TMP_OUT" 2> "${OUTPUT_DIR}/${BASENAME}_pal2nal.err"

    # Verify output is non-empty (PAL2NAL exits 0 even on mismatches)
    if [ -s "$TMP_OUT" ]; then
        mv "$TMP_OUT" "$OUT"
        # Remove empty error log
        [ -s "${OUTPUT_DIR}/${BASENAME}_pal2nal.err" ] || \
            rm -f "${OUTPUT_DIR}/${BASENAME}_pal2nal.err"
    else
        rm -f "$TMP_OUT"
        echo "[WARN] Empty output for $BASENAME — see ${BASENAME}_pal2nal.err"
        echo "       Common causes: sequence ID mismatch between protein and NT files,"
        echo "       or codon table mismatch (try --codontable 2 for mitochondrial)."
    fi
}

# ── Run: array mode or sequential fallback ────────────────────────────────────
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    IDX="$SLURM_ARRAY_TASK_ID"
    if (( IDX >= N )); then
        echo "Array task $IDX is out of range (0..$((N-1))). Nothing to do."
        exit 0
    fi
    run_pal2nal "${FILES[$IDX]}"
else
    echo "No SLURM_ARRAY_TASK_ID set — running all $N files sequentially."
    echo "To parallelise, resubmit as:"
    echo "  sbatch --array=0-$((N-1))%20 pal2nal.sbatch \\"
    echo "         --protein-aln-dir '$PROTEIN_ALN_DIR' \\"
    echo "         --nt-dir          '$NT_DIR' \\"
    echo "         --output-dir      '$OUTPUT_DIR'"
    echo ""
    n_ok=0; n_skip=0
    for FILE in "${FILES[@]}"; do
        run_pal2nal "$FILE" && n_ok=$((n_ok + 1)) || n_skip=$((n_skip + 1))
    done
    echo ""
    echo "Summary: $n_ok processed, $n_skip skipped/failed"
fi

echo ""
echo "======================================================"
echo "  Done  ($(date))"
echo "  Nucleotide alignments written to $OUTPUT_DIR/"
echo "======================================================"
