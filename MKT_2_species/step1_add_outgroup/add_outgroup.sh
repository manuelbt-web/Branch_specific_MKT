#!/bin/bash
# ==============================================================================
# Add outgroup sequences to per-gene FASTA files via an RBH ortholog table.
#
# For each gene FASTA in INPUT_DIR, this script looks up the corresponding
# outgroup sequence ID in the RBH table and appends it to a COPY of the gene
# FASTA written to OUTPUT_DIR.  Input files are never modified.
#
# Dependencies: seqkit ≥ 2.0 (must be in PATH or loaded via module).
#
# Module loading — adapt to your cluster before running:
#   module load bioinfo-ifb          # [CLUSTER] or bioinfo-cirad
#   module load seqkit/2.8.1         # [CLUSTER] adapt version
#
# ==============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT_DIR=""
OUTPUT_DIR=""
OUTGROUP_FASTA=""
RBH_TABLE=""
LOG_FILE=""

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
Usage: bash add_outgroup.sh -i INPUT_DIR -o OUTPUT_DIR -g OUTGROUP_FASTA -r RBH_TABLE [-l LOG]

Required:
  -i DIR    Directory with per-gene FASTA files (one .fasta per gene, all
            individuals of the focal species as separate records)
  -o DIR    Output directory — gene FASTAs with outgroup appended (created
            if absent; INPUT_DIR files are never modified)
  -g FILE   Outgroup FASTA (all outgroup CDS sequences in a single file)
  -r FILE   Reciprocal Best Hit table (TSV, column 1 = gene ID,
            column 2 = outgroup sequence ID; # lines are ignored)

Optional:
  -l FILE   Log file for skipped genes (default: OUTPUT_DIR/add_outgroup.log)
  -h        Show this message

Outputs:
  OUTPUT_DIR/<gene>.fasta   Original sequences + outgroup sequence at end
  OUTPUT_DIR/add_outgroup.log  Genes for which no ortholog was found

Example:
  # Load modules first (adapt to your cluster):
  #   module load bioinfo-ifb; module load seqkit/2.8.1

  bash add_outgroup.sh \
      -i cds_per_gene/ \
      -o cds_with_outgroup/ \
      -g H_vulgare_CDS.fasta \
      -r RBH_speltoides_vulgare.tab

Note:
  seqkit faidx creates an index file (<OUTGROUP_FASTA>.fai) in the same
  directory as the outgroup FASTA — ensure the directory is writable.
EOF
}

# ── Parse arguments ───────────────────────────────────────────────────────────
[ $# -eq 0 ] && { show_help; exit 0; }

while getopts ":i:o:g:r:l:h" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG"      ;;
        o) OUTPUT_DIR="$OPTARG"     ;;
        g) OUTGROUP_FASTA="$OPTARG" ;;
        r) RBH_TABLE="$OPTARG"      ;;
        l) LOG_FILE="$OPTARG"       ;;
        h) show_help; exit 0        ;;
        :) echo "ERROR: Option -$OPTARG requires an argument." >&2; exit 1 ;;
       \?) echo "ERROR: Unknown option: -$OPTARG." >&2; show_help >&2; exit 1 ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$INPUT_DIR" ]      || die "-i INPUT_DIR is required"
[ -n "$OUTPUT_DIR" ]     || die "-o OUTPUT_DIR is required"
[ -n "$OUTGROUP_FASTA" ] || die "-g OUTGROUP_FASTA is required"
[ -n "$RBH_TABLE" ]      || die "-r RBH_TABLE is required"
[ -d "$INPUT_DIR" ]      || die "Input directory not found: $INPUT_DIR"
[ -f "$OUTGROUP_FASTA" ] || die "Outgroup FASTA not found: $OUTGROUP_FASTA"
[ -f "$RBH_TABLE" ]      || die "RBH table not found: $RBH_TABLE"

command -v seqkit &>/dev/null || \
    die "seqkit not found — load the appropriate module or activate your conda env."

N_FASTA=$(find "$INPUT_DIR" -maxdepth 1 -name "*.fasta" | wc -l)
[ "$N_FASTA" -gt 0 ] || die "No .fasta files found in $INPUT_DIR"

mkdir -p "$OUTPUT_DIR"
[ -z "$LOG_FILE" ] && LOG_FILE="$OUTPUT_DIR/add_outgroup.log"
> "$LOG_FILE"

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  Add Outgroup Sequences"
echo "======================================================"
echo "  Input dir      : $INPUT_DIR  ($N_FASTA genes)"
echo "  Output dir     : $OUTPUT_DIR"
echo "  Outgroup FASTA : $OUTGROUP_FASTA"
echo "  RBH table      : $RBH_TABLE"
echo "  Log            : $LOG_FILE"
echo "  Start          : $(date)"
echo "======================================================"
echo ""

# Build seqkit index once so all lookups are O(1)
seqkit faidx "$OUTGROUP_FASTA" --quiet 2>/dev/null || true

TOTAL=0; FOUND=0; MISSING=0

for fasta in "$INPUT_DIR"/*.fasta; do
    [ -f "$fasta" ] || continue
    TOTAL=$((TOTAL + 1))
    gene_id=$(basename "$fasta" .fasta)

    # Look up outgroup ortholog (col 1 = gene_id, col 2 = outgroup_id)
    ortholog_id=$(grep -v '^#' "$RBH_TABLE" \
        | awk -v id="$gene_id" '$1 == id { print $2; exit }')

    if [ -z "$ortholog_id" ]; then
        echo "  [SKIP] $gene_id — no ortholog in RBH table" | tee -a "$LOG_FILE"
        MISSING=$((MISSING + 1))
        continue
    fi

    # Extract outgroup sequence (seqkit faidx exits 1 if ID not found)
    TMPFILE=$(mktemp)
    if ! seqkit faidx "$OUTGROUP_FASTA" "$ortholog_id" > "$TMPFILE" 2>/dev/null \
            || [ ! -s "$TMPFILE" ]; then
        echo "  [SKIP] $gene_id — ortholog '$ortholog_id' not in FASTA" | tee -a "$LOG_FILE"
        rm -f "$TMPFILE"
        MISSING=$((MISSING + 1))
        continue
    fi

    # Write focal sequences + outgroup to output (input file is unchanged)
    out_fasta="$OUTPUT_DIR/${gene_id}.fasta"
    cat "$fasta" "$TMPFILE" > "$out_fasta"
    rm -f "$TMPFILE"

    echo "  [OK]   $gene_id  →  $ortholog_id"
    FOUND=$((FOUND + 1))
done

echo ""
echo "======================================================"
echo "  Done  ($(date))"
echo "  Total genes   : $TOTAL"
echo "  With outgroup : $FOUND"
echo "  Skipped       : $MISSING  (details in $LOG_FILE)"
echo "======================================================"
