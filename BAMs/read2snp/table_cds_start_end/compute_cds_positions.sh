#!/bin/bash
# ==============================================================================
# compute_cds_positions.sh
# ==============================================================================
# Compute CDS (coding sequence) start and end positions for consensus
# transcript sequences, using only two inputs: the transcript FASTA and the
# GFF3 gene annotation.
#
# THREE STEPS RUN AUTOMATICALLY
# ------------------------------
#   Step 1  Measure the length of every transcript from the FASTA file.
#           → transcript_lengths.tab
#
#   Step 2  Extract 5' and 3' UTR lengths from the GFF3 annotation.
#           → utr_lengths.tab
#
#   Step 3  Combine lengths and UTRs to compute CDS boundaries.
#           → cds_positions.tab
#
# BACKGROUND — why not position 1?
# ---------------------------------
# Consensus transcript sequences in this project are assembled with a fixed
# flanking region added on BOTH the 5' and 3' ends of the annotated gene.
# This flanking provides sequence context for read mapping and SNP calling.
#
# As a result, the CDS does not start at position 1.  Its exact boundaries
# inside each consensus sequence are:
#
#   cds_start = FLANK + five_prime_UTR_length + 1   (1-based, inclusive)
#   cds_end   = transcript_length − FLANK − three_prime_UTR_length
#
# Default FLANK = 200 bp.  Change with --flank if your assembly used a
# different flanking window.
#
# NOTE ON STRAND
#   Consensus transcriptome assemblies are oriented 5'→3' regardless of
#   genomic strand, so the formula applies uniformly.  The strand column in
#   the output is extracted from the GFF3 for reference only.
#
# USAGE
# -----
#   bash compute_cds_positions.sh \
#       --fasta   sequences.fasta \
#       --gff3    annotation.gff3 \
#       --outdir  results/cds/
#
#   # Re-use an existing UTR table (skip step 2 — useful when re-running
#   # with a different --flank value):
#   bash compute_cds_positions.sh \
#       --fasta   sequences.fasta \
#       --utr     results/cds/utr_lengths.tab \
#       --outdir  results/cds/
#
# REQUIRED ARGUMENTS
#   --fasta  FILE    Transcript FASTA file (single- or multi-line format).
#                    Lengths are computed from this file; no pre-processing
#                    needed.
#   --outdir DIR     Output directory (created if absent).
#   And one of:
#   --gff3   FILE    GFF3 annotation (must contain mRNA, five_prime_UTR and
#                    three_prime_UTR features).  Used for UTR extraction.
#   --utr    FILE    Pre-existing UTR table (skips step 2).  Expected columns
#                    (tab-delimited, with header):
#                      mRNA_ID | five_prime_UTR_len | three_prime_UTR_len | strand
#
# OPTIONAL ARGUMENTS
#   --flank  N       Flanking region in bp on each side (default: 200).
#   --help           Print this help message and exit.
#
# OUTPUTS (all written to --outdir)
#   transcript_lengths.tab  Step 1: sequence ID → length in bp
#   utr_lengths.tab         Step 2: mRNA_ID | 5'UTR | 3'UTR | strand
#   cds_positions.tab       Step 3 (final): per-transcript CDS coordinates
#     Columns: full_transcript_id | transcript_id | cds_start | cds_end |
#              strand | five_prime_utr | three_prime_utr | original_length
#
# NEXT STEP
#   Feed cds_positions.tab to extraction_of_covered_contigs.py to extract
#   CDS sequences and filter by a sequence-completeness threshold.
#
# DEPENDENCY
#   GNU awk (gawk) — standard on Linux; needed for 3-arg match() in GFF3
#   parsing.  Verify: awk --version | head -1  (should print "GNU Awk …")
# ==============================================================================
set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
FASTA=""
GFF3=""
OUTPUT_DIR=""
FLANK=200
UTR_IN=""

# ── Helpers ───────────────────────────────────────────────────────────────────
die()  { echo "ERROR: $*" >&2; exit 1; }
info() { echo "[$(date '+%H:%M:%S')] $*"; }

show_help() {
    grep '^#' "$0" | sed 's/^# \{0,1\}//'
}

# ── Argument parsing ──────────────────────────────────────────────────────────
[ $# -eq 0 ] && { show_help; exit 0; }

while [ $# -gt 0 ]; do
    case "$1" in
        --fasta)  FASTA="$2";      shift 2 ;;
        --gff3)   GFF3="$2";       shift 2 ;;
        --outdir) OUTPUT_DIR="$2"; shift 2 ;;
        --flank)  FLANK="$2";      shift 2 ;;
        --utr)    UTR_IN="$2";     shift 2 ;;
        --help|-h) show_help; exit 0 ;;
        *) die "Unknown argument: '$1'  (run with --help for usage)" ;;
    esac
done

# ── Validate arguments ────────────────────────────────────────────────────────
[ -n "$FASTA" ]      || die "--fasta is required"
[ -f "$FASTA" ]      || die "FASTA file not found: $FASTA"
[ -n "$OUTPUT_DIR" ] || die "--outdir is required"

if [ -n "$UTR_IN" ]; then
    [ -f "$UTR_IN" ] || die "UTR file not found: $UTR_IN"
else
    [ -n "$GFF3" ] || die "--gff3 is required (or use --utr to supply an existing UTR table)"
    [ -f "$GFF3" ] || die "GFF3 file not found: $GFF3"
fi

[[ "$FLANK" =~ ^[0-9]+$ ]] || die "--flank must be a non-negative integer, got: $FLANK"

mkdir -p "$OUTPUT_DIR"
LEN_OUT="$OUTPUT_DIR/transcript_lengths.tab"
UTR_OUT="$OUTPUT_DIR/utr_lengths.tab"
CDS_OUT="$OUTPUT_DIR/cds_positions.tab"
T_START=$(date +%s)

# ── Run banner ─────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  compute_cds_positions.sh"
echo "======================================================"
echo "  FASTA file   : $FASTA"
[ -n "$GFF3" ]   && echo "  GFF3 file    : $GFF3"
[ -n "$UTR_IN" ] && echo "  UTR file     : $UTR_IN  (step 2 skipped)"
echo "  Flanking     : ${FLANK} bp"
echo "  Output dir   : $OUTPUT_DIR"
echo "======================================================"
echo ""

# ==============================================================================
# STEP 1 — Measure transcript lengths from the FASTA
# ==============================================================================
info "Step 1 — Computing transcript lengths from FASTA..."

awk '
/^>/ {
    if (seq_id != "") print seq_id "\t" seq_len
    seq_id  = substr($1, 2)   # strip ">" — use only the first whitespace-delimited token
    seq_len = 0
    next
}
{ seq_len += length($0) }    # accumulate over multi-line sequences
END {
    if (seq_id != "") print seq_id "\t" seq_len
}
' "$FASTA" > "$LEN_OUT"

N_SEQ=$(wc -l < "$LEN_OUT")
info "  Lengths table: $LEN_OUT  (${N_SEQ} sequences)"

# ==============================================================================
# STEP 2 — Extract UTR lengths from GFF3
# ==============================================================================
if [ -n "$UTR_IN" ]; then
    UTR_OUT="$UTR_IN"
    info "Step 2 skipped — using existing UTR file: $UTR_OUT"
else
    echo ""
    info "Step 2 — Extracting UTR lengths from GFF3..."

    awk '
    BEGIN {
        FS = "\t"; OFS = "\t";
        print "mRNA_ID", "five_prime_UTR_len", "three_prime_UTR_len", "strand";
    }
    $3 == "mRNA" {
        # Strip optional "transcript:" prefix from mRNA ID
        match($9, /ID=(transcript:)?([^;]+)/, idarr);
        transcript = idarr[2];
        strand     = $7;
        utr5[transcript] = 0;
        utr3[transcript] = 0;
        transcript_strand[transcript] = strand;
    }
    $3 == "five_prime_UTR" {
        match($9, /Parent=(transcript:)?([^;]+)/, parr);
        transcript = parr[2];
        utr5[transcript] += $5 - $4 + 1;
    }
    $3 == "three_prime_UTR" {
        match($9, /Parent=(transcript:)?([^;]+)/, parr);
        transcript = parr[2];
        utr3[transcript] += $5 - $4 + 1;
    }
    END {
        n = 0;
        for (t in utr5) {
            print t, utr5[t]+0, utr3[t]+0, transcript_strand[t];
            n++;
        }
        print "  " n " transcripts written to UTR table" > "/dev/stderr";
    }
    ' "$GFF3" > "$UTR_OUT"

    N_UTR=$(( $(wc -l < "$UTR_OUT") - 1 ))
    info "  UTR table    : $UTR_OUT  (${N_UTR} transcripts)"
fi

# ==============================================================================
# STEP 3 — Compute CDS positions
# ==============================================================================
echo ""
info "Step 3 — Computing CDS positions..."
info "  cds_start = ${FLANK} + five_prime_UTR_length + 1"
info "  cds_end   = transcript_length − ${FLANK} − three_prime_UTR_length"
echo ""

awk -F'\t' -v flank="$FLANK" '
# Strip "transcript:" prefix from an ID if present
function norm(id,    t) {
    t = id;
    gsub(/^transcript:/, "", t);
    return t;
}
BEGIN {
    OFS = "\t";
    print "full_transcript_id", "transcript_id", "cds_start", "cds_end",
          "strand", "five_prime_utr", "three_prime_utr", "original_length";
    warnings = 0;
}

# ── File 1: UTR table ─────────────────────────────────────────────────────────
NR == FNR {
    if (FNR == 1 && $2 ~ /[A-Za-z]/) next;   # skip header

    raw = $1;
    id  = norm(raw);

    # Register under raw, normalised, and each pipe-separated token so that
    # IDs like "GENE001|HOG12345" match against any of their parts
    utr_key[raw] = id;
    utr_key[id]  = id;
    n = split(id, parts, /\|/);
    for (i = 1; i <= n; i++) {
        if (parts[i] != "") utr_key[parts[i]] = id;
    }

    five_utr[id]   = ($2 == "" ? 0 : $2 + 0);
    three_utr[id]  = ($3 == "" ? 0 : $3 + 0);
    strand_col[id] = ($4 == "" ? "+" : $4);
    next;
}

# ── File 2: transcript lengths ─────────────────────────────────────────────────
{
    full_id  = $1;
    orig_len = ($2 == "" ? 0 : $2 + 0);
    idnorm   = norm(full_id);

    # Look up UTR info — try progressively more relaxed ID forms
    found = "";
    if      (idnorm  in utr_key) { found = utr_key[idnorm]; }
    else if (full_id in utr_key) { found = utr_key[full_id]; }
    else {
        m = split(idnorm, pparts, /\|/);
        for (j = 1; j <= m; j++) {
            if (pparts[j] in utr_key) { found = utr_key[pparts[j]]; break; }
        }
    }

    if (found == "") {
        print "Warning: no UTR annotation for " full_id " — skipping" > "/dev/stderr";
        warnings++;
        next;
    }

    f_utr = five_utr[found]  + 0;
    t_utr = three_utr[found] + 0;
    s     = strand_col[found];

    cds_start = flank + f_utr + 1;
    cds_end   = orig_len - flank - t_utr;

    if (cds_start >= cds_end) {
        print "Warning: invalid CDS for " full_id \
              " (start=" cds_start " >= end=" cds_end ") — skipping" > "/dev/stderr";
        warnings++;
        next;
    }

    print full_id, found, cds_start, cds_end, s, f_utr, t_utr, orig_len;
}
END {
    if (warnings > 0) {
        print warnings " transcript(s) skipped (see warnings above)" > "/dev/stderr";
    }
}
' "$UTR_OUT" "$LEN_OUT" > "$CDS_OUT"

N_CDS=$(( $(wc -l < "$CDS_OUT") - 1 ))
info "  CDS table    : $CDS_OUT  (${N_CDS} transcripts)"

ELAPSED=$(( $(date +%s) - T_START ))

echo ""
echo "======================================================"
echo "  Done  (${ELAPSED}s)"
echo ""
echo "  Outputs:"
printf "    %-32s  %s\n" "transcript_lengths.tab" "$LEN_OUT"
printf "    %-32s  %s\n" "utr_lengths.tab"        "$UTR_OUT"
printf "    %-32s  %s\n" "cds_positions.tab"      "$CDS_OUT"
echo ""
echo "  Next step:"
echo "    python extraction_of_covered_contigs.py \\"
echo "        --fasta            sequences.fasta \\"
echo "        --cds-table        $CDS_OUT \\"
echo "        --min-completeness 0.7 \\"
echo "        --min-fraction     0.5 \\"
echo "        --gene-list        retained_genes.txt \\"
echo "        --cds-dir          cds_sequences/"
echo "======================================================"
