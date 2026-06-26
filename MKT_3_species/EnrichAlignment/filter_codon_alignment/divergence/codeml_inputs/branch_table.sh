#!/usr/bin/env bash
# ==============================================================================
# branch_table.sh — Prepare codeml input directories after VESPA infer_genetree
# ==============================================================================
#
# For each gene directory produced by VESPA (Inferred_Genetree_*/), this script:
#   1. Creates a codeml_input/ subdirectory.
#   2. Copies the NT FASTA and cleans sequence headers:
#        - removes everything after '|' (species|gene_id → species)
#        - removes underscores (Aegilops_tauschii → Aegilopstauschii)
#   3. Copies the gene tree (.tre) and cleans its labels the same way.
#   4. Writes branch_table.txt listing the species for the VESPA branch model.
#        By default ALL species found in the FASTA headers are written
#        (one per line), so codeml runs on every branch of the phylogeny.
#        Pass --target-species to restrict to a subset.
#
# VESPA codeml_setup (next step: codeml_set_up.sh) reads codeml_input/ and
# branch_table.txt to generate the full codeml workspace — one
# cleaned_<species>/ directory per species listed in branch_table.txt.
#
# USAGE
#   # Default: run codeml on EVERY branch (species auto-detected from FASTA)
#   bash branch_table.sh --root-dir vespa_out/
#
#   # Restrict to a subset of species
#   bash branch_table.sh \
#       --root-dir       vespa_out/ \
#       --target-species "Aegilopsspeltoides,Aegilopsmutica"
#
#   # Custom gene-directory pattern (default: match all subdirs):
#   bash branch_table.sh \
#       --root-dir       vespa_out/ \
#       --gene-pattern   "Am*"
#
# NOTE ON SPECIES NAMES IN BRANCH TABLE
#   The names must match the cleaned headers (no pipes, no underscores).
#   Original header:  >Aegilops_speltoides|TraesCS1A02G000001
#   Cleaned name:      Aegilopsspeltoides
#   (The header is split at '|', then all '_' are removed.)
#   Auto-detection reads these cleaned names directly from cleaned.fasta.
#
# NOTE ON VESPA DIRECTORY STRUCTURE
#   After VESPA infer_genetree, directories are named:
#     Inferred_Genetree_<HOG_id>/<gene_dir>/
#   where <gene_dir> contains the original FASTA and the generated .tre file.
#   This script loops over <gene_dir>/ matching --gene-pattern within any
#   Inferred_Genetree_*/ directory under --root-dir.
# ==============================================================================
set -euo pipefail
shopt -s nullglob

# ── Defaults ──────────────────────────────────────────────────────────────────
ROOT_DIR=""
TARGET_SPECIES=""      # optional; if empty, all species from FASTA are used
GENE_PATTERN="*"       # glob pattern for the gene-level subdirectory
FASTA_EXT="fasta"
FORCE=false

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
branch_table.sh — Prepare VESPA codeml_input directories

USAGE
  # Default: run codeml on every branch (species auto-detected from FASTA headers)
  bash branch_table.sh --root-dir <dir>

  # Restrict to specific species only
  bash branch_table.sh \
      --root-dir       <dir> \
      --target-species "<Sp1>,<Sp2>,..."

REQUIRED
  -r | --root-dir       DIR     Root directory containing Inferred_Genetree_*/ dirs
                                (output from gene_tree.sh)

OPTIONAL
  -t | --target-species LIST    Comma-separated foreground species for the branch
                                model (cleaned names: no underscores, no pipes).
                                Example: "Aegilopsspeltoides,Aegilopsmutica"
                                Default: all species found in the FASTA headers
                                (one cleaned_<species>/ model per branch).
  -p | --gene-pattern   PAT     Glob pattern for gene subdirectory name inside
                                Inferred_Genetree_*/  (default: * = all)
                                Example: "Am*"  to match Amblyopyrum muticum genes
  -e | --fasta-ext      EXT     FASTA file extension (default: fasta)
  -f | --force                  Overwrite existing codeml_input/ directories
  -h | --help                   Show this message

SPECIES NAME FORMAT
  Names in --target-species must match the cleaned alignment headers:
    Original: >Aegilops_speltoides|TraesCS1A02G000001
    Cleaned:   Aegilopsspeltoides   (split at '|', remove '_')
  When --target-species is omitted, species are read from the cleaned FASTA
  headers after copying — no manual specification needed.

OUTPUT (per gene directory)
  <gene_dir>/branch_table.txt        — species list (one per line)
  <gene_dir>/codeml_input/cleaned.fasta  — cleaned NT alignment
  <gene_dir>/codeml_input/cleaned.tre    — cleaned gene tree
EOF
}

[ $# -eq 0 ] && { show_help; exit 0; }

# ── Argument parsing ──────────────────────────────────────────────────────────
while [ $# -gt 0 ]; do
    case "$1" in
        -r|--root-dir)        ROOT_DIR="$2";         shift 2 ;;
        -t|--target-species)  TARGET_SPECIES="$2";   shift 2 ;;
        -p|--gene-pattern)    GENE_PATTERN="$2";     shift 2 ;;
        -e|--fasta-ext)       FASTA_EXT="$2";        shift 2 ;;
        -f|--force)           FORCE=true;            shift ;;
        -h|--help)            show_help; exit 0 ;;
        *) echo "ERROR: Unknown argument: '$1'" >&2; show_help >&2; exit 1 ;;
    esac
done

die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$ROOT_DIR" ] || die "--root-dir is required"
[ -d "$ROOT_DIR" ] || die "Root directory not found: $ROOT_DIR"

ROOT_DIR="$(cd "$ROOT_DIR" && pwd)"

# Pre-compute branch_table content when --target-species is set (shared across all genes)
FIXED_BRANCH_TABLE_CONTENT=""
if [ -n "$TARGET_SPECIES" ]; then
    IFS=',' read -ra SPECIES_ARRAY <<< "$TARGET_SPECIES"
    for sp in "${SPECIES_ARRAY[@]}"; do
        sp="${sp// /}"
        FIXED_BRANCH_TABLE_CONTENT="${FIXED_BRANCH_TABLE_CONTENT}${sp}"$'\n'
    done
fi

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  branch_table.sh — Prepare codeml input directories"
echo "======================================================"
echo "  Root dir         : $ROOT_DIR"
if [ -n "$TARGET_SPECIES" ]; then
    echo "  Target species   : $TARGET_SPECIES"
else
    echo "  Target species   : [auto-detect from FASTA headers — all branches]"
fi
echo "  Gene pattern     : $GENE_PATTERN"
echo "  FASTA extension  : .$FASTA_EXT"
echo "======================================================"
echo ""

# ── Main loop ─────────────────────────────────────────────────────────────────
PROCESSED=0
SKIPPED=0
WARN=0

for gene_dir in "$ROOT_DIR"/Inferred_Genetree_*/${GENE_PATTERN}; do
    [ -d "$gene_dir" ] || continue
    gene_name="$(basename "$gene_dir")"
    codeml_dir="${gene_dir}/codeml_input"

    echo "[PROC] $gene_name"

    # Skip if outputs exist and --force not set
    if [ -d "$codeml_dir" ] && [ "$FORCE" = false ]; then
        if [ -f "${codeml_dir}/cleaned.fasta" ] && [ -f "${gene_dir}/branch_table.txt" ]; then
            echo "  [SKIP] codeml_input already exists (use --force to overwrite)"
            (( SKIPPED++ )) || true
            continue
        fi
    fi

    # 1. Create codeml_input/ directory
    mkdir -p "$codeml_dir"

    # 2. Find NT FASTA (prefer *_NT* or *_aligned_NT* naming)
    fasta_file=""
    for candidate in \
        "${gene_dir}"/*_aligned_NT."${FASTA_EXT}" \
        "${gene_dir}"/*_NT."${FASTA_EXT}" \
        "${gene_dir}"/*."${FASTA_EXT}"
    do
        if [ -f "$candidate" ]; then
            fasta_file="$candidate"
            break
        fi
    done

    if [ -n "$fasta_file" ] && [ -f "$fasta_file" ]; then
        cp "$fasta_file" "${codeml_dir}/cleaned.fasta"
        # Clean headers: remove everything after '|', then remove all underscores
        sed -E -i.bak '/^>/ { s/\|.*//; s/_//g }' "${codeml_dir}/cleaned.fasta" \
            && rm -f "${codeml_dir}/cleaned.fasta.bak"
        echo "  cleaned.fasta: $(basename "$fasta_file")"
    else
        echo "  [WARN] No .$FASTA_EXT file found in $gene_dir" >&2
        (( WARN++ )) || true
    fi

    # 3. Find gene tree (.tre)
    tree_file=""
    for candidate in \
        "${gene_dir}"/*.tre \
        "${gene_dir}"/*_genetree.tre \
        "${gene_dir}"/*_gene_tree.tre
    do
        if [ -f "$candidate" ]; then
            tree_file="$candidate"
            break
        fi
    done

    if [ -n "$tree_file" ] && [ -f "$tree_file" ]; then
        cp "$tree_file" "${codeml_dir}/cleaned.tre"
        # Clean tree labels: remove everything after '|' in node labels, remove '_'
        sed -E -i.bak 's/\|[^),;:]*//g; s/_//g' "${codeml_dir}/cleaned.tre" \
            && rm -f "${codeml_dir}/cleaned.tre.bak"
        echo "  cleaned.tre: $(basename "$tree_file")"
    else
        echo "  [WARN] No .tre file found in $gene_dir" >&2
        (( WARN++ )) || true
    fi

    # 4. Determine species list for branch_table.txt
    if [ -n "$TARGET_SPECIES" ]; then
        # Explicit list supplied by the user
        BRANCH_TABLE_CONTENT="$FIXED_BRANCH_TABLE_CONTENT"
        echo "  branch_table.txt: $TARGET_SPECIES"
    elif [ -f "${codeml_dir}/cleaned.fasta" ]; then
        # Auto-detect: collect every unique >Header from the cleaned FASTA.
        # Headers are already stripped of '|' and '_' by the sed step above.
        BRANCH_TABLE_CONTENT=""
        while IFS= read -r header; do
            sp="${header#>}"
            sp="${sp%% *}"   # first word only (before any space or tab)
            BRANCH_TABLE_CONTENT="${BRANCH_TABLE_CONTENT}${sp}"$'\n'
        done < <(grep '^>' "${codeml_dir}/cleaned.fasta" | sort -u)
        sp_list="$(printf '%s' "$BRANCH_TABLE_CONTENT" | tr '\n' ',' | sed 's/,$//')"
        echo "  branch_table.txt: [auto] $sp_list"
    else
        echo "  [WARN] Cannot write branch_table.txt — no FASTA and no --target-species" >&2
        (( WARN++ )) || true
        (( PROCESSED++ )) || true
        continue
    fi

    # 5. Write branch_table.txt
    printf '%s' "$BRANCH_TABLE_CONTENT" > "${gene_dir}/branch_table.txt"

    (( PROCESSED++ )) || true
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "======================================================"
echo "  Done  ($(date))"
echo "  Processed : $PROCESSED"
echo "  Skipped   : $SKIPPED"
echo "  Warnings  : $WARN"
echo "======================================================"
echo ""
echo "Next step:"
echo "  bash codeml_set_up.sh \\"
echo "      --root-dir       $ROOT_DIR \\"
echo "      --vespa-script   /path/to/VESPA/vespa.py"
