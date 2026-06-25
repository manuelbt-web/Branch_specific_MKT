#!/usr/bin/env bash
# ==============================================================================
# branch_table.sh — Prepare codeml input directories after VESPA infer_genetree
# ==============================================================================
#
# For each gene directory produced by VESPA (Inferred_Genetree_*/), this script:
#   1. Writes branch_table.txt listing the foreground (target) species for the
#      VESPA branch model (the species whose lineage is under selection).
#   2. Creates a codeml_input/ subdirectory.
#   3. Copies the NT FASTA and cleans sequence headers:
#        - removes everything after '|' (species|gene_id → species)
#        - removes underscores (Aegilops_tauschii → Aegilopstauschii)
#   4. Copies the gene tree (.tre) and cleans its labels the same way.
#
# VESPA codeml_setup (next step: codeml_set_up.sh) reads codeml_input/ and
# branch_table.txt to generate the full codeml workspace.
#
# USAGE
#   bash branch_table.sh \
#       --root-dir       vespa_out/ \
#       --target-species "Aegilopsspeltoides,Aegilopsmutica"
#
#   # Custom gene-directory pattern (default: match all subdirs):
#   bash branch_table.sh \
#       --root-dir       vespa_out/ \
#       --target-species "Aegilopsspeltoides" \
#       --gene-pattern   "Am*"
#
# NOTE ON SPECIES NAMES IN BRANCH TABLE
#   The names must match the cleaned headers (no pipes, no underscores).
#   Original header:  >Aegilops_speltoides|TraesCS1A02G000001
#   Cleaned name:      Aegilopsspeltoides
#   (The header is split at '|', then all '_' are removed.)
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
TARGET_SPECIES=""      # comma-separated, cleaned names (no underscores)
GENE_PATTERN="*"       # glob pattern for the gene-level subdirectory
FASTA_EXT="fasta"
FORCE=false

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
branch_table.sh — Prepare VESPA codeml_input directories

USAGE
  bash branch_table.sh \
      --root-dir       <dir> \
      --target-species "<Sp1>,<Sp2>,..."

REQUIRED
  -r | --root-dir       DIR     Root directory containing Inferred_Genetree_*/ dirs
                                (output from gene_tree.sh)
  -t | --target-species LIST    Comma-separated foreground species for the branch
                                model.  Names must be CLEANED (no underscores, no
                                pipe characters).
                                Example: "Aegilopsspeltoides,Aegilopsmutica"

OPTIONAL
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
  List the species in the same order as their appearance in the alignment.

OUTPUT (per gene directory)
  <gene_dir>/branch_table.txt        — foreground species list
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

[ -n "$ROOT_DIR"        ] || die "--root-dir is required"
[ -n "$TARGET_SPECIES"  ] || die "--target-species is required (cleaned names, comma-separated)"
[ -d "$ROOT_DIR"        ] || die "Root directory not found: $ROOT_DIR"

ROOT_DIR="$(cd "$ROOT_DIR" && pwd)"

# Parse species list into newline-separated block for branch_table.txt
BRANCH_TABLE_CONTENT=""
IFS=',' read -ra SPECIES_ARRAY <<< "$TARGET_SPECIES"
for sp in "${SPECIES_ARRAY[@]}"; do
    sp="${sp// /}"  # strip any spaces
    BRANCH_TABLE_CONTENT="${BRANCH_TABLE_CONTENT}${sp}"$'\n'
done

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  branch_table.sh — Prepare codeml input directories"
echo "======================================================"
echo "  Root dir         : $ROOT_DIR"
echo "  Target species   : $TARGET_SPECIES"
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

    # 1. Write branch_table.txt
    printf '%s' "$BRANCH_TABLE_CONTENT" > "${gene_dir}/branch_table.txt"
    echo "  branch_table.txt: $(echo "$TARGET_SPECIES" | tr ',' ' ')"

    # 2. Create codeml_input/ directory
    mkdir -p "$codeml_dir"

    # 3. Find NT FASTA (prefer *_NT* or *_aligned_NT* naming)
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

    # 4. Find gene tree (.tre)
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
