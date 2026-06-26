#!/usr/bin/env bash
# ==============================================================================
# codeml_set_up.sh — Run VESPA codeml_setup for each gene directory
# ==============================================================================
#
# Runs VESPA's `codeml_setup` command on every gene directory (Inferred_Genetree_*/
# subdirectories) that contains a codeml_input/ directory and a branch_table.txt
# (prepared by branch_table.sh).
#
# VESPA codeml_setup generates a complete codeml workspace per orthogroup:
#   Codeml_Setup_codeml_input/
#   ├── cleaned/                  (site models — no branch labelling)
#   │   ├── m0/Omega0/
#   │   ├── m1Neutral/Omega0/
#   │   ├── m2Selection/Omega0/
#   │   ├── m7/Omega0/
#   │   └── m8/Omega0/
#   └── cleaned_<species>/        (branch model — one dir per foreground species)
#       └── modelA/Omega0/
#
# The correction_codeml.sh script (next step) renames directories and rewrites
# codeml.ctl files with the correct initial omega value (0.5).
#
# VESPA requires Python 2.7. Supply via --python-exec or use --docker-image /
# --singularity-image (see gene_tree.sh and Dockerfile.vespa).
#
# USAGE
#   bash codeml_set_up.sh \
#       --root-dir     vespa_out/ \
#       --vespa-script /path/to/VESPA/vespa.py
#
#   # Docker:
#   bash codeml_set_up.sh \
#       --root-dir     vespa_out/ \
#       --docker-image vespa:latest
#
#   # Custom gene pattern:
#   bash codeml_set_up.sh \
#       --root-dir     vespa_out/ \
#       --docker-image vespa:latest \
#       --gene-pattern "Am*"
# ==============================================================================
set -euo pipefail
shopt -s nullglob

# ── Defaults ──────────────────────────────────────────────────────────────────
ROOT_DIR=""
VESPA_SCRIPT=""
PYTHON_EXEC=""
DOCKER_IMAGE=""
SINGULARITY_IMAGE=""
GENE_PATTERN="*"

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
codeml_set_up.sh — Run VESPA codeml_setup on each gene directory

USAGE
  bash codeml_set_up.sh \
      --root-dir     <dir> \
      { --docker-image <image> | --singularity-image <sif> |
        ( --vespa-script <path> [--python-exec <python2>] ) }

REQUIRED
  -r | --root-dir       DIR     Root directory with Inferred_Genetree_*/ dirs
                                (output from gene_tree.sh + branch_table.sh)

VESPA RUNTIME (choose one)
  --docker-image      IMAGE     Docker image (e.g. vespa:latest)
  --singularity-image FILE      Singularity/Apptainer .sif file
  --vespa-script      FILE      Path to vespa.py (requires Python 2.7)
  --python-exec       EXEC      Python 2.7 interpreter (default: python2)

OPTIONAL
  -p | --gene-pattern PAT       Gene subdirectory glob (default: * = all)
                                Example: "Am*"
  -h | --help                   Show this message

PREREQUISITES
  - gene_tree.sh must have completed (Inferred_Genetree_*/ dirs exist)
  - branch_table.sh must have completed (codeml_input/ and branch_table.txt exist)

OUTPUT (per gene directory)
  Inferred_Genetree_<id>/<gene>/Codeml_Setup_codeml_input/
    cleaned/                       site models (m0, m1Neutral, ...)
    cleaned_<foreground_species>/  branch model (modelA/Omega0/)
EOF
}

[ $# -eq 0 ] && { show_help; exit 0; }

# ── Argument parsing ──────────────────────────────────────────────────────────
while [ $# -gt 0 ]; do
    case "$1" in
        -r|--root-dir)          ROOT_DIR="$2";           shift 2 ;;
        --vespa-script)         VESPA_SCRIPT="$2";       shift 2 ;;
        --python-exec)          PYTHON_EXEC="$2";        shift 2 ;;
        --docker-image)         DOCKER_IMAGE="$2";       shift 2 ;;
        --singularity-image)    SINGULARITY_IMAGE="$2";  shift 2 ;;
        -p|--gene-pattern)      GENE_PATTERN="$2";       shift 2 ;;
        -h|--help)              show_help; exit 0 ;;
        *) echo "ERROR: Unknown argument: '$1'" >&2; show_help >&2; exit 1 ;;
    esac
done

die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$ROOT_DIR" ] || die "--root-dir is required"
[ -d "$ROOT_DIR" ] || die "Root directory not found: $ROOT_DIR"
ROOT_DIR="$(cd "$ROOT_DIR" && pwd)"

# Resolve VESPA runtime
if [ -n "$DOCKER_IMAGE" ] && [ -n "$SINGULARITY_IMAGE" ]; then
    die "Specify either --docker-image OR --singularity-image, not both."
fi

VESPA_MODE="local"
if [ -n "$DOCKER_IMAGE" ]; then
    command -v docker &>/dev/null || die "Docker not found."
    VESPA_MODE="docker"
elif [ -n "$SINGULARITY_IMAGE" ]; then
    [ -f "$SINGULARITY_IMAGE" ] || die "Singularity image not found: $SINGULARITY_IMAGE"
    SIF_CMD=""
    command -v singularity &>/dev/null && SIF_CMD="singularity"
    command -v apptainer  &>/dev/null && SIF_CMD="apptainer"
    [ -n "$SIF_CMD" ] || die "Neither singularity nor apptainer found."
    VESPA_MODE="singularity"
else
    [ -n "$VESPA_SCRIPT" ] || die "No VESPA runtime specified."
    [ -f "$VESPA_SCRIPT" ] || die "VESPA script not found: $VESPA_SCRIPT"
    VESPA_SCRIPT="$(cd "$(dirname "$VESPA_SCRIPT")" && pwd)/$(basename "$VESPA_SCRIPT")"
    if [ -z "$PYTHON_EXEC" ]; then
        if command -v python2 &>/dev/null; then
            PYTHON_EXEC=python2
        elif python -c 'import sys; sys.exit(0 if sys.version_info[0]==2 else 1)' 2>/dev/null; then
            PYTHON_EXEC=python
        else
            echo "WARNING: python2 not found — VESPA may fail with Python 3." >&2
            PYTHON_EXEC=python3
        fi
    fi
fi

# ── VESPA run helper (called from gene directory) ─────────────────────────────
run_vespa_codeml_setup() {
    local gene_dir="$1"
    local branch_file="${gene_dir}/branch_table.txt"

    case "$VESPA_MODE" in
        docker)
            docker run --rm \
                -v "${gene_dir}":/data \
                -w /data \
                "$DOCKER_IMAGE" \
                codeml_setup -input=codeml_input -branch_file=/data/branch_table.txt
            ;;
        singularity)
            "$SIF_CMD" exec \
                --bind "${gene_dir}":/data \
                "$SINGULARITY_IMAGE" \
                /usr/local/bin/vespa \
                codeml_setup -input=codeml_input -branch_file=branch_table.txt
            ;;
        local)
            (cd "$gene_dir" && "$PYTHON_EXEC" "$VESPA_SCRIPT" \
                codeml_setup -input=codeml_input -branch_file=branch_table.txt)
            ;;
    esac
}

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  VESPA codeml_setup"
echo "======================================================"
echo "  Root dir     : $ROOT_DIR"
echo "  Gene pattern : $GENE_PATTERN"
echo "  VESPA mode   : $VESPA_MODE"
echo "======================================================"
echo ""

# ── Main loop ─────────────────────────────────────────────────────────────────
PROCESSED=0; SKIPPED=0; FAILED=0

for gene_dir in "$ROOT_DIR"/Inferred_Genetree_*/${GENE_PATTERN}; do
    [ -d "$gene_dir" ] || continue
    gene_name="$(basename "$gene_dir")"

    # Check prerequisites
    if [ ! -f "${gene_dir}/branch_table.txt" ]; then
        echo "[WARN] $gene_name: missing branch_table.txt — run branch_table.sh first" >&2
        (( FAILED++ )) || true
        continue
    fi
    if [ ! -d "${gene_dir}/codeml_input" ]; then
        echo "[WARN] $gene_name: missing codeml_input/ — run branch_table.sh first" >&2
        (( FAILED++ )) || true
        continue
    fi

    # Skip if codeml_setup already ran
    if [ -d "${gene_dir}/Codeml_Setup_codeml_input" ]; then
        echo "[SKIP] $gene_name — Codeml_Setup_codeml_input already exists"
        (( SKIPPED++ )) || true
        continue
    fi

    echo "[RUN ] $gene_name"
    if run_vespa_codeml_setup "$gene_dir"; then
        (( PROCESSED++ )) || true
    else
        echo "  [WARN] codeml_setup failed for $gene_name" >&2
        (( FAILED++ )) || true
    fi
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "======================================================"
echo "  Done  ($(date))"
echo "  Processed : $PROCESSED"
echo "  Skipped   : $SKIPPED"
echo "  Failed    : $FAILED"
echo "======================================================"
echo ""
echo "Next step:"
echo "  bash correction_codeml.sh \\"
echo "      --root-dir $ROOT_DIR \\"
echo "      --mode all"
