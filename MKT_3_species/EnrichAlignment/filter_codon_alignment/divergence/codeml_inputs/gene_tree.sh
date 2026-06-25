#!/usr/bin/env bash
# ==============================================================================
# gene_tree.sh — Generate per-gene trees with VESPA infer_genetree
# ==============================================================================
#
# Runs VESPA's infer_genetree on every codon-alignment FASTA in <input-dir>.
# VESPA maps alignment sequence headers to taxa in the species tree topology,
# prunes absent species, and writes one gene tree (.tre) per alignment.
#
# VESPA requires Python 2.7.  Supply the interpreter via --python-exec, or
# run VESPA inside a container (see Dockerfile.vespa) via --docker-image or
# --singularity-image.
#
# SPECIES TREE FORMAT
#   Pure Newick topology — species names only, no branch lengths, no support:
#     ((Aegilops_tauschii,Triticum_urartu),Amblyopyrum_muticum);
#   Generate from OrthoFinder SpeciesTree_rooted.txt by stripping numbers:
#     sed 's/:[0-9.eE+-]*//g; s/[0-9]//g' SpeciesTree_rooted.txt \
#         > SpeciesTree_topology.txt
#
# VESPA INSTALLATION
#   See Dockerfile.vespa for the recommended container approach.
#   Local install: git clone https://github.com/aewebb80/VESPA.git
#                  conda activate vespa_py27
#                  # (Python 2.7.18 + dendropy==4.4.0)
#
# USAGE
#   # Local Python 2.7:
#   bash gene_tree.sh \
#       --input-dir    split_out/divergence/ \
#       --species-tree SpeciesTree_topology.txt \
#       --vespa-script /path/to/VESPA/vespa.py \
#       --python-exec  python2
#
#   # Docker (recommended):
#   bash gene_tree.sh \
#       --input-dir    split_out/divergence/ \
#       --species-tree SpeciesTree_topology.txt \
#       --docker-image vespa:latest
#
#   # Singularity / Apptainer (HPC without Docker):
#   bash gene_tree.sh \
#       --input-dir    split_out/divergence/ \
#       --species-tree SpeciesTree_topology.txt \
#       --singularity-image /path/to/vespa.sif
#
# OUTPUT
#   One gene tree per FASTA, written next to the input file as <gene>.tre.
#   Already-present clean tree files are skipped (idempotent).
# ==============================================================================
set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT_DIR=""
SPECIES_TREE=""
VESPA_SCRIPT=""
PYTHON_EXEC=""
DOCKER_IMAGE=""
SINGULARITY_IMAGE=""
EXT="fasta"
ALLOW_PARALOGS=false
FORCE=false

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
gene_tree.sh — Run VESPA infer_genetree on codon-alignment FASTAs

USAGE
  bash gene_tree.sh \
      --input-dir    <dir>  \
      --species-tree <file> \
      { --docker-image <image> | --singularity-image <sif> |
        ( --vespa-script <path> [--python-exec <python2>] ) }

REQUIRED
  -i | --input-dir       DIR    Directory containing *.<ext> FASTA files
  -s | --species-tree    FILE   Species topology Newick (names only, no lengths)
                                Example: ((Sp_A,Sp_B),Sp_C);

VESPA RUNTIME (choose one)
  --docker-image      IMAGE     Docker image (e.g. vespa:latest)
  --singularity-image FILE      Singularity/Apptainer .sif file
  --vespa-script      FILE      Path to vespa.py (requires Python 2.7 in PATH)
  --python-exec       EXEC      Python 2.7 executable (default: python2)

OPTIONAL
  -e | --ext          EXT       FASTA file extension (default: fasta)
  -a | --allow-paralogs         Pass -allow_paralogs to VESPA
  -f | --force                  Re-run even if gene tree already exists
  -h | --help                   Show this message

SPECIES TREE FORMAT
  Newick topology with species names that exactly match the alignment
  sequence headers (after '|' splitting and '_' removal by branch_table.sh).
  No branch lengths, no bootstrap values.
  Example: ((Aegilops_tauschii,Triticum_urartu),Amblyopyrum_muticum);

  Generate from OrthoFinder output:
    sed 's/:[0-9.eE+-]*//g; s/[0-9]//g' SpeciesTree_rooted.txt \
        > SpeciesTree_topology.txt

OUTPUT
  <input_dir>/<gene>.tre   — pruned gene tree (one per alignment)
EOF
}

[ $# -eq 0 ] && { show_help; exit 0; }

# ── Argument parsing ──────────────────────────────────────────────────────────
while [ $# -gt 0 ]; do
    case "$1" in
        -i|--input-dir)         INPUT_DIR="$2";          shift 2 ;;
        -s|--species-tree)      SPECIES_TREE="$2";       shift 2 ;;
        --vespa-script)         VESPA_SCRIPT="$2";       shift 2 ;;
        --python-exec)          PYTHON_EXEC="$2";        shift 2 ;;
        --docker-image)         DOCKER_IMAGE="$2";       shift 2 ;;
        --singularity-image)    SINGULARITY_IMAGE="$2";  shift 2 ;;
        -e|--ext)               EXT="$2";                shift 2 ;;
        -a|--allow-paralogs)    ALLOW_PARALOGS=true;     shift ;;
        -f|--force)             FORCE=true;              shift ;;
        -h|--help)              show_help; exit 0 ;;
        *) echo "ERROR: Unknown argument: '$1'" >&2; show_help >&2; exit 1 ;;
    esac
done

# ── Validation ────────────────────────────────────────────────────────────────
die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$INPUT_DIR"    ] || die "--input-dir is required"
[ -n "$SPECIES_TREE" ] || die "--species-tree is required"
[ -d "$INPUT_DIR"    ] || die "Input directory not found: $INPUT_DIR"
[ -f "$SPECIES_TREE" ] || die "Species tree file not found: $SPECIES_TREE"

# Use absolute paths so cd in the loop doesn't break relative references
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
SPECIES_TREE="$(cd "$(dirname "$SPECIES_TREE")" && pwd)/$(basename "$SPECIES_TREE")"

# Determine how VESPA will be called
if [ -n "$DOCKER_IMAGE" ] && [ -n "$SINGULARITY_IMAGE" ]; then
    die "Specify either --docker-image OR --singularity-image, not both."
fi

VESPA_MODE="local"
if [ -n "$DOCKER_IMAGE" ]; then
    command -v docker &>/dev/null || die "Docker not found. Install Docker or use --singularity-image."
    VESPA_MODE="docker"
elif [ -n "$SINGULARITY_IMAGE" ]; then
    [ -f "$SINGULARITY_IMAGE" ] || die "Singularity image not found: $SINGULARITY_IMAGE"
    SIF_CMD=""
    command -v singularity &>/dev/null && SIF_CMD="singularity"
    command -v apptainer  &>/dev/null && SIF_CMD="apptainer"
    [ -n "$SIF_CMD" ] || die "Neither singularity nor apptainer found in PATH."
    VESPA_MODE="singularity"
else
    [ -n "$VESPA_SCRIPT" ] || die "No VESPA runtime specified. Use --docker-image, --singularity-image, or --vespa-script."
    [ -f "$VESPA_SCRIPT" ] || die "VESPA script not found: $VESPA_SCRIPT"
    # Resolve Python 2.7 interpreter
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
    echo "Python interpreter: $PYTHON_EXEC ($(${PYTHON_EXEC} --version 2>&1))"
fi

# ── Build VESPA run function ──────────────────────────────────────────────────
# run_vespa CMD ARGS...  (must be called from the directory containing the input file)
run_vespa() {
    local work_dir="$1"; shift
    case "$VESPA_MODE" in
        docker)
            docker run --rm \
                -v "${work_dir}":/data \
                -v "${SPECIES_TREE%/*}":/speciestree_host:ro \
                -w /data \
                "$DOCKER_IMAGE" \
                "$@" -species_tree="/speciestree_host/$(basename "$SPECIES_TREE")"
            ;;
        singularity)
            "$SIF_CMD" exec \
                --bind "${work_dir}":/data \
                --bind "${SPECIES_TREE%/*}":/speciestree_host:ro \
                "$SINGULARITY_IMAGE" \
                /usr/local/bin/vespa "$@" -species_tree="/speciestree_host/$(basename "$SPECIES_TREE")"
            ;;
        local)
            (cd "$work_dir" && "$PYTHON_EXEC" "$VESPA_SCRIPT" "$@" -species_tree="$SPECIES_TREE")
            ;;
    esac
}

# ── Collect FASTA files ───────────────────────────────────────────────────────
shopt -s nullglob
FASTA_FILES=("$INPUT_DIR"/*."$EXT")
shopt -u nullglob
N=${#FASTA_FILES[@]}
[ "$N" -gt 0 ] || die "No .$EXT files found in $INPUT_DIR"

PARALOGS_ARG=""
[ "$ALLOW_PARALOGS" = true ] && PARALOGS_ARG="-allow_paralogs"

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  VESPA infer_genetree"
echo "======================================================"
echo "  Input dir    : $INPUT_DIR"
echo "  Species tree : $SPECIES_TREE"
echo "  Extension    : .$EXT"
echo "  Files found  : $N"
echo "  VESPA mode   : $VESPA_MODE"
[ "$ALLOW_PARALOGS" = true ] && echo "  Allow paralogs: yes"
echo "======================================================"
echo ""

# ── Main loop ─────────────────────────────────────────────────────────────────
declare -a RAN=() SKIPPED=() FAILED=()

for fasta in "${FASTA_FILES[@]}"; do
    base="${fasta%.*}"          # full path without extension
    stem="$(basename "$base")"  # filename without extension
    dir="$(dirname  "$fasta")"  # containing directory

    # Check for existing clean gene tree (no | or _ remaining after VESPA cleaning)
    existing_tree=""
    if [ "$FORCE" = false ]; then
        for candidate in \
            "${base}.tre" \
            "${base}_gene_tree.tre" \
            "${base}_genetree.tre"
        do
            if [ -s "$candidate" ] && ! grep -q '[|_]' "$candidate"; then
                existing_tree="$candidate"
                break
            fi
        done
    fi

    if [ -n "$existing_tree" ]; then
        echo "[SKIP] $stem — existing tree: $(basename "$existing_tree")"
        SKIPPED+=("$stem")
        continue
    fi

    echo "[RUN ] $stem"
    if run_vespa "$dir" infer_genetree \
            -input="${stem}.${EXT}" \
            $PARALOGS_ARG; then
        RAN+=("$stem")
    else
        echo "  [WARN] infer_genetree failed for $stem" >&2
        FAILED+=("$stem")
    fi
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "======================================================"
echo "  Done  ($(date))"
echo "  Ran     : ${#RAN[@]} / $N"
echo "  Skipped : ${#SKIPPED[@]} / $N"
echo "  Failed  : ${#FAILED[@]} / $N"
if [ ${#FAILED[@]} -gt 0 ]; then
    echo "  Failed genes:"
    printf '    %s\n' "${FAILED[@]}"
fi
echo "======================================================"
