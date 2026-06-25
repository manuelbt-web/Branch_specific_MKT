#!/usr/bin/env bash
# ==============================================================================
# run_codeml.sbatch — SLURM array job: run codeml for branch model + m0
# ==============================================================================
#
# By default this script submits ALL branch model jobs and ALL m0 site model
# jobs in a single SLURM array.  Additional site models (m1Neutral, m2Selection,
# m7, m8) can be added with --site-models.
#
# PREREQUISITES
#   correction_codeml.sh has run and created both taskfarm files:
#     <root>/codeml_taskfarm_fullpath.sh              (branch model)
#     <root>/codeml_taskfarm_fullpath_site_model.sh   (site models)
#
# SUBMISSION (from the run_codeml/ directory):
#
#   # Step 1 — count total tasks (branch + m0):
#   BRANCH=../codeml_inputs/codeml_taskfarm_fullpath.sh
#   SITE=../codeml_inputs/codeml_taskfarm_fullpath_site_model.sh
#   N=$(bash run_codeml.sbatch --count \
#         --taskfarm-branch "$BRANCH" --taskfarm-site "$SITE")
#
#   # Step 2 — submit:
#   sbatch --array=0-$((N-1))%15 run_codeml.sbatch \
#          --taskfarm-branch "$BRANCH" --taskfarm-site "$SITE"
#
#   # Branch model only (no site models):
#   N=$(wc -l < "$BRANCH")
#   sbatch --array=0-$((N-1))%15 run_codeml.sbatch \
#          --taskfarm-branch "$BRANCH" --site-models none
#
#   # Branch + multiple site models (m0, m7, m8):
#   N=$(bash run_codeml.sbatch --count \
#         --taskfarm-branch "$BRANCH" --taskfarm-site "$SITE" \
#         --site-models m0,m7,m8)
#   sbatch --array=0-$((N-1))%15 run_codeml.sbatch \
#          --taskfarm-branch "$BRANCH" --taskfarm-site "$SITE" \
#          --site-models m0,m7,m8
#
# ARGUMENTS (passed at submission and forwarded to each running task)
#   --taskfarm-branch FILE   Branch model taskfarm (default: codeml_taskfarm_fullpath.sh)
#   --taskfarm-site   FILE   Site model taskfarm   (default: codeml_taskfarm_fullpath_site_model.sh)
#   --site-models     LIST   Comma-separated site models to include (default: m0)
#                            Set to "none" to skip site models entirely.
#   --count                  Print total task count and exit (does NOT submit)
#   -h | --help              Show this message
#
# SLURM ARRAY LAYOUT
#   Tasks 0 to (N_branch-1)                → branch model jobs
#   Tasks N_branch to (N_branch+N_site-1)  → site model jobs (filtered by --site-models)
# ==============================================================================
#SBATCH --job-name=codeml
#SBATCH --mem=10G                    # [CLUSTER] adjust per gene / model complexity
#SBATCH --time=48:00:00             # [CLUSTER] branch model can be slow on long genes
#SBATCH --partition=cpu-dedicated   # [CLUSTER] set your partition name
##SBATCH --account=myaccount        # [CLUSTER] uncomment and set your account if required
#SBATCH --array=0-0%15              # override at submission: --array=0-$((N-1))%15
#SBATCH --output=codeml_%A_%a.out
#SBATCH --error=codeml_%A_%a.err
# ==============================================================================

set -euo pipefail

# ── Load cluster modules ───────────────────────────────────────────────────────
# [CLUSTER] Adjust to match your cluster's module system.
if command -v module &>/dev/null; then
    module load paml/4.10.7 2>/dev/null || true
fi

# ── Defaults ──────────────────────────────────────────────────────────────────
TASKFARM_BRANCH="codeml_taskfarm_fullpath.sh"
TASKFARM_SITE="codeml_taskfarm_fullpath_site_model.sh"
SITE_MODELS="m0"    # default: branch jobs + m0 site model
COUNT_ONLY=false

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
run_codeml.sbatch — SLURM array job: branch model + m0 (and optionally other site models)

USAGE (count tasks first, then submit):
  BRANCH=../codeml_inputs/codeml_taskfarm_fullpath.sh
  SITE=../codeml_inputs/codeml_taskfarm_fullpath_site_model.sh

  N=$(bash run_codeml.sbatch --count --taskfarm-branch "$BRANCH" --taskfarm-site "$SITE")
  sbatch --array=0-$((N-1))%15 run_codeml.sbatch \
         --taskfarm-branch "$BRANCH" --taskfarm-site "$SITE"

ARGUMENTS
  --taskfarm-branch FILE   Branch model taskfarm
                           (default: codeml_taskfarm_fullpath.sh)
  --taskfarm-site   FILE   Site model taskfarm
                           (default: codeml_taskfarm_fullpath_site_model.sh)
  --site-models     LIST   Comma-separated site models to include (default: m0)
                           Available: m0, m1Neutral, m2Selection, m7, m8
                           Use "none" to run branch model jobs only.
  --count                  Print total task count and exit (for --array sizing)
  -h | --help              Show this message

SLURM ARRAY LAYOUT
  Tasks 0 .. N_branch-1              → branch model (free-ratios, one ω per branch)
  Tasks N_branch .. N_branch+N_site-1 → site models (filtered by --site-models)
EOF
}

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --taskfarm-branch) TASKFARM_BRANCH="$2"; shift 2 ;;
        --taskfarm-site)   TASKFARM_SITE="$2";   shift 2 ;;
        --site-models)     SITE_MODELS="$2";     shift 2 ;;
        --count)           COUNT_ONLY=true;       shift ;;
        -h|--help)         show_help; exit 0 ;;
        *) echo "ERROR: Unknown argument: '$1'" >&2; show_help >&2; exit 1 ;;
    esac
done

# ── Load branch commands (all non-comment, non-blank lines) ───────────────────
if [ ! -f "$TASKFARM_BRANCH" ]; then
    echo "ERROR: branch taskfarm not found: $TASKFARM_BRANCH" >&2
    echo "       Run correction_codeml.sh --mode branch (or --mode all) first." >&2
    exit 1
fi
mapfile -t BRANCH_CMDS < <(grep -v '^#' "$TASKFARM_BRANCH" | grep -v '^[[:space:]]*$')
N_BRANCH=${#BRANCH_CMDS[@]}

# ── Load site commands filtered by --site-models ──────────────────────────────
N_SITE=0
SITE_CMDS=()

if [[ "$SITE_MODELS" != "none" && -n "$SITE_MODELS" ]]; then
    if [ ! -f "$TASKFARM_SITE" ]; then
        echo "ERROR: site taskfarm not found: $TASKFARM_SITE" >&2
        echo "       Run correction_codeml.sh --mode site (or --mode all) first." >&2
        echo "       Or use --site-models none to skip site models." >&2
        exit 1
    fi
    # Build an ERE pattern matching any of the requested model directory names.
    # The taskfarm lines look like:  cd /path/.../cleaned/m0/Omega0_5; codeml
    # so we match  /<model_name>/
    MODEL_PATTERN=$(echo "$SITE_MODELS" | sed 's/,/|/g')
    mapfile -t SITE_CMDS < <(
        grep -v '^#' "$TASKFARM_SITE" \
        | grep -v '^[[:space:]]*$' \
        | grep -E "/(${MODEL_PATTERN})/"
    )
    N_SITE=${#SITE_CMDS[@]}
fi

N_TOTAL=$(( N_BRANCH + N_SITE ))

# ── Count-only mode — print total and exit (used to size the SLURM array) ─────
if $COUNT_ONLY; then
    echo "$N_TOTAL"
    exit 0
fi

# ── Sanity check ──────────────────────────────────────────────────────────────
if (( N_TOTAL == 0 )); then
    echo "ERROR: No tasks found (branch=$N_BRANCH, site=$N_SITE)." >&2
    echo "       Check that taskfarm files exist and are non-empty." >&2
    exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

if (( TASK_ID >= N_TOTAL )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID=$TASK_ID is out of range (total=$N_TOTAL)." >&2
    exit 1
fi

# ── Select command for this task ──────────────────────────────────────────────
if (( TASK_ID < N_BRANCH )); then
    CMD="${BRANCH_CMDS[$TASK_ID]}"
    MODEL_TYPE="branch"
else
    IDX=$(( TASK_ID - N_BRANCH ))
    CMD="${SITE_CMDS[$IDX]}"
    MODEL_TYPE="site (${SITE_MODELS})"
fi

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  run_codeml.sbatch"
echo "  Job   : ${SLURM_JOB_ID:-local}  |  Task: $TASK_ID / $((N_TOTAL - 1))"
echo "  Type  : $MODEL_TYPE"
echo "  CMD   : $CMD"
echo "  Date  : $(date)"
echo "======================================================"

# ── Run codeml ────────────────────────────────────────────────────────────────
eval "$CMD"
