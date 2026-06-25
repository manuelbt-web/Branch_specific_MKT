#!/usr/bin/env bash
# ==============================================================================
# correction_codeml.sh — Fix VESPA codeml_setup directory structure and ctl files
# ==============================================================================
#
# VESPA's codeml_setup creates model directories with two issues that require
# correction before running codeml:
#
#   1. BRANCH MODEL (cleaned_<species>/modelA/Omega0/):
#        - VESPA names the branch model directory 'modelA'; rename → 'model_branch'
#        - VESPA creates an 'Omega0' directory (initial omega = 0); rename → 'Omega0_5'
#          and rewrite codeml.ctl with omega = 0.5 (better starting point for optimization)
#        - Remove all other model/omega directories (only keep model_branch/Omega0_5)
#        - model = 2, NSsites = 0  (free-ratios branch model in PAML)
#
#   2. SITE MODELS (cleaned/m0/Omega0/, cleaned/m1Neutral/Omega0/, ...):
#        - VESPA creates Omega0 (initial omega=0); rename → Omega0_5 and use omega=0.5
#        - Remove omega directories other than Omega0 (keep one starting value)
#        - Remove model directories not in --models list
#        - NOTE: For PAML site models, the correct parameters are:
#            model = 0  (applies to ALL site models)
#            NSsites = <N>  (0=M0, 1=M1a, 2=M2a, 7=M7, 8=M8)
#          The original VESPA codeml_setup may set wrong 'model' values;
#          this script enforces the correct PAML parameterization.
#
# After running this script, a taskfarm file is generated listing one
# 'cd <dir>; codeml' command per codeml workspace, for use with
# run_codeml.sbatch.
#
# REPLACES: correction_branch_model.sh  +  correction_model_m0.sh
#
# USAGE
#   bash correction_codeml.sh \
#       --root-dir    vespa_out/ \
#       --mode        all
#
#   # Branch model only:
#   bash correction_codeml.sh --root-dir vespa_out/ --mode branch
#
#   # Site models only, keep only m0 and m7:
#   bash correction_codeml.sh --root-dir vespa_out/ --mode site \
#       --models m0,m7
# ==============================================================================
set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
ROOT_DIR=""
MODE="all"                                           # branch | site | all
MODELS="m0,m1Neutral,m2Selection,m7,m8"             # site models to keep
GENE_PATTERN="*"
TASKFARM_BRANCH=""     # default: <root>/codeml_taskfarm_fullpath.sh
TASKFARM_SITE=""       # default: <root>/codeml_taskfarm_fullpath_site_model.sh

# ── Help ──────────────────────────────────────────────────────────────────────
show_help() {
cat <<'EOF'
correction_codeml.sh — Fix VESPA codeml_setup directories and codeml.ctl files

USAGE
  bash correction_codeml.sh \
      --root-dir <dir> \
      [--mode branch|site|all]

REQUIRED
  -r | --root-dir    DIR    Root directory with Inferred_Genetree_*/ dirs
                            (output from codeml_set_up.sh)

OPTIONAL
  -m | --mode        MODE   Which model type to correct (default: all)
                              branch  — correct branch model workspaces only
                              site    — correct site model workspaces only
                              all     — correct both (default)
  --models           LIST   Comma-separated site models to keep (default: m0,m1Neutral,m2Selection,m7,m8)
                            Models not in this list are deleted.
  -p | --gene-pattern PAT   Gene subdirectory pattern (default: *)
  --taskfarm-branch  FILE   Path for branch model taskfarm (default: <root>/codeml_taskfarm_fullpath.sh)
  --taskfarm-site    FILE   Path for site model taskfarm   (default: <root>/codeml_taskfarm_fullpath_site_model.sh)
  -h | --help               Show this message

CODEML PARAMETERS WRITTEN
  Branch model  : model = 2, NSsites = 0, fix_omega = 0, omega = 0.5
  Site models   : model = 0, NSsites = <N>, fix_omega = 0, omega = 0.5
                  where N = 0 (M0), 1 (M1a), 2 (M2a), 7 (M7), 8 (M8)

OUTPUT
  Corrected directories under Codeml_Setup_codeml_input/:
    cleaned_<species>/model_branch/Omega0_5/codeml.ctl  (branch model)
    cleaned/<model>/Omega0_5/codeml.ctl                 (site models)
  Taskfarm files:
    <root>/codeml_taskfarm_fullpath.sh
    <root>/codeml_taskfarm_fullpath_site_model.sh
EOF
}

[ $# -eq 0 ] && { show_help; exit 0; }

# ── Argument parsing ──────────────────────────────────────────────────────────
while [ $# -gt 0 ]; do
    case "$1" in
        -r|--root-dir)        ROOT_DIR="$2";          shift 2 ;;
        -m|--mode)            MODE="$2";              shift 2 ;;
        --models)             MODELS="$2";            shift 2 ;;
        -p|--gene-pattern)    GENE_PATTERN="$2";      shift 2 ;;
        --taskfarm-branch)    TASKFARM_BRANCH="$2";   shift 2 ;;
        --taskfarm-site)      TASKFARM_SITE="$2";     shift 2 ;;
        -h|--help)            show_help; exit 0 ;;
        *) echo "ERROR: Unknown argument: '$1'" >&2; show_help >&2; exit 1 ;;
    esac
done

die() { echo "ERROR: $*" >&2; exit 1; }

[ -n "$ROOT_DIR" ] || die "--root-dir is required"
[ -d "$ROOT_DIR" ] || die "Root directory not found: $ROOT_DIR"
[[ "$MODE" =~ ^(branch|site|all)$ ]] || die "--mode must be branch, site, or all"

ROOT_DIR="$(cd "$ROOT_DIR" && pwd)"

# Default taskfarm file paths
[ -z "$TASKFARM_BRANCH" ] && TASKFARM_BRANCH="${ROOT_DIR}/codeml_taskfarm_fullpath.sh"
[ -z "$TASKFARM_SITE"   ] && TASKFARM_SITE="${ROOT_DIR}/codeml_taskfarm_fullpath_site_model.sh"

# Parse model list
IFS=',' read -ra KEEP_MODELS <<< "$MODELS"

# NSsites mapping for PAML site models
# PAML: model=0, NSsites=N where N encodes the site model
declare -A NSSITE_MAP
NSSITE_MAP=( ["m0"]=0 ["m1Neutral"]=1 ["m2Selection"]=2 ["m7"]=7 ["m8"]=8 )

# ── Helper: write a standard codeml.ctl ──────────────────────────────────────
write_ctl_branch() {
    local ctl_file="$1"
    cat > "$ctl_file" <<'CTL'
seqfile  = align.phy
treefile = tree
outfile  = out

noisy    = 3
verbose  = 1
runmode  = 0

seqtype  = 1
ndata    = 1
icode    = 0
cleandata = 1

model    = 2
NSsites  = 0
CodonFreq = 7
estFreq  = 0
clock    = 0
fix_omega = 0
omega    = 0.5
CTL
}

write_ctl_site() {
    local ctl_file="$1"
    local nssite_val="$2"
    cat > "$ctl_file" <<CTL
seqfile  = align.phy
treefile = tree
outfile  = out

noisy    = 3
verbose  = 1
runmode  = 0

seqtype  = 1
ndata    = 1
icode    = 0
cleandata = 1

model    = 0
NSsites  = ${nssite_val}
CodonFreq = 7
estFreq  = 0
clock    = 0
fix_omega = 0
omega    = 0.5
CTL
}

# ── Banner ────────────────────────────────────────────────────────────────────
echo "======================================================"
echo "  correction_codeml.sh"
echo "======================================================"
echo "  Root dir     : $ROOT_DIR"
echo "  Mode         : $MODE"
[ "$MODE" != "branch" ] && echo "  Site models  : $MODELS"
echo "======================================================"
echo ""

# Clear taskfarm files
[ "$MODE" != "site"   ] && > "$TASKFARM_BRANCH"
[ "$MODE" != "branch" ] && > "$TASKFARM_SITE"

BRANCH_OK=0; BRANCH_WARN=0
SITE_OK=0;   SITE_WARN=0

# ── Process each Codeml_Setup_codeml_input directory ─────────────────────────
shopt -s nullglob

for codeml_setup_dir in \
    "$ROOT_DIR"/Inferred_Genetree_*/${GENE_PATTERN}/Codeml_Setup_codeml_input
do
    [ -d "$codeml_setup_dir" ] || continue

    # ── BRANCH MODEL (cleaned_<species>/ directories) ─────────────────────
    # VESPA creates one cleaned_<species>/ directory per focal species passed
    # to branch_table.sh (e.g. cleaned_Aegilopsspeltoides/, cleaned_Aegilopsmutica/).
    # The glob 'cleaned_*/' processes ALL of them, so codeml_taskfarm_fullpath.sh
    # will contain commands for EVERY focal species, not just one.
    if [[ "$MODE" == "branch" || "$MODE" == "all" ]]; then
        for branch_dir in "${codeml_setup_dir}"/cleaned_*/; do
            [ -d "$branch_dir" ] || continue

            # Check for VESPA error log before processing
            if grep -rqI "Error running GenerateCodemlWorkspace.pl" \
                    "${codeml_setup_dir}" 2>/dev/null; then
                echo "[WARN] Skipping $branch_dir — GenerateCodemlWorkspace.pl error found" >&2
                echo "# WARN: codeml_setup error — $branch_dir" >> "$TASKFARM_BRANCH"
                (( BRANCH_WARN++ )) || true
                continue
            fi

            # Rename modelA → model_branch
            if [ -d "${branch_dir}modelA" ] && [ ! -d "${branch_dir}model_branch" ]; then
                mv "${branch_dir}modelA" "${branch_dir}model_branch"
            fi

            if [ -d "${branch_dir}model_branch" ]; then
                # Rename Omega0 → Omega0_5
                if [ -d "${branch_dir}model_branch/Omega0" ] \
                   && [ ! -d "${branch_dir}model_branch/Omega0_5" ]; then
                    mv "${branch_dir}model_branch/Omega0" \
                       "${branch_dir}model_branch/Omega0_5"
                fi
                # Remove any other Omega* directories
                for d in "${branch_dir}model_branch"/Omega*/; do
                    [ "$(basename "${d%/}")" = "Omega0_5" ] && continue
                    rm -rf "$d"
                done

                if [ -d "${branch_dir}model_branch/Omega0_5" ]; then
                    write_ctl_branch "${branch_dir}model_branch/Omega0_5/codeml.ctl"
                    echo "cd ${branch_dir}model_branch/Omega0_5; codeml" \
                        >> "$TASKFARM_BRANCH"
                    (( BRANCH_OK++ )) || true
                fi
            fi

            # Remove all directories that are not model_branch
            for d in "${branch_dir}"*/; do
                local_name="$(basename "${d%/}")"
                [ "$local_name" = "model_branch" ] && continue
                rm -rf "$d"
            done
        done
    fi

    # ── SITE MODELS (cleaned/ directory) ──────────────────────────────────
    if [[ "$MODE" == "site" || "$MODE" == "all" ]]; then
        site_dir="${codeml_setup_dir}/cleaned"
        [ -d "$site_dir" ] || continue

        # Remove model directories not in the keep list
        for d in "${site_dir}"/*/; do
            model_name="$(basename "${d%/}")"
            keep=false
            for km in "${KEEP_MODELS[@]}"; do
                [ "$km" = "$model_name" ] && { keep=true; break; }
            done
            [ "$keep" = false ] && rm -rf "$d"
        done

        # Process each kept site model
        for model_name in "${KEEP_MODELS[@]}"; do
            model_dir="${site_dir}/${model_name}"
            [ -d "$model_dir" ] || continue

            # Get NSsites value for this model
            nssite_val="${NSSITE_MAP[$model_name]:-}"
            if [ -z "$nssite_val" ]; then
                echo "  [WARN] Unknown site model '$model_name' — skipping" >&2
                (( SITE_WARN++ )) || true
                continue
            fi

            # Remove Omega directories except Omega0
            for d in "${model_dir}"/Omega*/; do
                [ "$(basename "${d%/}")" = "Omega0" ] && continue
                rm -rf "$d"
            done

            # Rename Omega0 → Omega0_5
            if [ -d "${model_dir}/Omega0" ] \
               && [ ! -d "${model_dir}/Omega0_5" ]; then
                mv "${model_dir}/Omega0" "${model_dir}/Omega0_5"
            fi

            if [ -d "${model_dir}/Omega0_5" ]; then
                # Write correct codeml.ctl:
                # model = 0 (not the model number!) + NSsites = <nssite_val>
                write_ctl_site "${model_dir}/Omega0_5/codeml.ctl" "$nssite_val"
                echo "cd ${model_dir}/Omega0_5; codeml" >> "$TASKFARM_SITE"
                (( SITE_OK++ )) || true
            fi
        done
    fi
done

shopt -u nullglob

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "======================================================"
echo "  Done  ($(date))"
if [[ "$MODE" == "branch" || "$MODE" == "all" ]]; then
    echo "  Branch model workspaces : $BRANCH_OK  (warnings: $BRANCH_WARN)"
    echo "  Taskfarm (branch)       : $TASKFARM_BRANCH"
fi
if [[ "$MODE" == "site" || "$MODE" == "all" ]]; then
    echo "  Site model workspaces   : $SITE_OK  (warnings: $SITE_WARN)"
    echo "  Taskfarm (site)         : $TASKFARM_SITE"
fi
echo "======================================================"
echo ""
echo "Next step: submit codeml runs (branch + m0, from the run_codeml/ directory)"
echo "  N=\$(bash run_codeml.sbatch --count \\"
echo "         --taskfarm-branch $TASKFARM_BRANCH \\"
echo "         --taskfarm-site   $TASKFARM_SITE \\"
echo "         --site-models     m0)"
echo "  sbatch --array=0-\$((N-1))%15 run_codeml.sbatch \\"
echo "         --taskfarm-branch $TASKFARM_BRANCH \\"
echo "         --taskfarm-site   $TASKFARM_SITE \\"
echo "         --site-models     m0"
echo ""
echo "  # To run branch model only:"
echo "  N=\$(wc -l < $TASKFARM_BRANCH)"
echo "  sbatch --array=0-\$((N-1))%15 run_codeml.sbatch \\"
echo "         --taskfarm-branch $TASKFARM_BRANCH \\"
echo "         --site-models     none"
