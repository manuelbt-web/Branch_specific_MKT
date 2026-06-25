#!/bin/bash
# ==============================================================================
# COVERAGE ANALYSIS PIPELINE — Master Wrapper
# ==============================================================================
#
# Runs the full coverage analysis pipeline or any individual steps.
# Can be submitted as a SLURM job or run locally in bash.
#
# QUICK START
# ───────────
# Run all main steps (1–7) locally:
#   bash run_coverage_pipeline.sh --bam-dir /path/to/bams --output /path/to/out
#
# Run a single step (e.g. redo step 2 with custom thresholds):
#   bash run_coverage_pipeline.sh --steps 2 --output /path/to/out \
#        --min-depth 15 --min-ind 8
#
# Submit the full pipeline to SLURM (edit [CLUSTER] lines below first):
#   sbatch run_coverage_pipeline.sh --bam-dir /path/bams --output /path/out
#
# PIPELINE STEPS
# ──────────────
#  1  Generate depth matrix from BAM files            [needs: --bam-dir]
#  2  Filter positions by coverage threshold          [needs: step 1 output]
#  3  Bin contigs by % coverage                       [needs: step 2 output]
#  4  Calculate mean depth per position               [needs: steps 1 + 2]
#  5  Calculate mean depth per contig                 [needs: step 4 output]
#  6  Count positions at each depth threshold         [needs: steps 2 + 4]
#  7  Build summary tables (counts and percentages)   [needs: step 6 output]
#  8  Estimate whole-genome coverage  [optional]      [needs: --bam-dir --ref]
#  9  Count covered RBH pairs         [optional]      [needs: --sp1 --sp2 --rbh]
#
# Use --help to see all options.
# ==============================================================================

# ==============================================================================
# SLURM HEADERS — active only when submitted via sbatch; ignored in bash
# Edit lines marked [CLUSTER] before submitting to your cluster
# ==============================================================================
#SBATCH --job-name=coverage_pipeline
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --account=dedicated-cpu@cirad-normal   # [CLUSTER] your billing account
#SBATCH --partition=cpu-dedicated              # [CLUSTER] your queue/partition

# ==============================================================================
# MODULE LOADING — uncomment and edit for your cluster if samtools/seqkit are
# not already in PATH.  Required for steps 1 and 8 only.
# ==============================================================================
# module load bioinfo-ifb          # [CLUSTER] IFB module collection
# module load samtools/1.21        # [CLUSTER] required for steps 1, 8
# module load seqkit/2.9.0         # [CLUSTER] required for step 8
# ==============================================================================

# Resolve the directory containing this script (needed to call sub-scripts)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ==============================================================================
# COLOUR HELPERS  (disabled automatically when stdout is not a terminal)
# ==============================================================================
if [ -t 1 ]; then
    RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
    CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
else
    RED=''; GREEN=''; YELLOW=''; CYAN=''; BOLD=''; RESET=''
fi

_log()    { echo "$*" >> "$LOG_FILE"; }
info()    { printf "${CYAN}[INFO]${RESET}    %s\n" "$*";  _log "[INFO]    $*"; }
success() { printf "${GREEN}[OK]${RESET}      %s\n" "$*"; _log "[OK]      $*"; }
warn()    { printf "${YELLOW}[WARN]${RESET}    %s\n" "$*" >&2; _log "[WARN]    $*"; }
die()     { printf "${RED}[ERROR]${RESET}   %s\n" "$*" >&2; _log "[ERROR]   $*"; exit 1; }
hr()      { printf '%0.s─' {1..72}; printf '\n'; }

step_banner() {
    local num="$1" title="$2"
    {
        echo ""
        hr
        printf "${BOLD}  STEP %s — %s${RESET}\n" "$num" "$title"
        printf "  Started: %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
        hr
    } | tee -a "$LOG_FILE"
}

step_ok() {
    local t0="$1"
    local elapsed=$(( $(date +%s) - t0 ))
    success "Step completed in ${elapsed}s" | tee -a "$LOG_FILE"
}

check_tool() {
    command -v "$1" >/dev/null 2>&1 \
        || die "'$1' not found in PATH — load the module or add it to PATH"
}

check_input_file() {
    [ -f "$1" ] \
        || die "Required input file not found: $1
       (Make sure the prerequisite step has been run, or check --output path)"
}

check_input_dir() {
    [ -d "$1" ] \
        || die "Required directory not found: $1"
}

# ==============================================================================
# HELP
# ==============================================================================
show_help() {
cat << 'HELP'

Coverage Analysis Pipeline — Master Wrapper
════════════════════════════════════════════

USAGE
  bash run_coverage_pipeline.sh [OPTIONS]
  sbatch run_coverage_pipeline.sh [OPTIONS]   ← SLURM submission

ALWAYS REQUIRED
  --output DIR       Master output directory (created if it does not exist)

STEP SELECTION
  --steps LIST       Steps to run. Use a comma-separated list of numbers, or
                     'all' to run steps 1–7 (default: all).
                     Examples:
                       --steps 2              run only step 2
                       --steps 3,4,5          run steps 3, 4, and 5
                       --steps all            run steps 1 through 7
                       --steps 1,2,3,4,5,6,7,8  main pipeline + coverage estimate

STEP-SPECIFIC OPTIONS
  Step 1  --bam-dir DIR    Directory containing .bam files          [required]
  Step 2  --min-depth N    Min reads per individual per position     [default: 10]
          --min-ind   N    Min individuals that must meet min-depth  [default: 6]
  Step 8  --bam-dir DIR    Same BAM directory as step 1             [required]
          --ref     FILE   Reference genome/transcriptome FASTA      [required]
  Step 9  --sp1     FILE   List of covered contigs for species 1    [required]
          --sp2     FILE   List of covered contigs for species 2    [required]
          --rbh     FILE   Tab-separated reciprocal best-hit pairs  [required]
  Step 10 (no extra arguments needed)
          Renders Chromosome_coverage.Rmd into an interactive HTML report.
          Requires R and the rmarkdown package.
          Reads outputs from steps 2, 3, 5, 7, and optionally 8.

OTHER OPTIONS
  --dry-run          Print all commands without executing them
  --help             Show this message

OUTPUT DIRECTORY STRUCTURE
  OUTPUT/
  ├── step1_depth/                        depth matrix (one column per individual)
  ├── step2_covered/                      filtered positions + per-contig stats
  ├── step3_coverage_bins/                contigs binned by % coverage
  ├── step4_mean_depth/                   mean depth per position
  ├── step5_mean_depth_per_contig/        mean depth per contig
  ├── step6_positions_per_depth/          positions per depth threshold per contig
  ├── step7_summary_tables/               final coverage summary tables
  ├── step8_whole_coverage/               whole-genome coverage estimates per sample
  ├── step9_rbh_pairs/                    RBH pair coverage counts
  ├── step10_report/
  │   └── coverage_report.html            interactive HTML report (step 10)
  └── pipeline.log                        full log of this run

EXAMPLES
  # Full main pipeline:
  bash run_coverage_pipeline.sh \
      --bam-dir /data/bams --output /results/coverage

  # Full pipeline + HTML report:
  bash run_coverage_pipeline.sh \
      --steps 1,2,3,4,5,6,7,10 \
      --bam-dir /data/bams --output /results/coverage

  # Only step 2, with stricter thresholds:
  bash run_coverage_pipeline.sh \
      --steps 2 --output /results/coverage --min-depth 15 --min-ind 8

  # Re-run steps 6, 7 and regenerate the report:
  bash run_coverage_pipeline.sh --steps 6,7,10 --output /results/coverage

  # Full pipeline + whole-genome coverage estimate + report:
  bash run_coverage_pipeline.sh \
      --steps 1,2,3,4,5,6,7,8,10 \
      --bam-dir /data/bams --output /results/coverage \
      --ref /data/reference.fasta

  # Preview what would run without executing:
  bash run_coverage_pipeline.sh --dry-run \
      --bam-dir /data/bams --output /results/coverage

HELP
}

# ==============================================================================
# DEFAULTS
# ==============================================================================
STEPS="all"
OUTPUT_DIR=""
BAM_DIR=""
REF_FASTA=""
MIN_DEPTH=10
MIN_IND=6
SP1_FILE=""
SP2_FILE=""
RBH_FILE=""
DRY_RUN=false

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================
[ $# -eq 0 ] && { show_help; exit 0; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)    show_help; exit 0 ;;
        --steps)      STEPS="$2";     shift 2 ;;
        --output)     OUTPUT_DIR="$2"; shift 2 ;;
        --bam-dir)    BAM_DIR="$2";   shift 2 ;;
        --ref)        REF_FASTA="$2"; shift 2 ;;
        --min-depth)  MIN_DEPTH="$2"; shift 2 ;;
        --min-ind)    MIN_IND="$2";   shift 2 ;;
        --sp1)        SP1_FILE="$2";  shift 2 ;;
        --sp2)        SP2_FILE="$2";  shift 2 ;;
        --rbh)        RBH_FILE="$2";  shift 2 ;;
        --dry-run)    DRY_RUN=true;   shift ;;
        *) printf "Unknown option: %s\nRun with --help for usage.\n" "$1" >&2; exit 1 ;;
    esac
done

# Expand 'all' shortcut
[[ "$STEPS" == "all" ]] && STEPS="1,2,3,4,5,6,7"

# Parse steps into array
IFS=',' read -ra STEP_LIST <<< "$STEPS"

# ==============================================================================
# VALIDATE BASE REQUIREMENTS
# ==============================================================================
[ -n "$OUTPUT_DIR" ] || { echo "ERROR: --output is required"; exit 1; }

# Create output dir first so we can write the log
mkdir -p "$OUTPUT_DIR"
LOG_FILE="$OUTPUT_DIR/pipeline.log"

# Named output sub-directories
D1="$OUTPUT_DIR/step1_depth"
D2="$OUTPUT_DIR/step2_covered"
D3="$OUTPUT_DIR/step3_coverage_bins"
D4="$OUTPUT_DIR/step4_mean_depth"
D5="$OUTPUT_DIR/step5_mean_depth_per_contig"
D6="$OUTPUT_DIR/step6_positions_per_depth"
D7="$OUTPUT_DIR/step7_summary_tables"
D8="$OUTPUT_DIR/step8_whole_coverage"
D9="$OUTPUT_DIR/step9_rbh_pairs"
D10="$OUTPUT_DIR/step10_report"

# Key intermediate files referenced across steps
F1="$D1/depth_per_position_all_individuals_BAMs.depth"
F2_COVERED="$D2/depth_per_position_covered_all_individuals.depth"
F2_CONTIG="$D2/contig_percentage_covered_positions.txt"
F4_COVERED="$D4/mean_depth_per_position_covered.txt"
F4_ALL="$D4/mean_depth_per_position.txt"
F5="$D5/mean_depth_per_contig.txt"
F6="$D6/number_position_per_depth_per_contig_sorted_with_total_positions.txt"
F8="$D8/coverage_results.tsv"

# Validate step-specific required args up front (better UX than failing mid-run)
for step in "${STEP_LIST[@]}"; do
    case "$step" in
        1) [ -n "$BAM_DIR" ]   || die "Step 1 requires --bam-dir" ;;
        2|3|4|5|6|7) ;;  # checked just before each step runs
        8) [ -n "$BAM_DIR" ]   || die "Step 8 requires --bam-dir"
           [ -n "$REF_FASTA" ] || die "Step 8 requires --ref" ;;
        9) [ -n "$SP1_FILE" ]  || die "Step 9 requires --sp1"
           [ -n "$SP2_FILE" ]  || die "Step 9 requires --sp2"
           [ -n "$RBH_FILE" ]  || die "Step 9 requires --rbh" ;;
        10) ;;  # Rscript availability is checked just before running
        *) die "Unknown step: '$step'  (valid steps: 1–10)" ;;
    esac
done

# ==============================================================================
# RUN HELPER  (no-op in dry-run mode)
# ==============================================================================
run() {
    if $DRY_RUN; then
        printf "${YELLOW}[DRY-RUN]${RESET} %s\n" "$*"
        _log "[DRY-RUN] $*"
    else
        _log "Running: $*"
        "$@"
    fi
}

# Skip output check in dry-run (files won't exist)
check_output() { $DRY_RUN || check_input_file "$1"; }

# ==============================================================================
# PIPELINE HEADER
# ==============================================================================
{
    hr
    printf "${BOLD}  COVERAGE ANALYSIS PIPELINE${RESET}\n"
    hr
    printf "  Date        : %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
    printf "  Host        : %s\n" "$(hostname)"
    printf "  SLURM job   : %s\n" "${SLURM_JOB_ID:-none (local run)}"
    printf "  Steps       : %s\n" "$STEPS"
    printf "  Output dir  : %s\n" "$OUTPUT_DIR"
    [ -n "$BAM_DIR"   ] && printf "  BAM dir     : %s\n" "$BAM_DIR"
    [ -n "$REF_FASTA" ] && printf "  Reference   : %s\n" "$REF_FASTA"
    printf "  Min depth   : %d reads per individual\n" "$MIN_DEPTH"
    printf "  Min indiv.  : %d individuals must meet min depth\n" "$MIN_IND"
    $DRY_RUN && printf "  ${YELLOW}DRY-RUN — no files will be written${RESET}\n"
    hr
    echo ""
} | tee -a "$LOG_FILE"

# ==============================================================================
# STEP EXECUTION
# ==============================================================================
STEP_STATUS=()   # collects "N|OK|path" or "N|FAIL|reason"

for step in "${STEP_LIST[@]}"; do
    T0=$(date +%s)
    STEP_OK=true

    case "$step" in

    # ─────────────────────────────────────────────────────────────────────
    1)
        step_banner 1 "Generate depth matrix from BAM files"
        check_input_dir "$BAM_DIR"
        check_tool samtools
        mkdir -p "$D1"
        if ! run bash "$SCRIPT_DIR/depth_per_individual/depth_per_individual.sbatch" \
                "$BAM_DIR" "$D1"; then
            STEP_OK=false
        fi
        check_output "$F1"
        $STEP_OK && step_ok "$T0" || warn "Step 1 finished with errors"
        $STEP_OK && STEP_STATUS+=("1|OK|$F1") || STEP_STATUS+=("1|FAIL|see $LOG_FILE")
        ;;

    # ─────────────────────────────────────────────────────────────────────
    2)
        step_banner 2 "Filter positions by coverage threshold"
        info "Parameters: min_depth=$MIN_DEPTH  min_individuals=$MIN_IND"
        check_input_file "$F1"
        mkdir -p "$D2"
        if ! run bash "$SCRIPT_DIR/covered_positions/covered_positions.sh" \
                "$F1" "$D2" "$MIN_DEPTH" "$MIN_IND"; then
            STEP_OK=false
        fi
        check_output "$F2_COVERED"
        $STEP_OK && step_ok "$T0" || warn "Step 2 finished with errors"
        $STEP_OK && STEP_STATUS+=("2|OK|$D2") || STEP_STATUS+=("2|FAIL|see $LOG_FILE")
        ;;

    # ─────────────────────────────────────────────────────────────────────
    3)
        step_banner 3 "Bin contigs by percentage coverage"
        check_input_file "$F2_CONTIG"
        mkdir -p "$D3"
        if ! run bash "$SCRIPT_DIR/coverage_summary_contigs/percentage_coverage.sh" \
                "$F2_CONTIG" "$D3"; then
            STEP_OK=false
        fi
        $STEP_OK && step_ok "$T0" || warn "Step 3 finished with errors"
        $STEP_OK && STEP_STATUS+=("3|OK|$D3/contig_coverage_summary.txt") \
                 || STEP_STATUS+=("3|FAIL|see $LOG_FILE")
        ;;

    # ─────────────────────────────────────────────────────────────────────
    4)
        step_banner 4 "Calculate mean depth per position"
        check_input_file "$F2_COVERED"
        check_input_file "$F1"
        mkdir -p "$D4"
        if ! run bash "$SCRIPT_DIR/mean_depth/calculate_mean_depth.sh" \
                "$F2_COVERED" "$F1" "$D4"; then
            STEP_OK=false
        fi
        check_output "$F4_COVERED"
        $STEP_OK && step_ok "$T0" || warn "Step 4 finished with errors"
        $STEP_OK && STEP_STATUS+=("4|OK|$D4") || STEP_STATUS+=("4|FAIL|see $LOG_FILE")
        ;;

    # ─────────────────────────────────────────────────────────────────────
    5)
        step_banner 5 "Calculate mean depth per contig"
        check_input_file "$F4_COVERED"
        mkdir -p "$D5"
        if ! run bash "$SCRIPT_DIR/mean_depth_per_contig/mean_depth_per_contig.sh" \
                "$F4_COVERED" "$F5"; then
            STEP_OK=false
        fi
        check_output "$F5"
        $STEP_OK && step_ok "$T0" || warn "Step 5 finished with errors"
        $STEP_OK && STEP_STATUS+=("5|OK|$F5") || STEP_STATUS+=("5|FAIL|see $LOG_FILE")
        ;;

    # ─────────────────────────────────────────────────────────────────────
    6)
        step_banner 6 "Count positions at each depth threshold per contig"
        check_input_file "$F4_COVERED"
        check_input_file "$F2_CONTIG"
        mkdir -p "$D6"
        if ! run bash "$SCRIPT_DIR/number_position_per_depth/number_positions_per_depth_per_contig.sh" \
                "$F4_COVERED" "$F2_CONTIG" "$D6"; then
            STEP_OK=false
        fi
        check_output "$F6"
        $STEP_OK && step_ok "$T0" || warn "Step 6 finished with errors"
        $STEP_OK && STEP_STATUS+=("6|OK|$F6") || STEP_STATUS+=("6|FAIL|see $LOG_FILE")
        ;;

    # ─────────────────────────────────────────────────────────────────────
    7)
        step_banner 7 "Build summary tables (counts and percentages)"
        check_input_file "$F6"
        mkdir -p "$D7"
        if ! run bash "$SCRIPT_DIR/table_contigs_per_depth/table_contigs_per_depth_coverage.sh" \
                "$F6" "$D7"; then
            STEP_OK=false
        fi
        $STEP_OK && step_ok "$T0" || warn "Step 7 finished with errors"
        $STEP_OK && STEP_STATUS+=("7|OK|$D7") || STEP_STATUS+=("7|FAIL|see $LOG_FILE")
        ;;

    # ─────────────────────────────────────────────────────────────────────
    8)
        step_banner 8 "Estimate whole-genome coverage per sample"
        check_input_dir "$BAM_DIR"
        check_input_file "$REF_FASTA"
        check_tool samtools
        check_tool seqkit
        mkdir -p "$D8"
        if ! run bash "$SCRIPT_DIR/whole_coverage_estimation/coverage_estimation_samtool.sbatch" \
                "$BAM_DIR" "$F8" "$REF_FASTA"; then
            STEP_OK=false
        fi
        $STEP_OK && step_ok "$T0" || warn "Step 8 finished with errors"
        $STEP_OK && STEP_STATUS+=("8|OK|$F8") || STEP_STATUS+=("8|FAIL|see $LOG_FILE")
        ;;

    # ─────────────────────────────────────────────────────────────────────
    9)
        step_banner 9 "Count covered reciprocal best-hit pairs"
        check_input_file "$SP1_FILE"
        check_input_file "$SP2_FILE"
        check_input_file "$RBH_FILE"
        mkdir -p "$D9"
        RBH_OUT="$D9/rbh_pairs_coverage_result.txt"
        if ! run bash "$SCRIPT_DIR/count_paired_covered_contig_from_RBH/count_covered_pairs_contigs.sh" \
                "$SP1_FILE" "$SP2_FILE" "$RBH_FILE" > "$RBH_OUT"; then
            STEP_OK=false
        fi
        $STEP_OK && step_ok "$T0" || warn "Step 9 finished with errors"
        $STEP_OK && STEP_STATUS+=("9|OK|$RBH_OUT") || STEP_STATUS+=("9|FAIL|see $LOG_FILE")
        ;;

    # ─────────────────────────────────────────────────────────────────────
    10)
        step_banner 10 "Generate HTML coverage report (R Markdown)"
        check_tool Rscript
        # Verify rmarkdown is available in R
        if ! Rscript -e "stopifnot(requireNamespace('rmarkdown', quietly=TRUE))" 2>/dev/null; then
            die "R package 'rmarkdown' not found. Install it with:
       Rscript -e \"install.packages('rmarkdown')\""
        fi
        RMD_SRC="$SCRIPT_DIR/Chromosome_coverage.Rmd"
        REPORT_OUT="$D10/coverage_report.html"
        check_input_file "$RMD_SRC"
        mkdir -p "$D10"
        info "Rendering report from: $(basename "$RMD_SRC")"
        info "Output will be:        $REPORT_OUT"
        if ! run Rscript -e "rmarkdown::render(
                '$RMD_SRC',
                params           = list(output_dir = '$OUTPUT_DIR'),
                output_file      = '$REPORT_OUT',
                intermediates_dir = '$D10',
                quiet            = TRUE
            )"; then
            STEP_OK=false
        fi
        check_output "$REPORT_OUT"
        $STEP_OK && step_ok "$T0" || warn "Step 10 finished with errors"
        $STEP_OK && STEP_STATUS+=("10|OK|$REPORT_OUT") \
                 || STEP_STATUS+=("10|FAIL|see $LOG_FILE")
        ;;
    esac

done

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
{
    echo ""
    hr
    printf "${BOLD}  PIPELINE SUMMARY${RESET}\n"
    hr
    printf "  Finished : %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
    echo ""

    # Step status table
    printf "  %-10s  %-8s  %s\n" "STEP" "STATUS" "OUTPUT"
    printf "  %-10s  %-8s  %s\n" "──────────" "────────" "──────────────────────────────────────────"

    ALL_OK=true
    for entry in "${STEP_STATUS[@]}"; do
        IFS='|' read -r num status path <<< "$entry"
        if [ "$status" = "OK" ]; then
            printf "  Step %-5s  ${GREEN}%-8s${RESET}  %s\n" "$num" "$status" "$path"
        else
            printf "  Step %-5s  ${RED}%-8s${RESET}  %s\n" "$num" "$status" "$path"
            ALL_OK=false
        fi
    done

    echo ""

    # Output file inventory
    if ! $DRY_RUN && [ -d "$OUTPUT_DIR" ]; then
        printf "  %-52s  %7s  %s\n" "OUTPUT FILE" "SIZE" "DATA ROWS"
        printf "  %-52s  %7s  %s\n" "───────────────────────────────────────────────────" "───────" "─────────"
        while IFS= read -r f; do
            size=$(du -sh "$f" 2>/dev/null | cut -f1)
            relpath="${f#"$OUTPUT_DIR/"}"
            if [[ "$f" == *.html ]]; then
                printf "  %-52s  %7s  (report)\n" "$relpath" "$size"
            else
                rows=$(( $(wc -l < "$f") - 1 ))
                printf "  %-52s  %7s  %d rows\n" "$relpath" "$size" "$rows"
            fi
        done < <(find "$OUTPUT_DIR" \
                    \( -name "*.depth" -o -name "*.txt" -o -name "*.tsv" -o -name "*.html" \) \
                    ! -name "pipeline.log" | sort)
    fi

    echo ""

    if $ALL_OK; then
        success "All steps completed successfully."
    else
        warn "One or more steps finished with errors. Check: $LOG_FILE"
    fi

    printf "  Log file : %s\n" "$LOG_FILE"
    hr
    echo ""
} | tee -a "$LOG_FILE"
