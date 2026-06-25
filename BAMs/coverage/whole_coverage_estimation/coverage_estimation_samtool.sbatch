#!/bin/bash
# ==============================================================================
# CLUSTER CONFIGURATION — read this section before running on a new cluster
# ==============================================================================
#
# This script estimates mean sequencing coverage for each BAM file using the
# Lander-Waterman formula:  Coverage = (Reads * Mean_Read_Length) / Ref_Length
#
# Read length is auto-detected from each BAM via `samtools stats`, so no
# fixed read-length parameter is needed.
#
# Lines marked [CLUSTER] MUST be adapted before submission.
#
#   1. --account   : your project/billing account
#   2. --partition : the queue/partition name on your cluster
#   3. module load : replace with your cluster's module names
#
# Optional resource adjustments:
#   --mem   : 9G is sufficient for most cases (increase if BAMs are very large)
#   --time  : 10h is conservative for a typical number of samples
#
# Cluster examples:
#   IFB (ifb-core.ipsl.fr):
#     --partition=cpu-dedicated  --account=<your_project>@cpu
#     module load bioinfo-ifb; module load samtools/1.21; module load seqkit/2.9.0
#
#   CIRAD cluster:
#     --partition=cpu-dedicated  --account=dedicated-cpu@cirad-normal
#     module load bioinfo-cirad; module load samtools/1.14-bin; module load seqkit/2.8.1
#
#   Without a module system (tools already in PATH):
#     Comment out the module load lines below.
#
# Required software:  samtools >= 1.14, seqkit >= 2.8
# ==============================================================================

#SBATCH --job-name=coverage_bams
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=9G
#SBATCH --time=10:00:00
#SBATCH --account=dedicated-cpu@cirad-normal   # [CLUSTER] change to your account
#SBATCH --partition=cpu-dedicated              # [CLUSTER] change to your partition

# ==============================================================================
# MODULE LOADING — adapt to your cluster's module system
# ==============================================================================
module load bioinfo-ifb        # [CLUSTER] change or remove if not needed
module load samtools/1.21      # [CLUSTER] change to your available version
module load seqkit/2.9.0       # [CLUSTER] change to your available version
# ==============================================================================

# Check required arguments
if [ $# -ne 3 ]; then
    echo "Usage: sbatch $0 <bam_directory> <output_file> <reference_fasta>"
    echo ""
    echo "Arguments:"
    echo "  bam_directory    Directory containing .bam files"
    echo "  output_file      Path for the tab-delimited results file"
    echo "  reference_fasta  Reference genome/transcriptome FASTA file"
    echo ""
    echo "Output columns: Sample | Reads_mapped | Mean_Read_Length | Coverage"
    echo ""
    echo "Method: Coverage = (Reads_mapped * Mean_Read_Length) / Reference_length"
    echo "  - Paired-end: primary mapped reads counted, mean length per BAM"
    echo "  - Single-end: mapped reads counted, mean length per BAM"
    exit 1
fi

BAM_DIR="$1"
OUTPUT="$2"
REF_FASTA="$3"

{
    echo "=============================================="
    echo " Coverage Calculation Pipeline"
    echo " Job ID     : ${SLURM_JOB_ID:-N/A}"
    echo " Start Time : $(date)"
    echo "=============================================="
    echo " BAM Directory   : $BAM_DIR"
    echo " Output File     : $OUTPUT"
    echo " Reference FASTA : $REF_FASTA"
    echo "=============================================="
}

# Validate inputs
if [ ! -d "$BAM_DIR" ]; then
    echo "ERROR: BAM directory not found: $BAM_DIR"
    exit 1
fi
if [ ! -f "$REF_FASTA" ]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA"
    exit 1
fi

# Get total reference length using seqkit stats
echo "Calculating total reference length..."
REF_LENGTH=$(seqkit stats "$REF_FASTA" | awk 'NR==2 {gsub(",", "", $5); print $5}')

if [ -z "$REF_LENGTH" ] || [ "$REF_LENGTH" -eq 0 ]; then
    echo "ERROR: Could not determine reference length from: $REF_FASTA"
    exit 1
fi
echo "Total reference length: $REF_LENGTH bp"
echo "=============================================="

# Initialize output file
echo -e "Sample\tReads_mapped\tMean_Read_Length\tCoverage" > "$OUTPUT"

# Collect BAM files
mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 -name "*.bam" ! -name "*.bai" | sort)
BAM_COUNT="${#BAM_FILES[@]}"
if [ "$BAM_COUNT" -eq 0 ]; then
    echo "ERROR: No BAM files found in $BAM_DIR"
    exit 1
fi
echo "Found $BAM_COUNT BAM file(s) to process"
echo "=============================================="

COUNTER=1
for BAM in "${BAM_FILES[@]}"; do
    SAMPLE=$(basename "$BAM" .bam)
    echo "[$COUNTER/$BAM_COUNT] Processing: $SAMPLE"

    # Detect paired-end vs single-end
    PAIRED_READS=$(samtools view -c -f 1 "$BAM" 2>/dev/null)

    if [ "$PAIRED_READS" -gt 0 ]; then
        echo "  Type: Paired-end"
        # Count primary mapped reads (exclude unmapped and secondary/supplementary)
        READS=$(samtools view -c -F 260 "$BAM")
        READS=$((READS * 2))
    else
        echo "  Type: Single-end"
        READS=$(samtools view -c -F 260 "$BAM")
    fi

    if [ -z "$READS" ] || [ "$READS" -eq 0 ]; then
        echo "  WARNING: No mapped reads found — writing zeros"
        echo -e "$SAMPLE\t0\t0\t0.0000" >> "$OUTPUT"
        COUNTER=$((COUNTER + 1))
        continue
    fi

    # Auto-detect mean read length from samtools stats
    MEAN_LENGTH=$(samtools stats "$BAM" | awk '/^SN/ && /average length:/ {print int($4)}')

    # Calculate coverage
    COVERAGE=$(awk -v r="$READS" -v l="$MEAN_LENGTH" -v g="$REF_LENGTH" \
        'BEGIN { printf "%.4f", (r * l) / g }')

    {
        echo "  Reads mapped     : $READS"
        echo "  Mean read length : $MEAN_LENGTH bp"
        echo "  Coverage         : ${COVERAGE}x"
    }
    echo -e "$SAMPLE\t$READS\t$MEAN_LENGTH\t$COVERAGE" >> "$OUTPUT"
    echo "----------------------------------------------"

    COUNTER=$((COUNTER + 1))
done

{
    echo "=============================================="
    echo "Coverage calculation complete"
    echo "Results saved to: $OUTPUT"
    echo "End Time: $(date)"
    echo "=============================================="
}
