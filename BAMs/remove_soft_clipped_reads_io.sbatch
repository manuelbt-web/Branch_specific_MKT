#!/bin/bash
# ==============================================================================
# CLUSTER CONFIGURATION — read this section before running on a new cluster
# ==============================================================================
#
# This script filters soft-clipped reads from BAM files using samtools and awk.
# It is designed for SLURM-based HPC clusters.
#
# Lines marked [CLUSTER] below MUST be adapted before submission.
#
#   1. --account   : your project/billing account on the cluster
#   2. --partition : the queue or partition name for your cluster
#   3. module load : replace with your cluster's module names
#
# Optional resource adjustments:
#   --mem   : increase for very large BAM files (current: 9G)
#   --time  : increase for many or large BAM files (current: 20:00:00)
#
# Cluster examples:
#   IFB (ifb-core.ipsl.fr):
#     --partition=cpu-dedicated  --account=<your_project>@cpu
#     module load bioinfo-ifb; module load samtools/1.21
#
#   CIRAD cluster:
#     --partition=cpu-dedicated  --account=dedicated-cpu@cirad-normal
#     module load bioinfo-cirad; module load samtools/1.14-bin
#
#   Without a module system (samtools already in PATH):
#     Comment out the two module load lines below.
#
# Required software:  samtools >= 1.14
# ==============================================================================

#SBATCH --job-name=soft_clip_filter
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=9G
#SBATCH --time=20:00:00
#SBATCH --account=dedicated-cpu@cirad-normal   # [CLUSTER] change to your account
#SBATCH --partition=cpu-dedicated              # [CLUSTER] change to your partition

# ==============================================================================
# MODULE LOADING — adapt to your cluster's module system
# ==============================================================================
module load bioinfo-ifb        # [CLUSTER] change or remove if not needed
module load samtools/1.21      # [CLUSTER] change to your available version
# ==============================================================================

echo "=== Starting BAM Soft-Clip Filtering ==="
echo "Job started at: $(date)"
echo "Running on host: $(hostname)"

# Usage instructions
if [ $# -lt 2 ]; then
    echo "Usage: sbatch $0 <input_bam_directory> <output_directory> [threshold]"
    echo ""
    echo "Arguments:"
    echo "  input_bam_directory  Directory containing input .bam files"
    echo "  output_directory     Directory where filtered BAMs will be written"
    echo "  threshold            Maximum allowed soft-clipped bases per read (default: 1)"
    echo ""
    echo "Examples:"
    echo "  sbatch $0 /path/to/input_bams /path/to/output"
    echo "  sbatch $0 /path/to/input_bams /path/to/output 20"
    exit 1
fi

# Assign input parameters
INPUT_DIR="$1"
OUTPUT_DIR="$2"
THRESHOLD="${3:-1}"  # Default threshold: reads with >= 1 soft-clipped base are removed

echo "=== Job Parameters ==="
echo "Input BAM directory : $INPUT_DIR"
echo "Output directory    : $OUTPUT_DIR"
echo "Soft-clip threshold : $THRESHOLD (reads with >= $THRESHOLD soft-clipped bases are removed)"
echo ""

# Validate input directory
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Count input BAM files
BAM_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -name "*.bam" ! -name "*bai*" | wc -l)
if [ "$BAM_COUNT" -eq 0 ]; then
    echo "ERROR: No BAM files found in $INPUT_DIR"
    exit 1
fi
echo "Found $BAM_COUNT BAM file(s) to process"
echo ""

# Initialize counter
PROCESSED=0

# Process each BAM file
for INPUT_BAM in "$INPUT_DIR"/*.bam; do
    # Skip index files
    if [[ "$INPUT_BAM" == *"bai"* ]]; then
        continue
    fi

    PROCESSED=$((PROCESSED + 1))
    BASE_NAME=$(basename "$INPUT_BAM" .bam)
    OUTPUT_BAM="$OUTPUT_DIR/${BASE_NAME}_no_soft_clipped_reads.bam"

    echo "[$PROCESSED/$BAM_COUNT] Processing: $BASE_NAME"
    echo "  Filtering reads with >= $THRESHOLD soft-clipped bases..."

    # Filter soft-clipped reads using awk on the CIGAR string (field 6)
    # Reads with total soft-clipping >= threshold are discarded
    samtools view -h "$INPUT_BAM" | \
    awk -v threshold="$THRESHOLD" '
    BEGIN {
        FS = "\t";
        OFS = "\t";
        total_reads = 0;
        filtered_reads = 0;
    }
    {
        if ($1 ~ /^@/) {  # Pass through SAM header lines unchanged
            print $0;
            next;
        }

        total_reads++;
        if ($6 ~ /S/) {
            soft_start = 0;
            soft_end = 0;
            if (match($6, /^([0-9]+)S/, arr)) {
                soft_start = arr[1];
            }
            if (match($6, /([0-9]+)S$/, arr)) {
                soft_end = arr[1];
            }
            total_soft = soft_start + soft_end;

            if (total_soft >= threshold) {
                filtered_reads++;
                next;
            }
        }
        print $0;
    }
    END {
        printf "  Filtered %d of %d reads (%.2f%% removed)\n",
            filtered_reads, total_reads, (total_reads > 0 ? (filtered_reads/total_reads)*100 : 0) > "/dev/stderr";
    }' | samtools view -bS - > "$OUTPUT_BAM"

    # Index the output BAM file
    echo "  Indexing output BAM..."
    samtools index "$OUTPUT_BAM"

    echo "  Done. Output: $OUTPUT_BAM"
    echo ""
done

echo "=== Job Completed ==="
echo "Processed $PROCESSED BAM file(s)"
echo "Output directory: $OUTPUT_DIR"
echo "Job finished at: $(date)"
