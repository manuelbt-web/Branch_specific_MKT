#!/bin/bash
#SBATCH --job-name=read2snp
#SBATCH --output=./log_%j_%x_out.txt
#SBATCH --error=./log_%j_%x_err.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --time=300:00:00
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-long

echo "=== Starting read2snp analysis ==="
echo "Job started at: $(date)"
echo "Running on host: $(hostname)"


export PATH="/home/barrientosm/perl-5.16.3/perl-5.16.3:$PATH"
# Add local reads2snp to PATH
export PATH="$HOME/reads2snp_2.0:$PATH"

# Usage instructions
if [ $# -lt 3 ]; then
    echo "Usage: $0 <bam_list_path> <reference_fasta> <output_directory> [output_prefix]"
    echo "Example:"
    echo "  sbatch $0 /path/to/bam_list.txt /path/to/reference.fasta /path/to/output"
    echo "  sbatch $0 /path/to/bam_list.txt /path/to/reference.fasta /path/to/output tauschii"
    echo ""
    echo "output_prefix defaults to the basename of <output_directory> if omitted"
    echo "(e.g. output_directory=.../read2snp/tauschii -> prefix 'tauschii')."
    exit 1
fi

# Assign input parameters
BAM_LIST="$1"
REFERENCE="$2"
OUTPUT_DIR="$3"
OUTPUT_PREFIX="${4:-$(basename "$OUTPUT_DIR")}"

# Analysis parameters (customize as needed)
FIS_VALUE=0.18

echo "=== Analysis Parameters ==="
echo "BAM list: ${BAM_LIST}"
echo "Reference genome: ${REFERENCE}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "FIS threshold: ${FIS_VALUE}"
echo ""

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Verify input files exist
echo "Verifying input files..."
if [ ! -f "${BAM_LIST}" ]; then
    echo "ERROR: BAM list file not found at ${BAM_LIST}" >&2
    exit 1
fi

if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference file not found at ${REFERENCE}" >&2
    exit 1
fi

# Run reads2snp
echo "Starting reads2snp analysis..."
reads2snp_2.0.64.bin \
    -bamlist "${BAM_LIST}" \
    -bamref "${REFERENCE}" \
    -nbth 1 \
    -fis "${FIS_VALUE}" \
    -out "${OUTPUT_DIR}/${OUTPUT_PREFIX}"

# Check exit status
if [ $? -eq 0 ]; then
    echo "SNP extraction completed successfully."
    echo "Results saved to: ${OUTPUT_DIR}"
else
    echo "ERROR: reads2snp failed with exit code $?" >&2
    exit 1
fi

echo "Job finished at: $(date)"
