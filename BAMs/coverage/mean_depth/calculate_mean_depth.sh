#!/bin/bash

# Check if required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 /path/to/covered.depth /path/to/all.depth /path/to/output_dir"
    echo "Example: $0 ./covered_positions.depth ./all_positions.depth ./results"
    exit 1
fi

# Input files
COVERED_IN="$1"  # depth_per_position_covered_all_individuals.depth
ALL_IN="$2"      # depth_per_position_all_individuals_BAMs.depth
OUT_DIR="$3"

# Output files
MEAN_COVERED="${OUT_DIR}/mean_depth_per_position_covered.txt"
MEAN_ALL="${OUT_DIR}/mean_depth_per_position.txt"

# Create output directory
mkdir -p "$OUT_DIR" || {
    echo "Error: Failed to create output directory: $OUT_DIR" >&2
    exit 1
}

echo "=== Starting Mean Depth Calculation ==="
echo "Input files:"
echo "- Covered positions: $(basename "$COVERED_IN")"
echo "- All positions: $(basename "$ALL_IN")"
echo "Output directory: $OUT_DIR"
echo ""

# Function to calculate mean depth
calculate_mean() {
    local input="$1"
    local output="$2"
    local desc="$3"
    
    echo "Processing $desc..."
    echo -e "contig\tpos\tmean_depth\tsamples_used" > "$output"
    
    awk '
    BEGIN {
        FS = "\t";
        OFS = "\t";
        print "  0 positions processed..." > "/dev/stderr";
    }
    NR == 1 { next } # Skip header
    {
        sum = count = 0;
        for (i = 3; i <= NF; i++) {
            if ($i >= 0) { sum += $i; count++ }
        }
        if (count > 0) {
            print $1, $2, sum/count, count;
        }
        
        # Progress report
        if (NR % 100000 == 0) {
            printf "  %d positions processed...\r", NR-1 > "/dev/stderr";
        }
    }
    END {
        printf "  %d positions processed\n", NR-1 > "/dev/stderr";
    }' "$input" >> "$output" || {
        echo "Error processing $desc file" >&2
        exit 1
    }
}

# Calculate mean for covered positions
calculate_mean "$COVERED_IN" "$MEAN_COVERED" "covered positions"

# Calculate mean for all positions
calculate_mean "$ALL_IN" "$MEAN_ALL" "all positions"

echo -e "\n=== Analysis Complete ==="
echo "Output files created:"
echo "- $(basename "$MEAN_COVERED")"
echo "- $(basename "$MEAN_ALL")"
echo "Total positions processed:"
echo "- Covered: $(($(wc -l < "$MEAN_COVERED") - 1))"
echo "- All: $(($(wc -l < "$MEAN_ALL") - 1))"
