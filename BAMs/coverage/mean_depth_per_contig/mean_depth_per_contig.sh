#!/bin/bash
set -euo pipefail

# Check input parameters
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_depth_file.txt> <output_mean_depth.txt>"
    exit 1
fi

input="$1"
output="$2"

# Validate input file
if [[ ! -f "$input" ]]; then
    echo "Error: Input file not found: $input" >&2
    exit 1
fi

echo "Calculating mean depth per contig..."
echo "Input: $input"
echo "Output: $output"

# Create output file with header first
echo -e "contig\tmean_depth" > "$output"

# Process data and append to output file
awk '
NR == 1 { next }  # Skip input header
{
    contig = $1
    depth = $3
    
    sum[contig] += depth
    count[contig]++
}
END {
    for (contig in sum) {
        mean = sum[contig] / count[contig]
        printf "%s\t%.2f\n", contig, mean
    }
}' "$input" | sort >> "$output"  # Note the >> to append to file

echo "Calculation complete. Results saved to $output"
echo "First 5 lines with header:"
head -n 6 "$output"  # Show header + 5 contigs
