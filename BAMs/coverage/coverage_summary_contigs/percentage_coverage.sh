#!/bin/bash

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/input_file.txt /path/to/output_directory"
    exit 1
fi

input_file="$1"
output_dir="$2"
output_file="${output_dir}/contig_coverage_summary.txt"

# Validate input file
if [ ! -f "$input_file" ]; then
    echo "Error: Input file not found: $input_file" >&2
    exit 1
fi

mkdir -p "$output_dir" || {
    echo "Error: Could not create output directory: $output_dir" >&2
    exit 1
}

echo "=== Starting Coverage Analysis ==="
echo "Input file: $input_file"
echo "Output file: $output_file"
echo ""

# Initialize output
echo -e "Threshold_coverage\tCount_of_contigs" > "$output_file"

# Process file
awk -v output="$output_file" '
BEGIN {
    FS = "\t";
    OFS = "\t";

    # Initialize all bins (0-100% in 10% increments)
    for (i = 0; i <= 100; i += 10) counts[i] = 0;
}
NR > 1 {
    if ($3 > 0) {
        coverage = ($2 / $3) * 100;
    } else {
        coverage = 0;
    }

    # Assign to correct bin
    threshold = int(coverage / 10) * 10;

    # Special case for exactly 100%
    if (coverage == 100) threshold = 100;
    
    # For coverage > 100%, still count in 100% bin
    if (coverage > 100) threshold = 100;

    counts[threshold]++;
}
END {
    # Print all thresholds including 100%
    for (i = 0; i <= 100; i += 10) {
        print i"%", counts[i] >> output;
    }

    # Log total to stderr
    total = 0;
    for (i = 0; i <= 100; i += 10) total += counts[i];
    print "\nTotal contigs counted: " total > "/dev/stderr";
}' "$input_file" || {
    echo "Error: Failed to process input file" >&2
    exit 1
}

echo ""
echo "=== Analysis Complete ==="
echo "Coverage summary saved to: $output_file"
echo "Total contigs analyzed: $(($(wc -l < "$input_file") - 1))"
