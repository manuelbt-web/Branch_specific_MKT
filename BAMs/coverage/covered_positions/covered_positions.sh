#!/bin/bash

# Check for required arguments
if [ "$#" -lt 2 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <input_file.depth> <output_directory> [min_depth] [min_individuals]"
    echo ""
    echo "Arguments:"
    echo "  input_file.depth   Depth matrix from depth_per_individual.sbatch"
    echo "                     Columns: contig | pos | sample1 | sample2 | ..."
    echo "  output_directory   Where output files will be written"
    echo "  min_depth          Minimum read depth per individual to count a position"
    echo "                     as covered (default: 10)"
    echo "  min_individuals    Minimum number of individuals that must meet min_depth"
    echo "                     for a position to be retained (default: 6)"
    echo ""
    echo "Output files:"
    echo "  depth_per_position_covered_all_individuals.depth  Filtered positions only"
    echo "  contig_percentage_covered_positions.txt           Per-contig coverage stats"
    echo "  number_covered_position.txt                       Genome-wide summary"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_DIR="$2"
MIN_DEPTH="${3:-10}"
MIN_INDIVIDUALS="${4:-6}"

# Define output filenames
FILTERED_OUTPUT="${OUTPUT_DIR}/depth_per_position_covered_all_individuals.depth"
COUNT_OUTPUT="${OUTPUT_DIR}/contig_percentage_covered_positions.txt"
SUMMARY_OUTPUT="${OUTPUT_DIR}/number_covered_position.txt"

# Validate input
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR" || {
    echo "Error: Failed to create output directory: $OUTPUT_DIR" >&2
    exit 1
}

{
    echo "=== Starting Depth Coverage Analysis ==="
    echo "Input file      : $INPUT_FILE"
    echo "Output directory: $OUTPUT_DIR"
    echo "Min depth       : $MIN_DEPTH reads per individual"
    echo "Min individuals : $MIN_INDIVIDUALS individuals must meet min_depth"
    echo ""
    echo "Output files:"
    echo "  1. Filtered positions : $(basename "$FILTERED_OUTPUT")"
    echo "  2. Coverage per contig: $(basename "$COUNT_OUTPUT")"
    echo "  3. Genome-wide summary: $(basename "$SUMMARY_OUTPUT")"
    echo ""
}

# Initialize output files with headers
echo -e "contig\tpos\t$(head -1 "$INPUT_FILE" | cut -f3-)" > "$FILTERED_OUTPUT" || {
    echo "Error: Failed to initialize filtered output file" >&2
    exit 1
}

echo -e "Contig\tPositions_covered\tTotal_positions\tPercentage_covered" > "$COUNT_OUTPUT" || {
    echo "Error: Failed to initialize count output file" >&2
    exit 1
}

echo -e "total_positions\tcovered_positions\tpercentage_covered" > "$SUMMARY_OUTPUT" || {
    echo "Error: Failed to initialize summary output file" >&2
    exit 1
}

# Process depth data with awk
# A position is "covered" if at least MIN_INDIVIDUALS samples have depth >= MIN_DEPTH
echo "Processing depth data..."
awk -v filtered="$FILTERED_OUTPUT" \
    -v counts="$COUNT_OUTPUT" \
    -v min_depth="$MIN_DEPTH" \
    -v min_ind="$MIN_INDIVIDUALS" \
'
BEGIN {
    OFS = "\t";
    current_contig = "";
    contig_count = 0;
    total_positions_per_contig = 0;
    total_filtered_positions = 0;
    total_all_positions = 0;
}
NR == 1 { next }  # Skip header
{
    total_all_positions++;
    contig = $1;
    count = 0;

    # Count individuals meeting the minimum depth threshold
    for (i = 3; i <= NF; i++) {
        if ($i >= min_depth) count++;
    }

    # Flush stats when moving to a new contig
    if (contig != current_contig) {
        if (current_contig != "") {
            percentage = (contig_count / total_positions_per_contig) * 100;
            print current_contig, contig_count, total_positions_per_contig, percentage >> counts;
        }
        current_contig = contig;
        contig_count = 0;
        total_positions_per_contig = 0;
    }

    total_positions_per_contig++;
    if (count >= min_ind) {
        contig_count++;
        total_filtered_positions++;
        print $0 >> filtered;
    }
}
END {
    # Write stats for the last contig
    if (current_contig != "") {
        percentage = (contig_count / total_positions_per_contig) * 100;
        print current_contig, contig_count, total_positions_per_contig, percentage >> counts;
    }

    percentage_covered = (total_filtered_positions / total_all_positions) * 100;

    # Write summary to a temp file (shell picks it up below)
    print total_all_positions, total_filtered_positions, percentage_covered > "summary.tmp";

    print "Processed " total_all_positions " total positions" > "/dev/stderr";
    print "Retained " total_filtered_positions " covered positions (" percentage_covered "%)" > "/dev/stderr";
}
' "$INPUT_FILE" || {
    echo "Error: Failed to process input file" >&2
    exit 1
}

# Write genome-wide summary from temp file
read -r total_positions covered_positions percentage_covered < summary.tmp
printf "%s\t%s\t%.2f\n" "$total_positions" "$covered_positions" "$percentage_covered" >> "$SUMMARY_OUTPUT"
rm -f summary.tmp

{
    echo ""
    echo "=== Analysis Complete ==="
    echo "Results:"
    echo "  1. Filtered positions : $FILTERED_OUTPUT ($(( $(wc -l < "$FILTERED_OUTPUT") - 1 )) positions)"
    echo "  2. Coverage per contig: $COUNT_OUTPUT ($(( $(wc -l < "$COUNT_OUTPUT") - 1 )) contigs)"
    echo "  3. Genome-wide summary: $SUMMARY_OUTPUT"
    echo ""
    echo "  Total positions  : $total_positions"
    echo "  Covered positions: $covered_positions"
    echo "  Genome coverage  : ${percentage_covered}%"
}
