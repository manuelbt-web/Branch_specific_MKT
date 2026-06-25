#!/bin/bash
set -euo pipefail

# ======================
# VALIDATION FUNCTIONS
# ======================
validate_file() {
    if [ ! -f "$1" ]; then
        echo "ERROR: Missing file: $1" >&2
        exit 1
    fi
    if [ ! -s "$1" ]; then
        echo "ERROR: Empty file: $1" >&2
        exit 1
    fi
    if [ $(wc -l < "$1") -le 1 ]; then
        echo "ERROR: File contains only header: $1" >&2
        exit 1
    fi
}

validate_directory() {
    mkdir -p "$1" || { echo "ERROR: Cannot create directory: $1" >&2; exit 1; }
}

# ======================
# MAIN PROCESSING
# ======================
process_coverage() {
    local mean_depth="$1"
    local counts_file="$2"
    local output_dir="$3"
    local final_output="${output_dir}/number_position_per_depth_per_contig_sorted_with_total_positions.txt"
    local temp_file="${output_dir}/temp_processing.txt"
    local contig_order_file="${output_dir}/contig_order.txt"

    echo "=== Processing coverage data ==="
    echo "Input files:"
    echo "- Mean depth: $(realpath "$mean_depth")"
    echo "- Counts file: $(realpath "$counts_file")"
    echo "- Output directory: $(realpath "$output_dir")"

    # Create header
    {
        echo -ne "Contig\tTotal_Positions"
        echo -ne "\tPositions_10X\tPositions_15X\tPositions_20X\tPositions_25X"
        echo -ne "\tPositions_30X\tPositions_40X\tPositions_50X"
        echo -ne "\tPositions_60X\tPositions_70X\tPositions_80X"
        echo
    } > "$final_output"

    # First process counts file to get total positions per contig
    declare -A total_positions
    # Save contig order from counts file
    awk 'NR>1 {print $1}' "$counts_file" > "$contig_order_file"

    while IFS=$'\t' read -r contig covered total percentage; do
        # Skip header
        if [ "$contig" = "Contig" ]; then
            continue
        fi
        total_positions["$contig"]=$total
    done < "$counts_file"

    # Process mean depth file with AWK
    awk -v contig_order="$contig_order_file" '
    BEGIN {
        OFS = "\t"
        # Load all contigs from counts file first
        while ((getline contig < contig_order) > 0) {
            all_contigs[contig] = 1
            contig_list[++contig_count] = contig
        }
        close(contig_order)
    }
    {
        contig = $1
        depth = $3
        # Mark contig as seen in depth file
        seen_in_depth[contig] = 1
        # Count positions for each depth threshold
        if (depth >= 10) covered_10[contig]++
        if (depth >= 15) covered_15[contig]++
        if (depth >= 20) covered_20[contig]++
        if (depth >= 25) covered_25[contig]++
        if (depth >= 30) covered_30[contig]++
        if (depth >= 40) covered_40[contig]++
        if (depth >= 50) covered_50[contig]++
        if (depth >= 60) covered_60[contig]++
        if (depth >= 70) covered_70[contig]++
        if (depth >= 80) covered_80[contig]++
    }
    END {
        # Write results for all contigs in original order
        for (i = 1; i <= contig_count; i++) {
            contig = contig_list[i]
            # If contig wasn'\''t in depth file, set all counts to 0
            if (!(contig in seen_in_depth)) {
                print contig, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            } else {
                print contig, \
                    (contig in covered_10 ? covered_10[contig] : 0), \
                    (contig in covered_15 ? covered_15[contig] : 0), \
                    (contig in covered_20 ? covered_20[contig] : 0), \
                    (contig in covered_25 ? covered_25[contig] : 0), \
                    (contig in covered_30 ? covered_30[contig] : 0), \
                    (contig in covered_40 ? covered_40[contig] : 0), \
                    (contig in covered_50 ? covered_50[contig] : 0), \
                    (contig in covered_60 ? covered_60[contig] : 0), \
                    (contig in covered_70 ? covered_70[contig] : 0), \
                    (contig in covered_80 ? covered_80[contig] : 0)
            }
        }
    }' "$mean_depth" > "$temp_file"

    # Merge with total positions and sort
    while IFS=$'\t' read -r contig p10 p15 p20 p25 p30 p40 p50 p60 p70 p80; do
        total=${total_positions["$contig"]:-0}
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
               "$contig" "$total" "$p10" "$p15" "$p20" "$p25" "$p30" "$p40" "$p50" "$p60" "$p70" "$p80"
    done < "$temp_file" | sort -k1,1 >> "$final_output"

    rm -f "$temp_file" "$contig_order_file"

    echo "=== Results saved to ==="
    echo "$(realpath "$final_output")"
    echo "=== First few lines ==="
    head -n 5 "$final_output"
}

# ======================
# ARGUMENT HANDLING
# ======================
if [ $# -ne 3 ]; then
    echo "Usage: $0 <mean_depth_file> <counts_file> <output_directory>"
    echo "Example: $0 ./mean_depth.txt ./counts.txt ./results"
    exit 1
fi

# Validate inputs
validate_file "$1"
validate_file "$2"
validate_directory "$3"

# Execute processing
process_coverage "$1" "$2" "$3"
