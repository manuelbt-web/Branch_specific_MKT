#!/bin/bash
set -euo pipefail

# ======================
# ARGUMENT HANDLING
# ======================
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <sorted_depth_file> <output_directory>"
    echo "Example: $0 number_position_per_depth_per_contig_sorted.txt ./results"
    exit 1
fi

# ======================
# CONFIGURATION
# ======================
INPUT_FILE="$1"
OUTPUT_DIR="$2"
COUNT_OUTPUT="${OUTPUT_DIR}/coverage_depth_summary_counts.txt"
PERCENT_OUTPUT="${OUTPUT_DIR}/coverage_depth_summary_percentage.txt"
POSITION_COUNT_OUTPUT="${OUTPUT_DIR}/coverage_depth_min_positions_counts.txt"

# Coverage thresholds (percentage of positions covered at each depth)
PERCENT_THRESHOLDS=(5 10 20 30 40 50 60 70 80 90 100)

# Depth values (matching columns in input file)
DEPTHS=(10 15 20 25 30 40 50 60 70 80)

# ======================
# VALIDATION
# ======================
validate_file() {
    [[ -f "$1" ]] || { echo "ERROR: Missing file: $1" >&2; exit 1; }
    [[ -s "$1" ]] || { echo "ERROR: Empty file: $1" >&2; exit 1; }
    [[ $(wc -l < "$1") -ge 2 ]] || { echo "ERROR: File contains only header: $1" >&2; exit 1; }
}

validate_directory() {
    mkdir -p "$1" || { echo "ERROR: Cannot create directory: $1" >&2; exit 1; }
}

# ======================
# MAIN PROCESSING
# ======================
generate_tables() {
    local input="$1"
    local count_out="$2"
    local percent_out="$3"
    local position_out="$4"

    echo "=== Generating coverage summary tables ==="
    echo "Input file: $(realpath "$input")"
    echo "Output directory: $(realpath "$(dirname "$count_out")")"

    # Get total contigs (skip header)
    local total_contigs=$(awk 'NR>1 {count++} END {print count}' "$input")
    echo "Total contigs: $total_contigs"

    if [[ "$total_contigs" -eq 0 ]]; then
        echo "ERROR: No contigs found in input file (only header?)" >&2
        exit 1
    fi

    # -------- 1. PERCENTAGE THRESHOLDS TABLE ----------
    {
        # Create header
        echo -ne "%Coverage/Depth"
        for depth in "${DEPTHS[@]}"; do
            echo -ne "\t${depth}X"
        done
        echo
    } > "$count_out"
    
    # Create identical header for percentages
    cp "$count_out" "$percent_out"

    # Process each coverage percentage threshold
    for thr in "${PERCENT_THRESHOLDS[@]}"; do
        echo -n "$thr%" >> "$count_out"
        echo -n "$thr%" >> "$percent_out"

        # Process each depth column
        for col_idx in "${!DEPTHS[@]}"; do
            # Columns are: Contig, Total_Positions, Depth10, Depth15,...Depth80
            local column=$((col_idx + 3))  
            local depth=${DEPTHS[$col_idx]}

            # Count contigs meeting the threshold
            local count=$(awk -v col="$column" -v thr="$thr" '
                BEGIN { count = 0 }
                NR == 1 { next }  # Skip header
                {
                    total_pos = $2
                    covered_pos = $col
                    if (total_pos > 0) {
                        coverage_pct = (covered_pos / total_pos) * 100
                        if (coverage_pct >= thr) {
                            count++
                        }
                    }
                }
                END { print count }
            ' "$input")

            # Calculate percentage of total contigs
            local percent=$(awk -v count="$count" -v total="$total_contigs" '
                BEGIN { printf "%.2f", (total > 0 ? (count / total) * 100 : 0) }
            ')

            echo -ne "\t$count" >> "$count_out"
            echo -ne "\t$percent" >> "$percent_out"
        done
        
        echo >> "$count_out"
        echo >> "$percent_out"
    done

    # -------- 2. MINIMUM POSITIONS TABLE ----------
    {
        echo -ne "MinPositions/Depth"
        for depth in "${DEPTHS[@]}"; do
            echo -ne "\t${depth}X"
        done
        echo
    } > "$position_out"

    # Process each minimum position count threshold
    for pos_thr in 1 10 50 100 200 500 1000; do
        echo -n "$pos_thr" >> "$position_out"

        for col_idx in "${!DEPTHS[@]}"; do
            local column=$((col_idx + 3))
            local depth=${DEPTHS[$col_idx]}

            # Count contigs with at least X positions at this depth
            local count=$(awk -v col="$column" -v thr="$pos_thr" '
                BEGIN { count = 0 }
                NR == 1 { next }  # Skip header
                $col >= thr { count++ }
                END { print count }
            ' "$input")

            echo -ne "\t$count" >> "$position_out"
        done
        echo >> "$position_out"
    done

    echo "=== Results saved ==="
    echo "- Count of contigs meeting % coverage thresholds: $count_out"
    echo "- Percentage of contigs meeting % coverage thresholds: $percent_out"
    echo "- Count of contigs with minimum positions at depth: $position_out"
}

# ======================
# EXECUTION
# ======================
validate_file "$INPUT_FILE"
validate_directory "$OUTPUT_DIR"
generate_tables "$INPUT_FILE" "$COUNT_OUTPUT" "$PERCENT_OUTPUT" "$POSITION_COUNT_OUTPUT"
