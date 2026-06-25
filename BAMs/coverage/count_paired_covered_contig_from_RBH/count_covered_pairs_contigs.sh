#!/bin/bash

# ------------------- Script: Best-Hit Pair Coverage Analysis -------------------
# This script counts the number of best-hit pairs where both contigs are covered 
# by at least X% in two different species (coverage files provided as input).
# ------------------------------------------------------------------------------

# ----------- Argument Handling -----------

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <covered_species1.txt> <covered_species2.txt> <best_hits.tsv>"
    echo "Where:"
    echo "  <covered_species1.txt> - List of contigs meeting coverage threshold for species 1"
    echo "  <covered_species2.txt> - List of contigs meeting coverage threshold for species 2"
    echo "  <best_hits.tsv>        - Tab-separated file of best-hit pairs"
    exit 1
fi

file_covered_species1="$1"
file_covered_species2="$2"
file_best_hits="$3"

# ----------- Initialization -----------

echo "=== Best-Hit Pair Coverage Analyzer ==="
echo "Species 1 coverage file: $file_covered_species1"
echo "Species 2 coverage file: $file_covered_species2"
echo "Best-hit pairs file: $file_best_hits"
echo "--------------------------------------"

# Create timestamp for unique temp files
timestamp=$(date +%Y%m%d_%H%M%S)

# Temporary files with timestamp
temp_covered_species1="covered_species1_${timestamp}.txt"
temp_covered_species2="covered_species2_${timestamp}.txt"
temp_sorted_species1="sorted_species1_${timestamp}.txt"
temp_sorted_species2="sorted_species2_${timestamp}.txt"
temp_sorted_best_hits="sorted_best_hits_${timestamp}.txt"

# ----------- Step 1: Prepare Input Data -----------

echo -e "\n[1/3] Preparing input data..."
echo "  - Processing $file_covered_species1"
awk '{print $1}' "$file_covered_species1" | sort > "$temp_sorted_species1"

echo "  - Processing $file_covered_species2"
awk '{print $1}' "$file_covered_species2" | sort > "$temp_sorted_species2"

echo "  - Sorting best-hit pairs"
sort "$file_best_hits" > "$temp_sorted_best_hits"

# ----------- Step 2: Count Matching Pairs -----------

echo -e "\n[2/3] Analyzing best-hit pairs..."
matching_pairs=0
total_pairs=$(wc -l < "$temp_sorted_best_hits")
current_pair=0

while read -r line; do
    current_pair=$((current_pair + 1))
    if (( current_pair % 1000 == 0 )); then
        echo "  - Processed $current_pair/$total_pairs pairs..."
    fi
    
    contig_species1=$(echo "$line" | awk '{print $1}')
    contig_species2=$(echo "$line" | awk '{print $2}')
    
    if grep -Fxq "$contig_species1" "$temp_sorted_species1" && grep -Fxq "$contig_species2" "$temp_sorted_species2"; then
        matching_pairs=$((matching_pairs + 1))
    fi
done < "$temp_sorted_best_hits"

# ----------- Step 3: Display Results -----------

echo -e "\n[3/3] Results:"
echo "--------------------------------------"
echo "Total best-hit pairs analyzed: $total_pairs"
echo "Pairs where both contigs meet coverage threshold: $matching_pairs"
echo "Percentage of matching pairs: $(awk "BEGIN {printf \"%.2f\", ($matching_pairs/$total_pairs)*100}")%"
echo "--------------------------------------"

# ----------- Cleanup -----------

echo -e "\nCleaning up temporary files..."
rm -f "$temp_covered_species1" "$temp_covered_species2" \
      "$temp_sorted_species1" "$temp_sorted_species2" \
      "$temp_sorted_best_hits"

echo "=== Analysis complete ==="
