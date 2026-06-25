# Coverage Analysis Pipeline

## Description
Calculates position coverage statistics across multiple depth thresholds while maintaining exact contig order.

## Input Requirements

### Required Files
1. **Mean Depth File** (`mean_depth_per_position_covered.txt`):

contig pos mean_depth samples_used
chr1 100 25.4 7


2. **Counts File** (`countig_percentage_covered_positions.txt`):

Contig Positions_covered Total_positions Percentage
chr1 4500 5000 90.00

## Usage

```bash
./number_positions_per_depth_per_contig.sh <mean_depth_file> <counts_file> <output_directory>

./number_positions_per_depth_per_contig.sh \
    data/mean_depth_per_position_covered.txt \
    data/countig_percentage_covered_positions.txt \
    results/

#Creates in output directory:

number_position_per_depth_per_contig_sorted.txt

#Format:

Contig  10  15  20  25  30  40  50  60  70  80
chr1    500 480 450 400 380 350 300 250 200 100
chr2    800 780 700 650 600 550 500 400 300 200
