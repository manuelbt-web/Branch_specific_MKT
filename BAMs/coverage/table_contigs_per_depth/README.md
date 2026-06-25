# Coverage Summary Generator
## Description
Generates two comprehensive coverage tables (counts and percentages) from depth data.

# input sorted_depth_file

## Usage
```bash
./generate_coverage_tables.sh <sorted_depth_file> <output_directory>

#Input File

#number_position_per_depth_per_contig_sorted.txt format:

Contig  Total  10  15  20  25  30  40  50  60  70  80
chr1    5000   500 480 450 400 380 350 300 250 200 100
chr2    4000   800 780 700 650 600 550 500 400 300 200

 #   Counts Table (coverage_depth_summary_counts.txt):

    %/depth  Depth_10  Depth_15  ...  Depth_80
    5        1200      1150      ...  850
    10       1100      1050      ...  800
    ...

  #  Percentage Table (coverage_depth_summary_percentage.txt):

    %/depth  Depth_10  Depth_15  ...  Depth_80
    5        85.71     82.14     ...  60.71
    10       78.57     75.00     ...  57.14
    ...
