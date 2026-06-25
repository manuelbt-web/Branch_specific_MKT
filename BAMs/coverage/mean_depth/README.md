# Mean Depth Calculator
Bash script that calculates mean sequencing depth from position-level depth files.

## Usage
```bash
./calculate_mean_depth.sh covered.depth all.depth output_dir

#Input Files

 #   Covered Positions (covered.depth):

        Positions meeting coverage criteria

        Format: contig pos sample1_depth sample2_depth ...

  #  All Positions (all.depth):

        All genomic positions

        Same format as covered positions file

#Method

For each position in both files, the script:

    Sums depth values across all samples

    Counts non-negative values

    Calculates mean depth (sum/count)

    Reports contig, position, mean depth, and samples used

#Output Files

    mean_depth_per_position_covered.txt:
    Copy

    contig    pos    mean_depth    samples_used
    chr1      100    25.4          7

    mean_depth_per_position.txt:
    Copy

    contig    pos    mean_depth    samples_used
    chr1      100    18.2          10


# Run analysis
./calculate_mean_depth.sh \
    data/covered.depth \
    data/all.depth \
    results/

