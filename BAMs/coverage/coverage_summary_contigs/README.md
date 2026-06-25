## Description

This script processes contig coverage data to generate:
- Percentage coverage distribution (0-100% in 10% increments)
- Count of contigs in each coverage bin
- Summary statistics of genome-wide coverage

### Input

contig_percentage_covered_positions.txt

number of positions covered per contig

## Output

### Output File
`contig_coverage_summary.txt` (in specified output directory):


## Usage

```bash
./coverage_analyzer.sh /path/to/input.txt /path/to/output_dir


