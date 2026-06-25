# Depth Coverage Analyzer

A Bash script for analyzing sequencing depth coverage from BAM-derived depth files, 
calculating position coverage statistics.

### Input File Format
depth per position all individuals
- **File type**: Tab-delimited text file (.depth)
- **Structure**:
- **Required columns**:
1. Contig/chromosome name
2. Position
3..N: Depth values for each sample

## Output Files
The script generates three output files in the specified directory: 

1. **Filtered Positions** (`depth_per_position_covered_all_individuals.depth`)
   - Contains all positions passing coverage thresholds
   - Maintains original input format plus filtering results 

2. **Contig Coverage Statistics** (`contig_percentage_covered_positions.txt`)
   - Per-contig coverage percentages
   - Columns:
     ```
     Contig Positions_covered Total_positions Percentage_covered
     ``` 

3. **Summary Statistics** (`number_covered_position.txt`)
   - Genome-wide coverage summary
   - Columns:
     ```
     total_positions covered_positions percentage_covered
     ```
## Analysis Workflow
1. **Position Filtering**:
   - For each genomic position:
     - Counts samples meeting minimum depth threshold (default >10 reads)
     - Retains positions where N samples (default: 7) meet threshold 

2. **Coverage Calculation**:
   - Calculates coverage percentages for each contig
   - Computes genome-wide coverage statistics 

3. **Result Compilation**:
   - Generates three structured output files
   - Provides both detailed and summary views of coverage



## Usage
### Basic Command
```bash
./covered_positions.sh /path/to/input.depth /path/to/output_directory
