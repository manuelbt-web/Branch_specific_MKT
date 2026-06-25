# BAM Coverage Calculator

A SLURM-compatible script to calculate average coverage from BAM files.

##  Method

Calculates average coverage using the formula:  
`Coverage = (Mapped Reads * Read Length) / Reference Length`

Steps:
1. Computes total reference length using `seqkit stats`
2. For each BAM file:
   - Counts mapped reads (`samtools view -c -F 4`)
   - Calculates average coverage
3. Outputs tab-delimited results

##  Input Requirements

### Required Parameters (positional):
1. `BAM_DIR`: Path to directory containing `.bam` files
2. `OUTPUT_PATH`: Path for results file
3. `REF_FASTA`: Reference genome/transcriptome FASTA file
4. `READ_LENGTH`: Integer read length (e.g., 125)

##Output

#Tab-delimited file with columns:

    Sample: Sample name (from BAM filename)

    Reads_mapped: Number of mapped reads

    Average_Coverage: Calculated coverage (4 decimal places)

#Example:

Sample      Reads_mapped  Average_Coverage
sample1     1500000       12.3456
sample2     2000000       16.4608

# Usage
sbatch coverage_script.sh /path/to/bams /path/to/output.txt /path/to/reference.fasta 125
