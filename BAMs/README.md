# BAM Soft-Clipping Filter

## Key Features
- Filters reads exceeding user-defined soft-clip thresholds
- Generates indexed output BAM files

## Usage

### Basic Command
```bash
sbatch Bams_soft_clipped_reads_threshold_20.sbatch <INPUT_DIR> <OUTPUT_DIR> [THRESHOLD]

Example
bash


# With custom threshold (20 soft-clipped bases)
sbatch Bams_soft_clipped_reads_threshold_20.sbatch /data/input_bams /data/output 20

# Using default threshold (1 soft-clipped base)
sbatch Bams_soft_clipped_reads_threshold_20.sbatch /data/input_bams /data/output

#Output Files

    Filtered BAM files (suffix: _no_soft_clipped_reads.bam)

    Corresponding BAM indexes (.bai)

    Log files in SLURM output directory

