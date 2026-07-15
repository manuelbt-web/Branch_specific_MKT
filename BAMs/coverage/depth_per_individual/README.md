# BAM Depth Merger

A parallel processing pipeline that converts multiple BAM files into a unified depth matrix.

### BAM Files Directory
- Contains one or more `.bam` files
- Files should be properly aligned and indexed
- Naming convention: `[sample_name].bam`

## Output
### Final Output File 
`depth_per_position_all_individuals_BAMs.depth`: 
- Tab-delimited matrix 
- Header: `contig pos sample1 sample2 ... sampleN` 
- Each row shows depth at a genomic position

## Workflow
1. **Depth Calculation**
   - Runs `samtools depth` on each BAM file
   - Creates temporary `.depth` files 

2. **Format Standardization**
   - Ensures consistent tab-delimited format
   - Creates temporary reformatted versions

3. **Matrix Construction** 
   - Merges all samples' depth data
   - Preserves genomic coordinates
   - Adds sample-specific depth columns

4. **Cleanup**  
   - Removes all temporary files
   - Keeps only the final merged output

## Usage

**Cluster:** SLURM — edit the `[CLUSTER]` lines in `depth_per_individual.sbatch`
before submission (account, partition, module names).

### Basic Command
```bash
sbatch depth_per_individual.sbatch /path/to/bams /path/to/output_dir
```

Can also be run directly with bash (e.g. on a non-SLURM system, after
removing the `#SBATCH` header and `module load` lines):
```bash
bash depth_per_individual.sbatch /path/to/bams /path/to/output_dir
```

