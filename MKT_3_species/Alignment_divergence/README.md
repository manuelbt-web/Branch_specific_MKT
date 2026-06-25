# Alignment_divergence — MAFFT Protein Alignment

**Role in the pipeline:** Step 1 of the divergence track.
Takes single-copy ortholog protein sequences identified by OrthoFinder and
produces a multiple-sequence alignment (MSA) per orthogroup, which is then
cleaned in `clean_alignment/`.

---

## Contents

| Script | Type | Description |
|--------|------|-------------|
| `mafft_alignments.sbatch` | SLURM array | MAFFT protein MSA, one task per orthogroup |

---

## Dependencies

| Tool | Version | Installation |
|------|---------|--------------|
| **MAFFT** | ≥ 7.4 | `conda install -c bioconda mafft` or cluster module |

```bash
conda create -n mafft_env -c bioconda mafft
conda activate mafft_env
```

---

## Input

OrthoFinder `Single_Copy_Orthologue_Sequences/` directory.
Each `.fa` file contains one protein sequence per species for one single-copy
orthogroup (naming: `OG0000001.fa`).

Typical path after an OrthoFinder 2.x / 3.x run:
```
<orthofinder_outdir>/OrthoFinder/Results_DATE/Single_Copy_Orthologue_Sequences/
```

---

## Output

| File | Description |
|------|-------------|
| `<output_dir>/OG0000001_aln.fa` | MAFFT-aligned protein sequences (one per orthogroup) |
| `<output_dir>/<gene>_mafft.err` | MAFFT stderr, kept only if non-empty |

**Next step:** pass `<output_dir>/` to `clean_alignment/hmm_cleaner.sbatch`.

---

## Usage

```bash
# 1. Activate MAFFT
conda activate mafft_env

# 2. Count input files
N=$(ls -1 /path/to/Single_Copy_Orthologue_Sequences/*.fa | wc -l)

# 3. Submit as SLURM array (max 15 concurrent tasks)
sbatch --array=0-$((N-1))%15 mafft_alignments.sbatch \
       --input-dir  /path/to/Single_Copy_Orthologue_Sequences/ \
       --output-dir /path/to/mafft_out/

# 4. Or run sequentially without SLURM
bash mafft_alignments.sbatch \
     --input-dir  /path/to/Single_Copy_Orthologue_Sequences/ \
     --output-dir /path/to/mafft_out/
```

---

## Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-i` / `--input-dir DIR` | Yes | — | OrthoFinder Single_Copy_Orthologue_Sequences/ path |
| `-o` / `--output-dir DIR` | Yes | — | Output directory for aligned FASTA files |
| `-e` / `--ext EXT` | No | `fa` | Input file extension |
| `-t` / `--threads N` | No | `$SLURM_CPUS_PER_TASK` or 1 | MAFFT `--thread` value |
| `--mafft-opts OPTS` | No | `--auto` | Additional MAFFT algorithm options |
| `-h` / `--help` | No | — | Show usage and exit |

### MAFFT algorithm options

| Value | Algorithm | When to use |
|-------|-----------|-------------|
| `--auto` (default) | Automatic | General use; MAFFT selects the best method |
| `--linsi` | L-INS-i | High accuracy for short/medium alignments |
| `--ginsi` | G-INS-i | Globally similar sequences, few gaps |
| `--einsi` | E-INS-i | Conserved domains embedded in unconserved regions |

---

## Cluster configuration

Lines marked `[CLUSTER]` in the script must be adapted before submission:

| Line | What to change |
|------|----------------|
| `#SBATCH --partition=cpu-dedicated` | Your cluster partition |
| `#SBATCH --account=myaccount` | Uncomment and set if required |
| `#SBATCH --cpus-per-task=4` | Available CPUs (benefit saturates ~8 for MAFFT) |
| `#SBATCH --time=04:00:00` | Wall-time per task |
| `module load mafft/7.515` | Your cluster's MAFFT module name |

---

## Idempotency

Already-aligned files are skipped (the script checks for a non-empty `*_aln.fa`
before running MAFFT). Safe to re-submit the same command to resume after failures.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `mafft: command not found` | Module not loaded / conda env inactive | `conda activate mafft_env` or check `module load` |
| `No .fa files found` | Wrong `--input-dir` or extension | `ls Single_Copy_Orthologue_Sequences/*.fa` to verify |
| Array task out of range | N computed on wrong directory | Recount: `ls *.fa \| wc -l` |
| Some tasks failed | Timeout or memory | Re-submit unchanged; completed files are skipped |
| Non-empty `.err` files | Unusual sequences (very long, non-standard AA) | Inspect; consider `--einsi` for difficult genes |

---

## Navigation

← [MKT_3_species overview](../README.md)
→ [clean_alignment/ — HmmCleaner masking](clean_alignment/README.md)
