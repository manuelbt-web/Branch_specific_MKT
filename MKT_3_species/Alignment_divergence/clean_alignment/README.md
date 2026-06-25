# clean_alignment — HmmCleaner Alignment Masking

**Role in the pipeline:** Step 2 of the divergence track.
Takes MAFFT protein alignments from `../` and masks poorly aligned or
low-confidence columns using HmmCleaner, producing cleaned alignments ready
for coverage assessment and back-translation in `translate_back_in_nt/`.

---

## Contents

| Script | Type | Description |
|--------|------|-------------|
| `hmm_cleaner.sbatch` | SLURM array | Run HmmCleaner on each protein alignment |

---

## What HmmCleaner does

For each alignment, HmmCleaner:
1. Builds a profile HMM with `hmmbuild` (HMMER)
2. Re-aligns sequences with `hmmalign`
3. Scores each column against the profile
4. Replaces low-scoring (unreliable) columns with `?` characters

This removes poorly aligned, low-information columns that would bias dN/dS
estimates or MKT calculations.

Reference: Ranwez et al. 2018, doi:10.1371/journal.pone.0151312

---

## Dependencies

| Tool | Version | Installation |
|------|---------|--------------|
| **HmmCleaner** | ≥ 1.8 | See below |
| **HMMER** | ≥ 3.3 | `conda install -c bioconda hmmer` |
| **Perl** | ≥ 5.20 | system |

### Option A — conda / bioconda (recommended)

```bash
conda create -n hmmcleaner_env -c bioconda hmmcleaner hmmer
conda activate hmmcleaner_env
HMMC=$(which HmmCleaner.pl)
```

### Option B — source install

```bash
git clone https://github.com/ranwez/HmmCleaner.git /path/to/HmmCleaner
# Use paths below:
#   --hmmcleaner-bin /path/to/HmmCleaner/bin/HmmCleaner.pl
#   --hmmcleaner-lib /path/to/HmmCleaner/lib
```

---

## Input

MAFFT-aligned protein FASTA files from `../mafft_alignments.sbatch`.

```
<input_dir>/OG0000001_aln.fa
<input_dir>/OG0000002_aln.fa
...
```

---

## Output (per orthogroup, in `<output_dir>/`)

| File | Description |
|------|-------------|
| `<gene>_aln_hmm.fasta` | Cleaned alignment — unreliable columns replaced with `?` |
| `<gene>_aln_hmm.log` | Number of masked positions per sequence |
| `<gene>_aln_hmm.score` | Per-position HMM alignment score |

**Next step:** pass `<output_dir>/` to
`translate_back_in_nt/alignment_coverage.py`.

---

## Usage

```bash
conda activate hmmcleaner_env
HMMC=$(which HmmCleaner.pl)

# 1. Count input files
N=$(ls -1 /path/to/mafft_out/*_aln.fa | wc -l)

# 2. Submit as SLURM array (max 15 concurrent)
sbatch --array=0-$((N-1))%15 hmm_cleaner.sbatch \
       --input-dir      /path/to/mafft_out/ \
       --output-dir     /path/to/hmmcleaner_out/ \
       --hmmcleaner-bin "$HMMC"

# 3. Source install: also pass --hmmcleaner-lib
sbatch --array=0-$((N-1))%15 hmm_cleaner.sbatch \
       --input-dir      /path/to/mafft_out/ \
       --output-dir     /path/to/hmmcleaner_out/ \
       --hmmcleaner-bin /path/to/HmmCleaner/bin/HmmCleaner.pl \
       --hmmcleaner-lib /path/to/HmmCleaner/lib

# 4. Run sequentially without SLURM
bash hmm_cleaner.sbatch \
     --input-dir      /path/to/mafft_out/ \
     --output-dir     /path/to/hmmcleaner_out/ \
     --hmmcleaner-bin "$HMMC"
```

---

## Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-i` / `--input-dir DIR` | Yes | — | MAFFT output directory (`*_aln.fa` files) |
| `-o` / `--output-dir DIR` | Yes | — | Output directory for cleaned files |
| `-b` / `--hmmcleaner-bin PATH` | Yes | — | Full path to `HmmCleaner.pl` |
| `-l` / `--hmmcleaner-lib DIR` | No | — | HmmCleaner `lib/` (source install only) |
| `-e` / `--ext EXT` | No | `fa` | Input file extension |
| `--pattern PAT` | No | `*_aln` | Input filename stem pattern |
| `-p` / `--perl PATH` | No | `perl` | Perl interpreter (if not in PATH) |
| `-h` / `--help` | No | — | Show usage and exit |

---

## Cluster configuration

| Line | What to change |
|------|----------------|
| `#SBATCH --partition=cpu-dedicated` | Your cluster partition |
| `#SBATCH --account=myaccount` | Uncomment if required |
| `#SBATCH --mem=4G` | Memory per task (2–4 GB usually sufficient) |
| `#SBATCH --time=04:00:00` | Wall-time per task |
| `module load hmmer/3.3.2` | Your cluster's HMMER module name |

---

## Idempotency

Already-cleaned files (`*_aln_hmm.fasta`) are skipped. Safe to re-submit.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `hmmbuild: command not found` | HMMER not on PATH | `conda activate hmmcleaner_env` or load HMMER module |
| `HmmCleaner.pl not found` | Wrong `--hmmcleaner-bin` | `which HmmCleaner.pl` after activating env |
| Perl error / `Can't locate ...` | Missing `lib/` | Add `--hmmcleaner-lib /path/to/HmmCleaner/lib` |
| Empty `*_aln_hmm.fasta` | All columns below threshold | Inspect `.score`; alignment may be too divergent |
| Array task out of range | N computed incorrectly | Recount: `ls mafft_out/*_aln.fa \| wc -l` |

---

## Navigation

← [Alignment_divergence/ — MAFFT alignment](../README.md)
→ [translate_back_in_nt/ — coverage filter + PAL2NAL](translate_back_in_nt/README.md)
