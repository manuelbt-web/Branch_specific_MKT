# Step 3 — Codon-Aware Alignment (MACSE)

## Objective

Produce a codon-aware multiple sequence alignment for each gene using MACSE's
`alignSequences` program.

Unlike standard tools (MAFFT, MUSCLE), MACSE aligns nucleotide sequences
while explicitly modelling the reading frame and handling frameshifts.
Frameshift positions are marked with `!` in the NT alignment, preserving
full codon information for downstream dN/dS analysis.  The AA alignment
is also produced for quality inspection.

This step runs as a **SLURM array job** — one job per gene in parallel.

---

## Script

`MACSE_aligning_sequences.sbatch`

---

## Dependencies

| Tool | Version | Installation |
|---|---|---|
| **Java** | ≥ 8 | `conda install -c conda-forge openjdk` or cluster module |
| **MACSE** | v2.07 (inside MACSE_V2_PIPELINES v12.01) | see [Step 2 README](../step2_MACSE_trim/README.md) |
| **SLURM** | any | provided by the HPC cluster |

MACSE installation is shared with Step 2 — install once, reuse the same JAR.

---

## Input files

| Input | Description | Source |
|---|---|---|
| `-i INPUT_DIR` | Directory with trimmed nucleotide FASTA files | `step2_out/NT/` from Step 2 |
| `-j MACSE_JAR` | Path to the MACSE JAR file | Same JAR used in Step 2 |

### Expected FASTA format

Each file (`<gene>_NT.fasta`) must contain:
- Unaligned trimmed nucleotide CDS sequences
- One record per focal individual + one outgroup record

```
>GENE001|Ae_speltoides|ind1
ATGGCT---AGCCTG
>GENE001|Ae_speltoides|ind2
ATGGCTAGCAGCCTG
>HORVU.Morex.r3.1HG0000010.1
ATGGCTAGCAGCCTG
```

---

## Output files

| File | Description |
|---|---|
| `OUTPUT_DIR/<gene>_aligned_NT.fasta` | Codon-aware NT alignment — **input for Step 4** |
| `OUTPUT_DIR/<gene>_aligned_AA.fasta` | AA alignment (for quality inspection) |

### The `!` character in NT alignments

MACSE marks positions where a frameshift was introduced to maintain the
reading frame with `!`.  Downstream tools (dNdSpNpS) must be able to
handle this character.  Do not remove it — it carries alignment information.

```
>GENE001|Ae_speltoides|ind1
ATG---GCT!AGC
>HORVU.Morex.r3.1HG0000010.1
ATGCCCGCTAGC--
```

---

## Cluster configuration

Before submitting, edit the `[CLUSTER]` lines in the script:

```
#SBATCH --account=...     # [CLUSTER] your billing account
#SBATCH --partition=...   # [CLUSTER] your queue/partition name
module load bioinfo-ifb   # [CLUSTER] your cluster's bioinfo module
module load java-jdk      # [CLUSTER] your Java module name
```

---

## Usage

```bash
N=$(ls step2_out/NT/*.fasta | wc -l)
sbatch --array=0-$((N-1))%15 step3_MACSE_align/MACSE_aligning_sequences.sbatch \
    -i step2_out/NT/ \
    -o step3_out/ \
    -j "$MACSE_JAR"
```

### All options

| Option | Required | Default | Description |
|---|---|---|---|
| `-i DIR` | Yes | — | Input directory with trimmed FASTA files (`NT/` from Step 2) |
| `-o DIR` | Yes | — | Output directory for alignments (created if absent) |
| `-j FILE` | Yes | — | Path to MACSE JAR |
| `-h` | No | — | Show help |

### Monitor job progress

```bash
squeue -u $USER
ls step3_out/ | wc -l        # number of completed alignments (*_aligned_NT.fasta)
```

---

## Resource tuning

| Resource | Default | When to change |
|---|---|---|
| `--mem=30G` | 30 GB | Increase for genes with > 50 individuals or > 5 kb CDS |
| `--time=84:00:00` | 84 h | `alignSequences` is slower than trimming; reduce after benchmarking |
| `%15` | 15 simultaneous | Adjust to cluster policy |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Output file empty | Input sequences already aligned (e.g. wrong input directory) | Ensure `-i` points to `step2_out/NT/` (unaligned) |
| `OutOfMemoryError` | Java heap too small | Increase `--mem` and `JAVA_MEM` in the script |
| Alignment very slow | Large number of sequences | Consider splitting genes into batches or increasing `--cpus-per-task` |
| `!` characters cause downstream errors | Tool does not handle MACSE's frameshift notation | Pre-process with `sed 's/!/-/g'` or use a dNdSpNpS version that supports MACSE output |

---

## Next step

→ [Step 4: Prepare headers for dNdSpNpS](../step4_prep_dNdS/README.md)

```
step3_out/
├── GENE001_NT_aligned_NT.fasta   ← codon-aware NT alignment
├── GENE001_NT_aligned_AA.fasta
├── GENE002_NT_aligned_NT.fasta
└── ...
```

Pass `step3_out/` as `--input-dir` to
`step4_prep_dNdS/input_preparation_for_dNdSpNpS.py`.
