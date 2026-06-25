# Step 2 — Trim Non-Homologous Regions (MACSE)

## Objective

Remove sequence fragments that are not homologous across individuals before
alignment.  Non-homologous fragments (e.g. transposon insertions, assembly
artefacts, misassigned exons) cause spurious alignment columns and inflate
estimates of divergence.

MACSE's `trimNonHomologousFragments` program identifies and removes these
regions while respecting the codon reading frame.  It outputs:
- trimmed nucleotide sequences (used as input to Step 3)
- trimmed amino acid sequences (for quality inspection)
- per-gene trim statistics

This step runs as a **SLURM array job** — one job per gene in parallel.

---

## Script

`MACSE_trim_non_homologous.sbatch`

---

## Dependencies

| Tool | Version | Installation |
|---|---|---|
| **Java** | ≥ 8 | see below |
| **MACSE** | v2.07 (inside MACSE_V2_PIPELINES v12.01) | see below |
| **SLURM** | any | provided by the HPC cluster |

### Install Java

```bash
# Option A — conda:
conda install -c conda-forge openjdk

# Option B — cluster module (adapt version):
module load java-jdk   # [CLUSTER]
java -version          # verify
```

### Install MACSE (MACSE_V2_PIPELINES v12.01)

The MACSE JAR is bundled inside the MACSE_V2_PIPELINES repository archive:

```bash
wget https://github.com/ranwez/MACSE_V2_PIPELINES/archive/refs/tags/v12.01.tar.gz
tar xzf v12.01.tar.gz
MACSE_JAR=$(find MACSE_V2_PIPELINES-12.01/ -name "*.jar" | head -1)
echo "MACSE_JAR=$MACSE_JAR"   # keep this path for submission
```

---

## Input files

| Input | Description | Source |
|---|---|---|
| `-i INPUT_DIR` | Directory with per-gene FASTA files containing focal individuals + outgroup | Output of Step 1 (`step1_out/`) |
| `-j MACSE_JAR` | Path to the MACSE JAR file | Extracted from MACSE_V2_PIPELINES v12.01 |

### Expected FASTA format

Each FASTA file (`<gene>.fasta`) must contain:
- Unaligned nucleotide CDS sequences
- One sequence per focal individual + one outgroup sequence
- All sequences for the same gene (same reading frame)

```
>GENE001|Ae_speltoides|ind1
ATGGCTAGC...
>GENE001|Ae_speltoides|ind2
ATGGCTAGC...
>HORVU.Morex.r3.1HG0000010.1
ATGGCTAGC...
```

---

## Output files

All outputs are written to subdirectories of `OUTPUT_DIR`:

| File | Description |
|---|---|
| `OUTPUT_DIR/NT/<gene>_NT.fasta` | Trimmed nucleotide sequences — **input for Step 3** |
| `OUTPUT_DIR/AA/<gene>_AA.fasta` | Trimmed amino acid sequences |
| `OUTPUT_DIR/trim_stats/<gene>_trim_stats.csv` | Per-sequence homology statistics |
| `OUTPUT_DIR/mask_detail/<gene>_mask_detail_NT.fasta` | NT with per-position masking detail |

The `NT/` subdirectory is the critical output.  Pass it as `-i` to Step 3.

---

## Cluster configuration

Before submitting, edit the `[CLUSTER]` lines inside the script:

```
#SBATCH --account=...     # [CLUSTER] your billing account
#SBATCH --partition=...   # [CLUSTER] your queue/partition name
module load bioinfo-ifb   # [CLUSTER] your cluster's bioinfo module
module load java-jdk      # [CLUSTER] your Java module name
```

Common examples:

| Cluster | account | partition | Java module |
|---|---|---|---|
| IFB | `<project>@cpu` | `cpu-dedicated` | `java-jdk` |
| CIRAD | `dedicated-cpu@cirad-long` | `cpu-dedicated` | `java-jdk` |

---

## Usage

Calculate the number of genes, then submit with `--array`:

```bash
N=$(ls step1_out/*.fasta | wc -l)
sbatch --array=0-$((N-1))%15 step2_MACSE_trim/MACSE_trim_non_homologous.sbatch \
    -i step1_out/ \
    -o step2_out/ \
    -j "$MACSE_JAR"
```

The `%15` limits simultaneous tasks to 15 (adjust to match cluster policy).

### All options

| Option | Required | Default | Description |
|---|---|---|---|
| `-i DIR` | Yes | — | Input directory with per-gene FASTA files |
| `-o DIR` | Yes | — | Output directory (created if absent) |
| `-j FILE` | Yes | — | Path to MACSE JAR |
| `-m FLOAT` | No | `0.6` | Min fraction of homologous sites to retain a sequence |
| `-h` | No | — | Show help |

### Monitor job progress

```bash
squeue -u $USER              # list running jobs
sacct -j <JOBID> --format=JobID,State,ExitCode   # job accounting
ls step2_out/NT/ | wc -l    # count completed genes
```

---

## Resource tuning

| Resource | Default | When to change |
|---|---|---|
| `--mem=30G` | 30 GB | Increase to 60G+ for very long genes (> 5 kb) or > 50 individuals |
| `--time=24:00:00` | 24 h | Sufficient for typical transcriptomes; reduce if most jobs finish quickly |
| `%15` (array throttle) | 15 | Increase on clusters with many cores; decrease if I/O-limited |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `java not found` | Java module not loaded | Adapt the `module load java-jdk` line |
| `MACSE JAR not found` | Wrong path passed to `-j` | Re-check `MACSE_JAR=$(find ... -name "*.jar")` |
| Output NT file is empty | All sequences filtered (homology < `-m`) | Lower `-m` threshold (e.g. `-m 0.4`) or inspect the trim_stats CSV |
| Jobs fail with `OutOfMemoryError` | Java heap too small | The script allocates `-Xmx24g`; increase `--mem` in `#SBATCH` and the `JAVA_MEM` variable |

---

## Next step

→ [Step 3: MACSE alignment](../step3_MACSE_align/README.md)

```
step2_out/NT/
├── GENE001_NT.fasta   ← trimmed nucleotide sequences
├── GENE002_NT.fasta
└── ...
```

Pass `step2_out/NT/` as `-i` to `step3_MACSE_align/MACSE_aligning_sequences.sbatch`.
