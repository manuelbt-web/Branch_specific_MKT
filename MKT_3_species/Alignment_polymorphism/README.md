# Alignment_polymorphism — MACSE CDS Alignment (Polymorphism Track)

**Role in the pipeline:** Polymorphism track (runs in parallel with the
divergence track).
Takes trimmed, unaligned CDS sequences for the target species (typically
one FASTA per gene, multiple individuals or haplotypes) and produces a
codon-level multiple-sequence alignment using MACSE `alignSequences`.
The output serves as the fixed polymorphism alignment in the MACSE
`enrichAlignment` step (`../EnrichAlignment/`).

---

## Contents

| Script | Type | Description |
|--------|------|-------------|
| `MACSE_aligning_sequences.sbatch` | SLURM array | MACSE `alignSequences` on trimmed CDS FASTAs |

---

## Dependencies

| Tool | Version | Installation |
|------|---------|--------------|
| **MACSE** | ≥ 2.07 | Download `macse_vX.XX.jar` from [agap-ge2pop.org/macse](https://www.agap-ge2pop.org/macse/) |
| **Java** | ≥ 11 | `conda install -c conda-forge openjdk` or cluster module |

---

## Input

A directory of trimmed, unaligned nucleotide CDS FASTA files.
Typically the `NT/` subdirectory produced by a previous MACSE
`trimNonHomologousFragments` step.

```
<input_dir>/
├── gene_A.fasta    (one file per gene; each file = multiple CDS sequences)
├── gene_B.fasta
└── ...
```

Each file should contain CDS sequences from the intraspecific polymorphism
dataset (multiple individuals, accessions, or haplotypes of the target species).

---

## Output (per gene, in `<output_dir>/`)

| File | Description |
|------|-------------|
| `<gene>_aligned_AA.fasta` | MACSE amino acid alignment |
| `<gene>_aligned_NT.fasta` | MACSE nucleotide alignment (frameshifts coded as `!`) |

The `_aligned_NT.fasta` files are the primary output used in `EnrichAlignment/`.

**Next step:** pass `<output_dir>/` as `--poly-dir` to
`../EnrichAlignment/MACSE_enrichment.py`.

---

## Usage

```bash
# 1. Count input files
N=$(ls INPUT_DIR/NT/*.fasta | wc -l)

# 2. Submit as SLURM array
sbatch --array=0-$((N-1))%15 MACSE_aligning_sequences.sbatch \
    -i INPUT_DIR/NT/ \
    -o cds_aligned/ \
    -j /path/to/macse_v2.07.jar

# Or run sequentially (no SLURM; will default to task ID 0 — runs one file)
# For sequential processing of all files, wrap in a loop:
for F in INPUT_DIR/NT/*.fasta; do
    N=$(ls INPUT_DIR/NT/*.fasta | wc -l)
    for i in $(seq 0 $((N-1))); do
        SLURM_ARRAY_TASK_ID=$i bash MACSE_aligning_sequences.sbatch \
            -i INPUT_DIR/NT/ -o cds_aligned/ -j /path/to/macse_v2.07.jar
    done
done
```

---

## Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-i DIR` | Yes | — | Input directory with trimmed `.fasta` CDS files |
| `-o DIR` | Yes | — | Output directory (created if absent) |
| `-j FILE` | Yes | — | Path to MACSE JAR file |
| `-h` | No | — | Show usage and exit |

---

## Cluster configuration

Lines marked `[CLUSTER]` in the script must be adapted:

| Line | What to change |
|------|----------------|
| `#SBATCH --partition=cpu-dedicated` | Your cluster partition |
| `#SBATCH --account=dedicated-cpu@cirad-long` | Your account |
| `#SBATCH --mem=30G` | Memory (MACSE can be memory-intensive for large alignments) |
| `#SBATCH --time=84:00:00` | Wall-time (long genes may take several hours) |
| `module load bioinfo-ifb` | Your cluster's module environment |
| `module load java-jdk` | Your Java module name |

---

## Important notes

- MACSE represents frameshift-corrected positions as `!` in the NT output.
  These are handled downstream by `filter_codon_alignment.py`, which treats
  `!` as an invalid character and replaces affected codons with `NNN`.
- The `_aligned_NT.fasta` files are the **fixed reference alignment** (`-align`)
  in the `enrichAlignment` step. Their column structure must be preserved.
- If outputs already exist for a gene, the job exits silently (idempotent).

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `java not found` | Java module not loaded | `conda activate openjdk_env` or load Java module |
| `OutOfMemoryError` | Gene alignment too large | Increase `#SBATCH --mem` (try 60G) |
| Empty `_aligned_NT.fasta` | MACSE encountered irreconcilable frameshifts | Inspect MACSE stderr (`log_*_err.txt`); gene may be unsuitable |
| Task ID ≥ number of files | Array submitted with wrong N | Recount: `ls *.fasta \| wc -l` |

---

## Navigation

← [MKT_3_species overview](../README.md)
→ [EnrichAlignment/ — combine with divergence alignment](../EnrichAlignment/README.md)
