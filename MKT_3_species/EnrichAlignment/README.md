# EnrichAlignment — MACSE Enrichment and Codon Alignment

**Role in the pipeline:** Step 4–5.
Combines the divergence nucleotide alignments (from
`../Alignment_divergence/`) with intraspecific polymorphism CDS sequences
(from `../Alignment_polymorphism/`) into a single per-gene codon alignment
using MACSE `enrichAlignment`.
The enriched alignments then pass through `filter_codon_alignment/` before
codeml or MKT analysis.

---

## Contents

| Script | Type | Description |
|--------|------|-------------|
| `orthogroup_table.py` | Python | Build HOG → gene-ID mapping table (run once on login node) |
| `MACSE_enrichment.py` | Python | Run MACSE `enrichAlignment` on all orthogroups |
| `macse_enrichment.sbatch` | SLURM | SLURM wrapper for `MACSE_enrichment.py` |

Run `orthogroup_table.py` first, then `macse_enrichment.sbatch` (or
`MACSE_enrichment.py` directly).

---

## Dependencies

| Tool | Version | Installation |
|------|---------|--------------|
| **Python** | ≥ 3.8 | system / conda |
| **pandas** | ≥ 1.3 | `conda install pandas` |
| **MACSE** | ≥ 2.07 | [agap-ge2pop.org/macse](https://www.agap-ge2pop.org/macse/) |
| **Java** | ≥ 11 | `conda install -c conda-forge openjdk` |

---

## Step 1 — `orthogroup_table.py`

### What it does

Reads all `*HOG*.fasta` files in the nucleotide alignment directory and
extracts the mapping `orthogroup → gene_ID per species` from FASTA headers
(`>species|gene_id`). Produces the TSV used by `MACSE_enrichment.py`.

### Input

PAL2NAL output directory containing files like:
```
nt_aln/HOG0000001_aln_hmm_nt.fasta
nt_aln/HOG0000002_aln_hmm_nt.fasta
...
```

Each FASTA header must follow the format: `>species_name|gene_id`

### Output

`orthogroup_table.tsv`:
```
orthogroup    Aegilops_tauschii    Triticum_urartu    Amblyopyrum_muticum
HOG0000001    Aet_00001.1          Tu_1G000001.1      Am_g00001
HOG0000002    Aet_00002.1          Tu_1G000002.1      Am_g00002
```

### Usage

```bash
# Build the table (fast — run on login node, no SLURM needed)
python orthogroup_table.py \
    --input-dir  nt_aln/ \
    --output     orthogroup_table.tsv

# Verify which species columns are available
python MACSE_enrichment.py \
    --ortho-table orthogroup_table.tsv \
    --list-species
```

### Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-i` / `--input-dir DIR` | Yes | — | Directory with `*HOG*.fasta` files |
| `-o` / `--output FILE` | No | `orthogroup_table.tsv` | Output TSV path |
| `-e` / `--ext EXT` | No | fasta/fa/faa | FASTA file extension |
| `--sep CHAR` | No | `\|` | Header separator between species and gene ID |

---

## Step 2 — `MACSE_enrichment.py` (via `macse_enrichment.sbatch`)

### What it does

For each orthogroup:
1. Looks up the target-species gene ID in the orthogroup table
2. Finds the codon alignment file in `--ortho-dir`
3. Finds the CDS polymorphism FASTA for that gene in `--poly-dir`
4. Runs MACSE `enrichAlignment`:
   - `-align` = polymorphism CDS alignment (fixed reference, from `Alignment_polymorphism/`)
   - `-seq` = orthogroup divergence sequences (inserted to fit the fixed alignment)
5. Validates that sequences were successfully added

### What is MACSE `enrichAlignment`?

MACSE `enrichAlignment` inserts new sequences into an existing alignment
while keeping the existing alignment structure fixed (`-fixed_alignment_ON`).
The polymorphism sequences define the codon frame; the divergence sequences
are inserted without disrupting it.

### Input

| Source | Argument | Format |
|--------|----------|--------|
| `orthogroup_table.py` output | `--ortho-table` | TSV |
| PAL2NAL output | `--ortho-dir` | `*HOG*.fasta` per orthogroup |
| MACSE `alignSequences` output | `--poly-dir` | `*_aligned_NT.fasta` per gene |

### Output (in `<output-dir>/`)

| File | Description |
|------|-------------|
| `enriched/*_NT.fasta` | Combined divergence + polymorphism NT alignment |
| `enriched/*_AA.fasta` | Combined AA alignment |
| `logs/enrichment_summary.tsv` | Per-gene result: files found, MACSE status, sequences added |
| `logs/failed_alignments.txt` | Orthogroups that failed with reason |
| `logs/rerun_failed.sh` | Re-run script for failed genes |
| `logs/macse_commands.log` | All MACSE commands (reproducibility log) |

**Next step:** pass `<output-dir>/enriched/` to
`filter_codon_alignment/filter_codon_alignment.py`.

### Usage — SLURM (recommended)

```bash
# 1. Edit the USER SETTINGS block in macse_enrichment.sbatch:
#    ORTHO_TABLE, ORTHO_DIR, POLY_DIR, MACSE_JAR, SPECIES_COL, OUTPUT_DIR, THREADS
nano macse_enrichment.sbatch

# 2. Test paths without running MACSE
bash macse_enrichment.sbatch --dry-run

# 3. Submit
sbatch macse_enrichment.sbatch

# 4. Monitor
squeue -u $USER
tail -f logs/enrichment_summary.tsv
```

### Usage — Python directly

```bash
# Show which species columns are available
python MACSE_enrichment.py \
    --ortho-table orthogroup_table.tsv \
    --list-species

# Run enrichment
python MACSE_enrichment.py \
    --ortho-table   orthogroup_table.tsv \
    --ortho-dir     nt_aln/ \
    --poly-dir      cds_aligned/ \
    --output-dir    macse_enriched/ \
    --macse-jar     /path/to/macse_v2.07.jar \
    --species-col   Aegilops_speltoides \
    --threads       8

# Dry run
python MACSE_enrichment.py ... --dry-run

# Re-run after failures (already-done genes are skipped)
python MACSE_enrichment.py ...   # same command, no changes needed
```

### Options (`MACSE_enrichment.py`)

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-t` / `--ortho-table FILE` | Yes | — | Orthogroup table TSV |
| `-d` / `--ortho-dir DIR` | Yes | — | Orthogroup NT alignment directory |
| `-P` / `--poly-dir DIR` | Yes | — | CDS polymorphism FASTA directory |
| `-o` / `--output-dir DIR` | Yes | — | Output directory |
| `-j` / `--macse-jar FILE` | Yes | — | Path to macse JAR |
| `-s` / `--species-col COL` | Yes | — | Species column from orthogroup table |
| `--list-species` | No | — | Print available species and exit |
| `--threads N` | No | `4` | Parallel MACSE workers |
| `--dry-run` | No | off | Log commands without running MACSE |
| `--overwrite` | No | off | Re-run even if output already exists |
| `--java-bin PATH` | No | `java` | Java executable path |

### File-matching patterns

**Orthogroup alignment** (in `--ortho-dir`), first match wins:
1. `{HOG_ID}_aln_hmm_codon.fasta`
2. `{HOG_ID}_aln_hmm_nt.fasta`
3. `{HOG_ID}*.fasta`
4. `{HOG_ID}*.fa`

**Polymorphism CDS** (in `--poly-dir`), first match wins:
1. `{gene_id}_CDS_aligned_NT.fasta`
2. `{gene_id}*.fasta`
3. `{gene_id}*.fa`

---

## Cluster configuration (`macse_enrichment.sbatch`)

The sbatch script has a `USER SETTINGS` block at the top — edit those
variables instead of the `#SBATCH` lines.

| Variable | Description |
|----------|-------------|
| `ORTHO_TABLE` | Path to `orthogroup_table.tsv` |
| `ORTHO_DIR` | Path to nucleotide alignment directory |
| `POLY_DIR` | Path to MACSE `alignSequences` output |
| `MACSE_JAR` | Path to macse JAR |
| `SPECIES_COL` | Column name from the orthogroup table |
| `OUTPUT_DIR` | Output directory |
| `THREADS` | Number of parallel MACSE jobs inside Python |

Lines marked `[CLUSTER]` in the script must also be adapted (partition,
account, module load).

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `--species-col not found` | Column name mismatch | Run `--list-species` |
| `ortho_found=NO` for all genes | Wrong `--ortho-dir` or naming | `ls nt_aln/*HOG*.fasta` |
| `poly_found=NO` for all genes | Wrong `--poly-dir` or gene ID | Check `{gene_id}_CDS_aligned_NT.fasta` pattern |
| `all_added=NO` | MACSE codon/stop issue | Check `logs/{gene}.log` |
| `java not found` | Java not in PATH | Load Java module or `conda install openjdk` |
| Re-run doesn't resume | `--overwrite` not set | Add `--overwrite` or run `bash logs/rerun_failed.sh` |

---

## Navigation

← [translate_back_in_nt/ — PAL2NAL](../Alignment_divergence/clean_alignment/translate_back_in_nt/README.md)
← [Alignment_polymorphism/ — MACSE CDS alignment](../Alignment_polymorphism/README.md)
→ [filter_codon_alignment/ — per-codon quality filter](filter_codon_alignment/README.md)
