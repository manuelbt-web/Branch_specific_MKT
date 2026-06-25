# Step 5 — Run dNdSpiNpiS (Divergence Statistics)

## Objective

Compute divergence statistics (dN, dS, pN, pS, dN/dS, pN/pS, DoS, NI) for
each gene by comparing the focal ingroup sequences against the outgroup using
dNdSpiNpiS v1.0.

This step runs as a **SLURM array job** — one task per gene — on the
codon-aligned FASTA files produced by Step 4 (standardised headers, outgroup
sequence present).

---

## Script

`dNdSpiNpiS.sbatch`

---

## Dependencies

| Tool | Version | Installation |
|---|---|---|
| **dNdSpiNpiS** | v1.0 | Download binary from PopPhyl (see below) |
| **Bio++ libraries** | — | Bundled in the dNdSpiNpiS binary |
| **SLURM** | any | Provided by the HPC cluster |

### Install the dNdSpiNpiS binary

1. Go to the PopPhyl tools page:
   `https://kimura.univ-montp2.fr/PopPhyl/index.php?section=tools`
2. Download `dNdSpiNpiS_1.0` (Linux x86-64 static binary).
3. Make it executable:

```bash
chmod +x /path/to/dNdSpiNpiS_1.0
# Verify:
/path/to/dNdSpiNpiS_1.0 --help
```

The binary is statically linked and requires no library installation.

---

## Input files

| Input | Description | Source |
|---|---|---|
| `-i INPUT_DIR` | Directory with per-gene codon-aligned FASTA files (standardised headers) | `step4_out/` from Step 4 |
| `-b BINARY` | Full path to the `dNdSpiNpiS_1.0` binary | Downloaded above |
| `-g INGROUP` | Ingroup species name | Must match the `\|species\|` field in FASTA headers |
| `-k OUTGROUP` | Outgroup species name | Must match the `\|outgroup\|` field in FASTA headers |

### Expected FASTA header format (from Step 4)

```
>GENE001|Ae_speltoides|ind1
>GENE001|Ae_speltoides|ind2
>GENE001|H_vulgare|HORVU_...|outgroup
```

With this example, use `-g Ae_speltoides -k H_vulgare`.

---

## Output files

| File | Description |
|---|---|
| `OUTPUT_DIR/<gene>.out` | Raw dNdSpiNpiS result: dN, dS, pN, pS, dN/dS, pN/pS, DoS, NI |
| `log_<jobID>_<taskID>_out.txt` | SLURM stdout per task |
| `log_<jobID>_<taskID>_err.txt` | SLURM stderr per task |

### Merge all genes into one table

After all array tasks complete:

```bash
# header from any completed output file
head -1 "$(ls step5_out/*.out | head -1)" > all_genes_dNdSpNpS.tsv
# append data rows from all genes
for f in step5_out/*.out; do tail -n +2 "$f"; done >> all_genes_dNdSpNpS.tsv
```

---

## Cluster configuration

Before submitting, edit the `[CLUSTER]` lines in the script:

```bash
#SBATCH --account=...    # [CLUSTER] your billing account
#SBATCH --partition=...  # [CLUSTER] your partition/queue name
```

| Cluster | Account example | Partition example |
|---|---|---|
| IFB | `<project>@cpu` | `cpu-dedicated` |
| CIRAD | `dedicated-cpu@cirad-normal` | `cpu-dedicated` |

---

## Usage

```bash
# Count FASTA files to set the array size
N=$(ls step4_out/*.fasta | wc -l)

sbatch --array=0-$((N-1))%15 step5_run_dNdSpNpS/dNdSpiNpiS.sbatch \
    -i step4_out/ \
    -o step5_out/ \
    -b /path/to/dNdSpiNpiS_1.0 \
    -g Ae_speltoides \
    -k H_vulgare
```

### All options

| Option | Required | Default | Description |
|---|---|---|---|
| `-i DIR` | Yes | — | Input directory with per-gene FASTA files (Step 4 output) |
| `-o DIR` | Yes | — | Output directory (created if absent) |
| `-b FILE` | Yes | — | Path to the `dNdSpiNpiS_1.0` binary |
| `-g NAME` | Yes | — | Ingroup species name (must match FASTA headers) |
| `-k NAME` | Yes | — | Outgroup species name (must match FASTA headers) |
| `-c INT` | No | `8` | Number of GC bins |
| `-n INT` | No | `10` | Max % of N/gap sites per position |
| `-t FLOAT` | No | `1` | Transition/transversion ratio kappa |
| `-h` | No | — | Show help |

### Monitor progress

```bash
squeue -u $USER
ls step5_out/*.out | wc -l     # number of completed genes
```

---

## Resource tuning

| Resource | Default | When to change |
|---|---|---|
| `--mem=10G` | 10 GB | Sufficient for most genes |
| `--time=06:00:00` | 6 h | Reduce after benchmarking on your data |
| `%15` | 15 simultaneous | Adjust to cluster policy |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `binary is not executable` | `chmod +x` not run | `chmod +x /path/to/dNdSpiNpiS_1.0` |
| Empty `.out` files | Species name mismatch in headers | Verify `-g`/`-k` match the exact strings in FASTA headers from Step 4 |
| `No .fasta files found` | Wrong `-i` path | Check `ls INPUT_DIR/*.fasta` |
| All jobs finish instantly with empty output | Binary for wrong architecture | Download the Linux x86-64 binary or recompile from source |

---

## Previous step

← [Step 4: Prepare headers for dNdSpNpS](../step4_prep_dNdS/README.md)

## Next step

→ [Step 6: Polymorphism statistics](../step6_polymorphism/README.md)

Note: Steps 5 and 6 are **independent** — Step 6 does **not** use the Step 4
output directly; it starts from Step 4 output but first removes the outgroup
(see Step 6 README).
