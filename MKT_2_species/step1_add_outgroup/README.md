# Step 1 — Add Outgroup Sequences

## Objective

For each gene, append the outgroup ortholog sequence to the multi-individual
focal-species FASTA file.  The outgroup is needed by the MKT to polarise
mutations as ancestral or derived and to anchor the alignment for MACSE.

The script looks up the outgroup sequence ID in a Reciprocal Best Hit (RBH)
table produced by the `Orthology/` pipeline, extracts the matching sequence
from the outgroup CDS FASTA, and writes a new FASTA containing all focal
individuals followed by the outgroup sequence.

**Input files are never modified.**  All output is written to a dedicated
output directory.

---

## Script

`add_outgroup.sh`

---

## Dependencies

| Tool | Version | Installation |
|---|---|---|
| **seqkit** | ≥ 2.0 | `conda install -c bioconda seqkit` |
| **bash** | ≥ 4.0 | system |

### Install seqkit

```bash
# Option A — conda (recommended):
conda install -c bioconda seqkit

# Option B — pre-built binary:
wget https://github.com/shenwei356/seqkit/releases/download/v2.8.1/seqkit_linux_amd64.tar.gz
tar xzf seqkit_linux_amd64.tar.gz
mv seqkit /usr/local/bin/

# On a cluster — adapt to your module name:
module load seqkit/2.8.1   # [CLUSTER]
```

---

## Input files

| Input | Description | Format |
|---|---|---|
| `-i INPUT_DIR` | Directory with per-gene CDS FASTA files | One `.fasta` per gene; each file contains one sequence per focal individual |
| `-g OUTGROUP_FASTA` | All CDS sequences of the outgroup species in a single FASTA file | Standard multi-FASTA |
| `-r RBH_TABLE` | Reciprocal Best Hit table from the `Orthology/` pipeline | TSV, column 1 = focal gene ID, column 2 = outgroup sequence ID |

### Expected FASTA header format (focal individuals)

```
>GENE001|Ae_speltoides|ind1
ATGGCTAGC...
>GENE001|Ae_speltoides|ind2
ATGGCTAGC...
```

Gene FASTA filename (e.g. `GENE001.fasta`) must match the gene ID in column 1
of the RBH table.

### Expected RBH table format

```
# gene_id    outgroup_id    score_A→B    score_B→A
GENE001      HORVU.Morex.r3.1HG0000010.1    1245.0    1187.0
GENE002      HORVU.Morex.r3.1HG0000025.1     987.0     943.0
```

Lines starting with `#` are ignored.

---

## Output files

All outputs are written to `OUTPUT_DIR/` (never modifies the input directory).

| File | Description |
|---|---|
| `OUTPUT_DIR/<gene>.fasta` | Original focal sequences + outgroup sequence appended at end |
| `OUTPUT_DIR/add_outgroup.log` | Genes skipped (no ortholog found or ID absent from FASTA) |

### Output FASTA header for the outgroup sequence

The outgroup sequence is appended with its original header from the outgroup
FASTA.  The `input_preparation_for_dNdSpNpS.py` script (Step 4) will
later rename it to the standard format.

---

## Usage

```bash
# Adapt module load [CLUSTER] lines inside the script, then:
bash add_outgroup.sh \
    -i cds_per_gene/ \
    -o step1_out/ \
    -g H_vulgare_CDS.fasta \
    -r RBH_Ae_speltoides_H_vulgare.tab
```

### All options

| Option | Required | Default | Description |
|---|---|---|---|
| `-i DIR` | Yes | — | Input directory with per-gene FASTA files |
| `-o DIR` | Yes | — | Output directory (created if absent) |
| `-g FILE` | Yes | — | Outgroup CDS FASTA |
| `-r FILE` | Yes | — | RBH table |
| `-l FILE` | No | `OUTPUT_DIR/add_outgroup.log` | Log file |
| `-h` | No | — | Show help |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `seqkit not found` | Module not loaded | Add `module load seqkit/2.8.1` before the script or in `[CLUSTER]` lines |
| Many genes in log with "not in FASTA" | ID mismatch between RBH col 2 and outgroup FASTA headers | Check that the RBH table was produced with the same outgroup FASTA |
| `seqkit faidx` fails on read-only filesystem | Index file cannot be created next to the FASTA | Copy the outgroup FASTA to a writable location |

---

## Next step

→ [Step 2: MACSE trim non-homologous regions](../step2_MACSE_trim/README.md)

```
step1_out/
├── GENE001.fasta   ← focal sequences + outgroup
├── GENE002.fasta
└── ...
```

Pass `step1_out/` as `-i` to `step2_MACSE_trim/MACSE_trim_non_homologous.sbatch`.
