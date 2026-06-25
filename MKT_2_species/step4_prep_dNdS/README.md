# Step 4 — Prepare Alignments for dNdSpNpS

## Objective

Standardise the FASTA headers of the codon-aligned sequences so that the
dNdSpNpS analysis pipeline can correctly identify focal individuals and the
outgroup sequence in each alignment file.

For each aligned FASTA, the script:
1. Separates focal sequences (members of the species under study) from the
   outgroup sequence.
2. Writes focal sequences first, outgroup sequence last (required by
   dNdSpNpS to correctly compute polymorphism vs. divergence).
3. Renames sequence headers to the standard format:
   - **Focal**: `>GENE001|Ae_speltoides|ind1`
     (the `|sp|` placeholder in MACSE-output headers is replaced by the
     actual species name)
   - **Outgroup**: `>GENE001|H_vulgare|original_id|outgroup`

The script can process files in-place or write to a separate output directory,
and supports parallelism via `--threads`.

---

## Script

`input_preparation_for_dNdSpNpS.py`

---

## Dependencies

| Tool | Version | Installation |
|---|---|---|
| **Python** | ≥ 3.7 | system / conda |
| **BioPython** | ≥ 1.79 | `pip install biopython` |

### Install BioPython

```bash
# Option A — pip:
pip install biopython

# Option B — conda:
conda install -c conda-forge biopython

# Verify:
python -c "import Bio; print(Bio.__version__)"
```

---

## Input files

| Input | Description | Source |
|---|---|---|
| `--input-dir DIR` | Directory with codon-aligned FASTA files | `step3_out/` from Step 3 |

### Expected FASTA format

Each file is a codon-aligned FASTA from MACSE (`*_aligned_NT.fasta`).
Focal sequences must contain the token `|sp|` in their ID (added by the
MACSE pipeline conventions) or share a common gene-ID prefix.

```
>GENE001|sp|ind1
ATG---GCT!AGC
>GENE001|sp|ind2
ATGGCCCGT-AGC
>HORVU.Morex.r3.1HG0000010.1
ATGCCCGCTAGC-
```

If `|sp|` is absent from all headers, the script falls back to using the most
common first-field prefix (before `|`) as the focal group identifier.

---

## Output files

| File | Description |
|---|---|
| `OUTPUT_DIR/<gene>.fasta` | Reordered FASTA: focal sequences first, outgroup last, standardised headers |
| `OUTPUT_DIR/prep_dndspnps.log` | Processing log — one line per file (OK / WARN / ERROR) |

When `--output-dir` is omitted, files are modified **in-place** inside
`--input-dir`.

### Output header format

```
>GENE001|Ae_speltoides|ind1       ← focal: gene|species|individual
>GENE001|Ae_speltoides|ind2
>GENE001|H_vulgare|HORVU_...|outgroup   ← outgroup: gene|outgroup|original_id|outgroup
```

This format is read directly by the dNdSpNpS analysis scripts to distinguish
focal individuals (used to compute polymorphism, *pN* and *pS*) from the
outgroup (used to compute divergence, *dN* and *dS*).

---

## Usage

```bash
python step4_prep_dNdS/input_preparation_for_dNdSpNpS.py \
    --input-dir  step3_out/ \
    --species    Ae_speltoides \
    --outgroup   H_vulgare \
    --output-dir step4_out/ \
    --threads    4
```

### All options

| Option | Required | Default | Description |
|---|---|---|---|
| `--input-dir DIR` | Yes | — | Directory with aligned FASTA files |
| `--species NAME` | Yes | — | Focal species name; replaces `\|sp\|` token in headers |
| `--outgroup NAME` | Yes | — | Outgroup species name inserted into outgroup headers |
| `--output-dir DIR` | No | in-place | Write results here instead of modifying input files |
| `--log FILE` | No | `OUTPUT_DIR/prep_dndspnps.log` | Log file (appended) |
| `--threads N` | No | `1` | Number of parallel worker processes |

---

## Header naming conventions

| Field | Description | Example |
|---|---|---|
| Field 1 (gene prefix) | Gene identifier, shared by all sequences in the file | `GENE001` |
| Field 2 (species/outgroup) | Species name (focal) or outgroup name | `Ae_speltoides` / `H_vulgare` |
| Field 3 (individual/original ID) | Individual identifier or original outgroup sequence ID | `ind1` / `HORVU_Morex_r3_1HG0000010_1` |
| Field 4 (role) | `outgroup` only for outgroup records | `outgroup` |

Fields are separated by `|`.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `No .fasta files found` | Wrong `--input-dir` or files have different extension | Check `ls --input-dir/*.fasta` |
| `biopython` not found | BioPython not installed | `pip install biopython` |
| `No outgroup sequence found` | No non-focal sequences in the FASTA | Check Step 1 output — outgroup may be missing |
| `|sp|` not replaced | Focal headers do not contain `\|sp\|` | Script falls back to prefix matching; check log |
| Output log shows `[WARN]` for many genes | Mixed header formats | Standardise headers before running, or open an issue |

---

## Previous step

← [Step 3: MACSE alignment](../step3_MACSE_align/README.md)

The aligned FASTAs produced by Step 3 are the direct input to this script.
After Step 4, the standardised FASTA files in `step4_out/` are ready for
the dNdSpNpS analysis.
