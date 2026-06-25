# filter_codon_alignment — Per-Codon Quality Filter and Sequence Splitting

**Role in the pipeline:** Final preprocessing steps before codeml or MKT analysis.

This folder contains two scripts that are run sequentially:

1. **`filter_codon_alignment.py`** — removes codon positions that are
   uninformative (divergence sequences have gaps/ambiguities, or too few
   polymorphism sequences have valid codons). Failing codons are replaced with
   `NNN` so codeml treats them as missing data.

2. **`split_sequences_polymorphism_divergence.py`** — splits the filtered
   alignments into separate `divergence/` and `polymorphism/` FASTA files, one
   pair per orthogroup, ready for codeml and MKT analysis.

---

## Contents

| Script | Type | Description |
|--------|------|-------------|
| `filter_codon_alignment.py` | Python | Per-codon quality filter (replace failing codons with NNN) |
| `split_sequences_polymorphism_divergence.py` | Python | Split combined alignment into divergence and polymorphism files (3-species codeml track) |
| `divergence/divergence_two_species/trim_outgroup.py` | Python | Remove focal-species reference and second outgroup; keep 1 outgroup for standard MKT |

---

## Dependencies

| Tool | Version | Installation |
|------|---------|--------------|
| **Python** | ≥ 3.8 | system / conda (no third-party packages required) |

---

## Pipeline position

```
MACSE enrichAlignment output
  macse_enriched/enriched/*_NT.fasta
        │
        ▼
  filter_codon_alignment.py        ← Step 1: remove uninformative codons
        │
        ▼
  filtered_aln/*_filtered.fasta    (NNN-masked, codon-frame intact; poly + 3 div seqs)
        │
        ├──► divergence/divergence_two_species/trim_outgroup.py  ← Step 2a
        │         ↓
        │    alignment_for_standard_MKT/     → standard 2-species MKT (EggLib)
        │         (poly + 1 outgroup; focal-species ref and second outgroup removed)
        │
        └──► split_sequences_polymorphism_divergence.py          ← Step 2b
                  │
                  ├── split_out/divergence/*_divergence.fasta    → codeml dN/dS (3-species)
                  └── split_out/polymorphism/*_polymorphism.fasta → imputed MKT (EggLib)
```

---

## Script 1 — `filter_codon_alignment.py`

### Filtering criteria (per codon = every 3 columns)

1. **Divergence — strict**: ALL divergence sequences must have a fully valid
   codon (A/C/G/T only; no gaps `-`, no ambiguous bases, no MACSE frameshifts `!`)
2. **Polymorphism — minimum depth**: at least `--min-polym` sequences (default: 6)
   must have a valid codon

Codons failing either criterion are replaced with:
- `NNN` (default, **recommended**): codeml treats N as missing data; frame preserved
- `---` (optional): gap — may cause codeml to skip the column

### Why NNN and not gaps?

PAML / codeml handles `N` (IUPAC: any nucleotide) as missing data and
excludes `NNN` codons from likelihood calculations while preserving the
alignment frame. Gap codons `---` can cause column-removal issues in some
codeml settings and may confuse downstream parsers.

### Divergence sequence identification

| Option | Example | Notes |
|--------|---------|-------|
| `--divergence-pattern REGEX` | `"Aegilops_tauschii\|Triticum_urartu"` | Recommended — robust to variable orthogroup sizes |
| `--divergence LIST` | `"Aet_seq,Tu_seq,Am_seq"` | Exact headers, comma-separated |
| `--n-divergence N` | `3` (default) | Use last N sequences — works when MACSE output is ordered (polymorphism first, divergence last) |

### Input

MACSE `enrichAlignment` output:
```
macse_enriched/enriched/
├── HOG0000001_<species>_NT.fasta
├── HOG0000002_<species>_NT.fasta
└── ...
```

### Output (in `<output-dir>/`)

| File | Description |
|------|-------------|
| `<gene>_filtered.fasta` | Alignment with failing codons replaced by NNN |
| `<gene>_codon_pass.tsv` | Per-codon report: index, pass/fail, reason, polym depth |
| `<gene>_pass_positions.txt` | 1-based indices of passing codons |
| `filter_summary.tsv` | Across-gene summary table (directory mode only) |

### Usage — single file

```bash
python filter_codon_alignment.py \
    -i  macse_enriched/enriched/HOG0000001_NT.fasta \
    -o  filtered_aln/
```

### Usage — directory (recommended)

```bash
# Minimal (last 3 sequences as divergence)
python filter_codon_alignment.py \
    --input-dir   macse_enriched/enriched/ \
    --output-dir  filtered_aln/

# With explicit divergence pattern (robust for variable orthogroup sizes)
python filter_codon_alignment.py \
    --input-dir            macse_enriched/enriched/ \
    --output-dir           filtered_aln/ \
    --divergence-pattern   "Aegilops_tauschii|Triticum_urartu|Amblyopyrum_muticum"

# Parallel processing
python filter_codon_alignment.py \
    --input-dir            macse_enriched/enriched/ \
    --output-dir           filtered_aln/ \
    --divergence-pattern   "Aegilops_tauschii|Triticum_urartu|Amblyopyrum_muticum" \
    --threads              8
```

### Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-i` / `--input FILE` | Yes* | — | Single input FASTA file |
| `--input-dir DIR` | Yes* | — | Directory of FASTA files (batch mode) |
| `-o` / `--output-dir DIR` | Yes | — | Output directory (created if absent) |
| `--divergence-pattern REGEX` | No† | — | Regex matched against headers to identify divergence sequences |
| `--divergence LIST` | No† | — | Comma-separated exact divergence sequence headers |
| `--n-divergence N` | No | `3` | Use last N sequences as divergence (fallback) |
| `--min-polym N` | No | `6` | Minimum polymorphism sequences with valid codon per site |
| `--replace-with` | No | `N` | Replacement: `N` → NNN (recommended), `-` → gaps |
| `-e` / `--ext EXT` | No | fa/fasta | File extension filter for `--input-dir` |
| `-t` / `--threads N` | No | `1` | Parallel workers (directory mode) |
| `--verbose` | No | off | Print per-gene statistics |

*Exactly one of `--input` or `--input-dir` is required.
†`--divergence-pattern` and `--divergence` are mutually exclusive.

### Console output example

```
============================================================
  Codon Alignment Filter — Summary
============================================================
  Genes processed      : 234
  Replace failing with : NNN
  Divergence           : pattern 'Aegilops_tauschii|Triticum_urartu|...'
  Min polymorphism     : 6 valid codons per site

  Fraction of codons retained (per gene):
    Mean   : 0.783   Median : 0.801   Min : 0.124   Max : 0.999

  Genes with ≥ 50% codons retained : 210/234
  Genes with <  10% codons retained :   8/234

  Output directory     : filtered_aln/
  Summary table        : filtered_aln/filter_summary.tsv
============================================================
```

### `filter_summary.tsv` columns

| Column | Description |
|--------|-------------|
| `gene` | Gene / orthogroup ID |
| `n_codons` | Total codons |
| `n_pass` | Codons retained |
| `n_fail_div` | Codons removed — divergence invalid |
| `n_fail_poly` | Codons removed — polymorphism depth too low |
| `frac_pass` | Fraction retained |
| `n_div` | Number of divergence sequences |
| `n_poly` | Number of polymorphism sequences |
| `error` | Error message if any |

---

## Script 2 — `split_sequences_polymorphism_divergence.py`

### What it does

Takes the filtered FASTA files (output of `filter_codon_alignment.py`) and
writes two separate files per gene:
- `divergence/<gene>_divergence.fasta` — divergence sequences (multi-species orthologs)
- `polymorphism/<gene>_polymorphism.fasta` — polymorphism sequences (intraspecific)

These two tracks are then used independently:
- Divergence FASTAs → codeml to estimate dN/dS
- Polymorphism FASTAs → EggLib or similar for Pn/Ps estimation

Output subdirectories (`divergence/` and `polymorphism/`) are created
**automatically** inside the output directory.

### Divergence sequence identification

Same options as `filter_codon_alignment.py`:

| Option | Example | Notes |
|--------|---------|-------|
| `--divergence-pattern REGEX` | `"Aegilops_tauschii\|Triticum_urartu"` | Recommended |
| `--divergence LIST` | `"Aet_seq,Tu_seq,Am_seq"` | Exact headers |
| `--n-divergence N` | `3` (default) | Last N sequences |

### Input

Filtered FASTA files from `filter_codon_alignment.py`:
```
filtered_aln/
├── HOG0000001_NT_filtered.fasta
├── HOG0000002_NT_filtered.fasta
└── ...
```

### Output (auto-created inside `<output-dir>/`)

```
split_out/
├── divergence/
│   ├── HOG0000001_divergence.fasta
│   ├── HOG0000002_divergence.fasta
│   └── ...
├── polymorphism/
│   ├── HOG0000001_polymorphism.fasta
│   ├── HOG0000002_polymorphism.fasta
│   └── ...
└── split_summary.tsv
```

The `_NT_filtered` suffix is stripped automatically from input filenames when
generating output names. Use `--strip-suffix` to change or disable this.

### Usage — single file

```bash
python split_sequences_polymorphism_divergence.py \
    -i  filtered_aln/HOG0000001_NT_filtered.fasta \
    -o  split_out/
```

### Usage — directory (recommended)

```bash
# Minimal (last 3 sequences as divergence)
python split_sequences_polymorphism_divergence.py \
    --input-dir   filtered_aln/ \
    --output-dir  split_out/

# With explicit divergence pattern (recommended)
python split_sequences_polymorphism_divergence.py \
    --input-dir            filtered_aln/ \
    --output-dir           split_out/ \
    --divergence-pattern   "Aegilops_tauschii|Triticum_urartu|Amblyopyrum_muticum"
```

### Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-i` / `--input FILE` | Yes* | — | Single input FASTA file |
| `--input-dir DIR` | Yes* | — | Directory of FASTA files (batch mode) |
| `-o` / `--output-dir DIR` | Yes | — | Output directory; `divergence/` and `polymorphism/` created inside |
| `--divergence-pattern REGEX` | No† | — | Regex to identify divergence sequences |
| `--divergence LIST` | No† | — | Comma-separated exact divergence headers |
| `--n-divergence N` | No | `3` | Use last N sequences as divergence (fallback) |
| `-e` / `--ext EXT` | No | fa/fasta | File extension filter for `--input-dir` |
| `--strip-suffix SUFFIX` | No | `_NT_filtered` | Suffix stripped from input stem for output naming |
| `--overwrite` | No | off | Overwrite existing outputs (default: skip) |
| `--verbose` | No | off | Print one line per gene |

*Exactly one of `--input` or `--input-dir` is required.
†Mutually exclusive.

### Console output example

```
Found 234 file(s) in filtered_aln/
Divergence     : pattern 'Aegilops_tauschii|Triticum_urartu|Amblyopyrum_muticum'
Output         : split_out/  (divergence/ and polymorphism/ subdirectories)

============================================================
  Sequence Split — Summary
============================================================
  Files processed         : 234
  Files split             : 234
  Divergence spec         : pattern '...'
  Mean divergence seqs    : 3.0
  Mean polymorphism seqs  : 28.4

  divergence/   : split_out/divergence/
  polymorphism/ : split_out/polymorphism/
  Summary table : split_out/split_summary.tsv
============================================================
```

### `split_summary.tsv` columns

| Column | Description |
|--------|-------------|
| `gene` | Gene / orthogroup ID |
| `n_div` | Number of divergence sequences written |
| `n_poly` | Number of polymorphism sequences written |
| `n_total` | Total sequences in input |
| `warn` | Warning message (empty if clean) |
| `error` | Error message if any |

---

## Full two-step example

```bash
MACSE_OUT=macse_enriched/enriched/
DIV_PATTERN="Aegilops_tauschii|Triticum_urartu|Amblyopyrum_muticum"

# Step 1 — filter uninformative codons
python filter_codon_alignment.py \
    --input-dir            "$MACSE_OUT" \
    --output-dir           filtered_aln/ \
    --divergence-pattern   "$DIV_PATTERN" \
    --min-polym            6 \
    --threads              8

# Review: check filtered_aln/filter_summary.tsv for genes with low retention

# Step 2 — split into divergence and polymorphism tracks
python split_sequences_polymorphism_divergence.py \
    --input-dir            filtered_aln/ \
    --output-dir           split_out/ \
    --divergence-pattern   "$DIV_PATTERN"

# Result:
#   split_out/divergence/    → feed to codeml for dN/dS
#   split_out/polymorphism/  → feed to EggLib / MKT scripts for Pn/Ps
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `--divergence-pattern matched no headers` | Wrong species name or regex | `grep '^>' input.fasta \| head` to print headers; test regex |
| All codons fail — divergence | Divergence sequences have many `!`/`?` | Check MACSE output quality; inspect per-gene `.fasta` |
| All codons fail — polym depth | Too few polymorphism sequences | Lower `--min-polym` |
| 0 polymorphism sequences in split | All sequences matched as divergence | Check `--divergence-pattern` or `--n-divergence` |
| Alignment length not multiple of 3 | Upstream issue in PAL2NAL or MACSE | Re-run problematic genes through PAL2NAL |
| Outputs already exist warning | Previous partial run | Add `--overwrite` to reprocess |

---

## Navigation

← [EnrichAlignment/ — MACSE enrichAlignment](../README.md)
→ codeml (PAML) dN/dS using `split_out/divergence/`
→ MKT analysis using `split_out/polymorphism/`
