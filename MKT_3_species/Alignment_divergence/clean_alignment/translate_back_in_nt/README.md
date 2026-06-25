# translate_back_in_nt — Coverage Assessment and PAL2NAL Back-Translation

**Role in the pipeline:** Step 3 of the divergence track.
Takes HmmCleaner-cleaned protein alignments, filters out low-quality or
incomplete ones, and back-translates the passing alignments to codon-level
nucleotide alignments using PAL2NAL.
The output is a directory of per-orthogroup nucleotide alignments ready for
MACSE enrichAlignment (`../../EnrichAlignment/`).

---

## Contents

| Script | Type | Description |
|--------|------|-------------|
| `alignment_coverage.py` | Python | Assess completeness; filter alignments by quality |
| `pal2nal.sbatch` | SLURM array | Back-translate protein alignments → codon NT alignments |

Run `alignment_coverage.py` first, then `pal2nal.sbatch` on the passing files.

---

## Dependencies

| Tool | Version | Installation |
|------|---------|--------------|
| **Python** | ≥ 3.8 | system / conda |
| **Perl** | ≥ 5.20 | system |
| **PAL2NAL** | newest | **auto-downloaded from GitHub on first run** |

PAL2NAL is downloaded automatically to `~/.local/share/pal2nal/pal2nal.pl`
on the first run and reused in all subsequent jobs.
No manual installation is required — only Perl must be available.

---

## Step A — `alignment_coverage.py`

### What it does

For each alignment, computes per-sequence completeness
(fraction of valid, non-gap, non-masked characters) and reports whether the
alignment meets minimum depth and length thresholds.

A gene **passes** if ALL of the following hold:
- Number of sequences ≥ `--min-sequences` (default: 2)
- Alignment length ≥ `--min-aln-length` (default: 30)
- Fraction of sequences with completeness ≥ `--min-completeness` (default: 0.50)
  is at least `--min-fraction` (default: 0.80)

### Input

HmmCleaner output: `hmmcleaner_out/*.fasta`
(or any directory of alignment FASTA files)

### Output

| File | Description |
|------|-------------|
| `coverage_report.tsv` | Per-gene statistics and PASS/FAIL verdict |
| `passing_genes.txt` | List of gene IDs that passed all thresholds |
| `passing_aln/` | (optional) Copy of passing alignment files |

### Usage

```bash
# Assess protein alignments and copy passing ones
python alignment_coverage.py \
    --input-dir   hmmcleaner_out/ \
    --report      coverage_report.tsv \
    --gene-list   passing_genes.txt \
    --out-dir     passing_aln/ \
    --alphabet    aa

# Assess nucleotide alignments
python alignment_coverage.py \
    --input-dir   nt_aln/ \
    --report      nt_coverage_report.tsv \
    --alphabet    nt

# Parallel processing for large datasets
python alignment_coverage.py \
    --input-dir  hmmcleaner_out/ \
    --report     coverage_report.tsv \
    --threads    8
```

### Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-i` / `--input-dir DIR` | Yes* | — | Directory of alignment files |
| `--fasta FILE` | Yes* | — | Single alignment file (alternative) |
| `-r` / `--report FILE` | Yes | — | Output TSV report path |
| `-g` / `--gene-list FILE` | No | — | Write passing gene IDs here |
| `-o` / `--out-dir DIR` | No | — | Copy passing alignments here |
| `--filter-sequences` | No | off | When copying, remove sequences below `--min-completeness` |
| `-a` / `--alphabet` | No | `auto` | `aa`, `nt`, or `auto` |
| `-e` / `--ext EXT` | No | fa/fasta/faa | File extension filter |
| `--min-completeness` | No | `0.50` | Minimum valid-character fraction per sequence |
| `--min-fraction` | No | `0.80` | Minimum fraction of sequences passing completeness |
| `--min-sequences` | No | `2` | Minimum number of sequences per alignment |
| `--min-aln-length` | No | `30` | Minimum alignment length (aa or nt) |
| `-t` / `--threads` | No | `1` | Parallel worker processes |

*One of `--input-dir` or `--fasta` is required.

### Report columns

| Column | Description |
|--------|-------------|
| `gene` | Gene / orthogroup ID (filename stem) |
| `n_sequences` | Number of sequences |
| `alignment_length` | Alignment length (residues) |
| `mean_completeness` | Mean per-sequence completeness |
| `min_completeness` | Worst-sequence completeness |
| `n_pass_sequences` | Sequences meeting `--min-completeness` |
| `fraction_pass` | `n_pass / n_sequences` |
| `frac_cols_covered` | Fraction of columns with coverage ≥ `--min-fraction` |
| `alphabet` | `aa` or `nt` |
| `PASS` | `True` if all thresholds met |
| `FAIL_REASON` | Human-readable explanation (empty if PASS) |

---

## Step B — `pal2nal.sbatch`

### What it does

PAL2NAL back-translates a protein multiple-sequence alignment to a
codon-level nucleotide alignment using unaligned CDS sequences.
The protein alignment structure (gaps, column positions) is preserved exactly;
only the amino acids are replaced with the corresponding DNA triplets.

Reference: Suyama et al. 2006, Nucleic Acids Res. 34:W609-W612

### Input

| | Directory | File format |
|-|-----------|-------------|
| A | Protein alignments | `*_aln_hmm.fasta` (from coverage step `passing_aln/`) |
| B | Nucleotide CDS | `OG0000001.fa` or `OG0000001.fasta` (unaligned CDS, one per orthogroup) |

The script matches files by orthogroup ID (strips `_aln`, `_hmm`, `_aligned`
suffixes). Nucleotide files must be named `OG0000001.fa`, `OG0000001.fasta`,
`OG0000001_nt.fa`, or `OG0000001_nt.fasta`.

### Output

| File | Description |
|------|-------------|
| `<gene>_nt.fasta` | Codon-level nucleotide alignment |
| `<gene>_pal2nal.err` | PAL2NAL warnings (kept only if non-empty) |

**Next step:** pass `<output_dir>/` to
`../../EnrichAlignment/orthogroup_table.py`.

### Usage

```bash
# 1. Count protein alignments passing coverage filter
N=$(ls -1 passing_aln/*.fasta | wc -l)

# 2. Submit as SLURM array
sbatch --array=0-$((N-1))%20 pal2nal.sbatch \
       --protein-aln-dir passing_aln/ \
       --nt-dir           nt_cds/ \
       --output-dir       nt_aln/

# 3. Or run sequentially
bash pal2nal.sbatch \
     --protein-aln-dir passing_aln/ \
     --nt-dir           nt_cds/ \
     --output-dir       nt_aln/
```

### Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `-a` / `--protein-aln-dir DIR` | Yes | — | Protein alignment directory |
| `-n` / `--nt-dir DIR` | Yes | — | Nucleotide CDS directory |
| `-o` / `--output-dir DIR` | Yes | — | Output directory |
| `--pal2nal-bin PATH` | No | auto-download | Path to existing `pal2nal.pl` |
| `--pal2nal-cache DIR` | No | `~/.local/share/pal2nal` | Download cache |
| `--prot-ext EXT` | No | `fasta` | Protein file extension |
| `--nt-ext EXT` | No | `fa fasta` | NT file extensions to try |
| `--format FMT` | No | `fasta` | Output format: `fasta`, `paml`, `clustal`, `codon` |
| `--nogap` | No | off | Remove gap-containing columns |
| `--codontable N` | No | `1` | Codon table (1=universal, 2=vertebrate mt) |

---

## Cluster configuration (pal2nal.sbatch)

| Line | What to change |
|------|----------------|
| `#SBATCH --partition=cpu-dedicated` | Your cluster partition |
| `#SBATCH --account=myaccount` | Uncomment if required |
| `#SBATCH --time=02:00:00` | Wall-time per task |
| `module load perl/...` | Uncomment if Perl is a module on your cluster |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `[SKIP] — no nucleotide file found` | NT filename doesn't match | Check `nt_cds/` naming; expected `OG0000001.fa` |
| Empty output `.fasta` | Sequence ID mismatch between protein and NT | IDs must be identical in both files |
| Empty output `.fasta` | Wrong codon table | Try `--codontable 2` (mitochondrial) |
| `Download failed` | No internet on compute node | Run once on login node; cache is at `~/.local/share/pal2nal/` |
| Coverage: all genes fail | Thresholds too strict | Lower `--min-completeness` or `--min-fraction`; review distribution |

---

## Navigation

← [clean_alignment/ — HmmCleaner masking](../README.md)
→ [EnrichAlignment/ — combine with polymorphism](../../../EnrichAlignment/README.md)
