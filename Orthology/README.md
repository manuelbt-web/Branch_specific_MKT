# Ortholog Detection Pipeline

This folder contains two complementary approaches for detecting orthologs
between species, used as input to the McDonald-Kreitman test (MKT):

| Method | Script | When to use |
|---|---|---|
| **Reciprocal Best BLAST Hits (RBH)** | `reciprocal_blast.sbatch` | Pairwise (2 species); nucleotide or protein |
| **OrthoFinder** | `orthofinder.sbatch` | Multi-species orthogroup inference (≥ 2 species, protein only) |

---

## Method 1 — Reciprocal Best BLAST Hits (RBH)

Two genes are reciprocal best hits (and thus putative orthologs) if:
- Gene A in species 1 has gene B as its best BLAST match in species 2
- Gene B in species 2 has gene A as its best BLAST match in species 1

This is a fast pairwise approach, suitable for comparing exactly 2 species
at a time.

### How it works (one command)

`reciprocal_blast.sbatch` runs the complete pipeline automatically:

```
Species A FASTA + Species B FASTA
         │
         ▼  [1/3]  BLAST: A → B        →  A_vs_B.tab
         ▼  [2/3]  BLAST: B → A        →  B_vs_A.tab
         ▼  [3/3]  Reciprocal Best Hits →  RBH_A_B.tab
```

### Usage

Edit the `[CLUSTER]` lines in `reciprocal_blast.sbatch` before submission.

```bash
# Nucleotide BLAST (CDS sequences):
sbatch reciprocal_blast.sbatch \
    --mode    nucl \
    --query   speltoides_cds.fasta \
    --subject tauschii_cds.fasta \
    --outdir  results/blast/

# Protein BLAST (proteomes):
sbatch reciprocal_blast.sbatch \
    --mode    prot \
    --query   speltoides_proteome.fasta \
    --subject tauschii_proteome.fasta \
    --outdir  results/blast/
```

### Arguments

| Argument | Required | Description |
|---|---|---|
| `--mode nucl\|prot` | Yes | `blastn` (nucleotide) or `blastp` (protein) |
| `--query FILE` | Yes | FASTA file for species A |
| `--subject FILE` | Yes | FASTA file for species B |
| `--outdir DIR` | Yes | Output directory |
| `--evalue FLOAT` | No | E-value cutoff (default: 1e-20 for nucl, 1e-10 for prot) |
| `--threads N` | No | BLAST threads (default: SLURM_CPUS_PER_TASK or 8) |

### Outputs

| File | Description |
|---|---|
| `<A>_vs_<B>.tab` | BLAST hits: A as query, B as subject (outfmt 6) |
| `<B>_vs_<A>.tab` | BLAST hits: B as query, A as subject (outfmt 6) |
| `RBH_<A>_<B>.tab` | Reciprocal Best Hit pairs |

**BLAST outfmt 6 columns:**

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| qseqid | sseqid | pident | length | mismatch | gapopen | qstart | qend | sstart | send | evalue | bitscore |

**RBH output format (`RBH_<A>_<B>.tab`):**

```
#A_id           B_id            A_vs_B_score    B_vs_A_score
GENE001|HOG123  GENE789|HOG456  1245.0          1187.0
```

### Re-running the RBH step only

If BLAST already ran and you want to re-run just the RBH detection
(e.g. with a different score column or order), call `reciprocal_best_hits.py`
directly:

```bash
python reciprocal_best_hits.py \
    --a-vs-b  results/blast/speltoides_vs_tauschii.tab \
    --b-vs-a  results/blast/tauschii_vs_speltoides.tab \
    --output  results/RBH_speltoides_tauschii.tab \
    --score-col 12 \   # 12 = bitscore (default); use 11 for evalue
    --order high       # 'high' for bitscore, 'low' for evalue
```

### Dependencies

| Tool | Version | Notes |
|---|---|---|
| BLAST+ | ≥ 2.12 | `blastn` and/or `blastp` |
| Python | ≥ 3.7 | No extra packages needed |

---

## Method 2 — OrthoFinder

OrthoFinder is a comprehensive comparative genomics platform that:
- Finds orthogroups (groups of genes descended from a single ancestral gene)
- Infers rooted gene trees for every orthogroup
- Identifies all gene duplication events
- Infers a rooted species tree
- Provides pairwise ortholog tables per species pair

Use OrthoFinder when you have **3 or more species** or want full orthogroup
inference beyond simple pairwise RBH.

### Important — input format

> OrthoFinder requires **one protein FASTA file per species** in a single
> directory.  It does **NOT** accept a single combined file.

```
proteomes/
├── Ae_speltoides.fasta     ← one file per species
├── Ae_tauschii.fasta       ← filename (without extension) = species name
├── Ae_urartu.fasta
└── T_monococcum.fasta
```

Each FASTA must contain **protein sequences** (amino acids) — not CDS or mRNA.
If you only have CDS sequences, translate them first:

```bash
# Translate CDS to protein using seqkit:
seqkit translate cds.fasta > proteome.fasta

# Or with TransDecoder:
TransDecoder.LongOrfs -t transcripts.fasta
TransDecoder.Predict -t transcripts.fasta
```

### Installation

**Option A — conda (recommended):**

```bash
conda create -n orthofinder_env -c bioconda orthofinder diamond
conda activate orthofinder_env
orthofinder --help
```

**Option B — manual (OrthoFinder 3.0.1b1):**

```bash
wget https://github.com/davidemms/OrthoFinder/releases/download/3.0.1b1/OrthoFinder.tar.gz
tar xzf OrthoFinder.tar.gz
cd OrthoFinder_source
pip install scipy biopython
# Test:
python orthofinder.py --help
```

Dependencies installed by conda automatically: DIAMOND, MCL, FastME, MAFFT.
For the manual install, ensure these are in your `$PATH`.

### Usage

Edit the `[CLUSTER]` lines in `orthofinder.sbatch` before submission.

```bash
sbatch orthofinder.sbatch \
    --proteomes proteomes/ \
    --outdir    results/orthofinder/
```

### Arguments

| Argument | Default | Description |
|---|---|---|
| `--proteomes DIR` | required | Directory with one FASTA per species |
| `--outdir DIR` | required | Output directory for all OrthoFinder results |
| `--search TOOL` | `diamond` | Search tool: `diamond`, `diamond_ultra_sens`, `blast`, `mmseqs2_default` |
| `--threads N` | `SLURM_CPUS_PER_TASK` | Threads for search and analysis |

### Key outputs

| File/Directory | Description |
|---|---|
| `Orthogroups/Orthogroups.tsv` | All orthogroups with member genes per species |
| `Orthogroups/Orthogroups_SingleCopyOrthologues.txt` | Single-copy orthogroups (for phylogenetics) |
| `Orthologues/<sp1>__v__<sp2>/` | Pairwise ortholog tables for each species pair |
| `Species_Tree/SpeciesTree_rooted.txt` | Rooted species tree (Newick) |
| `Gene_Trees/` | One rooted gene tree per orthogroup |
| `Comparative_Genomics_Statistics/` | Summary statistics |

After running, pairwise ortholog tables are ready-made per species pair:

```bash
# List available pairwise tables for Ae_speltoides:
ls results/orthofinder/Orthologues/Orthologues_Ae_speltoides/
# → Ae_speltoides__v__Ae_tauschii.tsv
#   Ae_speltoides__v__Ae_urartu.tsv

# Columns: OrthoGroup | Ae_speltoides_genes | Ae_tauschii_genes
head results/orthofinder/Orthologues/Orthologues_Ae_speltoides/Ae_speltoides__v__Ae_tauschii.tsv
```

### Dependencies

| Tool | Version | Notes |
|---|---|---|
| OrthoFinder | 3.0.1b1 | Install via conda (see above) |
| DIAMOND | ≥ 2.0 | Included in conda install |
| Python | ≥ 3.7 | With `scipy` and `biopython` |
| MCL | any | Included in conda install |

---

## Choosing between RBH and OrthoFinder

| | RBH | OrthoFinder |
|---|---|---|
| Number of species | 2 only | 2 or more |
| Input sequences | Nucleotide or protein | Protein only |
| Speed | Fast | Slower (depends on dataset) |
| Output | Pairwise best-hit table | Full orthogroup + gene trees |
| Gene trees | No | Yes |
| Species tree | No | Yes |
| Suitable for MKT | Yes | Yes (use pairwise Orthologues/ tables) |

For the Branch-specific MKT in this project, **both methods are valid**:
- RBH gives a simple, reproducible pairwise list for 2-species comparisons.
- OrthoFinder gives richer context (gene trees, duplication events) when
  analysing 3+ species simultaneously.

---

## Folder structure

```
Blast/
├── reciprocal_blast.sbatch     Complete RBH pipeline: BLAST + RBH in one command
├── reciprocal_best_hits.py     RBH helper script (called internally; also usable standalone)
├── orthofinder.sbatch          OrthoFinder multi-species orthogroup inference
└── README.md                   This file
```
