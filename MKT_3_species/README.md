# MKT 3-Species — Alignment, Cleaning and Enrichment Pipeline

This directory contains the full pipeline for preparing codon-level alignments
that combine divergence (multi-species orthologs) and polymorphism
(intraspecific CDS sequences) for McDonald-Kreitman Test (MKT) or
dN/dS (codeml) analyses.

---

## Directory structure

```
MKT_3_species/
├── Alignment_divergence/              — Step 1–2: align and clean orthologs
│   ├── mafft_alignments.sbatch        SLURM array — MAFFT protein alignment
│   └── clean_alignment/               — Step 2: quality filtering
│       ├── hmm_cleaner.sbatch         SLURM array — HmmCleaner masking
│       └── translate_back_in_nt/      — Step 3: coverage + back-translation
│           ├── alignment_coverage.py  assess completeness; filter alignments
│           └── pal2nal.sbatch         SLURM array — protein → codon NT alignment
│
├── Alignment_polymorphism/            — Step 1b (parallel track): align CDS
│   └── MACSE_aligning_sequences.sbatch  SLURM array — MACSE alignSequences
│
└── EnrichAlignment/                   — Step 4–5: enrich + filter
    ├── orthogroup_table.py            build HOG → gene-ID mapping table
    ├── MACSE_enrichment.py            add polymorphism sequences to divergence aln
    ├── macse_enrichment.sbatch        SLURM wrapper for MACSE_enrichment.py
    └── filter_codon_alignment/        — Step 5: filter codons for codeml / MKT
        └── filter_codon_alignment.py  per-codon quality filter; outputs NNN-masked aln
```

---

## Pipeline overview

```
OrthoFinder results
└── Single_Copy_Orthologue_Sequences/    (one protein FASTA per orthogroup)

  ─── DIVERGENCE TRACK ────────────────────────────────────────────────────
  Alignment_divergence/mafft_alignments.sbatch          [SLURM array]
  └── mafft_out/*.fa              (MAFFT --auto protein alignments)

  Alignment_divergence/clean_alignment/hmm_cleaner.sbatch   [SLURM array]
  └── hmmcleaner_out/*.fasta      (HmmCleaner-masked protein alignments)

  translate_back_in_nt/alignment_coverage.py            [Python, login node]
  └── passing_aln/*.fasta         (alignments passing coverage thresholds)
      coverage_report.tsv         (per-gene statistics)

  translate_back_in_nt/pal2nal.sbatch                   [SLURM array]
  └── nt_aln/*.fasta              (codon-level nucleotide alignments)

  ─── POLYMORPHISM TRACK ──────────────────────────────────────────────────
  Alignment_polymorphism/MACSE_aligning_sequences.sbatch [SLURM array]
  └── cds_aligned/*.fasta         (per-gene MACSE CDS alignments)

  ─── ENRICHMENT + FILTERING ──────────────────────────────────────────────
  EnrichAlignment/orthogroup_table.py                    [Python, login node]
  └── orthogroup_table.tsv        (HOG → species gene-ID mapping)

  EnrichAlignment/macse_enrichment.sbatch                [SLURM single job]
  └── macse_enriched/enriched/*.fasta  (combined divergence + polymorphism)

  EnrichAlignment/filter_codon_alignment/filter_codon_alignment.py  [Python]
  └── filtered_aln/*.fasta        (NNN-masked alignments ready for codeml/MKT)
```

---

## Dependencies at a glance

| Tool | Version | Folder | Installation |
|------|---------|--------|--------------|
| **MAFFT** | ≥ 7.4 | `Alignment_divergence/` | `conda install -c bioconda mafft` |
| **HMMER** | ≥ 3.3 | `clean_alignment/` | `conda install -c bioconda hmmer` |
| **HmmCleaner** | ≥ 1.8 | `clean_alignment/` | `conda install -c bioconda hmmcleaner` |
| **Perl** | ≥ 5.20 | `translate_back_in_nt/` | system |
| **PAL2NAL** | newest | `translate_back_in_nt/` | **auto-downloaded** from GitHub |
| **Python** | ≥ 3.8 | multiple | system / conda |
| **pandas** | ≥ 1.3 | `EnrichAlignment/` | `conda install pandas` |
| **MACSE** | ≥ 2.07 | `Alignment_polymorphism/`, `EnrichAlignment/` | [agap-ge2pop.org/macse](https://www.agap-ge2pop.org/macse/) |
| **Java** | ≥ 11 | `Alignment_polymorphism/`, `EnrichAlignment/` | `conda install -c conda-forge openjdk` |

---

## Subfolder documentation

Detailed usage, options, and troubleshooting for every script:

| Folder | README | Contents |
|--------|--------|----------|
| `Alignment_divergence/` | [README](Alignment_divergence/README.md) | MAFFT alignment |
| `Alignment_divergence/clean_alignment/` | [README](Alignment_divergence/clean_alignment/README.md) | HmmCleaner masking |
| `Alignment_divergence/clean_alignment/translate_back_in_nt/` | [README](Alignment_divergence/clean_alignment/translate_back_in_nt/README.md) | Coverage filter + PAL2NAL |
| `Alignment_polymorphism/` | [README](Alignment_polymorphism/README.md) | MACSE CDS alignment |
| `EnrichAlignment/` | [README](EnrichAlignment/README.md) | Orthogroup table + MACSE enrichAlignment |
| `EnrichAlignment/filter_codon_alignment/` | [README](EnrichAlignment/filter_codon_alignment/README.md) | Per-codon quality filter |

---

## Quick start (full pipeline)

```bash
# 0. Paths — adapt to your dataset
SCO=results/orthofinder/OrthoFinder/Results_DATE/Single_Copy_Orthologue_Sequences
HMMC=$(which HmmCleaner.pl)           # after conda activate hmmcleaner_env
MACSE_JAR=/path/to/macse_v2.07.jar

# 1. MAFFT protein alignment
N=$(ls "$SCO"/*.fa | wc -l)
sbatch --array=0-$((N-1))%15 Alignment_divergence/mafft_alignments.sbatch \
       --input-dir  "$SCO" --output-dir results/mafft_out/

# 2. HmmCleaner masking
N=$(ls results/mafft_out/*_aln.fa | wc -l)
sbatch --array=0-$((N-1))%15 \
       Alignment_divergence/clean_alignment/hmm_cleaner.sbatch \
       --input-dir      results/mafft_out/ \
       --output-dir     results/hmmcleaner_out/ \
       --hmmcleaner-bin "$HMMC"

# 3a. Coverage assessment
python Alignment_divergence/clean_alignment/translate_back_in_nt/alignment_coverage.py \
    --input-dir  results/hmmcleaner_out/ \
    --report     results/coverage_report.tsv \
    --gene-list  results/passing_genes.txt \
    --out-dir    results/passing_aln/ \
    --alphabet   aa

# 3b. PAL2NAL (protein → codon NT alignment)
N=$(ls results/passing_aln/*.fasta | wc -l)
sbatch --array=0-$((N-1))%20 \
       Alignment_divergence/clean_alignment/translate_back_in_nt/pal2nal.sbatch \
       --protein-aln-dir results/passing_aln/ \
       --nt-dir           nt_cds/ \
       --output-dir       results/nt_aln/

# [PARALLEL] MACSE CDS alignment for polymorphism track
N=$(ls cds_trimmed/*.fasta | wc -l)
sbatch --array=0-$((N-1))%15 \
       Alignment_polymorphism/MACSE_aligning_sequences.sbatch \
       -i cds_trimmed/ -o results/cds_aligned/ -j "$MACSE_JAR"

# 4. Build orthogroup table
python EnrichAlignment/orthogroup_table.py \
    --input-dir results/nt_aln/ \
    --output    results/orthogroup_table.tsv

# 5. MACSE enrichAlignment (divergence + polymorphism)
python EnrichAlignment/MACSE_enrichment.py \
    --ortho-table   results/orthogroup_table.tsv \
    --ortho-dir     results/nt_aln/ \
    --poly-dir      results/cds_aligned/ \
    --output-dir    results/macse_enriched/ \
    --macse-jar     "$MACSE_JAR" \
    --species-col   Aegilops_speltoides \
    --threads       8

# 6. Per-codon filter (replace failing codons with NNN for codeml)
python EnrichAlignment/filter_codon_alignment/filter_codon_alignment.py \
    --input-dir           results/macse_enriched/enriched/ \
    --output-dir          results/filtered_aln/ \
    --divergence-pattern  "Aegilops_tauschii|Triticum_urartu|Amblyopyrum_muticum" \
    --min-polym           6 \
    --threads             8
```

---

## Previous step

← [Orthology: OrthoFinder](../Orthology/README.md)

## Next step

→ codeml (PAML) dN/dS estimation using `filtered_aln/*.fasta`
