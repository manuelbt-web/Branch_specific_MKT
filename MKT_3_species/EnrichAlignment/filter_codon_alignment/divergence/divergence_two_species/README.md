# divergence_two_species — Standard 2-species MKT preparation

This directory contains scripts for the standard (2-species) McDonald-Kreitman
Test pipeline: reducing the enriched alignment to a single outgroup, computing
dN/dS + pN/pS statistics, and merging all results into one analysis table.

---

## Directory contents

| Script | Description |
|--------|-------------|
| `trim_outgroup.py` | Remove focal-species reference and second outgroup; keep 1 outgroup for standard MKT |
| `dNdSpiNpiS.sbatch` | SLURM array: run dNdSpiNpiS v1.0 per gene to get fixN, fixS, dN, dS, pN, pS |
| `merge_mkt_results.py` | Merge polymorphism stats + dNdSpiNpiS counts + codeml results into one table |

---

## Pipeline

```
filter_codon_alignment.py output
  filtered_aln/*_NT_filtered.fasta    (poly seqs + 3 div seqs: focal ref + OG1 + OG2)
        │
        ▼  Step 1 — trim_outgroup.py
  alignment_for_standard_MKT/         (poly seqs + OG1 only)
        │
        ├──► Step 2 — dNdSpiNpiS.sbatch
        │    dNdSpiNpiS_output/*.out   (fixN, fixS, dN, dS, pN, pS, DoS, NI per gene)
        │
        └──► [polymorphism_stats.py is run separately on the polymorphism sequences]
             poly_stats/*.tsv           (Pi, piS, piNS, PnMinus, PnGreater, … per gene)

  Step 3 — merge_mkt_results.py
    + focal_species_results.tsv        (from parse_results.py codeml branch model)
        │
        ▼
  merged_mkt_results.tsv              ← final analysis table
```

---

## Step 1 — `trim_outgroup.py`

Removes the first divergence sequence (focal-species reference) and the last
divergence sequence (second outgroup), keeping only the middle outgroup for
standard MKT.

```bash
python trim_outgroup.py \
    --input-dir          filtered_aln/ \
    --output-dir         mkt_ready/ \
    --divergence-pattern "Aegilops_speltoides|Aegilops_mutica|Aegilops_tauschii"
```

Output: `mkt_ready/alignment_for_standard_MKT/*.fasta`

See [trim_outgroup.py](trim_outgroup.py) options for `--n-divergence`, `--threads`, etc.

---

## Step 2 — `dNdSpiNpiS.sbatch`

Run on the `alignment_for_standard_MKT/` directory.

```bash
# Download binary from:
# https://kimura.univ-montp2.fr/PopPhyl/index.php?section=tools

N=$(ls mkt_ready/alignment_for_standard_MKT/*.fasta | wc -l)
sbatch --array=0-$((N-1))%15 dNdSpiNpiS.sbatch \
    -i mkt_ready/alignment_for_standard_MKT/ \
    -o dNdSpiNpiS_output/ \
    -b /path/to/dNdSpiNpiS_1.0 \
    -g Ae_speltoides \
    -k Ae_mutica
```

Output: `dNdSpiNpiS_output/<gene>.out` with columns including `fixN`, `fixS`.

---

## Step 3 — `merge_mkt_results.py`

Merges three result sources into one table, joined by `gene_id`.

```bash
python merge_mkt_results.py \
    --polymorphism-dir  poly_stats/ \
    --divergence-dir    dNdSpiNpiS_output/ \
    --codeml-file       codeml_results/focal_species_results.tsv \
    --output            merged_mkt_results.tsv \
    --focal-species     Aegilopsspeltoides
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-p / --polymorphism-dir DIR` | required | Directory of polymorphism TSV files |
| `-d / --divergence-dir DIR` | required | Directory of dNdSpiNpiS `.out` files |
| `-c / --codeml-file FILE` | required | `focal_species_results.tsv` from parse_results.py |
| `-o / --output FILE` | required | Output TSV file path |
| `--focal-species NAME` | all species | Filter codeml file to one species |
| `--strip-suffix SUFFIX` | `_NT_filtered` | Strip this suffix from filenames for gene_id matching |
| `--verbose` | off | Print merged column list |

### Join key

All three sources are joined on `gene_id`:
- Polymorphism TSVs: `gene_id` column
- dNdSpiNpiS: filename stem (e.g. `EVM0000002.1_NT_filtered.out` → `EVM0000002.1` after stripping `_NT_filtered`)
- codeml: `id_gene` column (renamed to `gene_id` before joining)

If `gene_id` values don't match (< 80% match rate), check `--strip-suffix`.

---

## Output columns — `merged_mkt_results.tsv`

### From polymorphism_stats.py

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier (join key) |
| `number_sites` | Effective nucleotide sites (lseff) |
| `num_codons_eff` | Effective codon count |
| `Pi_per_site` | Nucleotide diversity π per site |
| `pi_S` | Synonymous π per site |
| `pi_NS` | Non-synonymous π per site |
| `thetaW_per_site` | Watterson's θ per site |
| `Tajimas_D` | Tajima's D |
| `num_sites_S_egglib` | Effective synonymous sites |
| `num_sites_NS_egglib` | Effective non-synonymous sites |
| `polymorphism_count` | Total polymorphic sites |
| `num_pol_S` | Synonymous polymorphic sites |
| `num_pol_NS` | Non-synonymous polymorphic sites |
| `sample_size` | Number of alleles analysed |
| `PnMinus` | P_N^(<T) — low-frequency NS polymorphisms |
| `PnGreater` | P_N^(>T) — high-frequency NS polymorphisms |
| `PsMinus` | P_S^(<T) — low-frequency S polymorphisms |
| `PsGreater` | P_S^(>T) — high-frequency S polymorphisms |
| `ratioPs` | Neutral expectation P_N^(>T) / P_S^(>T) |
| `PnNeutral_estimation` | Neutral Pn (float) |
| `PnNeutral` | Neutral Pn (integer, for impMKT) |
| `deleterious` | Estimated SDM count |
| `downsampled` | `yes` if haploid downsampling was applied |
| `nt_engine` | `egglib` or `builtin` |
| `S` | Segregating nucleotide sites |

### From dNdSpiNpiS

| Column | Description |
|--------|-------------|
| `dn_counts` | Fixed non-synonymous differences (fixN) |
| `ds_counts` | Fixed synonymous differences (fixS) |
| `pN` | Polymorphic non-synonymous rate |
| `pS` | Polymorphic synonymous rate |
| `dN` | Divergence non-synonymous rate |
| `dS` | Divergence synonymous rate |
| `pN/pS` | Polymorphism ratio |
| `dN/dS` | Divergence ratio |
| `DoS` | Direction of Selection |
| `NI` | Neutrality Index |

### From parse_results.py (codeml)

| Column | Description |
|--------|-------------|
| `ortholog` | HOG identifier |
| `N_codeml` | Expected non-synonymous sites (PAML) |
| `S_codeml` | Expected synonymous sites (PAML) |
| `dN` | Branch dN (focal species) |
| `dS` | Branch dS (focal species) |
| `omega` | Branch dN/dS (focal species) |
| `lnL_branch_model` | Log-likelihood of branch model |
| `lnL_m0_model` | Log-likelihood of M0 model |
| `omega_m0` | Global dN/dS from M0 |
| `omega_diff` | omega − omega_m0 |
| `species` | Focal species label |

---

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `< 80% genes matched dNdSpiNpiS` | Filename suffix not stripped | Adjust `--strip-suffix` |
| `focal-species not found` | Wrong cleaned species name | Check `species` column in `focal_species_results.tsv` |
| `fixN not found in dNdSpiNpiS output` | Different column naming in tool version | Check .out file header; the tool may use `Dn` instead of `fixN` |
| `No .tsv files found in polymorphism directory` | Wrong directory | Pass the directory with per-gene TSV files, not a single summary file |

---

← [divergence/ — codeml 3-species pipeline](../README.md)
← [filter_codon_alignment/ — codon filter](../../README.md)
