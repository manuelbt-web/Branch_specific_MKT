# Step 6 — Polymorphism Statistics and SDM Estimation (EggLib)

## Objective

Compute per-gene polymorphism statistics and estimate strongly deleterious
mutations (SDMs) from codon-aligned FASTA files.  The user runs **two scripts**:

1. **`remove_outgroup_sequence.py`** (Step 6a) — strip the outgroup sequence so
   that polymorphism statistics are computed on focal individuals only.
2. **`polymorphisms_stats.py`** (Step 6b) — compute nucleotide and codon-level
   polymorphism statistics using EggLib 3.5.2 `CodingDiversity`, build the
   per-gene folded site frequency spectrum (SFS), and automatically estimate
   strongly deleterious mutations (SDMs) using the impMKT frequency-threshold
   approach (Fay et al. 2001; Messer & Petrov 2013).

SDM estimation is performed automatically at the end of `polymorphisms_stats.py`.
A standalone `estimation_of_SDM.py` is also provided for re-running with a
different frequency cutoff without reprocessing the alignments.

---

## Scripts

| Script | Role |
|---|---|
| `remove_outgroup_sequence.py` | Step 6a — remove outgroup by header tag |
| `polymorphisms_stats.py` | Step 6b — compute statistics + SFS + SDM estimation |
| `estimation_of_SDM.py` | Standalone — re-estimate SDMs from the saved SFS with a different cutoff |

---

## Dependencies

| Tool | Version | Installation |
|---|---|---|
| **Python** | ≥ 3.9 | system / conda |
| **EggLib** | 3.5.2 | Pixi (see below) — **required** |
| **pandas** | any | `pip install pandas` (required for standalone `estimation_of_SDM.py`) |

### Install EggLib 3.5.2 with Pixi

`egglib` on bioconda only ships `linux-64` builds:

```bash
# 1. Install Pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | bash

# 2. Create a Pixi environment with EggLib
pixi init egglib_env
cd egglib_env
# Add bioconda channel (required for EggLib)
sed -i 's/channels = \["conda-forge"\]/channels = ["conda-forge", "bioconda"]/' pixi.toml
pixi add "python=3.11.*" "egglib==3.5.2"

# 3. Run the polymorphism script through this environment
pixi run python polymorphisms_stats.py --input-dir ... --output ...
```

On native Windows, `pixi add egglib` cannot resolve (no `win-64` build exists
on bioconda) — install WSL (`wsl --install`), then set up the Pixi environment
above inside the WSL distribution. The Windows filesystem is reachable from
WSL under `/mnt/c/...`, so `--input-dir`/`--output` can point at paths on the
Windows drive directly.

**EggLib is a hard requirement, not an optional accelerator.** Its
`CodingDiversity` codon classification systematically differs from the
built-in pure-Python fallback — validated on a real dataset (Ae. speltoides,
3,560 genes), EggLib matched a hand-curated reference table on 94% of genes
vs. only 61% for the fallback, which was enough to change branch-specific MKT
candidate counts (46 vs. the correct 73). `polymorphisms_stats.py` therefore
**exits with an error if EggLib is not importable**, unless you explicitly
pass `--allow-builtin-fallback` (not recommended — only for quick,
non-publication smoke tests). The SFS produced by the pure-Python fallback is
also codon-based rather than site-based, though the file format is identical.

---

## Input files

### Step 6a — `remove_outgroup_sequence.py`

| Input | Description | Source |
|---|---|---|
| `--input-dir DIR` | Per-gene aligned FASTA files with outgroup | `step4_out/` from Step 4 |

### Step 6b — `polymorphisms_stats.py`

| Input | Description | Source |
|---|---|---|
| `--input-dir DIR` | Per-gene FASTA files without outgroup | `step6a_no_outgroup/` from Step 6a |

---

## Output files

### Step 6a — `remove_outgroup_sequence.py`

| File | Description |
|---|---|
| `step6a_no_outgroup/<gene>.fasta` | FASTA with only focal sequences |

### Step 6b — `polymorphisms_stats.py`

| File | Description |
|---|---|
| `step6_stats/polymorphism_stats.tsv` | Main output: one row per gene (see columns below) |
| `step6_stats/folded_sfs.tsv` | Folded SFS (optional; required only if re-running `estimation_of_SDM.py`) |

#### Columns in `polymorphism_stats.tsv`

| Column | Description |
|---|---|
| `gene_id` | Gene identifier (FASTA filename stem) |
| `number_sites` | Effective number of nucleotide sites analysed (lseff) |
| `num_codons_eff` | Effective number of codons (nS + nN) / 3 from EggLib, else lseff / 3 |
| `Pi_per_site` | Nucleotide diversity π per site |
| `pi_S` | Synonymous nucleotide diversity per site (EggLib CodingDiversity Pi_S) |
| `pi_NS` | Non-synonymous nucleotide diversity per site (EggLib Pi_NS) |
| `thetaW_per_site` | Watterson's θ per site |
| `Tajimas_D` | Tajima's D neutrality test statistic |
| `num_sites_S_egglib` | Effective synonymous sites (EggLib nseff_S; NA if not exposed) |
| `num_sites_NS_egglib` | Effective non-synonymous sites (EggLib nseff; NA if not exposed) |
| `polymorphism_count` | Total polymorphic sites (num_pol_S + num_pol_NS) |
| `num_pol_S` | Synonymous polymorphic sites |
| `num_pol_NS` | Non-synonymous polymorphic sites |
| `sample_size` | Number of alleles analysed (n; haploid-adjusted when `--haploid`) |
| `PnMinus` | P_N^(<T) — NS polymorphisms at frequency ≤ T |
| `PnGreater` | P_N^(>T) — NS polymorphisms at frequency > T |
| `PsMinus` | P_S^(<T) — S polymorphisms at frequency ≤ T |
| `PsGreater` | P_S^(>T) — S polymorphisms at frequency > T |
| `ratioPs` | Neutral expectation P_N^(>T) / P_S^(>T) |
| `PnNeutral_estimation` | Neutral Pn (float): P_N − deleterious, capped to [0, P_N] |
| `PnNeutral` | Neutral Pn rounded to nearest integer (used in impMKT) |
| `deleterious` | Estimated SDM count = P_N^(<T) − ratioPs × P_S^(<T) |
| `downsampled` | `yes` if haploid downsampling was applied (`--haploid`) |
| `nt_engine` | `egglib` or `builtin` (which engine computed the statistics) |
| `S` | Number of segregating (polymorphic) nucleotide sites |

#### Columns in `folded_sfs.tsv`

| Column | Description |
|---|---|
| `gene` | Gene identifier |
| `sample_size` | Number of alleles analysed (n, haploid-aware) |
| `bin` | Minor allele count (0 = monomorphic, 1 … n/2 = polymorphic) |
| `frequency` | Folded allele frequency = bin / sample_size |
| `syn_count` | Number of synonymous sites at this frequency class |
| `nonsyn_count` | Number of non-synonymous sites at this frequency class |

---

## SDM formula

```
P_SDM ≈ P_N^(<T) − (P_N^(>T) / P_S^(>T)) × P_S^(<T)
```

The neutral expectation **P_N^(>T) / P_S^(>T)** is estimated from the
high-frequency fraction of the SFS (above threshold T = 0.15 by default).
Low-frequency NS variants in excess of this expectation are attributed to
slightly deleterious mutations under purifying selection.

For selfing species (Ae. tauschii, T. urartu), all statistics are computed on
n/2 alleles by randomly drawing one haploid sequence per individual (`--haploid`).

---

## Usage

### Step 6a — Remove the outgroup

```bash
python step6_polymorphism/remove_outgroup_sequence.py \
    --input-dir  step4_out/ \
    --output-dir step6a_no_outgroup/
```

| Option | Required | Default | Description |
|---|---|---|---|
| `--input-dir DIR` | Yes | — | Aligned FASTAs with outgroup (Step 4 output) |
| `--output-dir DIR` | Yes | — | Output directory for focal-only FASTAs |
| `--outgroup-tag TAG` | No | `\|outgroup` | String in FASTA header identifying the outgroup |
| `--min-sequences N` | No | `2` | Minimum focal sequences required |

### Step 6b — Compute statistics + SDM estimation

```bash
python step6_polymorphism/polymorphisms_stats.py \
    --input-dir   step6a_no_outgroup/ \
    --output      step6_stats/polymorphism_stats.tsv \
    --sfs-output  step6_stats/folded_sfs.tsv \
    --freq-cutoff 0.15 \
    --max-missing 0.30 \
    --threads     4
```

For selfing species (Ae. tauschii, T. urartu), add `--haploid`:

```bash
python step6_polymorphism/polymorphisms_stats.py \
    --input-dir   step6a_no_outgroup/ \
    --output      step6_stats/polymorphism_stats.tsv \
    --sfs-output  step6_stats/folded_sfs.tsv \
    --freq-cutoff 0.15 \
    --max-missing 0.30 \
    --haploid \
    --seed        42 \
    --threads     4
```

| Option | Required | Default | Description |
|---|---|---|---|
| `--input-dir DIR` | Yes | — | Focal-only FASTAs (Step 6a output) |
| `--output FILE` | Yes | — | Output TSV with all statistics + SDM columns |
| `--sfs-output FILE` | No | none | Optional: also write the intermediate folded SFS |
| `--max-missing FLOAT` | No | `0.30` | Max fraction of missing/gap bases per site |
| `--freq-cutoff FLOAT` | No | `0.15` | Frequency threshold T for SDM estimation |
| `--haploid` | No | off | Draw one sequence per individual (selfing species) |
| `--seed INT` | No | `42` | Random seed for haploid downsampling |
| `--threads N` | No | `1` | Number of parallel worker processes |

### (Optional) Re-run SDM estimation with a different cutoff

If you already have `folded_sfs.tsv` and want to test a different frequency
threshold without reprocessing all alignments:

```bash
python step6_polymorphism/estimation_of_SDM.py \
    --sfs     step6_stats/folded_sfs.tsv \
    --cutoff  0.10 \
    --output  step6_stats/sdm_cutoff010.tsv
```

| Option | Required | Default | Description |
|---|---|---|---|
| `--sfs FILE` | Yes | — | Folded SFS TSV (`--sfs-output` from Step 6b) |
| `--cutoff FLOAT` | No | `0.15` | Frequency threshold T |
| `--output FILE` | Yes | — | Output TSV with SDM estimates |

---

## Full two-step example

```bash
# Step 6a — strip outgroup
python step6_polymorphism/remove_outgroup_sequence.py \
    --input-dir  step4_out/ \
    --output-dir step6a_no_outgroup/

# Step 6b — compute statistics, SFS, and SDM estimates in one run
python step6_polymorphism/polymorphisms_stats.py \
    --input-dir   step6a_no_outgroup/ \
    --output      step6_stats/polymorphism_stats.tsv \
    --sfs-output  step6_stats/folded_sfs.tsv \
    --freq-cutoff 0.15 \
    --threads     4
```

For selfing species, add `--haploid` to Step 6b only; Step 6a is unchanged.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `No .fasta files found` | Wrong `--input-dir` path | Check `ls INPUT_DIR/*.fasta` |
| All genes skipped | Wrong `--outgroup-tag` | Verify tag matches FASTA headers from Step 4 (default: `\|outgroup`) |
| `ERROR: EggLib is not installed` | EggLib not importable in this environment | Activate the Pixi environment (`pixi shell` / `pixi run python ...`) — see Dependencies above. Windows: run it from WSL. |
| `nt_engine: builtin` in all rows | Ran with `--allow-builtin-fallback` | Only for smoke tests — re-run without this flag through the EggLib environment for real results |
| `NaN` for Tajima's D | Fewer than 2 segregating sites | Expected for monomorphic genes; not an error |
| `NA` for `num_sites_S_egglib` | EggLib 3.5.2 does not expose `nseff_S` attribute | Check EggLib version; the pipeline remains functional |
| `PnNeutral = NA` for many genes | No high-frequency S polymorphisms | Use a lower `--freq-cutoff` or check alignment quality |
| Very slow | Large number of genes | Increase `--threads` |

---

## Previous step

← [Step 4: Prepare headers for dNdSpNpS](../step4_prep_dNdS/README.md)
← [Step 5: Divergence statistics (dNdSpiNpiS)](../step5_run_dNdSpNpS/README.md)

Steps 5 and 6 are independent — both start from `step4_out/` and can be run
in any order or simultaneously.
