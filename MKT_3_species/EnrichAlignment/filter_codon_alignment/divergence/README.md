# Divergence — codeml (PAML) dN/dS Pipeline

This directory contains the complete pipeline for estimating dN/dS ratios
from the divergence codon alignments using PAML's codeml program via VESPA.

**Input:** `split_out/divergence/*.fasta` (from `../split_sequences_polymorphism_divergence.py`)
**Output:** per-orthogroup TSV tables of dN, dS, omega (dN/dS) per branch

---

## Directory structure

```
divergence/
├── Dockerfile.vespa                  Build VESPA + Python 2.7 container
├── codeml_inputs/
│   ├── gene_tree.sh                  Step 1: VESPA infer_genetree
│   ├── branch_table.sh               Step 2: Prepare codeml_input/ directories
│   ├── codeml_set_up.sh              Step 3: VESPA codeml_setup
│   └── correction_codeml.sh          Step 4: Fix directories; generate taskfarm
└── run_codeml/
    ├── run_codeml.sbatch             Step 5: SLURM array — run codeml
    └── parse_results/
        └── parse_results.py          Step 6: Parse codeml output → TSV
```

---

## Pipeline overview

```
split_out/divergence/*.fasta      (codon alignments, divergence sequences only)
        │
        ▼  Step 1 — gene_tree.sh
Inferred_Genetree_<HOG>/
└── <gene>/
    ├── <gene>.fasta
    └── <gene>.tre                (gene tree, topology mapped from species tree)
        │
        ▼  Step 2 — branch_table.sh
└── <gene>/
    ├── branch_table.txt          (foreground species for branch model)
    └── codeml_input/
        ├── cleaned.fasta         (cleaned headers: no '|', no '_')
        └── cleaned.tre           (cleaned tree labels)
        │
        ▼  Step 3 — codeml_set_up.sh (VESPA codeml_setup)
└── <gene>/
    └── Codeml_Setup_codeml_input/
        ├── cleaned/              (site model workspaces: m0, m1Neutral, ...)
        └── cleaned_<sp>/         (branch model workspace: modelA/)
        │
        ▼  Step 4 — correction_codeml.sh
        ├── cleaned/<model>/Omega0_5/codeml.ctl    (fixed codeml.ctl)
        └── cleaned_<sp>/model_branch/Omega0_5/codeml.ctl
        codeml_taskfarm_fullpath.sh                  (branch model, one cd;codeml per gene)
        codeml_taskfarm_fullpath_site_model.sh       (site models, one cd;codeml per gene/model)
        │
        ▼  Step 5 — run_codeml.sbatch (SLURM)
        cleaned/<model>/Omega0_5/out    (codeml output files)
        │
        ▼  Step 6 — parse_results.py
codeml_results/
├── orthogroup_HOG*_cleaned_m0_codeml.tsv
└── orthogroup_HOG*_cleaned_<sp>_branch_codeml.tsv
```

---

## VESPA and Python 2.7

VESPA requires **Python 2.7**, which is end-of-life. The recommended approach
is to run VESPA inside a container. Build the container once, then use it for
steps 1–3.

### Build the Docker image

```bash
docker build -f Dockerfile.vespa -t vespa:latest .
```

### Convert to Singularity (HPC without Docker)

```bash
# On a machine with Docker:
singularity build vespa.sif docker-daemon://vespa:latest

# Or pull directly if image is pushed to a registry:
singularity build vespa.sif docker://your-registry/vespa:latest
```

### Test the image

```bash
docker run --rm vespa:latest --help
# or
singularity exec vespa.sif vespa --help
```

---

## Species tree format

The species tree must be a **pure Newick topology** — species names only,
no branch lengths, no bootstrap values.

```
((Aegilops_tauschii,Triticum_urartu),Amblyopyrum_muticum);
```

Generate from OrthoFinder output by stripping numeric annotations:

```bash
sed 's/:[0-9.eE+-]*//g; s/[0-9]//g' \
    OrthoFinder/Results_*/Species_Tree/SpeciesTree_rooted.txt \
    > SpeciesTree_topology.txt
```

Verify: the file should contain only species names and Newick punctuation `(,);`.

---

## Dependencies

| Tool | Version | Used in | Installation |
|------|---------|---------|--------------|
| **VESPA** | latest | Steps 1–3 | See Dockerfile.vespa |
| **Python 2.7** | 2.7.18 | Steps 1–3 (via container) | See Dockerfile.vespa |
| **dendropy** | 4.4.0 | Steps 1–3 (in container) | See Dockerfile.vespa |
| **PAML / codeml** | ≥ 4.9 | Step 5 | `conda install -c bioconda paml` |
| **Python 3** | ≥ 3.8 | Step 6 | system |
| **pandas** | ≥ 1.3 | Step 6 | `pip install pandas` |
| **biopython** | ≥ 1.79 | Step 6 | `pip install biopython` |

---

## Step-by-step usage

### Step 1 — Infer gene trees

```bash
bash codeml_inputs/gene_tree.sh \
    --input-dir    split_out/divergence/ \
    --species-tree SpeciesTree_topology.txt \
    --docker-image vespa:latest

# or Singularity:
bash codeml_inputs/gene_tree.sh \
    --input-dir           split_out/divergence/ \
    --species-tree        SpeciesTree_topology.txt \
    --singularity-image   /path/to/vespa.sif
```

### Step 2 — Prepare codeml input directories

```bash
bash codeml_inputs/branch_table.sh \
    --root-dir       split_out/divergence/ \
    --target-species "Aegilopsspeltoides,Aegilopsmutica"
```

The `--target-species` names must be **cleaned** (no underscores):
- Original header: `>Aegilops_speltoides|gene001`
- Cleaned name: `Aegilopsspeltoides`

List all species that should be foreground branches in the branch model,
comma-separated.

### Step 3 — VESPA codeml_setup

```bash
bash codeml_inputs/codeml_set_up.sh \
    --root-dir     split_out/divergence/ \
    --docker-image vespa:latest
```

### Step 4 — Correct directories and generate taskfarm files

```bash
bash codeml_inputs/correction_codeml.sh \
    --root-dir  split_out/divergence/ \
    --mode      all

# Two taskfarm files are written:
#   split_out/divergence/codeml_taskfarm_fullpath.sh             (branch model)
#   split_out/divergence/codeml_taskfarm_fullpath_site_model.sh  (all site models)
```

### Step 5 — Run codeml (SLURM)

Both branch model jobs and m0 site model jobs are submitted in a **single SLURM
array**.  The array is split internally: tasks `0..N_branch-1` run the branch
model; tasks `N_branch..N_total-1` run the selected site model(s).

```bash
BRANCH=split_out/divergence/codeml_taskfarm_fullpath.sh
SITE=split_out/divergence/codeml_taskfarm_fullpath_site_model.sh

# Count total tasks (branch + m0):
N=$(bash run_codeml/run_codeml.sbatch --count \
      --taskfarm-branch "$BRANCH" \
      --taskfarm-site   "$SITE" \
      --site-models     m0)

# Submit:
sbatch --array=0-$((N-1))%15 run_codeml/run_codeml.sbatch \
       --taskfarm-branch "$BRANCH" \
       --taskfarm-site   "$SITE" \
       --site-models     m0
```

**Branch model only** (skip site models):

```bash
N=$(wc -l < "$BRANCH")
sbatch --array=0-$((N-1))%15 run_codeml/run_codeml.sbatch \
       --taskfarm-branch "$BRANCH" \
       --site-models     none
```

**Adding more site models** (e.g., branch + m0 + m7 + m8):

```bash
N=$(bash run_codeml/run_codeml.sbatch --count \
      --taskfarm-branch "$BRANCH" --taskfarm-site "$SITE" \
      --site-models m0,m7,m8)
sbatch --array=0-$((N-1))%15 run_codeml/run_codeml.sbatch \
       --taskfarm-branch "$BRANCH" --taskfarm-site "$SITE" \
       --site-models m0,m7,m8
```

Wait for completion: `squeue -u $USER`

### Step 6 — Parse results

```bash
python run_codeml/parse_results/parse_results.py \
    --input-dir   split_out/divergence/ \
    --output-dir  codeml_results/

# Restrict to specific focal species only:
python run_codeml/parse_results/parse_results.py \
    --input-dir       split_out/divergence/ \
    --output-dir      codeml_results/ \
    --focal-species   "Aegilopsspeltoides,Aegilopsmutica"
```

Two TSV files are written to `codeml_results/`:

**`focal_species_results.tsv`** — one row per (focal species × gene):

```
species         id_gene        ortholog     N       S      dN      dS      omega   lnL_branch_model  lnL_m0_model  omega_m0  omega_diff
Aegilopsspeltoides  EVM0020492.1  HOG0005558  1616.5  666.5  0.0108  0.0563  0.1911  -3469.207368      -3477.012    0.1450    0.0461
```

**`all_branches_results.tsv`** — one row per branch per run (includes all species branches, not just the focal one).

---

## codeml model parameters

| Model | `model` | `NSsites` | Description |
|-------|---------|-----------|-------------|
| M0 | 0 | 0 | One ratio — single ω for all branches |
| M1a | 0 | 1 | Nearly neutral — two ω classes |
| M2a | 0 | 2 | Positive selection — three ω classes |
| M7 | 0 | 7 | Beta distribution (null) |
| M8 | 0 | 8 | Beta + positive selection |
| Branch (free ratios) | 2 | 0 | One ω per branch |

All models use:
- `omega = 0.5` initial value (set by correction_codeml.sh)
- `CodonFreq = 7` (F3×4 codon frequencies)
- `cleandata = 1` (skip codons with gaps/ambiguities)
- `icode = 0` (universal genetic code)

---

## Output columns

**`focal_species_results.tsv`** — one row per (focal species × gene):

| Column | Description |
|--------|-------------|
| `species` | Focal species (cleaned name from `cleaned_<species>/` dir) |
| `id_gene` | Gene identifier (directory name below `Inferred_Genetree_*/`) |
| `ortholog` | HOG identifier |
| `N` | Expected non-synonymous sites (focal branch) |
| `S` | Expected synonymous sites (focal branch) |
| `dN` | Non-synonymous substitution rate (focal branch) |
| `dS` | Synonymous substitution rate (focal branch) |
| `omega` | Foreground dN/dS from branch model (model=2) |
| `lnL_branch_model` | Log-likelihood of the branch model |
| `lnL_m0_model` | Log-likelihood of the M0 model (joined by ortholog + gene) |
| `omega_m0` | Global dN/dS from M0 |
| `omega_diff` | `omega` − `omega_m0` (focal species departure from genome-wide mean) |

**`all_branches_results.tsv`** — one row per branch per model=2 run:

| Column | Description |
|--------|-------------|
| `species_focal` | Focal species of this model=2 run |
| `id_gene` | Gene identifier |
| `ortholog` | HOG identifier |
| `branch` | PAML branch label (e.g. `3..1`) |
| `species_branch` | Species at the child node (or `node_N` if name not resolved) |
| `N`, `S` | Expected sites for this branch |
| `dN`, `dS` | Substitution rates for this branch |
| `omega` | dN/dS for this branch |
| `t` | Branch length |
| `lnL_branch_model` | Log-likelihood of this model=2 run |
| `is_foreground` | `True` if this branch is the focal species foreground branch |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `No .fasta files found` in gene_tree.sh | Wrong `--ext` | Check extension in `split_out/divergence/` |
| Gene tree contains `\|` or `_` still | Skipped (FORCE=false) | Add `--force` to gene_tree.sh |
| `Missing branch_table.txt` in codeml_set_up.sh | branch_table.sh not run | Run Step 2 first |
| `Codeml_Setup_codeml_input` already exists, skipped | Previous partial run | It is idempotent — re-run safely |
| codeml exits immediately with error | Wrong `codeml.ctl` | Check that correction_codeml.sh ran successfully |
| Duplicate species in branch TSV | Branch lengths identical | Pass `--species-order` to parse_results.py |
| `No 'out' file` in parse_results.py | codeml failed or not finished | Check `codeml_<jobID>_<taskID>.err` SLURM logs |
| `ERROR: branch taskfarm not found` | Wrong path to taskfarm | Verify `--taskfarm-branch` path |
| `ERROR: task ID N is out of range` | `--array` size set before `--count` | Always run `--count` first, then submit |

---

## Navigation

← [filter_codon_alignment/ — codon filter + split](../README.md)
→ MKT analysis using `polymorphism/` track
