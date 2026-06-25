#!/usr/bin/env python3
"""
parse_results.py
================
Parse PAML codeml output from a VESPA branch-model pipeline and write two TSV
result tables combining branch-model (model=2) and M0 results.

DIRECTORY STRUCTURE EXPECTED
  <input-dir>/
    Inferred_Genetree_<HOG>/
      <gene_id>/
        Codeml_Setup_codeml_input/
          cleaned/m0/Omega0_5/out          ← M0 (global omega, null model)
          cleaned_<species>/model_branch/Omega0_5/out   ← branch model per focal species

OUTPUTS (in <output-dir>/)

  focal_species_results.tsv
    One row per (focal species × gene).  Combines branch-model foreground branch
    values with M0 values joined on (ortholog, gene_id).
    Columns:
      species         cleaned name of focal species (from cleaned_<species>/ dir)
      id_gene         gene directory name (e.g. EVM0020492.1)
      ortholog        HOG identifier (e.g. HOG0005558)
      N               expected non-synonymous sites (focal branch)
      S               expected synonymous sites (focal branch)
      dN              non-synonymous substitution rate (focal branch)
      dS              synonymous substitution rate (focal branch)
      omega           dN/dS foreground omega (from model=2)
      lnL_branch_model  log-likelihood of branch model
      lnL_m0_model    log-likelihood of M0 model
      omega_m0        global dN/dS from M0
      omega_diff      omega (branch) − omega_m0

  all_branches_results.tsv
    One row per branch per (focal species × gene).  Contains ALL branches from
    each model=2 run, allowing inspection of non-focal branches as well.
    Columns:
      species_focal   focal species of this model=2 run
      id_gene         gene directory name
      ortholog        HOG identifier
      branch          PAML branch label (e.g. 3..1)
      species_branch  species at the child node of this branch (or node_N)
      N               expected non-synonymous sites
      S               expected synonymous sites
      dN              non-synonymous substitution rate
      dS              synonymous substitution rate
      omega           dN/dS for this branch
      t               branch length
      lnL_branch_model  log-likelihood of this model=2 run
      is_foreground   True if this branch is the focal (foreground) species branch

USAGE
  python parse_results.py \\
      --input-dir   split_out/divergence/ \\
      --output-dir  codeml_results/

  # Restrict to specific focal species:
  python parse_results.py \\
      --input-dir       split_out/divergence/ \\
      --output-dir      codeml_results/ \\
      --focal-species   "Aegilopsspeltoides,Aegilopsmutica"

DEPENDENCIES
  pip install pandas
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import pandas as pd


# =============================================================================
# Path info extraction
# =============================================================================

def extract_path_info(out_file: str) -> Dict[str, Optional[str]]:
    """
    Extract HOG id, gene_id, and focal species from a codeml 'out' file path.

    Expected path segments:
      .../<root>/Inferred_Genetree_<HOG>/<gene_id>/Codeml_Setup_codeml_input/
          cleaned_<species>/model_branch/Omega0_5/out      (branch model)
          cleaned/m0/Omega0_5/out                          (M0 model)
    """
    parts = os.path.normpath(out_file).split(os.sep)
    hog: Optional[str] = None
    gene_id: Optional[str] = None
    species: Optional[str] = None

    for i, part in enumerate(parts):
        # Match Inferred_Genetree_HOG<N> or Inferred_Genetree_<N>
        m = re.match(r"Inferred_Genetree_((?:HOG)?\d+)", part, re.IGNORECASE)
        if m:
            raw = m.group(1)
            hog = raw if raw.upper().startswith("HOG") else f"HOG{raw}"
            # The directory immediately after Inferred_Genetree_*/ is the gene id
            if i + 1 < len(parts):
                gene_id = parts[i + 1]

        # cleaned_<species>/ → focal species for branch model
        if part.startswith("cleaned_") and len(part) > 8:
            species = part[8:]   # strip "cleaned_"

    return {"hog": hog, "gene_id": gene_id, "species": species}


# =============================================================================
# PAML output parsing — low-level helpers
# =============================================================================

def _extract_lnL(txt: str) -> Optional[float]:
    """Extract log-likelihood from codeml output text."""
    for pat in [
        r"lnL\([^)]*\)\s*:\s*([\-0-9\.eE+]+)",
        r"lnL\s*=\s*([\-0-9\.eE+]+)",
    ]:
        m = re.search(pat, txt)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                pass
    return None


def _safe_float(s: str) -> Optional[float]:
    try:
        v = float(s)
        return v
    except (ValueError, TypeError):
        return None


def _extract_node_species_map(txt: str) -> Dict[int, str]:
    """
    Parse PAML's labeled tree section to map node numbers → species names.

    PAML prints terminal nodes as 'SpeciesName_NodeNum' in the labeled tree:
      ((Aegilopsspeltoides_1: 0.015, Aegilopstauschii_2: 0.045)_3: ...)
    Since species names have been cleaned (no underscores), the trailing _N
    is the only underscore, making the regex unambiguous.
    """
    node_map: Dict[int, str] = {}

    # Try "Tree with node labels for Rod Page's TreeView"
    m = re.search(
        r"(?:Tree with node labels[^\n]*|TreeView[^\n]*)\s*\n\s*([^\n]+)",
        txt, re.IGNORECASE
    )
    if not m:
        # Fallback: any Newick-like line with _N suffixes
        m = re.search(r"\([^\n]*_\d+[^\n]*\)\s*;", txt)

    if m:
        tree_str = m.group(0) if not m.lastindex else m.group(1)
        for mm in re.finditer(r"([A-Za-z][A-Za-z0-9]+)_(\d+)", tree_str):
            node_map[int(mm.group(2))] = mm.group(1)

    return node_map


def _parse_branch_table(txt: str, node_map: Dict[int, str]) -> List[Dict]:
    """
    Parse the per-branch dN/dS table in PAML codeml output.

    PAML branch table format:
      branch     t       N       S    dN/dS     dN      dS
      3..1   0.0000  1616.5   666.5   0.1911   0.0108   0.0563
    """
    rows: List[Dict] = []

    # Find the table section by its header
    header_m = re.search(
        r"^\s*branch\s+t\s+N\s+S\s+dN/dS\s+dN\s+dS\s*$(.+?)(?:\n\s*\n|\Z)",
        txt, re.MULTILINE | re.DOTALL | re.IGNORECASE
    )
    if not header_m:
        return rows

    for line in header_m.group(1).splitlines():
        line = line.strip()
        if not line:
            continue
        # Match: parent..child  t  N  S  dN/dS  dN  dS
        m = re.match(
            r"(\d+)\.\.(\d+)\s+"
            r"([\-0-9\.eE+]+)\s+"           # t
            r"([0-9\.eE+NnAa\-]+)\s+"       # N (can be NaN)
            r"([0-9\.eE+NnAa\-]+)\s+"       # S
            r"([0-9\.eE+NnAa\-]+)\s+"       # dN/dS
            r"([0-9\.eE+NnAa\-]+)\s+"       # dN
            r"([0-9\.eE+NnAa\-]+)",         # dS
            line
        )
        if not m:
            continue

        parent_node = int(m.group(1))
        child_node  = int(m.group(2))
        rows.append({
            "branch":         f"{parent_node}..{child_node}",
            "parent_node":    parent_node,
            "child_node":     child_node,
            "species_branch": node_map.get(child_node, f"node_{child_node}"),
            "t":    _safe_float(m.group(3)),
            "N":    _safe_float(m.group(4)),
            "S":    _safe_float(m.group(5)),
            "dNdS": _safe_float(m.group(6)),
            "dN":   _safe_float(m.group(7)),
            "dS":   _safe_float(m.group(8)),
        })

    return rows


# =============================================================================
# Per-file parsers
# =============================================================================

def parse_m0_out(out_file: str) -> Dict:
    """
    Parse M0 codeml output. Returns {lnL_m0, omega_m0}.
    """
    try:
        with open(out_file, encoding="utf-8", errors="replace") as f:
            txt = f.read()
    except OSError:
        return {}

    lnL = _extract_lnL(txt)

    omega = None
    for pat in [
        r"omega\s*\(dN/dS\)\s*=\s*([\-0-9\.eE+]+)",
        r"(?:^|\s)omega\s*=\s*([\-0-9\.eE+]+)",
        r"dN/dS\s*=\s*([\-0-9\.eE+]+)",
    ]:
        mm = re.search(pat, txt, re.MULTILINE | re.IGNORECASE)
        if mm:
            omega = _safe_float(mm.group(1))
            if omega is not None:
                break

    return {"lnL_m0": lnL, "omega_m0": omega}


def parse_branch_out(out_file: str, focal_species: str) -> Dict:
    """
    Parse a branch-model (model=2, NSsites=0) codeml output file.

    Returns:
        lnL_branch    log-likelihood
        fg_omega      foreground omega (first value in 'w (dN/dS) for branches')
        focal_branch  dict for the focal (foreground) branch row, or None
        all_branches  list of dicts for every branch in the tree
    """
    try:
        with open(out_file, encoding="utf-8", errors="replace") as f:
            txt = f.read()
    except OSError:
        return {}

    lnL = _extract_lnL(txt)

    # Foreground omega is the FIRST value printed for model=2
    fg_omega: Optional[float] = None
    m = re.search(r"w \(dN/dS\) for branches:\s+([\-0-9\.eE+]+)", txt)
    if m:
        fg_omega = _safe_float(m.group(1))

    # Node → species name
    node_map = _extract_node_species_map(txt)

    # Per-branch estimates
    branches = _parse_branch_table(txt, node_map)

    # Find the focal (foreground) branch
    # Strategy 1: match by species name from node_map
    focal: Optional[Dict] = None
    for b in branches:
        if b["species_branch"] == focal_species:
            focal = b
            break

    # Strategy 2: the foreground terminal branch has dNdS == fg_omega
    if focal is None and fg_omega is not None:
        candidates = [
            b for b in branches
            if b["dNdS"] is not None and abs(b["dNdS"] - fg_omega) < 1e-4
        ]
        # Prefer terminal branches (those that appear in node_map)
        terminal = [b for b in candidates if b["child_node"] in node_map]
        if terminal:
            focal = terminal[0]
        elif candidates:
            focal = candidates[0]

    return {
        "lnL_branch": lnL,
        "fg_omega":    fg_omega,
        "focal_branch": focal,
        "all_branches": branches,
    }


# =============================================================================
# Argument parsing
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-i", "--input-dir", required=True, metavar="DIR",
        help="Root directory containing Inferred_Genetree_*/ subdirectories.",
    )
    p.add_argument(
        "-o", "--output-dir", required=True, metavar="DIR",
        help="Directory where focal_species_results.tsv and "
             "all_branches_results.tsv will be written.",
    )
    p.add_argument(
        "--focal-species", default="", metavar="LIST",
        help="Comma-separated focal species names (cleaned, no underscores) to "
             "include. Default: all species found in cleaned_*/ directories. "
             "Example: \"Aegilopsspeltoides,Aegilopsmutica\"",
    )
    p.add_argument(
        "--verbose", action="store_true",
        help="Print one line per parsed file.",
    )
    return p.parse_args()


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()

    base_dir = os.path.abspath(args.input_dir)
    if not os.path.isdir(base_dir):
        sys.exit(f"ERROR: Input directory not found: {base_dir}")

    out_dir = os.path.abspath(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)

    focal_filter = set(
        s.strip() for s in args.focal_species.split(",") if s.strip()
    )

    # ── Parse M0 results keyed by (hog, gene_id) ─────────────────────────────
    m0_pattern = os.path.join(
        base_dir, "**", "Inferred_Genetree_*", "*",
        "Codeml_Setup_codeml_input", "cleaned", "m0", "Omega0_5", "out"
    )
    m0_files = glob.glob(m0_pattern, recursive=True)
    print(f"Found {len(m0_files)} M0 out file(s)")

    m0_data: Dict[Tuple[str, str], Dict] = {}
    for f in m0_files:
        info = extract_path_info(f)
        if not info["hog"] or not info["gene_id"]:
            continue
        key = (info["hog"], info["gene_id"])
        m0_data[key] = parse_m0_out(f)

    # ── Parse branch model results ────────────────────────────────────────────
    branch_pattern = os.path.join(
        base_dir, "**", "Inferred_Genetree_*", "*",
        "Codeml_Setup_codeml_input", "cleaned_*",
        "model_branch", "Omega0_5", "out"
    )
    branch_files = glob.glob(branch_pattern, recursive=True)
    print(f"Found {len(branch_files)} branch model out file(s)")

    focal_rows:      List[Dict] = []
    all_branch_rows: List[Dict] = []

    skipped = 0

    for f in branch_files:
        info = extract_path_info(f)
        hog     = info["hog"]
        gene_id = info["gene_id"]
        species = info["species"]

        if not hog or not gene_id or not species:
            if args.verbose:
                print(f"  [SKIP] Cannot parse path: {f}")
            skipped += 1
            continue

        if focal_filter and species not in focal_filter:
            continue

        result = parse_branch_out(f, species)
        if not result:
            skipped += 1
            continue

        lnL_branch  = result.get("lnL_branch")
        fg_omega    = result.get("fg_omega")
        focal_br    = result.get("focal_branch")
        all_branches = result.get("all_branches", [])

        # Join with M0 on (hog, gene_id)
        key = (hog, gene_id)
        m0  = m0_data.get(key, {})
        lnL_m0   = m0.get("lnL_m0")
        omega_m0 = m0.get("omega_m0")

        omega_diff: Optional[float] = None
        if fg_omega is not None and omega_m0 is not None:
            omega_diff = round(fg_omega - omega_m0, 6)

        # ── Focal species row ──────────────────────────────────────────────
        if focal_br is not None:
            focal_rows.append({
                "species":           species,
                "id_gene":           gene_id,
                "ortholog":          hog,
                "N":                 focal_br.get("N"),
                "S":                 focal_br.get("S"),
                "dN":                focal_br.get("dN"),
                "dS":                focal_br.get("dS"),
                "omega":             fg_omega,
                "lnL_branch_model":  lnL_branch,
                "lnL_m0_model":      lnL_m0,
                "omega_m0":          omega_m0,
                "omega_diff":        omega_diff,
            })
        else:
            print(
                f"  [WARN] {hog}/{gene_id}/{species}: foreground branch not found "
                f"in output — row omitted from focal_species_results.tsv",
                file=sys.stderr,
            )
            skipped += 1

        # ── All-branches rows ──────────────────────────────────────────────
        for br in all_branches:
            # A branch is foreground if: species name matches, OR dN/dS matches fg_omega
            is_fg = br["species_branch"] == species or (
                fg_omega is not None
                and br.get("dNdS") is not None
                and abs(br["dNdS"] - fg_omega) < 1e-4
            )
            all_branch_rows.append({
                "species_focal":     species,
                "id_gene":           gene_id,
                "ortholog":          hog,
                "branch":            br["branch"],
                "species_branch":    br["species_branch"],
                "N":                 br.get("N"),
                "S":                 br.get("S"),
                "dN":                br.get("dN"),
                "dS":                br.get("dS"),
                "omega":             br.get("dNdS"),
                "t":                 br.get("t"),
                "lnL_branch_model":  lnL_branch,
                "is_foreground":     is_fg,
            })

        if args.verbose:
            status = "OK  " if focal_br else "WARN"
            print(f"  [{status}] {hog} / {gene_id} / {species}")

    # ── Write output files ────────────────────────────────────────────────────
    FOCAL_COLS = [
        "species", "id_gene", "ortholog",
        "N", "S", "dN", "dS", "omega",
        "lnL_branch_model", "lnL_m0_model", "omega_m0", "omega_diff",
    ]
    ALL_BRANCH_COLS = [
        "species_focal", "id_gene", "ortholog",
        "branch", "species_branch",
        "N", "S", "dN", "dS", "omega", "t",
        "lnL_branch_model", "is_foreground",
    ]

    focal_tsv  = os.path.join(out_dir, "focal_species_results.tsv")
    allbr_tsv  = os.path.join(out_dir, "all_branches_results.tsv")

    pd.DataFrame(focal_rows, columns=FOCAL_COLS).to_csv(
        focal_tsv, sep="\t", index=False
    )
    pd.DataFrame(all_branch_rows, columns=ALL_BRANCH_COLS).to_csv(
        allbr_tsv, sep="\t", index=False
    )

    print(f"\nDone.")
    print(f"  {focal_tsv}  ({len(focal_rows)} rows)")
    print(f"  {allbr_tsv} ({len(all_branch_rows)} rows)")
    if skipped:
        print(f"  {skipped} file(s) / branch(es) skipped — check warnings above.")


if __name__ == "__main__":
    main()
