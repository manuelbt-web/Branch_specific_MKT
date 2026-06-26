#!/usr/bin/env python3
"""
merge_mkt_results.py
====================
Merge three result sources into one MKT-ready analysis table:

  1. Polymorphism statistics     (from polymorphism_stats.py output directory)
  2. dNdSpiNpiS divergence counts (from dNdSpiNpiS.sbatch output directory)
       fixN → renamed to dn_counts
       fixS → renamed to ds_counts
  3. codeml branch-model results  (focal_species_results.tsv from parse_results.py)
       id_gene → joined as gene_id
       N → N_codeml   (avoid collision with polymorphism columns)
       S → S_codeml

JOIN KEY: gene_id
  • polymorphism stats       : 'gene_id' column
  • dNdSpiNpiS output        : filename stem (e.g. EVM0000002.1.out → EVM0000002.1)
                               or first 'alignment' column if present
  • codeml focal_results.tsv : 'id_gene' column (renamed to gene_id before joining)

USAGE
  python merge_mkt_results.py \\
      --polymorphism-dir  poly_stats/ \\
      --divergence-dir    dNdSpiNpiS_output/ \\
      --codeml-file       codeml_results/focal_species_results.tsv \\
      --output            merged_mkt_results.tsv

  # If codeml file contains multiple species, restrict to one:
  python merge_mkt_results.py \\
      --polymorphism-dir  poly_stats/ \\
      --divergence-dir    dNdSpiNpiS_output/ \\
      --codeml-file       codeml_results/focal_species_results.tsv \\
      --focal-species     Aegilopsspeltoides \\
      --output            merged_speltoides.tsv

  # Strip filename suffix from dNdSpiNpiS output filenames before joining:
  python merge_mkt_results.py \\
      ... \\
      --strip-suffix _NT_filtered

OUTPUT COLUMNS (in order)
  gene_id          | ortholog (HOG)    | [all polymorphism columns]
  dn_counts (fixN) | ds_counts (fixS)  | [other dNdSpiNpiS columns]
  N_codeml  | S_codeml | dN | dS | omega | lnL_branch_model | lnL_m0_model
  omega_m0  | omega_diff | species (focal) | lnL (codeml)

DEPENDENCIES
  pip install pandas
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import List, Optional

import pandas as pd


# =============================================================================
# Loading helpers
# =============================================================================

def _strip(name: str, suffix: str) -> str:
    """Strip suffix from a string (case-insensitive)."""
    if suffix and name.lower().endswith(suffix.lower()):
        return name[: len(name) - len(suffix)]
    return name


# ── 1. Polymorphism stats ─────────────────────────────────────────────────────

def load_polymorphism(poly_dir: str) -> pd.DataFrame:
    """
    Load all TSV files from the polymorphism directory and concatenate.
    Expects a 'gene_id' column in each file.
    """
    tsv_files = sorted(
        f for f in os.listdir(poly_dir)
        if f.endswith(".tsv") and not f.startswith(".")
    )
    if not tsv_files:
        sys.exit(f"ERROR: No .tsv files found in polymorphism directory: {poly_dir}")

    dfs: List[pd.DataFrame] = []
    for fname in tsv_files:
        path = os.path.join(poly_dir, fname)
        try:
            df = pd.read_csv(path, sep="\t", dtype=str)
            dfs.append(df)
        except Exception as exc:
            print(f"  [WARN] Cannot read {fname}: {exc}", file=sys.stderr)

    if not dfs:
        sys.exit("ERROR: No polymorphism TSV files could be read.")

    result = pd.concat(dfs, ignore_index=True)

    if "gene_id" not in result.columns:
        sys.exit(
            "ERROR: Polymorphism TSV files must contain a 'gene_id' column.\n"
            f"       Columns found: {list(result.columns)}"
        )

    # Convert numeric columns to float where possible
    result = _coerce_numeric(result, exclude=["gene_id", "downsampled", "nt_engine"])

    print(f"  Polymorphism : {len(result)} genes  ({len(tsv_files)} files)")
    return result


# ── 2. dNdSpiNpiS output ─────────────────────────────────────────────────────

_FIXN_ALIASES = {"fixn", "dn_fix", "fixeddifferencesn", "fixationn", "dn"}
_FIXS_ALIASES = {"fixs", "ds_fix", "fixeddifferencess", "fixations", "ds"}


def load_dnds(dnds_dir: str, strip_suffix: str) -> pd.DataFrame:
    """
    Load all .out files from the dNdSpiNpiS output directory.
    Gene ID is taken from the first 'alignment'-like column if present,
    otherwise from the filename stem (with optional suffix stripping).
    fixN → dn_counts, fixS → ds_counts.
    """
    out_files = sorted(
        f for f in os.listdir(dnds_dir)
        if f.endswith(".out") and not f.startswith("log_") and not f.startswith(".")
    )
    if not out_files:
        sys.exit(f"ERROR: No .out files found in dNdSpiNpiS directory: {dnds_dir}")

    dfs: List[pd.DataFrame] = []
    for fname in out_files:
        # gene_id from filename stem, with optional suffix strip
        stem = re.sub(r"\.out$", "", fname, flags=re.IGNORECASE)
        gene_id_from_filename = _strip(stem, strip_suffix)
        path = os.path.join(dnds_dir, fname)
        try:
            df = pd.read_csv(path, sep=r"\s+", engine="python", dtype=str)
            df.columns = [c.strip() for c in df.columns]

            # Use the filename stem as gene_id for reliable joining.
            # dNdSpiNpiS Contig_name column holds only the gene name within the
            # alignment (e.g. "EVM0000002.1"), not the full HOG-enriched filename
            # stem needed for joining (e.g. "EVM0000002_1_HOG0018301_NT_filtered").
            # Always prefer the filename stem; drop first column if it looks like
            # a name column to avoid duplication.
            gene_id = gene_id_from_filename
            first_col = df.columns[0] if len(df.columns) > 0 else ""
            first_col_lower = first_col.lower()
            is_name_col = first_col_lower in (
                "alignment", "gene", "gene_id", "locus", "name", "id",
                "contig_name", "contig",
            ) or (len(df) > 0 and not str(df.iloc[0, 0]).replace(".", "").isnumeric())
            if is_name_col:
                # Drop the name column — gene_id is added from filename stem
                df = df.iloc[:, 1:]

            df["gene_id"] = gene_id
            dfs.append(df)
        except Exception as exc:
            print(f"  [WARN] Cannot read {fname}: {exc}", file=sys.stderr)

    if not dfs:
        sys.exit("ERROR: No dNdSpiNpiS .out files could be read.")

    result = pd.concat(dfs, ignore_index=True)
    result.columns = [c.strip() for c in result.columns]

    # Rename fixN / fixS to dn_counts / ds_counts
    rename: dict[str, str] = {}
    for col in result.columns:
        col_lower = col.lower()
        if col_lower in _FIXN_ALIASES and "dn_counts" not in rename.values():
            rename[col] = "dn_counts"
        elif col_lower in _FIXS_ALIASES and "ds_counts" not in rename.values():
            rename[col] = "ds_counts"
    # Prefer exact 'fixN' / 'fixS' names over aliases
    for exact, new_name in [("fixN", "dn_counts"), ("fixS", "ds_counts")]:
        if exact in result.columns:
            rename[exact] = new_name
    if rename:
        result = result.rename(columns=rename)

    if "dn_counts" not in result.columns:
        print(
            "  [WARN] Column 'fixN' not found in dNdSpiNpiS output — "
            "dn_counts will be NA.\n"
            f"         Columns found: {list(result.columns)}",
            file=sys.stderr,
        )
    if "ds_counts" not in result.columns:
        print(
            "  [WARN] Column 'fixS' not found in dNdSpiNpiS output — "
            "ds_counts will be NA.",
            file=sys.stderr,
        )

    result = _coerce_numeric(result, exclude=["gene_id"])
    print(f"  dNdSpiNpiS   : {len(result)} genes  ({len(out_files)} files)")
    return result


# ── 3. codeml focal species results ──────────────────────────────────────────

def load_codeml(
    codeml_file: str,
    focal_species: Optional[str],
    strip_suffix: str,
) -> pd.DataFrame:
    """
    Load focal_species_results.tsv from parse_results.py.
    Filters to --focal-species if given.
    Renames id_gene → gene_id, N → N_codeml, S → S_codeml.
    """
    df = pd.read_csv(codeml_file, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]

    if focal_species:
        if "species" not in df.columns:
            print(
                "[WARN] codeml file has no 'species' column — "
                "cannot filter by --focal-species",
                file=sys.stderr,
            )
        else:
            before = len(df)
            df = df[df["species"] == focal_species]
            if df.empty:
                avail = pd.read_csv(codeml_file, sep="\t")["species"].unique().tolist()
                sys.exit(
                    f"ERROR: --focal-species '{focal_species}' not found.\n"
                    f"       Available: {avail}"
                )
            print(
                f"  codeml       : filtered to '{focal_species}' "
                f"({len(df)}/{before} rows kept)"
            )
    else:
        if "species" in df.columns and df["species"].nunique() > 1:
            print(
                f"  [WARN] codeml file contains "
                f"{df['species'].nunique()} species: "
                f"{df['species'].unique().tolist()}\n"
                f"         Use --focal-species to restrict.",
                file=sys.stderr,
            )
        print(f"  codeml       : {len(df)} rows")

    # Rename columns for the merge
    rename = {}
    if "id_gene" in df.columns:
        rename["id_gene"] = "gene_id"
    if "N" in df.columns:
        rename["N"] = "N_codeml"
    if "S" in df.columns:
        rename["S"] = "S_codeml"
    if rename:
        df = df.rename(columns=rename)

    # Strip suffix from gene_id if requested
    if "gene_id" in df.columns and strip_suffix:
        df["gene_id"] = df["gene_id"].apply(lambda x: _strip(str(x), strip_suffix))

    df = _coerce_numeric(
        df, exclude=["gene_id", "species", "ortholog"]
    )
    return df


# =============================================================================
# Merge
# =============================================================================

def merge_all(
    poly: pd.DataFrame,
    dnds: pd.DataFrame,
    codeml: pd.DataFrame,
) -> pd.DataFrame:
    """
    Left-join:  polymorphism  ←→  dNdSpiNpiS  ←→  codeml
    All joins are on gene_id.
    """
    # Select dNdSpiNpiS columns to carry forward
    # Always include dn_counts / ds_counts; also bring other dNdSpiNpiS metrics
    # that don't duplicate polymorphism columns.
    dnds_keep = ["gene_id"] + [
        c for c in dnds.columns
        if c != "gene_id" and c not in poly.columns
    ]
    # Make sure dn_counts and ds_counts are included even if above filter catches them
    for col in ("dn_counts", "ds_counts"):
        if col in dnds.columns and col not in dnds_keep:
            dnds_keep.insert(1, col)
    dnds_sub = dnds[[c for c in dnds_keep if c in dnds.columns]].copy()

    # Step 1: polymorphism ← dNdSpiNpiS
    merged = pd.merge(poly, dnds_sub, on="gene_id", how="left")
    n_div = merged["dn_counts"].notna().sum() if "dn_counts" in merged.columns else 0
    print(f"  poly ← dNdSpiNpiS : {len(merged)} rows, {n_div} matched")

    # Select codeml columns to carry forward
    codeml_keep = ["gene_id"] + [
        c for c in codeml.columns
        if c != "gene_id" and c not in merged.columns
    ]
    # Always include the key codeml columns
    for col in ("ortholog", "N_codeml", "S_codeml", "dN", "dS", "omega",
                "lnL_branch_model", "lnL_m0_model", "omega_m0", "omega_diff",
                "species"):
        if col in codeml.columns and col not in codeml_keep:
            codeml_keep.append(col)
    codeml_sub = codeml[[c for c in codeml_keep if c in codeml.columns]].copy()

    # Step 2: merged ← codeml
    merged = pd.merge(merged, codeml_sub, on="gene_id", how="left")
    n_codeml = merged["ortholog"].notna().sum() if "ortholog" in merged.columns else 0
    print(f"  poly ← codeml     : {len(merged)} rows, {n_codeml} matched")

    return merged


# =============================================================================
# Utilities
# =============================================================================

def _coerce_numeric(df: pd.DataFrame, exclude: List[str]) -> pd.DataFrame:
    """Try to convert each column (not in exclude) to numeric."""
    for col in df.columns:
        if col in exclude:
            continue
        try:
            df[col] = pd.to_numeric(df[col], errors="ignore")
        except Exception:
            pass
    return df


# =============================================================================
# Argument parsing
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-p", "--polymorphism-dir", required=True, metavar="DIR",
        help="Directory of per-gene polymorphism statistics TSV files "
             "(output of polymorphism_stats.py).",
    )
    p.add_argument(
        "-d", "--divergence-dir", required=True, metavar="DIR",
        help="Directory of per-gene dNdSpiNpiS .out files "
             "(output of dNdSpiNpiS.sbatch). "
             "Columns fixN and fixS are renamed to dn_counts and ds_counts.",
    )
    p.add_argument(
        "-c", "--codeml-file", required=True, metavar="FILE",
        help="focal_species_results.tsv produced by parse_results.py. "
             "Contains columns: species, id_gene, ortholog, N, S, dN, dS, omega, ...",
    )
    p.add_argument(
        "-o", "--output", required=True, metavar="FILE",
        help="Output TSV file path for the merged results table.",
    )
    p.add_argument(
        "--focal-species", default="", metavar="NAME",
        help="Filter codeml results to this species (cleaned name, no underscores). "
             "Required when focal_species_results.tsv contains multiple species. "
             "Example: Aegilopsspeltoides",
    )
    p.add_argument(
        "--strip-suffix", default="_NT_filtered", metavar="SUFFIX",
        help="Suffix to strip from filenames when building the gene_id join key. "
             "Applied to dNdSpiNpiS filename stems "
             "(default: '_NT_filtered'). Use '' to disable.",
    )
    p.add_argument(
        "--verbose", action="store_true",
        help="Print full column list of the merged output.",
    )
    return p.parse_args()


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()

    # Validate inputs
    for path, label in [
        (args.polymorphism_dir, "--polymorphism-dir"),
        (args.divergence_dir, "--divergence-dir"),
    ]:
        if not os.path.isdir(path):
            sys.exit(f"ERROR: Directory not found ({label}): {path}")
    if not os.path.isfile(args.codeml_file):
        sys.exit(f"ERROR: File not found (--codeml-file): {args.codeml_file}")

    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print("=" * 54)
    print("  merge_mkt_results.py")
    print("=" * 54)
    print(f"  Polymorphism dir  : {args.polymorphism_dir}")
    print(f"  dNdSpiNpiS dir    : {args.divergence_dir}")
    print(f"  codeml file       : {args.codeml_file}")
    if args.focal_species:
        print(f"  Focal species     : {args.focal_species}")
    print(f"  Strip suffix      : '{args.strip_suffix}'")
    print(f"  Output            : {args.output}")
    print("=" * 54)
    print()

    print("Loading inputs:")
    poly   = load_polymorphism(args.polymorphism_dir)
    dnds   = load_dnds(args.divergence_dir, args.strip_suffix)
    codeml = load_codeml(
        args.codeml_file,
        args.focal_species if args.focal_species else None,
        args.strip_suffix,
    )

    print("\nMerging:")
    merged = merge_all(poly, dnds, codeml)

    # Write output
    merged.to_csv(args.output, sep="\t", index=False)

    print(f"\nDone.")
    print(f"  Output : {args.output}")
    print(f"  Rows   : {len(merged)}")
    print(f"  Columns: {len(merged.columns)}")
    if args.verbose:
        print(f"  Column list: {list(merged.columns)}")

    # Warn about low match rates
    n_poly = len(merged)
    for col, label in [("dn_counts", "dNdSpiNpiS"), ("ortholog", "codeml")]:
        if col in merged.columns:
            n_matched = merged[col].notna().sum()
            pct = 100 * n_matched / n_poly if n_poly > 0 else 0
            if pct < 80:
                print(
                    f"  [WARN] Only {pct:.0f}% of genes matched {label} data "
                    f"({n_matched}/{n_poly}). Check --strip-suffix or gene ID format."
                )


if __name__ == "__main__":
    main()
