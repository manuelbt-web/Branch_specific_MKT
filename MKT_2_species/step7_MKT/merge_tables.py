#!/usr/bin/env python3
"""
merge_tables.py — Step 7a: Merge polymorphism statistics + dNdSpiNpiS divergence.

Takes the output of polymorphisms_stats.py (Step 6b, which already contains all
polymorphism stats and SDM estimates) and the merged dNdSpiNpiS TSV (Step 5),
and produces a single unified table ready for the MKT analysis in Step 7b.

From polymorphisms_stats.py (Step 6b):
  gene_id, number_sites, num_codons_eff, Pi_per_site, pi_S, pi_NS,
  thetaW_per_site, Tajimas_D, num_sites_S_egglib, num_sites_NS_egglib,
  polymorphism_count, num_pol_S, num_pol_NS, sample_size,
  PnMinus, PnGreater, PsMinus, PsGreater,
  ratioPs, PnNeutral_estimation, PnNeutral, deleterious

From dNdSpiNpiS (Step 5):
  Dn_counts (fixN), Ds_counts (fixS)

Output columns:
  gene_id, number_sites, num_codons_eff, Pi_per_site, pi_S, pi_NS,
  thetaW_per_site, Tajimas_D, num_sites_S_egglib, num_sites_NS_egglib,
  polymorphism_count, num_pol_S, num_pol_NS, sample_size,
  PnMinus, PnGreater, PsMinus, PsGreater,
  ratioPs, PnNeutral_estimation, PnNeutral, deleterious,
  Dn_counts, Ds_counts

Dependencies: pandas  (pip install pandas)

Usage:
  python step7_MKT/merge_tables.py \\
      --dndspnps     all_genes_dNdSpNpS.tsv \\
      --polymorphism step6_stats/polymorphism_stats.tsv \\
      --output       step7_merged.tsv \\
      [--species     Ae_speltoides]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: pandas is required.  Install with: pip install pandas")


# ── Column order in output ────────────────────────────────────────────────────
OUTPUT_COLUMNS = [
    "gene_id",
    "number_sites", "num_codons_eff",
    "Pi_per_site", "pi_S", "pi_NS",
    "thetaW_per_site", "Tajimas_D",
    "num_sites_S_egglib", "num_sites_NS_egglib",
    "polymorphism_count", "num_pol_S", "num_pol_NS",
    "sample_size",
    "PnMinus", "PnGreater",
    "PsMinus", "PsGreater",
    "ratioPs", "PnNeutral_estimation", "PnNeutral", "deleterious",
    "Dn_counts", "Ds_counts",
]


# ── Helpers ───────────────────────────────────────────────────────────────────

def _num(df: pd.DataFrame, col: str) -> "pd.Series":
    return pd.to_numeric(df[col], errors="coerce")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--dndspnps", required=True, metavar="FILE",
        help="Merged dNdSpiNpiS output TSV (all_genes_dNdSpNpS.tsv from Step 5)",
    )
    p.add_argument(
        "--polymorphism", required=True, metavar="FILE",
        help="Polymorphism statistics TSV (polymorphism_stats.tsv from Step 6b). "
             "Must include SDM columns produced by polymorphisms_stats.py.",
    )
    p.add_argument(
        "--output", required=True, metavar="FILE",
        help="Output merged TSV (input for MKT_analysis.Rmd)",
    )
    p.add_argument(
        "--species", default=None, metavar="NAME",
        help="Focal species name (e.g. 'Ae_speltoides'). "
             "Used only for informational output; not required for the merge.",
    )
    return p.parse_args()


# ── Main logic ────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    for path, label in [
        (args.dndspnps,    "--dndspnps"),
        (args.polymorphism, "--polymorphism"),
    ]:
        if not Path(path).is_file():
            sys.exit(f"ERROR: {label} file not found: {path}")

    print("=" * 60)
    print("  Step 7a — Merge tables for MKT analysis")
    print("=" * 60)

    dnd  = pd.read_csv(args.dndspnps,     sep="\t", dtype=str)
    poly = pd.read_csv(args.polymorphism, sep="\t", dtype=str)
    dnd.columns  = [c.strip() for c in dnd.columns]
    poly.columns = [c.strip() for c in poly.columns]

    print(f"  dNdSpiNpiS   : {len(dnd)} genes, {len(dnd.columns)} columns")
    print(f"  Polymorphism : {len(poly)} genes, {len(poly.columns)} columns")
    if args.species:
        print(f"  Species      : {args.species!r}")

    # Verify the polymorphism table has the expected SDM columns
    required_poly = ["gene_id", "PnMinus", "PnGreater", "PsMinus", "PsGreater",
                     "PnNeutral_estimation", "PnNeutral", "deleterious"]
    missing = [c for c in required_poly if c not in poly.columns]
    if missing:
        sys.exit(
            f"ERROR: polymorphism table is missing required columns: {missing}\n"
            f"       Ensure polymorphisms_stats.py was run (which includes SDM estimation)."
        )

    # ── Extract divergence counts from dNdSpiNpiS ─────────────────────────────
    dnd_sub = pd.DataFrame()
    dnd_sub["gene_id"]   = dnd["Contig_name"].str.strip()
    dnd_sub["Dn_counts"] = _num(dnd, "fixN")
    dnd_sub["Ds_counts"] = _num(dnd, "fixS")

    # ── Merge polymorphism stats (full) + divergence ───────────────────────────
    poly["gene_id"] = poly["gene_id"].str.strip()
    out = poly.merge(dnd_sub, on="gene_id", how="left")

    # ── Merge quality report ──────────────────────────────────────────────────
    n_matched = out["Dn_counts"].notna().sum()
    n_miss    = out["Dn_counts"].isna().sum()
    print(f"  Merge result : {n_matched} genes with divergence data | "
          f"{n_miss} without dNdSpiNpiS match")
    if n_miss > 0:
        print("  Unmatched genes (first 10):",
              out.loc[out["Dn_counts"].isna(), "gene_id"].tolist()[:10])

    # ── Reorder columns ───────────────────────────────────────────────────────
    present = [c for c in OUTPUT_COLUMNS if c in out.columns]
    extra   = [c for c in out.columns if c not in OUTPUT_COLUMNS]
    out = out[present + extra]

    # ── Write output ──────────────────────────────────────────────────────────
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False, na_rep="NA")

    print()
    print("=" * 60)
    print(f"  Output : {args.output}")
    print(f"  Rows   : {len(out)}    Columns : {len(out.columns)}")
    print("=" * 60)


if __name__ == "__main__":
    main()
