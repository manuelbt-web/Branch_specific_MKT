#!/usr/bin/env python3
"""
estimation_of_SDM.py — Step 6c: Estimate strongly deleterious mutations (SDMs)
from the per-gene folded SFS produced by polymorphisms_stats.py.

For each gene, the folded allele-frequency spectrum is split at a cutoff T
into a low-frequency fraction (<= T, candidate SDMs) and a high-frequency
fraction (> T, effectively neutral).  The SDM count is estimated as:

  P_SDM = P_N^(<T) − (P_N^(>T) / P_S^(>T)) × P_S^(<T)

where P_N and P_S denote non-synonymous and synonymous polymorphism counts
and the neutral expectation ratio P_N^(>T) / P_S^(>T) is used to rescale the
low-frequency synonymous fraction.

Reference: Fay et al. 2001 (Genetics); Messer & Petrov 2013 (PNAS).

Input:
  --sfs FILE    folded SFS TSV from polymorphisms_stats.py --sfs-output
                Expected columns: gene, sample_size, bin, frequency,
                                  syn_count, nonsyn_count

Output:
  --output FILE per-gene TSV with SDM estimates:
                gene, sample_size, syn_count, nonsyn_count,
                PnMinus, PnGreater, PsMinus, PsGreater,
                ratioPs, deleterious, PnNeutral

Usage:
  python estimation_of_SDM.py \\
      --sfs     step6_stats/folded_sfs.tsv \\
      --cutoff  0.15 \\
      --output  step6_stats/sdm_per_gene.tsv

Dependencies: pandas  (pip install pandas)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: pandas is required.  Install with: pip install pandas")

import math


OUTPUT_COLUMNS = [
    "gene", "sample_size",
    "syn_count", "nonsyn_count",
    "PnMinus", "PnGreater",
    "PsMinus", "PsGreater",
    "ratioPs", "deleterious", "PnNeutral",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--sfs", required=True, metavar="FILE",
        help="Per-gene folded SFS TSV (--sfs-output from polymorphisms_stats.py)",
    )
    p.add_argument(
        "--cutoff", type=float, default=0.15, metavar="FLOAT",
        help="Folded allele-frequency cutoff T. Sites with freq <= T are treated "
             "as low-frequency (candidate SDMs); sites with freq > T are treated "
             "as high-frequency (effectively neutral). Default: 0.15",
    )
    p.add_argument(
        "--output", required=True, metavar="FILE",
        help="Output TSV with per-gene SDM estimates",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not Path(args.sfs).is_file():
        sys.exit(f"ERROR: SFS file not found: {args.sfs}")

    print("=" * 60)
    print("  SDM estimation from folded SFS")
    print("=" * 60)
    print(f"  SFS file : {args.sfs}")
    print(f"  Cutoff T : {args.cutoff}")
    print(f"  Output   : {args.output}")
    print("=" * 60)

    df = pd.read_csv(args.sfs, sep="\t", dtype={"gene": str})

    # Compute frequency from bin/sample_size when the column is missing or empty
    if "frequency" not in df.columns or df["frequency"].isna().all():
        df["frequency"] = pd.to_numeric(df["bin"], errors="coerce") / pd.to_numeric(df["sample_size"], errors="coerce")
    else:
        df["frequency"] = pd.to_numeric(df["frequency"], errors="coerce")
        mask = df["frequency"].isna() & df["sample_size"].notna()
        df.loc[mask, "frequency"] = (
            pd.to_numeric(df.loc[mask, "bin"],         errors="coerce") /
            pd.to_numeric(df.loc[mask, "sample_size"], errors="coerce")
        )

    df["syn_count"]    = pd.to_numeric(df["syn_count"],    errors="coerce").fillna(0)
    df["nonsyn_count"] = pd.to_numeric(df["nonsyn_count"], errors="coerce").fillna(0)

    cutoff = args.cutoff
    rows: list[dict] = []

    for gene, gdf in df.groupby("gene", sort=False):
        # sample_size: take first non-NA value
        ss_vals = pd.to_numeric(gdf["sample_size"], errors="coerce").dropna()
        sample_size = int(ss_vals.iloc[0]) if not ss_vals.empty else None

        Ps = gdf["syn_count"].sum()     # total synonymous polymorphisms P_S
        Pn = gdf["nonsyn_count"].sum()  # total NS polymorphisms P_N

        freq      = gdf["frequency"]
        mask_low  = freq <= cutoff
        mask_high = freq >  cutoff

        PnMinus   = gdf.loc[mask_low,  "nonsyn_count"].sum()   # P_N^(<T)
        PnGreater = gdf.loc[mask_high, "nonsyn_count"].sum()   # P_N^(>T)
        PsMinus   = gdf.loc[mask_low,  "syn_count"].sum()      # P_S^(<T)
        PsGreater = gdf.loc[mask_high, "syn_count"].sum()      # P_S^(>T)

        # Neutral Pn/Ps ratio from high-frequency category
        ratioPs = (PnGreater / PsGreater) if PsGreater > 0 else float("nan")

        # Estimated SDMs: P_SDM = P_N^(<T) − (P_N^(>T)/P_S^(>T)) × P_S^(<T)
        if math.isnan(ratioPs):
            deleterious = float("nan")
            PnNeutral   = float("nan")
        else:
            deleterious = PnMinus - ratioPs * PsMinus
            # Correct neutral Pn; cap to [0, Pn]
            PnNeutral = float(Pn) - deleterious
            PnNeutral = max(0.0, min(PnNeutral, float(Pn)))

        rows.append({
            "gene":         gene,
            "sample_size":  sample_size,
            "syn_count":    int(Ps),
            "nonsyn_count": int(Pn),
            "PnMinus":      int(PnMinus),
            "PnGreater":    int(PnGreater),
            "PsMinus":      int(PsMinus),
            "PsGreater":    int(PsGreater),
            "ratioPs":      ratioPs,
            "deleterious":  deleterious,
            "PnNeutral":    PnNeutral,
        })

    out = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False, na_rep="NA", float_format="%.6g")

    n_total      = len(out)
    n_computable = out["PnNeutral"].notna().sum()
    print(f"\n  Done — {n_total} genes written to {args.output}")
    print(f"  SDM computable : {n_computable}/{n_total} genes (PsGreater > 0)")
    if n_total > n_computable:
        print(f"  Skipped        : {n_total - n_computable} genes had no high-freq S polymorphisms")
    print("=" * 60)


if __name__ == "__main__":
    main()
