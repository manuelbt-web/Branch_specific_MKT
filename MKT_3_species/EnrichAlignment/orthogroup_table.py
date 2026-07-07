#!/usr/bin/env python3
"""
orthogroup_table.py
===================
Build an orthogroup-to-species gene-ID mapping table from per-HOG FASTA files.

Each input FASTA file must:
  - Have a filename containing the HOG identifier (e.g., HOG0000001_aln_hmm.fasta)
  - Have sequence headers formatted as:  >species_name|gene_id

OUTPUT
  A TSV with one row per orthogroup:
    orthogroup  Species_A  Species_B  Species_C  ...

This table is the required input for MACSE_enrichment.py.

USAGE
  python orthogroup_table.py \\
      --input-dir  nt_alignments/ \\
      --output     orthogroup_table.tsv

  # Custom separator (default is '|'):
  python orthogroup_table.py \\
      --input-dir  nt_alignments/ \\
      --output     orthogroup_table.tsv \\
      --sep        '|'
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit(
        "ERROR: pandas is required.\n"
        "  Install: pip install pandas   or   conda install pandas"
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-i", "--input-dir", required=True, metavar="DIR",
        help="Directory containing per-HOG FASTA files. "
             "Filenames must include 'HOGxxxx' (e.g. HOG0000001_aln_hmm.fasta).",
    )
    p.add_argument(
        "-o", "--output", default="orthogroup_table.tsv", metavar="FILE",
        help="Output TSV file. Default: orthogroup_table.tsv",
    )
    p.add_argument(
        "-e", "--ext", default=None, metavar="EXT",
        help="FASTA file extension to search for. "
             "Default: tries .fasta, .fa, .faa",
    )
    p.add_argument(
        "--sep", default="|", metavar="CHAR",
        help="Separator between species name and gene ID in FASTA headers. "
             "Default: '|'  (e.g. '>Aegilops_speltoides|TraesCS1A02G000001')",
    )
    return p.parse_args()


def collect_files(input_dir: Path, ext: str | None) -> list[Path]:
    """Return sorted list of HOG FASTA files in input_dir."""
    exts = [ext.lstrip(".")] if ext else ["fasta", "fa", "faa"]
    files: list[Path] = []
    for e in exts:
        files.extend(
            f for f in sorted(input_dir.glob(f"*.{e}"))
            if re.search(r"HOG", f.name, re.IGNORECASE)
        )
    return sorted(set(files), key=lambda f: f.name)


def extract_hog_id(filename: str) -> str | None:
    """Extract the HOG identifier from a filename (e.g. 'HOG0000001')."""
    m = re.search(r"(HOG\d+)", filename, re.IGNORECASE)
    return m.group(1).upper() if m else None


def parse_fasta_headers(fasta_path: Path, sep: str) -> dict[str, str]:
    """
    Return {species: gene_id} from FASTA headers.
    Skips headers that do not contain the separator (warns).
    If a species appears twice, keeps the first occurrence.
    """
    row: dict[str, str] = {}
    bad_headers: list[str] = []
    with open(fasta_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if sep not in header:
                bad_headers.append(header[:80])
                continue
            species, gene = header.split(sep, 1)
            species = species.strip()
            gene    = gene.strip()
            if species not in row:
                row[species] = gene
    return row, bad_headers


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)

    if not input_dir.is_dir():
        sys.exit(f"ERROR: Directory not found: {input_dir}")

    files = collect_files(input_dir, args.ext)
    if not files:
        ext_tried = args.ext or "fasta / fa / faa"
        sys.exit(
            f"ERROR: No HOG FASTA files found in {input_dir}\n"
            f"  Extension tried : {ext_tried}\n"
            f"  File names must contain 'HOG' (e.g. HOG0000001.fasta).\n"
            f"  Use --ext to specify a different extension."
        )

    print(f"Found {len(files)} FASTA file(s) in {input_dir}")

    data: dict[str, dict[str, str]] = {}
    skipped:  list[str] = []
    warnings: list[str] = []

    for i, fasta in enumerate(files, 1):
        hog = extract_hog_id(fasta.name)
        if hog is None:
            skipped.append(f"{fasta.name}: no HOG identifier found in filename")
            continue

        row, bad = parse_fasta_headers(fasta, args.sep)
        for h in bad:
            warnings.append(
                f"{fasta.name}: header without '{args.sep}' separator skipped: '{h}'"
            )

        if not row:
            skipped.append(f"{fasta.name}: no valid sequence headers found")
            continue

        if hog in data:
            warnings.append(f"Duplicate HOG '{hog}' — overwriting with {fasta.name}")
        data[hog] = row

        if len(files) > 20 and i % max(1, len(files) // 10) == 0:
            print(f"  {i}/{len(files)} files processed…")

    # ── Report issues ──────────────────────────────────────────────────────────
    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for w in warnings[:20]:
            print(f"  {w}")
        if len(warnings) > 20:
            print(f"  … and {len(warnings) - 20} more")

    if skipped:
        print(f"\nSkipped ({len(skipped)}):")
        for s in skipped:
            print(f"  {s}")

    if not data:
        sys.exit("ERROR: No valid orthogroup data found. Check the file format.")

    # ── Drop spurious species columns from contaminated headers ─────────────────
    # A "species" appearing in only a handful of rows is almost always not a
    # real species name but a gene/individual ID that leaked into a divergence
    # FASTA (e.g. a mis-split file whose header was ">EVM0001702.1|sp|..."
    # instead of ">Aegilops_speltoides|..."). Real species names appear in most
    # or all orthogroups. Keep only species present in at least 1% of files
    # (minimum 2), and report anything dropped so the underlying data issue is
    # visible rather than silently producing bogus one-off columns.
    all_species: dict[str, int] = {}
    for row in data.values():
        for sp in row:
            all_species[sp] = all_species.get(sp, 0) + 1

    min_files = max(2, round(0.01 * len(files)))
    real_species = {sp for sp, n in all_species.items() if n >= min_files}
    dropped_species = {sp: n for sp, n in all_species.items() if sp not in real_species}

    if dropped_species:
        print(
            f"\nDropped {len(dropped_species)} spurious 'species' column(s) "
            f"(present in fewer than {min_files} file(s) — likely contaminated "
            f"headers, e.g. a gene ID where a species name was expected):"
        )
        for sp, n in sorted(dropped_species.items(), key=lambda kv: -kv[1])[:20]:
            print(f"  {sp}  ({n} file(s))")
        if len(dropped_species) > 20:
            print(f"  … and {len(dropped_species) - 20} more")
        data = {
            hog: {sp: gene for sp, gene in row.items() if sp in real_species}
            for hog, row in data.items()
        }

    # ── Build DataFrame ────────────────────────────────────────────────────────
    df = pd.DataFrame.from_dict(data, orient="index")
    df.index.name = "orthogroup"
    df = df.reindex(sorted(df.columns), axis=1)
    df = df.reset_index()

    # ── Write output ───────────────────────────────────────────────────────────
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, sep="\t", index=False)

    n_species = len(df.columns) - 1
    species_names = sorted(df.columns[1:].tolist())

    print(f"\nOrthogroup table written to: {output}")
    print(f"  Orthogroups : {len(df)}")
    print(f"  Species     : {n_species}")
    for sp in species_names:
        n_present = df[sp].notna().sum()
        print(f"    {sp:<40s}  ({n_present}/{len(df)} orthogroups)")

    print(f"\nNext step — run MACSE enrichment:")
    print(f"  python MACSE_enrichment.py \\")
    print(f"      --ortho-table {output} \\")
    print(f"      --ortho-dir   <nt_alignment_dir> \\")
    print(f"      --poly-dir    <cds_polymorphism_dir> \\")
    print(f"      --output-dir  macse_enriched/ \\")
    print(f"      --macse-jar   macse_v2.07.jar \\")
    print(f"      --species-col <SPECIES>   # choose from: {', '.join(species_names)}")


if __name__ == "__main__":
    main()
