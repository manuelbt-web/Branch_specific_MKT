#!/usr/bin/env python3
"""
prepare_gene_coords.py
======================
Extract gene genomic coordinates from a GFF3 annotation file and write
them as a 4-column BED file (no header) for use with recombinaison.Rmd.

OUTPUT FORMAT (tab-separated, no header):
    chr    start    end    gene_id

  Coordinates are 0-based half-open (BED convention):
    start = GFF column 4  − 1
    end   = GFF column 5  (unchanged)

  gene_id is taken from the GFF attribute specified by --id-attr (default: ID).

USAGE
  # Basic: extract all "gene" features
  python prepare_gene_coords.py \\
      --gff  annotation.gff3 \\
      --out  data/speltoides_gene_coords.bed

  # Only mRNA features, using 'Name' attribute as gene_id
  python prepare_gene_coords.py \\
      --gff  annotation.gff3 \\
      --feature mRNA \\
      --id-attr Name \\
      --out  data/speltoides_gene_coords.bed

  # Strip version suffix (.1) from gene IDs: EVM.model.contig.1.1 → EVM.model.contig.1
  python prepare_gene_coords.py \\
      --gff  annotation.gff3 \\
      --strip-version \\
      --out  data/speltoides_gene_coords.bed

  # Preview only — print first 10 lines and unique chromosomes, do not write
  python prepare_gene_coords.py \\
      --gff  annotation.gff3 \\
      --preview

COMMON GFF3 ID ATTRIBUTE NAMES
  MAKER / EVidenceModeler  : ID=gene:EVM0000001  or  ID=EVM.model.contig.1.1
  Ensembl plants           : ID=gene:TRAES3B000001MC
  BRAKER / AUGUSTUS        : ID=gene1
  Use --id-attr to adjust if needed.

DEPENDENCIES
  Standard library only — no external packages required.
"""

from __future__ import annotations

import argparse
import os
import re
import sys


# =============================================================================
# GFF3 parsing
# =============================================================================

def parse_attributes(attr_string: str) -> dict[str, str]:
    """Parse a GFF3 attribute string into a dict."""
    attrs: dict[str, str] = {}
    for item in attr_string.strip().split(";"):
        item = item.strip()
        if "=" in item:
            key, _, val = item.partition("=")
            attrs[key.strip()] = val.strip()
    return attrs


def extract_id(attr_dict: dict[str, str], id_attr: str, strip_prefix: bool) -> str | None:
    """
    Extract the gene identifier from the parsed attribute dict.

    Tries id_attr first, then falls back to 'Name', then 'gene_id'.
    If strip_prefix is True, removes everything up to the last colon
    (e.g. 'gene:TRAES001' → 'TRAES001').
    """
    val = attr_dict.get(id_attr)
    if val is None:
        for fallback in ("Name", "gene_id", "transcript_id"):
            val = attr_dict.get(fallback)
            if val is not None:
                break
    if val is None:
        return None
    if strip_prefix and ":" in val:
        val = val.rsplit(":", 1)[1]
    return val


def read_gff(
    gff_path: str,
    feature_type: str,
    id_attr: str,
    strip_prefix: bool,
    strip_version: bool,
    chr_prefix: str | None,
    chr_strip: str | None,
) -> list[tuple[str, int, int, str]]:
    """
    Parse the GFF3 file and return a list of (chr, start0, end, gene_id) tuples.
    start0 is 0-based (BED convention).
    """
    records: list[tuple[str, int, int, str]] = []
    seen_ids: set[str] = set()
    skipped_no_id = 0
    skipped_dup    = 0

    if gff_path.endswith(".gz"):
        import gzip
        opener = lambda p: gzip.open(p, "rt", encoding="utf-8")
    else:
        opener = lambda p: open(p, encoding="utf-8")

    with opener(gff_path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            feat = parts[2]
            if feat != feature_type:
                continue

            try:
                start1 = int(parts[3])
                end    = int(parts[4])
            except ValueError:
                print(f"  [WARN] Line {lineno}: non-integer coordinates — skipped.",
                      file=sys.stderr)
                continue

            chrom = parts[0]
            # Chromosome name normalisation
            if chr_strip and chrom.startswith(chr_strip):
                chrom = chrom[len(chr_strip):]
            if chr_prefix:
                chrom = chr_prefix + chrom

            attrs = parse_attributes(parts[8])
            gene_id = extract_id(attrs, id_attr, strip_prefix)

            if gene_id is None:
                skipped_no_id += 1
                continue

            if strip_version:
                # Remove trailing .version (e.g. .1, .2) from gene IDs
                gene_id = re.sub(r"\.\d+$", "", gene_id)

            if gene_id in seen_ids:
                skipped_dup += 1
                continue
            seen_ids.add(gene_id)

            records.append((chrom, start1 - 1, end, gene_id))

    if skipped_no_id:
        print(f"  [WARN] {skipped_no_id} features skipped: no '{id_attr}' attribute found.",
              file=sys.stderr)
    if skipped_dup:
        print(f"  [INFO] {skipped_dup} duplicate gene_ids deduplicated (first occurrence kept).",
              file=sys.stderr)

    return records


# =============================================================================
# Preview mode
# =============================================================================

def preview(gff_path: str, feature_type: str, n: int = 10) -> None:
    """Print a summary of the GFF3 file without writing any output."""
    feature_counts: dict[str, int] = {}
    chromosomes: set[str] = set()
    sample_lines: list[str] = []
    n_sample = 0

    with open(gff_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            feat = parts[2]
            feature_counts[feat] = feature_counts.get(feat, 0) + 1
            chromosomes.add(parts[0])
            if feat == feature_type and n_sample < n:
                sample_lines.append(line.rstrip())
                n_sample += 1

    print("\n── GFF3 PREVIEW ─────────────────────────────────────────")
    print(f"  File    : {gff_path}")
    print(f"\n  Feature type counts:")
    for k, v in sorted(feature_counts.items(), key=lambda x: -x[1])[:15]:
        marker = "  ← selected" if k == feature_type else ""
        print(f"    {k:30s}  {v:>8d}{marker}")
    print(f"\n  Unique chromosomes ({len(chromosomes)}):")
    for c in sorted(chromosomes)[:20]:
        print(f"    {c}")
    if len(chromosomes) > 20:
        print(f"    ... ({len(chromosomes) - 20} more)")
    print(f"\n  First {n_sample} '{feature_type}' feature lines:")
    for l in sample_lines:
        print(f"    {l[:120]}")
    print("─────────────────────────────────────────────────────────\n")
    print("To extract gene coordinates, run without --preview.")


# =============================================================================
# Argument parsing
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-g", "--gff", required=True, metavar="FILE",
        help="Input GFF3 annotation file (may be gzipped with .gz extension "
             "if gzip is available).",
    )
    p.add_argument(
        "-o", "--out", metavar="FILE", default=None,
        help="Output BED file. Required unless --preview is used.",
    )
    p.add_argument(
        "--feature", default="gene", metavar="TYPE",
        help="GFF3 feature type to extract (column 3). Default: gene.",
    )
    p.add_argument(
        "--id-attr", default="ID", metavar="ATTR",
        help="GFF3 attribute key containing the gene identifier. "
             "Default: ID. Common alternatives: Name, gene_id, Parent.",
    )
    p.add_argument(
        "--strip-prefix", action="store_true",
        help="Strip prefix up to the last colon in gene IDs "
             "(e.g. 'gene:TRAES001' → 'TRAES001'). "
             "Useful for Ensembl-formatted annotations.",
    )
    p.add_argument(
        "--strip-version", action="store_true",
        help="Remove trailing version suffix (.1, .2, ...) from gene IDs "
             "(e.g. 'EVM.model.Chr1A.1.10' → 'EVM.model.Chr1A.1'). "
             "Use only if the MKT pipeline uses un-versioned IDs.",
    )
    p.add_argument(
        "--chr-prefix", metavar="PREFIX", default=None,
        help="Add this prefix to all chromosome names in the output "
             "(e.g. 'Chr' to convert '1A' → 'Chr1A').",
    )
    p.add_argument(
        "--chr-strip", metavar="PREFIX", default=None,
        help="Strip this prefix from all chromosome names "
             "(e.g. 'chr' to convert 'chr1A' → '1A').",
    )
    p.add_argument(
        "--preview", action="store_true",
        help="Print GFF3 summary and sample lines; do not write output. "
             "Use this to identify the correct --feature and --id-attr values.",
    )
    return p.parse_args()


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()

    if not os.path.isfile(args.gff):
        sys.exit(f"ERROR: GFF file not found: {args.gff}")

    if args.preview:
        preview(args.gff, args.feature)
        return

    if args.out is None:
        sys.exit("ERROR: --out is required (unless --preview is used).")

    print("=" * 54)
    print("  prepare_gene_coords.py")
    print("=" * 54)
    gz = args.gff.endswith(".gz")
    print(f"  Input GFF    : {args.gff}{' (gzip)' if gz else ''}")
    print(f"  Feature type : {args.feature}")
    print(f"  ID attribute : {args.id_attr}")
    print(f"  Strip prefix : {args.strip_prefix}")
    print(f"  Strip version: {args.strip_version}")
    if args.chr_prefix:
        print(f"  Add chr prefix: {args.chr_prefix}")
    if args.chr_strip:
        print(f"  Strip chr prefix: {args.chr_strip}")
    print(f"  Output BED   : {args.out}")
    print("=" * 54)
    print()

    print("Parsing GFF3...")
    records = read_gff(
        args.gff,
        feature_type  = args.feature,
        id_attr       = args.id_attr,
        strip_prefix  = args.strip_prefix,
        strip_version = args.strip_version,
        chr_prefix    = args.chr_prefix,
        chr_strip     = args.chr_strip,
    )

    if not records:
        sys.exit(
            f"ERROR: No '{args.feature}' features extracted.\n"
            "       Check --feature and --id-attr; run --preview to inspect the file."
        )

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.out, "w", encoding="utf-8") as fh:
        for chrom, start, end, gene_id in records:
            fh.write(f"{chrom}\t{start}\t{end}\t{gene_id}\n")

    n_chrs = len({r[0] for r in records})
    print(f"Done.")
    print(f"  Genes extracted : {len(records)}")
    print(f"  Chromosomes     : {n_chrs}")
    print(f"  Output          : {args.out}")

    # Quick sanity check
    print("\nSample output (first 5 lines):")
    for rec in records[:5]:
        print(f"  {rec[0]}\t{rec[1]}\t{rec[2]}\t{rec[3]}")

    print(
        "\nNEXT STEP: verify that gene_id values match what is in "
        "branch_specific_MKT_results.tsv.\n"
        "  head -1 branch_specific_MKT_results.tsv  # check gene_id column\n"
        "  head -5 " + args.out
    )


if __name__ == "__main__":
    main()
