#!/usr/bin/env python3
"""
prepare_recmap.py
=================
Convert a recombination map text file (e.g. from Glemin et al.
ms-rec-triticeae, https://github.com/sylvainglemin/ms-rec-triticeae)
to a 4-column BED file for use with recombinaison.Rmd.

OUTPUT FORMAT (tab-separated, no header):
    chr    start    end    rec_rate

  rec_rate is in the SAME units as the source file (cM/Mb by default
  in Glemin et al. maps). The GLM parameter rec_scale in
  render_recombinaison.R must match; default is 1e6 (converts
  cM/Mb → cM/bp equivalent used internally).

USAGE
  # Preview column names in the file before choosing --rec-col
  python prepare_recmap.py \\
      --input  Ae_speltoides_genomeB_1cM_remapped.txt \\
      --preview

  # Convert a Glemin et al. file (defaults match the actual format):
  python prepare_recmap.py \\
      --input  Ae_speltoides_genomeB_1cM_remapped.txt \\
      --out    data/Ae_speltoides_recmap.bed

  # Equivalent explicit form:
  python prepare_recmap.py \\
      --input     Ae_speltoides_genomeB_1cM_remapped.txt \\
      --chr-col   Chromosome \\
      --start-col Start \\
      --end-col   End \\
      --rec-col   recRate \\
      --out       data/Ae_speltoides_recmap.bed

  # Convert using column indices (0-based) if the file has no header
  python prepare_recmap.py \\
      --input  recmap.txt \\
      --no-header \\
      --col-indices 0,1,2,3 \\
      --out  data/recmap.bed

  # Scale rec_rate (divide by 1e6 to convert cM/Mb → per-bp fraction)
  python prepare_recmap.py \\
      --input  recmap.txt --rec-col mean_rho \\
      --scale  1e-6 \\
      --out    data/recmap.bed

  # Add or strip chromosome name prefix
  python prepare_recmap.py \\
      --input recmap.txt --rec-col rho \\
      --chr-prefix Chr \\
      --out data/recmap.bed

GLEMIN ET AL. FORMAT (ms-rec-triticeae)
  The recombination maps at
  https://github.com/sylvainglemin/ms-rec-triticeae/tree/main/outputs/recombination
  are space-separated text files with a header row. The format is:

    poscM  Chromosome  Start  End  recRate  nb_complete_site  piSyn  f0

  Example rows:
    0.02  1B  1418981.8    1425703.2      NA                   1261  0.00251  0.292
    0.06  1B  4302527.33   4307091    7.69e-07                  440  0        NA

  Notes:
  - Chromosome names do NOT include a prefix (e.g. "1B", not "Chr1B")
  - Start / End can be floating-point — they are converted to integers automatically
  - recRate can be NA — those intervals are skipped
  - The default column names in this script match this format exactly

DEPENDENCIES
  Standard library only — no external packages required.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from typing import Optional


# =============================================================================
# Parsing
# =============================================================================

def detect_sep(line: str) -> str:
    """Guess field separator from a single line."""
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return None  # will use split() (handles multiple spaces)


def read_recmap(
    input_path: str,
    chr_col: str,
    start_col: str,
    end_col: str,
    rec_col: str,
    no_header: bool,
    col_indices: Optional[list[int]],
    scale: float,
    filter_negative: bool,
    chr_prefix: Optional[str],
    chr_strip: Optional[str],
) -> list[tuple[str, int, int, float]]:

    records: list[tuple[str, int, int, float]] = []
    n_skipped_neg  = 0
    n_skipped_na   = 0
    n_skipped_coord = 0

    with open(input_path, encoding="utf-8") as fh:
        first_data_line = None
        header_dict: dict[str, int] = {}

        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            sep = detect_sep(line)
            fields = line.split(sep) if sep else line.split()
            fields = [f.strip() for f in fields]

            if no_header or col_indices:
                # No-header mode: use column indices directly
                if col_indices is None:
                    sys.exit("ERROR: --no-header requires --col-indices.")
                if lineno == 1:
                    first_data_line = fields
                try:
                    ci, si, ei, ri = col_indices
                    chrom = fields[ci]
                    start = int(fields[si])
                    end   = int(fields[ei])
                    rec   = float(fields[ri])
                except (IndexError, ValueError) as exc:
                    print(f"  [WARN] Line {lineno}: parse error ({exc}) — skipped.",
                          file=sys.stderr)
                    n_skipped_na += 1
                    continue
            else:
                if not header_dict:
                    # First non-comment line is the header
                    header_dict = {name: i for i, name in enumerate(fields)}
                    continue

                if first_data_line is None:
                    first_data_line = fields

                def get_col(name: str, lineno: int) -> Optional[str]:
                    idx = header_dict.get(name)
                    if idx is None:
                        return None
                    try:
                        return fields[idx]
                    except IndexError:
                        return None

                chrom_val = get_col(chr_col, lineno)
                start_val = get_col(start_col, lineno)
                end_val   = get_col(end_col, lineno)
                rec_val   = get_col(rec_col, lineno)

                if any(v is None for v in (chrom_val, start_val, end_val, rec_val)):
                    n_skipped_na += 1
                    continue

                chrom = chrom_val
                try:
                    start = int(float(start_val))
                    end   = int(float(end_val))
                except ValueError:
                    n_skipped_coord += 1
                    continue

                try:
                    rec = float(rec_val)
                except ValueError:
                    if rec_val.lower() in ("na", "nan", ".", ""):
                        n_skipped_na += 1
                    else:
                        print(f"  [WARN] Line {lineno}: cannot parse rec_rate "
                              f"'{rec_val}' — skipped.", file=sys.stderr)
                        n_skipped_na += 1
                    continue

            if math.isnan(rec):
                n_skipped_na += 1
                continue

            if filter_negative and rec < 0:
                n_skipped_neg += 1
                continue

            rec_scaled = rec * scale

            # Chromosome normalisation
            if chr_strip and chrom.startswith(chr_strip):
                chrom = chrom[len(chr_strip):]
            if chr_prefix:
                chrom = chr_prefix + chrom

            records.append((chrom, start, end, rec_scaled))

    if n_skipped_na:
        print(f"  [INFO] {n_skipped_na} rows skipped: missing or non-numeric value.",
              file=sys.stderr)
    if n_skipped_neg:
        print(f"  [INFO] {n_skipped_neg} rows skipped: negative rec_rate "
              "(--filter-negative).", file=sys.stderr)
    if n_skipped_coord:
        print(f"  [WARN] {n_skipped_coord} rows skipped: non-integer coordinates.",
              file=sys.stderr)

    return records


# =============================================================================
# Preview
# =============================================================================

def preview(input_path: str, n: int = 5) -> None:
    """Show header and first few data lines of the recombination map."""
    with open(input_path, encoding="utf-8") as fh:
        lines = []
        for line in fh:
            line = line.rstrip("\n")
            if line and not line.startswith("#"):
                lines.append(line)
            if len(lines) >= n + 1:
                break

    if not lines:
        print("File appears empty or contains only comments.")
        return

    sep = detect_sep(lines[0])
    header_fields = lines[0].split(sep) if sep else lines[0].split()

    print("\n── RECMAP PREVIEW ───────────────────────────────────────")
    print(f"  File     : {input_path}")
    print(f"  Separator: {'TAB' if sep == chr(9) else 'comma' if sep == ',' else 'whitespace'}")
    print(f"\n  Columns ({len(header_fields)}):")
    defaults = {"Chromosome": "--chr-col", "Start": "--start-col",
                "End": "--end-col", "recRate": "--rec-col"}
    for i, col in enumerate(header_fields):
        col = col.strip()
        tag = f"  ← default {defaults[col]}" if col in defaults else ""
        print(f"    [{i}] {col}{tag}")
    print(f"\n  First {min(n, len(lines) - 1)} data rows:")
    for line in lines[1:n + 1]:
        print(f"    {line[:120]}")
    print("─────────────────────────────────────────────────────────")
    print("\nDefaults (Chromosome / Start / End / recRate) match the")
    print("Glemin et al. ms-rec-triticeae format. If your file uses")
    print("different names, pass --chr-col, --start-col, --end-col,")
    print("--rec-col when running without --preview.\n")


# =============================================================================
# Argument parsing
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-i", "--input", required=True, metavar="FILE",
        help="Input recombination map text file (tab- or comma-separated).",
    )
    p.add_argument(
        "-o", "--out", metavar="FILE", default=None,
        help="Output BED file. Required unless --preview is used.",
    )
    # Column specification
    col_group = p.add_argument_group("column specification (header-based, default)")
    col_group.add_argument(
        "--chr-col", default="Chromosome", metavar="NAME",
        help="Column name for chromosome. Default: Chromosome "
             "(matches Glemin et al. ms-rec-triticeae format).",
    )
    col_group.add_argument(
        "--start-col", default="Start", metavar="NAME",
        help="Column name for interval start (bp). Default: Start. "
             "Float values are automatically rounded to the nearest integer.",
    )
    col_group.add_argument(
        "--end-col", default="End", metavar="NAME",
        help="Column name for interval end (bp). Default: End. "
             "Float values are automatically rounded to the nearest integer.",
    )
    col_group.add_argument(
        "--rec-col", default="recRate", metavar="NAME",
        help="Column name for recombination rate. Default: recRate "
             "(matches Glemin et al. ms-rec-triticeae format). "
             "NA values are skipped. Run --preview to see available column names.",
    )
    # No-header mode
    no_hdr = p.add_argument_group("no-header mode (column indices)")
    no_hdr.add_argument(
        "--no-header", action="store_true",
        help="Use this flag if the file has no header row.",
    )
    no_hdr.add_argument(
        "--col-indices", metavar="CHR,START,END,REC", default=None,
        help="Comma-separated 0-based column indices for chr, start, end, rec_rate. "
             "Requires --no-header. Example: --col-indices 0,1,2,3",
    )
    # Transformations
    p.add_argument(
        "--scale", type=float, default=1.0, metavar="FACTOR",
        help="Multiply rec_rate by this factor before writing. "
             "Default: 1.0 (no scaling). "
             "Use 1e-6 to convert cM/Mb → cM/bp if needed.",
    )
    p.add_argument(
        "--filter-negative", action="store_true",
        help="Discard intervals with negative rec_rate "
             "(some LDhat outputs include these as artefacts).",
    )
    p.add_argument(
        "--chr-prefix", metavar="PREFIX", default=None,
        help="Add prefix to chromosome names "
             "(e.g. 'Chr' to convert '1A' → 'Chr1A').",
    )
    p.add_argument(
        "--chr-strip", metavar="PREFIX", default=None,
        help="Strip prefix from chromosome names "
             "(e.g. 'Chr' to convert 'Chr1A' → '1A').",
    )
    p.add_argument(
        "--preview", action="store_true",
        help="Print column names and sample rows; do not write output.",
    )
    return p.parse_args()


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f"ERROR: Input file not found: {args.input}")

    if args.preview:
        preview(args.input)
        return

    if args.out is None:
        sys.exit("ERROR: --out is required (unless --preview is used).")

    col_indices: Optional[list[int]] = None
    if args.col_indices:
        try:
            col_indices = [int(x) for x in args.col_indices.split(",")]
            if len(col_indices) != 4:
                raise ValueError
        except ValueError:
            sys.exit("ERROR: --col-indices must be exactly 4 comma-separated integers "
                     "(e.g. 0,1,2,3).")

    print("=" * 54)
    print("  prepare_recmap.py")
    print("=" * 54)
    print(f"  Input        : {args.input}")
    if args.no_header and col_indices:
        print(f"  Col indices  : chr={col_indices[0]} start={col_indices[1]} "
              f"end={col_indices[2]} rec={col_indices[3]}")
    else:
        print(f"  chr col      : {args.chr_col}")
        print(f"  start col    : {args.start_col}")
        print(f"  end col      : {args.end_col}")
        print(f"  rec_rate col : {args.rec_col}")
    if args.scale != 1.0:
        print(f"  Scale        : {args.scale}")
    if args.chr_prefix:
        print(f"  Add chr pfx  : {args.chr_prefix}")
    if args.chr_strip:
        print(f"  Strip chr pfx: {args.chr_strip}")
    print(f"  Output       : {args.out}")
    print("=" * 54)
    print()

    print("Parsing recombination map...")
    records = read_recmap(
        input_path     = args.input,
        chr_col        = args.chr_col,
        start_col      = args.start_col,
        end_col        = args.end_col,
        rec_col        = args.rec_col,
        no_header      = args.no_header,
        col_indices    = col_indices,
        scale          = args.scale,
        filter_negative = args.filter_negative,
        chr_prefix     = args.chr_prefix,
        chr_strip      = args.chr_strip,
    )

    if not records:
        sys.exit(
            "ERROR: No records extracted.\n"
            "       Run --preview to check column names, then adjust --rec-col etc."
        )

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.out, "w", encoding="utf-8") as fh:
        for chrom, start, end, rec in records:
            fh.write(f"{chrom}\t{start}\t{end}\t{rec:.10g}\n")

    n_chrs = len({r[0] for r in records})
    rec_vals = [r[3] for r in records]
    print(f"Done.")
    print(f"  Intervals written : {len(records)}")
    print(f"  Chromosomes       : {n_chrs}")
    print(f"  rec_rate range    : [{min(rec_vals):.6g}, {max(rec_vals):.6g}]")
    print(f"  rec_rate median   : {sorted(rec_vals)[len(rec_vals)//2]:.6g}")
    print(f"  Output            : {args.out}")

    print("\nSample output (first 5 lines):")
    for rec in records[:5]:
        print(f"  {rec[0]}\t{rec[1]}\t{rec[2]}\t{rec[3]:.10g}")

    print(
        "\nNEXT STEP: verify that chromosome names match gene_coords_bed.\n"
        f"  cut -f1 {args.out} | sort -u\n"
        "  cut -f1 data/gene_coords.bed | sort -u"
    )


if __name__ == "__main__":
    main()
