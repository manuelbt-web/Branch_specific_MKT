#!/usr/bin/env python3
"""
trim_outgroup.py
================
Prepare codon alignments for standard 2-species MKT by retaining only ONE
outgroup sequence from the divergence block.

After filter_codon_alignment.py, each FASTA contains polymorphism sequences
followed by D divergence sequences in this order:

  [0]       focal-species reference  → REMOVED  (redundant with polymorphism)
  [1..D-2]  first outgroup(s)        → KEPT     (divergence used in MKT)
  [D-1]     second outgroup          → REMOVED  (not needed for standard MKT)

Example (D=3, focal = Aegilops speltoides):
  Input:
    >EVM0000002.1|sp|Tr156_...   ← polymorphism (kept)
    ...
    >Aegilops_speltoides|EVM0000002.1            ← div[0]  REMOVED
    >Aegilops_mutica|Ammut_EIv1.0_0527770.1      ← div[1]  KEPT
    >Aegilops_tauschii|transcript_AET4G...        ← div[2]  REMOVED

  Output:
    >EVM0000002.1|sp|Tr156_...   ← polymorphism (kept)
    ...
    >Aegilops_mutica|Ammut_EIv1.0_0527770.1      ← single outgroup (kept)

Trimmed files are written to <output-dir>/alignment_for_standard_MKT/.
A summary TSV is written to the same directory.

USAGE
  # Explicit divergence identification (recommended):
  python trim_outgroup.py \\
      --input-dir          filtered_aln/ \\
      --output-dir         mkt_ready/ \\
      --divergence-pattern "Aegilops_speltoides|Aegilops_mutica|Aegilops_tauschii"

  # Single file:
  python trim_outgroup.py \\
      --input              filtered_aln/HOG001_NT_filtered.fasta \\
      --output-dir         mkt_ready/ \\
      --divergence-pattern "Aegilops_speltoides|Aegilops_mutica|Aegilops_tauschii"

  # Exact headers (comma-separated):
  python trim_outgroup.py \\
      --input-dir  filtered_aln/ \\
      --output-dir mkt_ready/ \\
      --divergence "Aegilops_speltoides|EVM0000002.1,Aegilops_mutica|...,..."

  # Last-N heuristic (default N=3, no pattern needed):
  python trim_outgroup.py \\
      --input-dir  filtered_aln/ \\
      --output-dir mkt_ready/

OPTIONS
  --input FILE               Single FASTA file to process
  --input-dir DIR            Directory of FASTA files to process (batch mode)
  -o / --output-dir DIR      Root output directory; files go to
                             <output-dir>/alignment_for_standard_MKT/
  --divergence-pattern REGEX Python regex matched against sequence headers
                             to identify divergence sequences (recommended)
  --divergence LIST          Comma-separated exact headers for divergence seqs
  --n-divergence N           Fallback: last N sequences are divergence (default: 3)
                             Used only when no --divergence-pattern or --divergence
  --ext EXT                  File extension in --input-dir mode (default: fasta)
  --threads N                Parallel workers in --input-dir mode (default: 1)
  --verbose                  Print one line per file

DEPENDENCIES
  Python 3.8+, no external packages required.
"""

import argparse
import csv
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple


# =============================================================================
# FASTA I/O
# =============================================================================

def read_fasta(path: str) -> List[Tuple[str, str]]:
    """Read FASTA; return list of (header, sequence) without the leading '>'."""
    seqs: List[Tuple[str, str]] = []
    header = ""
    parts: List[str] = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    seqs.append((header, "".join(parts)))
                header = line[1:]
                parts = []
            else:
                parts.append(line)
    if header:
        seqs.append((header, "".join(parts)))
    return seqs


def write_fasta(seqs: List[Tuple[str, str]], path: str, width: int = 60) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for header, seq in seqs:
            f.write(f">{header}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i : i + width] + "\n")


# =============================================================================
# Classification and trimming
# =============================================================================

def classify(
    seqs: List[Tuple[str, str]],
    pattern: Optional[re.Pattern],
    n_div: int,
) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    """Return (divergence_seqs, polymorphism_seqs)."""
    if pattern:
        div  = [(h, s) for h, s in seqs if pattern.search(h)]
        poly = [(h, s) for h, s in seqs if not pattern.search(h)]
    else:
        n    = min(n_div, len(seqs))
        div  = seqs[-n:]
        poly = seqs[:-n]
    return div, poly


def trim(
    seqs: List[Tuple[str, str]],
    pattern: Optional[re.Pattern],
    n_div: int,
) -> Tuple[List[Tuple[str, str]], Dict]:
    """
    Remove first and last divergence sequences; keep all polymorphism and
    middle divergence sequences.  Returns (trimmed_seqs, report_dict).
    """
    div, poly = classify(seqs, pattern, n_div)

    report: Dict = {
        "n_poly":        len(poly),
        "n_div_in":      len(div),
        "removed_first": "",
        "removed_last":  "",
        "n_div_kept":    0,
        "warn":          "",
        "error":         "",
    }

    if len(div) == 0:
        report["error"] = "no divergence sequences found"
        return poly, report

    report["removed_first"] = div[0][0]

    if len(div) == 1:
        report["warn"] = (
            "only 1 divergence sequence found — it is the focal-species "
            "reference and will be removed; output will have no outgroup"
        )
        report["removed_last"] = "(same as first)"
        return poly, report

    report["removed_last"] = div[-1][0]

    if len(div) == 2:
        report["warn"] = (
            "only 2 divergence sequences; after removing first (focal ref) "
            "and last (second outgroup) no divergence sequences remain"
        )
        return poly, report

    # Keep indices 1 … D-2 (middle outgroup(s))
    kept_div = div[1:-1]
    report["n_div_kept"] = len(kept_div)

    return poly + kept_div, report


# =============================================================================
# Per-file worker
# =============================================================================

def process_file(
    in_path: str,
    out_path: str,
    pattern: Optional[re.Pattern],
    n_div: int,
) -> Dict:
    """Read, trim, write one file. Returns a report dict."""
    seqs = read_fasta(in_path)
    report: Dict = {
        "file":          os.path.basename(in_path),
        "n_poly":        0,
        "n_div_in":      0,
        "n_div_kept":    0,
        "removed_first": "",
        "removed_last":  "",
        "warn":          "",
        "error":         "",
    }

    if not seqs:
        report["error"] = "empty or unreadable file"
        return report

    trimmed, r = trim(seqs, pattern, n_div)
    report.update(r)

    if trimmed and not r.get("error"):
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        write_fasta(trimmed, out_path)
    elif not trimmed:
        report["error"] = (report.get("error") or "") + " — nothing written"

    return report


# =============================================================================
# Argument parsing
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    inp = p.add_mutually_exclusive_group(required=True)
    inp.add_argument("--input", metavar="FILE",
                     help="Single FASTA file to process.")
    inp.add_argument("--input-dir", metavar="DIR",
                     help="Directory of FASTA files (batch mode).")

    p.add_argument("-o", "--output-dir", required=True, metavar="DIR",
                   help="Output root; trimmed files go to "
                        "<output-dir>/alignment_for_standard_MKT/")

    div_grp = p.add_mutually_exclusive_group()
    div_grp.add_argument(
        "--divergence-pattern", metavar="REGEX",
        help="Regex matched against headers to flag divergence sequences. "
             "Example: \"Aegilops_speltoides|Aegilops_mutica|Aegilops_tauschii\"",
    )
    div_grp.add_argument(
        "--divergence", metavar="LIST",
        help="Comma-separated exact sequence headers for divergence sequences.",
    )

    p.add_argument(
        "--n-divergence", type=int, default=3, metavar="N",
        help="Fallback: treat the last N sequences as divergence (default: 3). "
             "Ignored when --divergence-pattern or --divergence is given.",
    )
    p.add_argument(
        "--ext", default="fasta", metavar="EXT",
        help="File extension to scan in --input-dir mode (default: fasta).",
    )
    p.add_argument(
        "--threads", type=int, default=1, metavar="N",
        help="Parallel workers for --input-dir mode (default: 1).",
    )
    p.add_argument("--verbose", action="store_true",
                   help="Print one status line per file.")
    return p.parse_args()


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()

    # Build the pattern used to classify divergence sequences
    pattern: Optional[re.Pattern] = None
    if args.divergence_pattern:
        try:
            pattern = re.compile(args.divergence_pattern)
        except re.error as exc:
            sys.exit(f"ERROR: Invalid --divergence-pattern regex: {exc}")
    elif args.divergence:
        headers = [h.strip() for h in args.divergence.split(",") if h.strip()]
        pattern = re.compile("|".join(re.escape(h) for h in headers))

    # Output subdirectory
    out_subdir = os.path.join(
        os.path.abspath(args.output_dir), "alignment_for_standard_MKT"
    )
    os.makedirs(out_subdir, exist_ok=True)

    # Collect input files
    if args.input:
        if not os.path.isfile(args.input):
            sys.exit(f"ERROR: File not found: {args.input}")
        in_files = [os.path.abspath(args.input)]
    else:
        if not os.path.isdir(args.input_dir):
            sys.exit(f"ERROR: Directory not found: {args.input_dir}")
        ext = args.ext.lstrip(".")
        in_files = sorted(
            os.path.join(args.input_dir, f)
            for f in os.listdir(args.input_dir)
            if f.endswith(f".{ext}")
        )
        if not in_files:
            sys.exit(f"ERROR: No .{ext} files found in {args.input_dir}")

    print(f"Input : {len(in_files)} file(s)")
    print(f"Output: {out_subdir}/")
    if pattern:
        print(f"Divergence pattern: {args.divergence_pattern or '(from --divergence list)'}")
    else:
        print(f"Divergence detection: last {args.n_divergence} sequences")
    print()

    # Build job list
    jobs = [
        (f, os.path.join(out_subdir, os.path.basename(f)), pattern, args.n_divergence)
        for f in in_files
    ]

    # Execute
    reports: List[Dict] = []
    if args.threads > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            futures = {ex.submit(process_file, *j): j[0] for j in jobs}
            for fut in as_completed(futures):
                r = fut.result()
                reports.append(r)
                if args.verbose:
                    _print_status(r)
    else:
        for j in jobs:
            r = process_file(*j)
            reports.append(r)
            if args.verbose:
                _print_status(r)

    # Counters
    n_ok   = sum(1 for r in reports if not r["error"] and not r["warn"])
    n_warn = sum(1 for r in reports if r["warn"]  and not r["error"])
    n_err  = sum(1 for r in reports if r["error"])

    print(f"Done.  {n_ok} OK  |  {n_warn} warnings  |  {n_err} errors")

    for r in reports:
        if r["warn"]:
            print(f"  [WARN] {r['file']}: {r['warn']}")
        if r["error"]:
            print(f"  [ERR ] {r['file']}: {r['error']}")

    # Summary TSV
    summary_path = os.path.join(out_subdir, "trim_outgroup_summary.tsv")
    _write_summary(reports, summary_path)
    print(f"\nSummary: {summary_path}")


def _print_status(r: Dict) -> None:
    tag = "ERR " if r["error"] else ("WARN" if r["warn"] else "OK  ")
    print(
        f"  [{tag}] {r['file']}"
        f"  poly={r['n_poly']}"
        f"  div_in={r['n_div_in']}"
        f"  div_kept={r['n_div_kept']}"
        + (f"  removed: [{r['removed_first']}] … [{r['removed_last']}]"
           if r["removed_first"] else "")
    )


def _write_summary(reports: List[Dict], path: str) -> None:
    fields = [
        "file", "n_poly", "n_div_in", "n_div_kept",
        "removed_first", "removed_last", "warn", "error",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(reports)


if __name__ == "__main__":
    main()
