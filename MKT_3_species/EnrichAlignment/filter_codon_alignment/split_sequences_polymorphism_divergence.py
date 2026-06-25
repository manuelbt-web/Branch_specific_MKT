#!/usr/bin/env python3
"""
split_sequences_polymorphism_divergence.py
==========================================
Split combined divergence-polymorphism FASTA alignments (produced by
filter_codon_alignment.py) into separate per-category files.

Output directories are created automatically inside <output-dir>:
  <output-dir>/
  ├── divergence/     — one FASTA per gene (divergence sequences only)
  └── polymorphism/   — one FASTA per gene (polymorphism sequences only)

DIVERGENCE SEQUENCE IDENTIFICATION (choose one — checked in priority order)
  --divergence-pattern REGEX   match headers by regular expression  [recommended]
  --divergence         LIST    comma-separated exact sequence headers
  --n-divergence       N       use the last N sequences (default: 3)

The last-N default works because MACSE enrichAlignment appends divergence
sequences after the fixed polymorphism block. Use --divergence-pattern when
orthogroups have variable numbers of divergence species.

USAGE — minimal (single file):
  python split_sequences_polymorphism_divergence.py \\
      -i  filtered_aln/HOG0000001_NT_filtered.fasta \\
      -o  split_out/

USAGE — directory of files (recommended):
  python split_sequences_polymorphism_divergence.py \\
      --input-dir   filtered_aln/ \\
      --output-dir  split_out/

USAGE — with explicit divergence pattern:
  python split_sequences_polymorphism_divergence.py \\
      --input-dir            filtered_aln/ \\
      --output-dir           split_out/ \\
      --divergence-pattern   "Aegilops_tauschii|Triticum_urartu|Amblyopyrum_muticum"

OUTPUTS (in <output-dir>/):
  divergence/<gene>_divergence.fasta     — divergence sequences
  polymorphism/<gene>_polymorphism.fasta — polymorphism sequences
  split_summary.tsv                      — per-gene counts (directory mode)
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ── FASTA IO ──────────────────────────────────────────────────────────────────

def read_fasta_ordered(path: Path) -> List[Tuple[str, str]]:
    """Parse FASTA; return list of (header, sequence) in file order."""
    seqs: List[Tuple[str, str]] = []
    header: Optional[str] = None
    parts: List[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(parts)))
                header = line[1:]
                parts = []
            else:
                parts.append(line.strip())
        if header is not None:
            seqs.append((header, "".join(parts)))
    return seqs


def write_fasta(path: Path, entries: List[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for h, s in entries:
            fh.write(f">{h}\n{s}\n")


# ── Divergence/polymorphism split ─────────────────────────────────────────────

def split_sequences(
    seqs:         List[Tuple[str, str]],
    div_list:     Optional[List[str]],
    div_pattern:  Optional["re.Pattern"],
    n_divergence: int,
) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]], str]:
    """
    Split (header, seq) pairs into (divergence_seqs, polymorphism_seqs, warning).
    warning is an empty string when there are no issues.
    """
    warn = ""
    if div_list:
        div_set   = set(div_list)
        div_seqs  = [(h, s) for h, s in seqs if h in div_set]
        poly_seqs = [(h, s) for h, s in seqs if h not in div_set]
        missing   = [h for h in div_list if h not in {h for h, _ in div_seqs}]
        if missing:
            warn = f"{len(missing)} divergence header(s) not found in this file"

    elif div_pattern:
        div_seqs  = [(h, s) for h, s in seqs if div_pattern.search(h)]
        poly_seqs = [(h, s) for h, s in seqs if not div_pattern.search(h)]
        if not div_seqs:
            warn = "pattern matched no sequences — check --divergence-pattern"

    else:
        n = min(n_divergence, len(seqs))
        if n < n_divergence:
            warn = (f"only {len(seqs)} sequence(s) in file; "
                    f"used all {n} as divergence (fewer than --n-divergence {n_divergence})")
        div_seqs  = seqs[-n:]
        poly_seqs = seqs[:-n]

    return div_seqs, poly_seqs, warn


# ── Output filename generation ────────────────────────────────────────────────

def make_base(stem: str, strip_suffix: str) -> str:
    """Strip the trailing suffix from the stem (e.g. '_NT_filtered')."""
    if strip_suffix and stem.endswith(strip_suffix):
        return stem[: -len(strip_suffix)]
    return stem


# ── Summary helpers ───────────────────────────────────────────────────────────

def _mean(vals: List[int]) -> float:
    return sum(vals) / len(vals) if vals else 0.0


def print_summary(
    results: List[Dict],
    args:    argparse.Namespace,
    out_dir: Path,
) -> None:
    ok   = [r for r in results if not r["error"]]
    warn = [r for r in results if r["warn"]]
    err  = [r for r in results if r["error"]]

    div_spec = (
        f"pattern '{args.divergence_pattern}'" if args.divergence_pattern
        else f"explicit list ({len(args.divergence.split(','))} headers)" if args.divergence
        else f"last {args.n_divergence} sequences"
    )

    print()
    print("=" * 60)
    print("  Sequence Split — Summary")
    print("=" * 60)
    print(f"  Files processed    : {len(results)}")
    if ok:
        print(f"  Files split        : {len(ok)}")
    if warn:
        print(f"  Warnings           : {len(warn)}")
    if err:
        print(f"  Errors             : {len(err)}")
    print(f"  Divergence spec    : {div_spec}")
    if ok:
        div_counts  = [r["n_div"]  for r in ok]
        poly_counts = [r["n_poly"] for r in ok]
        n_empty_poly = sum(1 for c in poly_counts if c == 0)
        print(f"  Mean divergence seqs   : {_mean(div_counts):.1f}")
        print(f"  Mean polymorphism seqs : {_mean(poly_counts):.1f}")
        if n_empty_poly:
            print(f"  [WARN] {n_empty_poly} gene(s) have 0 polymorphism sequences")
    print()
    print(f"  divergence/   : {out_dir / 'divergence'}/")
    print(f"  polymorphism/ : {out_dir / 'polymorphism'}/")
    if len(results) > 1:
        print(f"  Summary table : {out_dir / 'split_summary.tsv'}")
    print("=" * 60)
    print()


def write_summary(results: List[Dict], out_dir: Path) -> None:
    fields = ["gene", "n_div", "n_poly", "n_total", "warn", "error"]
    with open(out_dir / "split_summary.tsv", "w", encoding="utf-8") as fh:
        fh.write("\t".join(fields) + "\n")
        for r in sorted(results, key=lambda x: x["gene"]):
            fh.write("\t".join(str(r.get(f, "")) for f in fields) + "\n")


# ── Core processing (one file) ────────────────────────────────────────────────

def process_file(
    fpath:        Path,
    div_out:      Path,
    poly_out:     Path,
    div_list:     Optional[List[str]],
    div_pattern:  Optional["re.Pattern"],
    n_divergence: int,
    strip_suffix: str,
    overwrite:    bool,
    verbose:      bool,
) -> Dict:
    stem = make_base(fpath.stem, strip_suffix)
    d_path = div_out  / f"{stem}_divergence.fasta"
    p_path = poly_out / f"{stem}_polymorphism.fasta"

    if not overwrite and (d_path.exists() or p_path.exists()):
        if verbose:
            print(f"[SKIP] {fpath.name} — outputs already exist")
        return {"gene": stem, "n_div": 0, "n_poly": 0, "n_total": 0,
                "warn": "", "error": "skipped (outputs exist)"}

    seqs = read_fasta_ordered(fpath)
    if not seqs:
        print(f"[WARN] {fpath.name}: empty file — skipped", file=sys.stderr)
        return {"gene": stem, "n_div": 0, "n_poly": 0, "n_total": 0,
                "warn": "empty file", "error": "empty"}

    div_seqs, poly_seqs, warn = split_sequences(seqs, div_list, div_pattern, n_divergence)

    if warn:
        print(f"[WARN] {fpath.name}: {warn}", file=sys.stderr)

    write_fasta(d_path, div_seqs)
    write_fasta(p_path, poly_seqs)

    msg = (f"[OK  ] {fpath.name}  →  "
           f"div={len(div_seqs)}, poly={len(poly_seqs)}")
    if warn:
        msg += f"  [!] {warn}"
    if verbose or warn:
        print(msg)
    elif not verbose:
        pass  # silent in non-verbose mode (summary printed at end)

    return {
        "gene":    stem,
        "n_div":   len(div_seqs),
        "n_poly":  len(poly_seqs),
        "n_total": len(seqs),
        "warn":    warn,
        "error":   "",
    }


# ── Argument parsing ──────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Input — mutually exclusive: single file OR directory
    inp = p.add_mutually_exclusive_group(required=True)
    inp.add_argument(
        "-i", "--input", metavar="FILE",
        help="Single input FASTA file.",
    )
    inp.add_argument(
        "--input-dir", metavar="DIR",
        help="Directory of filtered FASTA files to process in batch.",
    )

    # Output — single path; subdirectories created automatically
    p.add_argument(
        "-o", "--output-dir", required=True, metavar="DIR",
        help=("Output directory. Two subdirectories are created automatically: "
              "divergence/ and polymorphism/."),
    )

    # Divergence identification — mutually exclusive; last-N is the fallback
    div = p.add_mutually_exclusive_group()
    div.add_argument(
        "--divergence-pattern", metavar="REGEX",
        help=("Regular expression matched against sequence headers to identify "
              "divergence sequences (recommended). "
              "Example: \"Aegilops_tauschii|Triticum_urartu\""),
    )
    div.add_argument(
        "--divergence", metavar="LIST",
        help="Comma-separated list of exact divergence sequence headers.",
    )
    p.add_argument(
        "--n-divergence", type=int, default=3, metavar="N",
        help=("When neither --divergence-pattern nor --divergence is given, "
              "treat the last N sequences as divergence (default: 3)."),
    )

    # File handling
    p.add_argument(
        "-e", "--ext", default=None, metavar="EXT",
        help="File extension filter for --input-dir (e.g. 'fasta'). "
             "Default: .fa, .fasta.",
    )
    p.add_argument(
        "--strip-suffix", default="_NT_filtered", metavar="SUFFIX",
        help=("Suffix to strip from the input filename stem when generating "
              "output names. Default: '_NT_filtered'. "
              "Pass '' to disable stripping."),
    )
    p.add_argument(
        "--overwrite", action="store_true",
        help="Overwrite existing output files (default: skip).",
    )
    p.add_argument(
        "--verbose", action="store_true",
        help="Print one line per gene.",
    )

    return p.parse_args()


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    out_dir  = Path(args.output_dir)
    div_out  = out_dir / "divergence"
    poly_out = out_dir / "polymorphism"
    div_out.mkdir(parents=True, exist_ok=True)
    poly_out.mkdir(parents=True, exist_ok=True)

    # Parse divergence options once
    div_list: Optional[List[str]] = (
        [h.strip() for h in args.divergence.split(",") if h.strip()]
        if args.divergence else None
    )
    div_pattern: Optional["re.Pattern"] = None
    if args.divergence_pattern:
        try:
            div_pattern = re.compile(args.divergence_pattern)
        except re.error as e:
            sys.exit(f"ERROR: Invalid --divergence-pattern regex: {e}")

    # ── Single-file mode ──────────────────────────────────────────────────────
    if args.input:
        fpath = Path(args.input)
        if not fpath.is_file():
            sys.exit(f"ERROR: Input file not found: {fpath}")

        div_spec = (
            f"pattern '{args.divergence_pattern}'" if args.divergence_pattern
            else f"explicit list" if div_list
            else f"last {args.n_divergence} sequences"
        )
        print(f"Input        : {fpath}")
        print(f"Divergence   : {div_spec}")
        print(f"Output dir   : {out_dir}/")
        print(f"  divergence/   and  polymorphism/  created")
        print()

        result = process_file(
            fpath, div_out, poly_out, div_list, div_pattern,
            args.n_divergence, args.strip_suffix, args.overwrite, True,
        )

        stem   = result["gene"]
        n_div  = result["n_div"]
        n_poly = result["n_poly"]

        if not result["error"]:
            print()
            print(f"Divergence seqs   : {n_div}  → {div_out / (stem + '_divergence.fasta')}")
            print(f"Polymorphism seqs : {n_poly}  → {poly_out / (stem + '_polymorphism.fasta')}")
            if n_poly == 0:
                print("[WARN] No polymorphism sequences — polymorphism FASTA is empty.")
        return

    # ── Directory mode ────────────────────────────────────────────────────────
    in_dir = Path(args.input_dir)
    if not in_dir.is_dir():
        sys.exit(f"ERROR: Directory not found: {in_dir}")

    exts = [args.ext.lstrip(".")] if args.ext else ["fa", "fasta"]
    files: List[Path] = []
    for ext in exts:
        files.extend(sorted(in_dir.glob(f"*.{ext}")))
    files = sorted(set(files), key=lambda f: f.name)

    if not files:
        sys.exit(
            f"ERROR: No FASTA files found in {in_dir} "
            f"(tried: {', '.join('.' + e for e in exts)}). "
            "Use --ext to specify the file extension."
        )

    div_spec = (
        f"pattern '{args.divergence_pattern}'" if args.divergence_pattern
        else f"explicit list ({len(div_list)} sequences)" if div_list
        else f"last {args.n_divergence} sequences"
    )
    print(f"Found {len(files)} file(s) in {in_dir}")
    print(f"Divergence     : {div_spec}")
    print(f"Output         : {out_dir}/  (divergence/ and polymorphism/ subdirectories)")
    print()

    results: List[Dict] = []
    n_total = len(files)
    for i, fpath in enumerate(files):
        result = process_file(
            fpath, div_out, poly_out, div_list, div_pattern,
            args.n_divergence, args.strip_suffix, args.overwrite, args.verbose,
        )
        results.append(result)
        if n_total > 20 and (i + 1) % max(1, n_total // 10) == 0:
            print(f"  {i+1}/{n_total} processed…")

    write_summary(results, out_dir)
    print_summary(results, args, out_dir)


if __name__ == "__main__":
    main()
