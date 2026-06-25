#!/usr/bin/env python3
"""
alignment_coverage.py
=====================
Assess and filter protein or nucleotide multiple-sequence alignments by
per-sequence completeness and per-column coverage.

Designed for the MKT_3_species pipeline: processes a directory of
per-orthogroup alignment files (e.g., from MAFFT or HmmCleaner).
Also accepts a single alignment file.

A gene alignment PASSES if ALL of the following hold:
  1. n_sequences       ≥ --min-sequences     (default: 2)
  2. alignment_length  ≥ --min-aln-length    (default: 30 aa / 90 nt)
  3. fraction_pass     ≥ --min-fraction      (default: 0.80)
     where fraction_pass = (sequences with completeness ≥ --min-completeness)
                           / n_sequences

Completeness of a sequence = (valid non-gap non-masked characters) /
                              (total alignment length)

For protein (aa): valid chars = 20 standard amino acids;
                  masked chars (HmmCleaner '?') are treated as gaps.
For nucleotide (nt): valid chars = A C G T.

USAGE — directory of alignments (primary mode):
  python alignment_coverage.py \\
      --input-dir   hmmcleaner_out/ \\
      --report      coverage_report.tsv \\
      [--gene-list  passing_genes.txt] \\
      [--out-dir    passing_alignments/] \\
      [--alphabet   aa|nt|auto] \\
      [--ext        fasta]

USAGE — single alignment file:
  python alignment_coverage.py \\
      --fasta    my_alignment.fasta \\
      --report   coverage_report.tsv
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── Character sets ────────────────────────────────────────────────────────────

NT_VALID   = frozenset("ACGTacgt")
AA_VALID   = frozenset("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")
# Characters that are protein-specific (never found in strict nucleotide data)
AA_ONLY    = frozenset("DEFHIKLMPQRSVWYdefhiklmpqrsvwy")
GAP_CHARS  = frozenset("-.")
MASK_CHARS = frozenset("?Xx")   # HmmCleaner '?', unknown AA 'X', unknown nt 'x'


# ── FASTA parser (no BioPython dependency) ────────────────────────────────────

def read_fasta(path: Path) -> List[Tuple[str, str]]:
    """Parse a FASTA file; returns list of (sequence_id, sequence) tuples."""
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    chunks: List[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].split()[0]   # first word = sequence ID
                chunks = []
            elif header is not None and line:
                chunks.append(line)
        if header is not None:
            records.append((header, "".join(chunks)))
    return records


# ── Alphabet detection ────────────────────────────────────────────────────────

def detect_alphabet(seqs: List[str], sample: int = 2000) -> str:
    """
    Heuristic alphabet detection from alignment sequences.
    Returns 'aa' if any protein-only letter is encountered, else 'nt'.
    """
    seen = 0
    for seq in seqs:
        for c in seq:
            if c in GAP_CHARS or c in MASK_CHARS:
                continue
            if c in AA_ONLY:
                return "aa"
            seen += 1
            if seen >= sample:
                break
        if seen >= sample:
            break
    return "nt"


# ── Per-sequence statistics ───────────────────────────────────────────────────

def sequence_completeness(seq: str, valid_set: frozenset) -> float:
    """Fraction of alignment positions that contain a valid (non-gap, non-masked) character."""
    if not seq:
        return 0.0
    n_valid = sum(1 for c in seq if c in valid_set)
    return n_valid / len(seq)


# ── Per-alignment assessment ──────────────────────────────────────────────────

def assess_alignment(
    records: List[Tuple[str, str]],
    alphabet: str,
    min_completeness: float,
    min_fraction: float,
    min_sequences: int,
    min_aln_length: int,
) -> Dict:
    """
    Assess a single alignment and return a statistics dict.

    Returns keys: n_sequences, alignment_length, mean_completeness,
    min_completeness, max_completeness, n_pass_sequences, fraction_pass,
    mean_col_coverage, PASS, FAIL_REASON.
    """
    valid_set = AA_VALID if alphabet == "aa" else NT_VALID
    n_seqs = len(records)
    aln_len = len(records[0][1]) if records else 0

    if n_seqs == 0 or aln_len == 0:
        return _empty_stats(fail="empty file")

    # Per-sequence completeness
    completeness = [sequence_completeness(seq, valid_set) for _, seq in records]
    n_pass = sum(1 for c in completeness if c >= min_completeness)
    frac_pass = n_pass / n_seqs
    mean_comp = sum(completeness) / n_seqs
    min_comp  = min(completeness)
    max_comp  = max(completeness)

    # Per-column coverage: fraction of sequences with a valid char per position
    col_valid = [0] * aln_len
    for _, seq in records:
        for i, c in enumerate(seq[:aln_len]):
            if c in valid_set:
                col_valid[i] += 1
    col_cov = [v / n_seqs for v in col_valid]
    mean_col_cov = sum(col_cov) / aln_len
    frac_cols = sum(1 for cv in col_cov if cv >= min_fraction) / aln_len

    # PASS / FAIL
    fails: List[str] = []
    if n_seqs < min_sequences:
        fails.append(f"n_sequences={n_seqs} < {min_sequences}")
    if aln_len < min_aln_length:
        fails.append(f"aln_length={aln_len} < {min_aln_length}")
    if frac_pass < min_fraction:
        fails.append(f"fraction_pass={frac_pass:.2f} < {min_fraction:.2f}")

    return {
        "n_sequences":       n_seqs,
        "alignment_length":  aln_len,
        "mean_completeness": round(mean_comp, 4),
        "min_completeness":  round(min_comp,  4),
        "max_completeness":  round(max_comp,  4),
        "n_pass_sequences":  n_pass,
        "fraction_pass":     round(frac_pass,     4),
        "mean_col_coverage": round(mean_col_cov,  4),
        "frac_cols_covered": round(frac_cols,     4),
        "PASS":              len(fails) == 0,
        "FAIL_REASON":       "; ".join(fails),
    }


def _empty_stats(fail: str = "empty") -> Dict:
    return {
        "n_sequences": 0, "alignment_length": 0,
        "mean_completeness": 0.0, "min_completeness": 0.0, "max_completeness": 0.0,
        "n_pass_sequences": 0, "fraction_pass": 0.0,
        "mean_col_coverage": 0.0, "frac_cols_covered": 0.0,
        "PASS": False, "FAIL_REASON": fail,
    }


# ── Worker (called in parallel) ───────────────────────────────────────────────

def _process_file(args_tuple) -> Dict:
    fpath, alphabet, min_comp, min_frac, min_seq, min_len = args_tuple
    try:
        records = read_fasta(fpath)
        alph = detect_alphabet([s for _, s in records]) if alphabet == "auto" else alphabet
        stats = assess_alignment(records, alph, min_comp, min_frac, min_seq, min_len)
        stats["gene"]     = fpath.stem
        stats["file"]     = str(fpath)
        stats["alphabet"] = alph
        return stats
    except Exception as exc:
        return {**_empty_stats(fail=f"ERROR: {exc}"),
                "gene": fpath.stem, "file": str(fpath), "alphabet": "?"}


# ── Summary printing ──────────────────────────────────────────────────────────

def _pct(n: int, total: int) -> str:
    return f"{100 * n / total:.1f}%" if total else "N/A"


def _median(vals: List[float]) -> float:
    if not vals:
        return float("nan")
    s = sorted(vals)
    m = len(s) // 2
    return (s[m] + s[m - 1]) / 2 if len(s) % 2 == 0 else s[m]


def print_assessment(rows: List[Dict], args: argparse.Namespace, n_files: int) -> None:
    """Print a human-readable coverage assessment to stdout."""
    total = len(rows)
    n_pass = sum(1 for r in rows if r["PASS"])
    n_fail = total - n_pass

    print()
    print("=" * 62)
    print("  Alignment Coverage Assessment")
    print("=" * 62)
    print(f"  Input        : {args.input_dir or args.fasta}  ({n_files} file(s))")
    alphabets = Counter(r["alphabet"] for r in rows)
    print(f"  Alphabet     : {dict(alphabets)}")
    print(f"  Thresholds   :")
    print(f"    --min-completeness  {args.min_completeness}  "
          f"(non-gap valid chars per sequence)")
    print(f"    --min-fraction      {args.min_fraction}  "
          f"(fraction of sequences passing completeness)")
    print(f"    --min-sequences     {args.min_sequences}  "
          f"(minimum sequences per alignment)")
    print(f"    --min-aln-length    {args.min_aln_length}  "
          f"(minimum alignment length in residues)")
    print("=" * 62)
    print()
    print(f"  PASS  : {n_pass:5d} / {total}  ({_pct(n_pass, total)})")
    print(f"  FAIL  : {n_fail:5d} / {total}  ({_pct(n_fail, total)})")
    print()

    # Coverage distribution
    frac_pass_vals = [r["fraction_pass"] for r in rows]
    bins = [(0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0), (1.0, 1.001)]
    labels = ["0.0 – 0.2", "0.2 – 0.4", "0.4 – 0.6", "0.6 – 0.8", "0.8 – 1.0", "= 1.00"]
    print("  Fraction-pass distribution  (sequences meeting --min-completeness):")
    for (lo, hi), label in zip(bins, labels):
        cnt = sum(1 for v in frac_pass_vals if lo <= v < hi)
        bar = "#" * min(30, int(cnt / max(total, 1) * 30 * 5))
        marker = "  ← passing" if lo >= args.min_fraction else ""
        print(f"    {label} : {cnt:5d} ({_pct(cnt, total):5s})  {bar}{marker}")
    print()

    # Mean completeness statistics
    comp_vals = [r["mean_completeness"] for r in rows if r["n_sequences"] > 0]
    if comp_vals:
        print("  Per-sequence completeness  (mean across sequences per gene):")
        print(f"    Median : {_median(comp_vals):.3f}")
        print(f"    Mean   : {sum(comp_vals)/len(comp_vals):.3f}")
        print(f"    Min    : {min(comp_vals):.3f}")
        print(f"    Max    : {max(comp_vals):.3f}")
        print()

    # Alignment length statistics
    lens = [r["alignment_length"] for r in rows if r["alignment_length"] > 0]
    if lens:
        alph_label = "nt" if alphabets.get("nt", 0) > alphabets.get("aa", 0) else "aa"
        print(f"  Alignment length  ({alph_label}):")
        print(f"    Median : {int(_median(lens))}")
        print(f"    Mean   : {int(sum(lens)/len(lens))}")
        print(f"    Min    : {min(lens)}")
        print(f"    Max    : {max(lens)}")
        print()

    # FAIL reason breakdown
    fail_rows = [r for r in rows if not r["PASS"] and r["FAIL_REASON"]]
    if fail_rows:
        reason_counts: Dict[str, int] = {}
        for r in fail_rows:
            for reason in r["FAIL_REASON"].split("; "):
                key = reason.split("=")[0].replace(" ", "_")
                reason_counts[key] = reason_counts.get(key, 0) + 1
        print("  FAIL reasons:")
        for reason, cnt in sorted(reason_counts.items(), key=lambda x: -x[1]):
            print(f"    {reason:<30s}: {cnt}")
        print()

    # Recommendations
    pct_pass = 100 * n_pass / total if total else 0
    print("  Recommendation:")
    if pct_pass >= 80:
        print(f"    ✓ {pct_pass:.1f}% of alignments passed — good coverage overall.")
    elif pct_pass >= 50:
        print(f"    ⚠ {pct_pass:.1f}% of alignments passed — moderate coverage.")
        print("      Consider relaxing --min-completeness or --min-fraction.")
    else:
        print(f"    ✗ {pct_pass:.1f}% of alignments passed — low coverage.")
        print("      Check alignment quality; consider lower thresholds or")
        print("      different HmmCleaner settings.")
    if args.out_dir is None:
        print(f"    → Use --out-dir to copy passing alignments for downstream analysis.")
    else:
        print(f"    ✓ Passing alignments copied to:  {args.out_dir}")
    print("=" * 62)
    print()


# ── CLI ───────────────────────────────────────────────────────────────────────

REPORT_FIELDS = [
    "gene", "file", "alphabet",
    "n_sequences", "alignment_length",
    "mean_completeness", "min_completeness", "max_completeness",
    "n_pass_sequences", "fraction_pass",
    "mean_col_coverage", "frac_cols_covered",
    "PASS", "FAIL_REASON",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Input — one of --input-dir or --fasta
    inp = p.add_mutually_exclusive_group(required=True)
    inp.add_argument(
        "-i", "--input-dir", metavar="DIR",
        help="Directory of alignment files (one per gene/orthogroup). "
             "Primary mode for the MKT_3_species pipeline.",
    )
    inp.add_argument(
        "--fasta", metavar="FILE",
        help="Single alignment FASTA file (alternative to --input-dir).",
    )

    # Output
    p.add_argument(
        "-r", "--report", required=True, metavar="FILE",
        help="Output TSV report with per-gene statistics and PASS/FAIL verdict.",
    )
    p.add_argument(
        "-g", "--gene-list", metavar="FILE",
        help="Output plain-text list of passing gene IDs (one per line).",
    )
    p.add_argument(
        "-o", "--out-dir", metavar="DIR",
        help="If provided, copy passing alignment files to this directory.",
    )
    p.add_argument(
        "--filter-sequences", action="store_true",
        help="When copying to --out-dir, remove sequences that fail "
             "--min-completeness (write only passing sequences).",
    )

    # Alphabet
    p.add_argument(
        "-a", "--alphabet", choices=["auto", "nt", "aa"], default="auto",
        help="Sequence alphabet: 'aa' (protein), 'nt' (nucleotide), "
             "or 'auto' to detect from each file. Default: auto.",
    )

    # File selection
    p.add_argument(
        "-e", "--ext", default=None, metavar="EXT",
        help="File extension filter when using --input-dir "
             "(e.g. 'fasta', 'fa'). If omitted, accepts .fa, .fasta, .faa.",
    )

    # Thresholds
    p.add_argument(
        "--min-completeness", type=float, default=0.50, metavar="FLOAT",
        help="Minimum fraction of valid (non-gap, non-masked) characters "
             "per sequence [0–1]. Default: 0.50.",
    )
    p.add_argument(
        "--min-fraction", type=float, default=0.80, metavar="FLOAT",
        help="Minimum fraction of sequences per alignment that must meet "
             "--min-completeness [0–1]. Default: 0.80.",
    )
    p.add_argument(
        "--min-sequences", type=int, default=2, metavar="N",
        help="Minimum number of sequences per alignment. Default: 2.",
    )
    p.add_argument(
        "--min-aln-length", type=int, default=30, metavar="N",
        help="Minimum alignment length in residues (aa or nt). Default: 30.",
    )

    # Performance
    p.add_argument(
        "-t", "--threads", type=int, default=1, metavar="N",
        help="Number of parallel worker processes. Default: 1.",
    )

    return p.parse_args()


# ── Main ──────────────────────────────────────────────────────────────────────

def collect_files(args: argparse.Namespace) -> List[Path]:
    """Return list of alignment files to process."""
    if args.fasta:
        p = Path(args.fasta)
        if not p.is_file():
            sys.exit(f"ERROR: File not found: {p}")
        return [p]

    d = Path(args.input_dir)
    if not d.is_dir():
        sys.exit(f"ERROR: Directory not found: {d}")

    if args.ext:
        exts = [args.ext.lstrip(".")]
    else:
        exts = ["fa", "fasta", "faa"]

    files: List[Path] = []
    for ext in exts:
        files.extend(sorted(d.glob(f"*.{ext}")))

    if not files:
        sys.exit(
            f"ERROR: No alignment files found in {d} "
            f"(tried extensions: {', '.join(exts)}). "
            "Use --ext to specify a different extension."
        )
    return files


def write_filtered_fasta(
    fpath: Path,
    out_path: Path,
    alphabet: str,
    min_completeness: float,
) -> None:
    """Write FASTA with only sequences passing the completeness threshold."""
    valid_set = AA_VALID if alphabet == "aa" else NT_VALID
    records = read_fasta(fpath)
    with open(out_path, "w", encoding="utf-8") as fh:
        for header, seq in records:
            comp = sequence_completeness(seq, valid_set)
            if comp >= min_completeness:
                fh.write(f">{header}\n{seq}\n")


def main() -> None:
    args = parse_args()

    # Validate thresholds
    for name, val in [
        ("--min-completeness", args.min_completeness),
        ("--min-fraction",     args.min_fraction),
    ]:
        if not 0.0 <= val <= 1.0:
            sys.exit(f"ERROR: {name} must be between 0 and 1, got {val}")

    files = collect_files(args)
    n_files = len(files)

    print(f"Found {n_files} alignment file(s) to assess.")

    # Build task list
    task_args = [
        (f, args.alphabet, args.min_completeness, args.min_fraction,
         args.min_sequences, args.min_aln_length)
        for f in files
    ]

    # Process (parallel or sequential)
    rows: List[Dict] = []
    if args.threads > 1 and n_files > 1:
        print(f"Processing with {args.threads} threads...")
        with ProcessPoolExecutor(max_workers=args.threads) as pool:
            futures = {pool.submit(_process_file, t): t[0] for t in task_args}
            done = 0
            for fut in as_completed(futures):
                rows.append(fut.result())
                done += 1
                if done % max(1, n_files // 10) == 0:
                    print(f"  {done}/{n_files} processed...")
    else:
        for i, t in enumerate(task_args):
            rows.append(_process_file(t))
            if n_files > 20 and (i + 1) % max(1, n_files // 10) == 0:
                print(f"  {i+1}/{n_files} processed...")

    # Sort by gene name
    rows.sort(key=lambda r: r["gene"])

    # Write TSV report
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=REPORT_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Report written to: {args.report}")

    # Write gene list
    passing = [r for r in rows if r["PASS"]]
    if args.gene_list:
        Path(args.gene_list).parent.mkdir(parents=True, exist_ok=True)
        with open(args.gene_list, "w", encoding="utf-8") as fh:
            for r in passing:
                fh.write(r["gene"] + "\n")
        print(f"Gene list written to: {args.gene_list}  ({len(passing)} genes)")

    # Copy passing alignments to out_dir
    if args.out_dir:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        copied = 0
        for r in passing:
            src = Path(r["file"])
            dst = out_dir / src.name
            if args.filter_sequences:
                write_filtered_fasta(src, dst, r["alphabet"], args.min_completeness)
            else:
                shutil.copy2(src, dst)
            copied += 1
        print(f"Copied {copied} passing alignments to: {args.out_dir}")

    # Print human-readable assessment
    print_assessment(rows, args, n_files)


if __name__ == "__main__":
    main()
