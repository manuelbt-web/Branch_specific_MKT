#!/usr/bin/env python3
"""
filter_codon_alignment.py
=========================
Filter per-codon sites from combined divergence-polymorphism alignments
produced by MACSE enrichAlignment, keeping only codons informative for
MKT or codeml (dN/dS) analysis.

FILTERING CRITERIA (applied to every codon position, i.e. every 3-nt column)
  1. Divergence — strict: ALL divergence sequences must have a fully valid
     codon (A/C/G/T only; no gaps, no ambiguous bases, no MACSE '!' frameshift)
  2. Polymorphism — minimum depth: at least --min-polym sequences must have
     a valid codon (default: 6)

Codons that fail either criterion are replaced with NNN (default) so that
codeml can read the alignment with N treated as missing data, preserving the
codon frame and alignment column structure.

DIVERGENCE SEQUENCE IDENTIFICATION (choose one — options are checked in order)
  --divergence-pattern REGEX   match headers by regular expression  [recommended]
  --divergence         LIST    comma-separated exact sequence headers
  --n-divergence       N       use the last N sequences (default: 3)

USAGE — single file (minimal):
  python filter_codon_alignment.py \\
      -i  combined_NT.fasta \\
      -o  filtered_out/

USAGE — directory of files:
  python filter_codon_alignment.py \\
      --input-dir   macse_enriched/enriched/ \\
      --output-dir  filtered_out/

USAGE — with explicit divergence pattern (recommended for robustness):
  python filter_codon_alignment.py \\
      --input-dir            macse_enriched/enriched/ \\
      --output-dir           filtered_out/ \\
      --divergence-pattern   "Aegilops_tauschii|Triticum_urartu|Amblyopyrum_muticum"

OUTPUTS (per gene, written in <output-dir>/):
  <gene>_filtered.fasta       — alignment with failing codons replaced by NNN
  <gene>_codon_pass.tsv       — per-codon report with reasons
  <gene>_pass_positions.txt   — passing codon indices (1-based)
  filter_summary.tsv          — across-gene summary table (directory mode)
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── Character sets ────────────────────────────────────────────────────────────

VALID_BASES = frozenset("ACGTacgt")
GAP_CHARS   = frozenset("-.")
AMBIG_CHARS = frozenset("NRYSWKMBDHVnryswkmbdhv!?")  # includes MACSE frameshifts


# ── FASTA IO ──────────────────────────────────────────────────────────────────

def read_fasta(path: Path) -> OrderedDict:
    """Parse a FASTA file; preserves insertion order."""
    seqs: OrderedDict = OrderedDict()
    header: Optional[str] = None
    parts: List[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(parts)
                header = line[1:]
                parts = []
            else:
                parts.append(line.strip())
        if header is not None:
            seqs[header] = "".join(parts)
    return seqs


def write_fasta(seqs: Dict[str, str], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for header, seq in seqs.items():
            fh.write(f">{header}\n{seq}\n")


# ── Codon validation ──────────────────────────────────────────────────────────

def is_valid_codon(codon: str) -> bool:
    """Return True if the codon contains only unambiguous ACGT bases (no gaps, no N)."""
    if len(codon) != 3:
        return False
    return all(c in VALID_BASES for c in codon)


# ── Divergence-header resolution ──────────────────────────────────────────────

def resolve_divergence(
    headers: List[str],
    div_list:    Optional[List[str]],
    div_pattern: Optional[str],
    n_divergence: int,
) -> Tuple[List[str], List[str]]:
    """
    Return (divergence_headers, polymorphism_headers).
    Priority: explicit list > regex pattern > last-N heuristic.
    """
    if div_list:
        missing = [h for h in div_list if h not in headers]
        if missing:
            sys.exit(
                f"ERROR: The following --divergence headers were not found in the alignment:\n"
                + "\n".join(f"  '{h}'" for h in missing)
                + f"\nAvailable headers:\n"
                + "\n".join(f"  '{h}'" for h in headers)
            )
        div = div_list
    elif div_pattern:
        try:
            pat = re.compile(div_pattern)
        except re.error as e:
            sys.exit(f"ERROR: Invalid --divergence-pattern regex: {e}")
        div = [h for h in headers if pat.search(h)]
        if not div:
            sys.exit(
                f"ERROR: --divergence-pattern '{div_pattern}' matched no headers.\n"
                f"Available headers:\n" + "\n".join(f"  '{h}'" for h in headers)
            )
    else:
        if len(headers) < n_divergence:
            sys.exit(
                f"ERROR: Alignment has only {len(headers)} sequences but "
                f"--n-divergence is {n_divergence}."
            )
        div = headers[-n_divergence:]

    poly = [h for h in headers if h not in set(div)]
    return div, poly


# ── Core filter ───────────────────────────────────────────────────────────────

def filter_alignment(
    seqs:        OrderedDict,
    div_headers: List[str],
    poly_headers: List[str],
    min_polym:   int,
    replace_char: str,
) -> Tuple[Dict, List[Tuple], List[int]]:
    """
    Apply the two-criterion codon filter.

    Returns:
      out_seqs    : dict header → filtered sequence string
      report_rows : list of (codon_idx, nt_start, passed, reason, n_valid_poly, n_poly_total)
      pass_indices: 1-based list of passing codon indices
    """
    headers  = list(seqs.keys())
    aln_len  = len(next(iter(seqs.values())))
    n_codons = aln_len // 3
    repl     = replace_char * 3

    out_seqs    = {h: [] for h in headers}
    report_rows: List[Tuple] = []
    pass_indices: List[int] = []

    for ci in range(n_codons):
        s, e = ci * 3, ci * 3 + 3

        div_codons  = [seqs[h][s:e] for h in div_headers]
        poly_codons = [seqs[h][s:e] for h in poly_headers]

        if not all(is_valid_codon(c) for c in div_codons):
            passed = False
            reason = "divergence_invalid"
        else:
            n_valid = sum(1 for c in poly_codons if is_valid_codon(c))
            if n_valid >= min_polym:
                passed = True
                reason = "ok"
            else:
                passed = False
                reason = f"poly_depth={n_valid}<{min_polym}"

        if passed:
            pass_indices.append(ci + 1)
            for h in headers:
                out_seqs[h].append(seqs[h][s:e])
        else:
            for h in headers:
                out_seqs[h].append(repl)

        n_valid_poly = sum(1 for c in poly_codons if is_valid_codon(c))
        report_rows.append((ci + 1, s + 1, int(passed), reason,
                            n_valid_poly, len(poly_codons)))

    return {h: "".join(v) for h, v in out_seqs.items()}, report_rows, pass_indices


# ── Write outputs ─────────────────────────────────────────────────────────────

def write_outputs(
    stem:         str,
    out_dir:      Path,
    out_seqs:     Dict[str, str],
    report_rows:  List[Tuple],
    pass_indices: List[int],
) -> None:
    write_fasta(out_seqs, out_dir / f"{stem}_filtered.fasta")

    with open(out_dir / f"{stem}_codon_pass.tsv", "w", encoding="utf-8") as tf:
        tf.write("codon_index\tnt_start\tpass\treason\tn_valid_polym\tn_polym_total\n")
        for row in report_rows:
            tf.write("\t".join(str(x) for x in row) + "\n")

    with open(out_dir / f"{stem}_pass_positions.txt", "w", encoding="utf-8") as pf:
        pf.write("\n".join(str(i) for i in pass_indices) + "\n")


# ── Worker (for parallel directory mode) ──────────────────────────────────────

def _process_file(args_tuple) -> Dict:
    fpath, out_dir, div_list, div_pattern, n_divergence, min_polym, replace_char = args_tuple
    stem = fpath.stem
    try:
        seqs = read_fasta(fpath)
        if not seqs:
            return {"gene": stem, "error": "empty file",
                    "n_codons": 0, "n_pass": 0, "n_fail_div": 0, "n_fail_poly": 0,
                    "frac_pass": 0.0, "n_div": 0, "n_poly": 0}

        headers = list(seqs.keys())
        aln_len = len(next(iter(seqs.values())))

        if aln_len % 3 != 0:
            print(f"[WARN] {stem}: alignment length {aln_len} is not a multiple of 3 "
                  f"— last {aln_len % 3} nt ignored.", file=sys.stderr)

        div_hdrs, poly_hdrs = resolve_divergence(
            headers, div_list, div_pattern, n_divergence
        )
        out_seqs, report_rows, pass_indices = filter_alignment(
            seqs, div_hdrs, poly_hdrs, min_polym, replace_char
        )
        write_outputs(stem, out_dir, out_seqs, report_rows, pass_indices)

        n_codons   = len(report_rows)
        n_pass     = len(pass_indices)
        n_fail_div  = sum(1 for r in report_rows if "divergence" in r[3])
        n_fail_poly = sum(1 for r in report_rows if "poly_depth" in r[3])
        frac_pass  = n_pass / n_codons if n_codons else 0.0
        return {
            "gene": stem, "error": "",
            "n_codons": n_codons, "n_pass": n_pass,
            "n_fail_div": n_fail_div, "n_fail_poly": n_fail_poly,
            "frac_pass": round(frac_pass, 4),
            "n_div": len(div_hdrs), "n_poly": len(poly_hdrs),
        }
    except SystemExit:
        raise
    except Exception as exc:
        return {"gene": stem, "error": str(exc),
                "n_codons": 0, "n_pass": 0, "n_fail_div": 0, "n_fail_poly": 0,
                "frac_pass": 0.0, "n_div": 0, "n_poly": 0}


# ── Summary output ────────────────────────────────────────────────────────────

def _median(vals: List[float]) -> float:
    if not vals:
        return float("nan")
    s = sorted(vals)
    m = len(s) // 2
    return (s[m] + s[m - 1]) / 2 if len(s) % 2 == 0 else s[m]


def print_summary(results: List[Dict], args: argparse.Namespace, out_dir: Path) -> None:
    ok = [r for r in results if not r["error"]]
    err = [r for r in results if r["error"]]
    frac_vals = [r["frac_pass"] for r in ok]

    div_spec = (
        f"pattern '{args.divergence_pattern}'" if args.divergence_pattern
        else f"headers: {args.divergence}" if args.divergence
        else f"last {args.n_divergence} sequences"
    )
    replace_label = "NNN" if args.replace_with == "N" else "---"

    print()
    print("=" * 60)
    print("  Codon Alignment Filter — Summary")
    print("=" * 60)
    print(f"  Genes processed      : {len(results)}")
    if err:
        print(f"  Errors               : {len(err)}")
    print(f"  Replace failing with : {replace_label}")
    print(f"  Divergence           : {div_spec}")
    print(f"  Min polymorphism     : {args.min_polym} valid codons per site")
    print()
    if ok:
        mean_frac = sum(frac_vals) / len(frac_vals)
        n_gt50  = sum(1 for v in frac_vals if v >= 0.50)
        n_lt10  = sum(1 for v in frac_vals if v < 0.10)
        print(f"  Fraction of codons retained (per gene):")
        print(f"    Mean   : {mean_frac:.3f}")
        print(f"    Median : {_median(frac_vals):.3f}")
        print(f"    Min    : {min(frac_vals):.3f}")
        print(f"    Max    : {max(frac_vals):.3f}")
        print()
        print(f"  Genes with ≥ 50% codons retained : {n_gt50}/{len(ok)}")
        print(f"  Genes with <  10% codons retained :  {n_lt10}/{len(ok)}")
    print()
    print(f"  Output directory     : {out_dir}/")
    if len(results) > 1:
        print(f"  Summary table        : {out_dir / 'filter_summary.tsv'}")
    print("=" * 60)
    print()


def write_dir_summary(results: List[Dict], out_dir: Path) -> None:
    fields = ["gene", "n_codons", "n_pass", "n_fail_div", "n_fail_poly",
              "frac_pass", "n_div", "n_poly", "error"]
    with open(out_dir / "filter_summary.tsv", "w", encoding="utf-8") as fh:
        fh.write("\t".join(fields) + "\n")
        for r in sorted(results, key=lambda x: x["gene"]):
            fh.write("\t".join(str(r.get(f, "")) for f in fields) + "\n")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Input — mutually exclusive: single file OR directory
    inp = p.add_mutually_exclusive_group(required=True)
    inp.add_argument(
        "-i", "--input", metavar="FILE",
        help="Single nucleotide alignment FASTA file (MACSE enrichAlignment output).",
    )
    inp.add_argument(
        "--input-dir", metavar="DIR",
        help="Directory of FASTA files to process in batch.",
    )

    # Output
    p.add_argument(
        "-o", "--output-dir", required=True, metavar="DIR",
        help="Output directory. Created if it does not exist.",
    )

    # Divergence specification (mutually exclusive; last-N is the fallback)
    div = p.add_mutually_exclusive_group()
    div.add_argument(
        "--divergence-pattern", metavar="REGEX",
        help="Regular expression matched against sequence headers to identify "
             "divergence sequences (recommended). "
             "Example: \"Aegilops_tauschii|Triticum_urartu\"",
    )
    div.add_argument(
        "--divergence", metavar="LIST",
        help="Comma-separated list of exact divergence sequence headers.",
    )
    p.add_argument(
        "--n-divergence", type=int, default=3, metavar="N",
        help="When neither --divergence-pattern nor --divergence is given, "
             "treat the last N sequences as divergence (default: 3).",
    )

    # Filtering thresholds
    p.add_argument(
        "--min-polym", type=int, default=6, metavar="N",
        help="Minimum number of polymorphism sequences with a valid codon at "
             "a site for the site to pass (default: 6).",
    )
    p.add_argument(
        "--replace-with", choices=["N", "-"], default="N",
        help="Character used to replace failing codons: "
             "'N' → NNN (default, recommended for codeml / PAML), "
             "'-' → --- (gap).",
    )

    # Performance
    p.add_argument(
        "-e", "--ext", default=None, metavar="EXT",
        help="File extension filter for --input-dir "
             "(e.g. 'fasta'). Default: .fa, .fasta.",
    )
    p.add_argument(
        "-t", "--threads", type=int, default=1, metavar="N",
        help="Parallel workers for --input-dir mode (default: 1).",
    )
    p.add_argument(
        "--verbose", action="store_true",
        help="Print per-gene divergence/polymorphism sequence counts.",
    )

    return p.parse_args()


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Parse explicit divergence list once (shared across files)
    div_list: Optional[List[str]] = (
        [h.strip() for h in args.divergence.split(",") if h.strip()]
        if args.divergence else None
    )

    # ── Single-file mode ──────────────────────────────────────────────────────
    if args.input:
        fpath = Path(args.input)
        if not fpath.is_file():
            sys.exit(f"ERROR: Input file not found: {fpath}")

        seqs = read_fasta(fpath)
        if not seqs:
            sys.exit(f"ERROR: No sequences found in {fpath}")

        headers  = list(seqs.keys())
        aln_len  = len(next(iter(seqs.values())))

        if aln_len % 3 != 0:
            print(f"Warning: alignment length {aln_len} is not a multiple of 3 — "
                  f"last {aln_len % 3} nt ignored.", file=sys.stderr)

        div_hdrs, poly_hdrs = resolve_divergence(
            headers, div_list, args.divergence_pattern, args.n_divergence
        )

        if args.verbose or True:  # always print for single-file mode
            replace_label = "NNN" if args.replace_with == "N" else "---"
            print(f"Input            : {fpath}  ({len(headers)} sequences, {aln_len} nt)")
            print(f"Divergence seqs  : {len(div_hdrs)}  ({', '.join(div_hdrs)})")
            print(f"Polymorphism seqs: {len(poly_hdrs)}")
            print(f"Min polym depth  : {args.min_polym}")
            print(f"Replace failing  : {replace_label}")
            print()

        out_seqs, report_rows, pass_indices = filter_alignment(
            seqs, div_hdrs, poly_hdrs, args.min_polym, args.replace_with
        )

        stem = fpath.stem
        write_outputs(stem, out_dir, out_seqs, report_rows, pass_indices)

        n_codons = len(report_rows)
        n_pass   = len(pass_indices)
        n_fail_div  = sum(1 for r in report_rows if "divergence" in r[3])
        n_fail_poly = sum(1 for r in report_rows if "poly_depth" in r[3])
        frac = n_pass / n_codons if n_codons else 0.0

        print(f"Codons total     : {n_codons}")
        print(f"Codons retained  : {n_pass}  ({100*frac:.1f}%)")
        print(f"  Failed — divergence invalid : {n_fail_div}")
        print(f"  Failed — polym depth low    : {n_fail_poly}")
        print()
        print(f"Filtered FASTA   : {out_dir / (stem + '_filtered.fasta')}")
        print(f"Codon report     : {out_dir / (stem + '_codon_pass.tsv')}")
        print(f"Pass positions   : {out_dir / (stem + '_pass_positions.txt')}")
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

    print(f"Found {len(files)} file(s) in {in_dir}")
    replace_label = "NNN" if args.replace_with == "N" else "---"
    div_spec = (
        f"pattern '{args.divergence_pattern}'" if args.divergence_pattern
        else f"explicit list ({len(div_list)} sequences)" if div_list
        else f"last {args.n_divergence} sequences"
    )
    print(f"Divergence       : {div_spec}")
    print(f"Replace failing  : {replace_label}")
    print(f"Min polym depth  : {args.min_polym}")
    print()

    task_args = [
        (f, out_dir, div_list, args.divergence_pattern,
         args.n_divergence, args.min_polym, args.replace_with)
        for f in files
    ]

    results: List[Dict] = []
    n_total = len(files)

    if args.threads > 1 and n_total > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as pool:
            futures = {pool.submit(_process_file, t): t[0].name for t in task_args}
            done = 0
            for fut in as_completed(futures):
                res = fut.result()
                results.append(res)
                done += 1
                if res["error"]:
                    print(f"[WARN] {res['gene']}: {res['error']}", file=sys.stderr)
                elif args.verbose:
                    print(f"[OK  ] {res['gene']}  "
                          f"{res['n_pass']}/{res['n_codons']} codons  "
                          f"({100*res['frac_pass']:.1f}%)")
                if done % max(1, n_total // 10) == 0:
                    print(f"  {done}/{n_total} processed…")
    else:
        for i, t in enumerate(task_args):
            res = _process_file(t)
            results.append(res)
            if res["error"]:
                print(f"[WARN] {res['gene']}: {res['error']}", file=sys.stderr)
            elif args.verbose:
                print(f"[OK  ] {res['gene']}  "
                      f"{res['n_pass']}/{res['n_codons']} codons  "
                      f"({100*res['frac_pass']:.1f}%)")
            if n_total > 20 and (i + 1) % max(1, n_total // 10) == 0:
                print(f"  {i+1}/{n_total} processed…")

    write_dir_summary(results, out_dir)
    print_summary(results, args, out_dir)


if __name__ == "__main__":
    main()
