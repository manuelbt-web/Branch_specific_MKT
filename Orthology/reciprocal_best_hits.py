#!/usr/bin/env python3
"""
reciprocal_best_hits.py
========================
Find Reciprocal Best Hits (RBH) from two pairwise BLAST tabular output files.

Two genes are reciprocal best hits if:
  - Gene A in species 1 has gene B as its best BLAST match in species 2, AND
  - Gene B in species 2 has gene A as its best BLAST match in species 1.

RBH is a widely used proxy for ortholog detection between two species.

INPUTS
------
Two BLAST tabular files (format 6, 12 columns), produced by reciprocal_blast.sbatch:

    blastn/blastp -outfmt 6 -query A.fa -subject B.fa -out A_vs_B.tab
    blastn/blastp -outfmt 6 -query B.fa -subject A.fa -out B_vs_A.tab

Column indices (1-based):
    1  = query ID       6  = q.start    11 = bitscore
    2  = subject ID     7  = q.end      12 = evalue
    3  = pident         8  = s.start
    4  = length         9  = s.end
    5  = mismatch       10 = s.end

By default, column 11 (bitscore) is used for scoring, keeping the highest.

USAGE
-----
  python reciprocal_best_hits.py --help
"""

import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--a-vs-b",
        required=True,
        dest="a_vs_b",
        metavar="FILE",
        help="BLAST tabular output (outfmt 6): species A as query, B as subject",
    )
    p.add_argument(
        "--b-vs-a",
        required=True,
        dest="b_vs_a",
        metavar="FILE",
        help="BLAST tabular output (outfmt 6): species B as query, A as subject",
    )
    p.add_argument(
        "--score-col",
        dest="score_col",
        type=int,
        default=12,
        metavar="N",
        help="1-based column number to use as score (default: 12 = bitscore in outfmt 6). "
             "Use 11 for evalue (pair with --order low).",
    )
    p.add_argument(
        "--order",
        choices=["high", "low"],
        default="high",
        help="Keep hit with highest (bitscore) or lowest (evalue) score per query "
             "(default: high)",
    )
    p.add_argument(
        "--output",
        required=True,
        metavar="FILE",
        help="Output file: tab-delimited table of RBH pairs",
    )
    return p.parse_args()


def load_best_hits(filepath, score_col_0based, want_highest):
    """
    Parse a BLAST tabular (outfmt 6) file and return the best hit per query.

    Returns: {query_id: (subject_id, float_score, str_score)}
    """
    best = {}
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) <= score_col_0based:
                continue
            query   = parts[0]
            subject = parts[1]
            try:
                score = float(parts[score_col_0based])
            except ValueError:
                continue

            if query not in best or \
               (want_highest and score > best[query][1]) or \
               (not want_highest and score < best[query][1]):
                best[query] = (subject, score, parts[score_col_0based])
    return best


def main():
    args = parse_args()

    for path, flag in [(args.a_vs_b, "--a-vs-b"), (args.b_vs_a, "--b-vs-a")]:
        if not os.path.isfile(path):
            print(f"ERROR: {flag} file not found: {path}", file=sys.stderr)
            sys.exit(1)

    if args.output in [args.a_vs_b, args.b_vs_a]:
        print("ERROR: --output would overwrite one of the input files", file=sys.stderr)
        sys.exit(1)

    score_col    = args.score_col - 1     # convert to 0-based
    want_highest = (args.order == "high")

    print(f"Loading A→B hits : {args.a_vs_b}")
    best_a_vs_b = load_best_hits(args.a_vs_b, score_col, want_highest)
    print(f"  {len(best_a_vs_b)} unique query sequences")

    print(f"Loading B→A hits : {args.b_vs_a}")
    best_b_vs_a = load_best_hits(args.b_vs_a, score_col, want_highest)
    print(f"  {len(best_b_vs_a)} unique query sequences")

    # A and B are reciprocal best hits if:
    #   best hit of A → B   is   B
    #   best hit of B → A   is   A
    count = 0
    with open(args.output, "w") as out:
        out.write("#A_id\tB_id\tA_vs_B_score\tB_vs_A_score\n")
        for a, (b, _, score_str_ab) in best_a_vs_b.items():
            if b in best_b_vs_a:
                a_prime, _, score_str_ba = best_b_vs_a[b]
                if a == a_prime:
                    out.write(f"{a}\t{b}\t{score_str_ab}\t{score_str_ba}\n")
                    count += 1

    print(f"\nDone: {count} reciprocal best hits written to {args.output}")


if __name__ == "__main__":
    main()
