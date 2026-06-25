#!/usr/bin/env python3
"""
Remove outgroup sequences from per-gene FASTA alignment files.

Sequences are identified as outgroup by the presence of a tag in their
header (default: '|outgroup'), as set by step 4.  Unlike position-based
removal, this is robust to changes in sequence ordering.

Input files are never modified — cleaned FASTAs are written to OUTPUT_DIR.

Usage:
  python remove_outgroup_sequence.py \
      --input-dir  step4_out/ \
      --output-dir step6_no_outgroup/ \
      [--outgroup-tag "|outgroup"] \
      [--min-sequences 2]

No external dependencies (pure Python).
"""

import argparse
import sys
from pathlib import Path
from typing import List, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input-dir", required=True, metavar="DIR",
        help="Directory with aligned FASTA files (one per gene, focal + outgroup)",
    )
    p.add_argument(
        "--output-dir", required=True, metavar="DIR",
        help="Output directory for FASTA files without the outgroup "
             "(created if absent; input files are never modified)",
    )
    p.add_argument(
        "--outgroup-tag", default="|outgroup", metavar="TAG",
        help="String that identifies the outgroup sequence in the FASTA header "
             "(default: '|outgroup', as set by step 4)",
    )
    p.add_argument(
        "--min-sequences", type=int, default=2, metavar="N",
        help="Minimum number of sequences required after outgroup removal; "
             "genes with fewer are skipped (default: 2)",
    )
    return p.parse_args()


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    """Return list of (header_line, sequence) for every record in a FASTA file."""
    records: List[Tuple[str, str]] = []
    header: str = ""
    seq_parts: List[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq_parts)))
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
    if header:
        records.append((header, "".join(seq_parts)))
    return records


def write_fasta(path: Path, records: List[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(f"{header}\n")
            # write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def remove_outgroup(
    fasta_path: Path,
    output_dir: Path,
    outgroup_tag: str,
    min_sequences: int,
) -> Tuple[str, int, int]:
    """Return (status, n_kept, n_removed)."""
    records = read_fasta(fasta_path)
    if not records:
        return "empty", 0, 0

    kept   = [(h, s) for h, s in records if outgroup_tag not in h]
    removed = len(records) - len(kept)

    if len(kept) < min_sequences:
        return (
            f"skipped — only {len(kept)} focal sequence(s) remain "
            f"(min required: {min_sequences})",
            len(kept), removed,
        )

    write_fasta(output_dir / fasta_path.name, kept)
    return "ok", len(kept), removed


def main() -> None:
    args = parse_args()

    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    if not input_dir.is_dir():
        sys.exit(f"ERROR: Input directory not found: {input_dir}")

    fasta_files = sorted(input_dir.glob("*.fasta"))
    if not fasta_files:
        sys.exit(f"ERROR: No .fasta files found in {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("  Remove Outgroup Sequences")
    print("=" * 60)
    print(f"  Input dir     : {input_dir}  ({len(fasta_files)} files)")
    print(f"  Output dir    : {output_dir}")
    print(f"  Outgroup tag  : '{args.outgroup_tag}'")
    print(f"  Min sequences : {args.min_sequences}")
    print("=" * 60)

    ok = skipped = errors = 0

    for fasta in fasta_files:
        try:
            status, n_kept, n_removed = remove_outgroup(
                fasta, output_dir, args.outgroup_tag, args.min_sequences
            )
        except Exception as exc:
            print(f"  [ERROR]  {fasta.stem}: {exc}")
            errors += 1
            continue

        if status == "ok":
            print(f"  [OK]     {fasta.stem} — {n_kept} focal kept, {n_removed} outgroup removed")
            ok += 1
        else:
            print(f"  [SKIP]   {fasta.stem} — {status}")
            skipped += 1

    print()
    print("=" * 60)
    print(f"  Done  —  {ok} written  |  {skipped} skipped  |  {errors} errors")
    print("=" * 60)

    if errors > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
