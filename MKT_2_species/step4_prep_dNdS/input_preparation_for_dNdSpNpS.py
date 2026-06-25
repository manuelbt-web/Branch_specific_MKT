#!/usr/bin/env python3
"""
Reorder outgroup sequences to the end of each FASTA file and standardise
sequence headers in the format expected by the dNdSpNpS analysis pipeline.

For each FASTA in INPUT_DIR:
  1. Identify focal sequences (those with '|sp|' in their ID, or sharing the
     most common gene-ID prefix) and outgroup sequences (all others).
  2. Select one representative outgroup record (preferring '|outgroup|-labelled
     or name-matching records; falling back to the longest).
  3. Rename headers:
       Focal    : replace the '|sp|' token with '|<species>|'
       Outgroup : <gene_prefix>|<outgroup_name>|<original_id>|outgroup
  4. Write focal records first, outgroup last, to OUTPUT_DIR (or in-place if
     --output-dir is omitted).

Dependencies: BioPython  (pip install biopython)

Usage:
  python input_preparation_for_dNdSpNpS.py \\
      --input-dir  aligned/ \\
      --species    Ae_speltoides \\
      --outgroup   H_vulgare \\
      [--output-dir prepared/] \\
      [--log       prepared/prep.log] \\
      [--threads   4]
"""

import os
import sys
import glob
import argparse
import logging
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ── Core logic ────────────────────────────────────────────────────────────────

def reorder_and_rename(fasta_path: str, species: str, outgroup_name: str,
                       output_path: str) -> str:
    """
    Reorder one FASTA file so that focal sequences come first and the outgroup
    last, and rename headers for dNdSpNpS.

    Returns a status string for logging.
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return f"[WARN]  {fasta_path} — no sequences found, skipped"

    # ── Identify focal prefix ─────────────────────────────────────────────────
    # Prefer a record already labelled with the focal marker '|sp|'
    focal_prefix = None
    for rec in records:
        if "|sp|" in rec.id:
            focal_prefix = rec.id.split("|")[0]
            break

    if focal_prefix is None:
        # Fallback: most common first-field prefix
        prefixes = [r.id.split("|")[0] for r in records if r.id]
        if prefixes:
            focal_prefix = Counter(prefixes).most_common(1)[0][0]

    if focal_prefix is None:
        focal_prefix = records[0].id.split("|")[0]

    # ── Classify records ──────────────────────────────────────────────────────
    focal_records = []
    outgroup_candidates = []

    for rec in records:
        is_focal = ("|sp|" in rec.id) or (rec.id.split("|")[0] == focal_prefix)
        if is_focal:
            focal_records.append(rec)
        else:
            outgroup_candidates.append(rec)

    # ── Select one representative outgroup record ─────────────────────────────
    if not outgroup_candidates:
        # Check whether an already-labelled outgroup exists among focal records
        existing = [r for r in focal_records if "|outgroup" in r.id]
        if existing:
            return f"[INFO]  {fasta_path} — outgroup already labelled, no changes needed"
        return f"[INFO]  {fasta_path} — no outgroup sequence found, skipped"

    labelled = [r for r in outgroup_candidates if "|outgroup" in r.id]
    if labelled:
        outgroup_rec = labelled[0]
    else:
        by_name = [r for r in outgroup_candidates if outgroup_name in r.id]
        outgroup_rec = by_name[0] if by_name else max(outgroup_candidates,
                                                      key=lambda r: len(r.seq))

    # ── Rename focal headers: replace |sp| with |<species>| ──────────────────
    renamed_focal = []
    for rec in focal_records:
        if "|sp|" in rec.id:
            new_id = rec.id.replace("|sp|", f"|{species}|")
            rec = SeqRecord(rec.seq, id=new_id, name=new_id, description="")
        renamed_focal.append(rec)

    # ── Rename outgroup header ────────────────────────────────────────────────
    orig_id = outgroup_rec.id.replace("|", "_")
    new_out_id = f"{focal_prefix}|{outgroup_name}|{orig_id}|outgroup"
    renamed_out = SeqRecord(outgroup_rec.seq, id=new_out_id,
                            name=new_out_id, description="")

    # ── Write output ──────────────────────────────────────────────────────────
    tmp_path = output_path + ".tmp"
    with open(tmp_path, "w") as fh:
        SeqIO.write(renamed_focal + [renamed_out], fh, "fasta")
    os.replace(tmp_path, output_path)

    n_focal = len(renamed_focal)
    return f"[OK]    {os.path.basename(fasta_path)} — {n_focal} focal + 1 outgroup ({new_out_id})"


def _worker(args):
    return reorder_and_rename(*args)


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input-dir", required=True, metavar="DIR",
        help="Directory containing aligned FASTA files to process",
    )
    p.add_argument(
        "--species", required=True, metavar="NAME",
        help="Species name that replaces the '|sp|' token in focal sequence headers "
             "(e.g. 'Ae_speltoides')",
    )
    p.add_argument(
        "--outgroup", required=True, metavar="NAME",
        help="Outgroup species name inserted into the outgroup sequence header "
             "(e.g. 'H_vulgare')",
    )
    p.add_argument(
        "--output-dir", metavar="DIR", default=None,
        help="Directory for processed FASTA files. If omitted, files are "
             "modified in-place inside --input-dir.",
    )
    p.add_argument(
        "--log", metavar="FILE", default=None,
        help="Log file (appended). Defaults to OUTPUT_DIR/prep_dndspnps.log "
             "or INPUT_DIR/prep_dndspnps.log when --output-dir is omitted.",
    )
    p.add_argument(
        "--threads", type=int, default=1, metavar="N",
        help="Number of parallel worker processes (default: 1)",
    )
    return p.parse_args()


def main():
    args = parse_args()

    # Validate inputs
    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        sys.exit(f"ERROR: Input directory not found: {input_dir}")

    fasta_files = sorted(input_dir.glob("*.fasta"))
    if not fasta_files:
        sys.exit(f"ERROR: No .fasta files found in {input_dir}")

    # Output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = input_dir  # in-place

    # Log file
    log_path = args.log or str(output_dir / "prep_dndspnps.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_path, mode="a"),
        ],
    )
    log = logging.getLogger()

    log.info("=" * 62)
    log.info("  input_preparation_for_dNdSpNpS")
    log.info("=" * 62)
    log.info(f"  Input dir   : {input_dir}  ({len(fasta_files)} files)")
    log.info(f"  Output dir  : {output_dir}")
    log.info(f"  Species     : {args.species}")
    log.info(f"  Outgroup    : {args.outgroup}")
    log.info(f"  Threads     : {args.threads}")
    log.info(f"  Log         : {log_path}")
    log.info("=" * 62)

    # Build task list: (fasta_path, species, outgroup_name, output_path)
    tasks = []
    for fasta in fasta_files:
        out_path = str(output_dir / fasta.name)
        # When processing in-place, source == dest is fine (atomic os.replace)
        tasks.append((str(fasta), args.species, args.outgroup, out_path))

    # Process files (sequential or parallel)
    ok = 0; warn = 0; skipped = 0

    if args.threads > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as pool:
            futures = {pool.submit(_worker, t): t[0] for t in tasks}
            for fut in as_completed(futures):
                try:
                    msg = fut.result()
                except Exception as exc:
                    msg = f"[ERROR] {futures[fut]}: {exc}"
                log.info(msg)
                if "[OK]" in msg:
                    ok += 1
                elif "[WARN]" in msg or "[ERROR]" in msg:
                    warn += 1
                else:
                    skipped += 1
    else:
        for t in tasks:
            try:
                msg = reorder_and_rename(*t)
            except Exception as exc:
                msg = f"[ERROR] {t[0]}: {exc}"
            log.info(msg)
            if "[OK]" in msg:
                ok += 1
            elif "[WARN]" in msg or "[ERROR]" in msg:
                warn += 1
            else:
                skipped += 1

    log.info("")
    log.info("=" * 62)
    log.info(f"  Done  —  {ok} processed  |  {skipped} skipped  |  {warn} warnings/errors")
    log.info("=" * 62)


if __name__ == "__main__":
    main()
