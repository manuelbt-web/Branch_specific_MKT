#!/usr/bin/env python3
"""
extraction_of_covered_contigs.py
=================================
Extract CDS sequences from a multi-FASTA consensus file and retain only those
genes whose sequences are sufficiently complete across individuals.

DEFINITIONS
-----------
  sequence completeness : fraction of non-N bases in the extracted CDS region
                          (N = missing/uncertain base in the consensus sequence).
                          This is NOT the same as sequencing depth (coverage).

WORKFLOW
--------
  1. Read the CDS position table produced by compute_cds_positions.sh.
  2. For each gene, extract the CDS region from every consensus sequence
     (i.e., every individual in the multi-FASTA).
  3. Compute sequence completeness = 1 - (count_N / cds_length).
  4. Count how many sequences per gene meet --min-completeness.
  5. Retain a gene if at least --min-fraction of its sequences are complete.
  6. Write:
       - a plain-text list of retained gene IDs
       - one FASTA file per retained gene containing its complete CDS sequences

NOTE ON FASTA FORMAT
--------------------
  The multi-FASTA is expected to contain multiple sequences per gene — one
  per individual.  Sequences are grouped by the first '|'-delimited token of
  their FASTA header (e.g. ">GENE001|sample_1" → gene "GENE001").

DEPENDENCIES
------------
  Python >= 3.7
  Biopython >= 1.79

  Install with:
    pip install biopython
  or:
    conda install -c conda-forge biopython

USAGE
-----
  python extraction_of_covered_contigs.py --help
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ==============================================================================
# Argument parsing
# ==============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        prog="extraction_of_covered_contigs.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    p.add_argument(
        "--fasta",
        required=True,
        metavar="FILE",
        help="Multi-FASTA file containing consensus sequences for all individuals. "
             "FASTA headers must use '|' to separate gene ID from sample ID "
             "(e.g. '>GENE001|sample_1').",
    )
    p.add_argument(
        "--cds-table",
        required=True,
        dest="cds_table",
        metavar="FILE",
        help="CDS position table produced by compute_cds_positions.sh "
             "(tab-delimited, with header: full_transcript_id | transcript_id | "
             "cds_start | cds_end | ...).",
    )
    p.add_argument(
        "--min-completeness",
        required=True,
        dest="min_completeness",
        type=float,
        metavar="FLOAT",
        help="Minimum fraction of non-N bases in the CDS for a single sequence "
             "to be counted as 'complete'.  Range 0–1 (e.g. 0.7 = 70%% non-N).",
    )
    p.add_argument(
        "--min-fraction",
        required=True,
        dest="min_fraction",
        type=float,
        metavar="FLOAT",
        help="Minimum fraction of sequences per gene that must meet "
             "--min-completeness for the gene to be retained.  "
             "Range 0–1 (e.g. 0.5 = at least 50%% of individuals).",
    )
    p.add_argument(
        "--gene-list",
        required=True,
        dest="gene_list",
        metavar="FILE",
        help="Output file: list of retained gene IDs, one per line.",
    )
    p.add_argument(
        "--cds-dir",
        required=True,
        dest="cds_dir",
        metavar="DIR",
        help="Output directory: one FASTA file per retained gene "
             "(<GENE_ID>_CDS.fasta), containing only the complete CDS sequences.",
    )

    return p.parse_args()


# ==============================================================================
# Helper functions
# ==============================================================================

def load_cds_table(path):
    """Return {transcript_id: (cds_start_1based, cds_end_1based)}."""
    positions = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene_id   = row["transcript_id"]
            cds_start = int(row["cds_start"])
            cds_end   = int(row["cds_end"])
            positions[gene_id] = (cds_start, cds_end)
    return positions


def sequence_completeness(seq):
    """Return the fraction of non-N bases (0.0–1.0).  Returns 0.0 for empty sequences."""
    length = len(seq)
    if length == 0:
        return 0.0
    n_count = seq.upper().count("N")
    return 1.0 - (n_count / length)


# ==============================================================================
# Main
# ==============================================================================

def main():
    args = parse_args()

    # ── Validate inputs ────────────────────────────────────────────────────────
    errors = []
    for path, flag in [(args.fasta, "--fasta"), (args.cds_table, "--cds-table")]:
        if not os.path.isfile(path):
            errors.append(f"{flag}: file not found: {path}")

    for val, flag in [
        (args.min_completeness, "--min-completeness"),
        (args.min_fraction,     "--min-fraction"),
    ]:
        if not 0.0 <= val <= 1.0:
            errors.append(f"{flag}: must be between 0 and 1, got {val}")

    if errors:
        for err in errors:
            print(f"ERROR: {err}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.cds_dir, exist_ok=True)

    # ── Load CDS positions ─────────────────────────────────────────────────────
    print(f"Loading CDS table : {args.cds_table}")
    cds_pos = load_cds_table(args.cds_table)
    print(f"  {len(cds_pos)} genes with CDS annotations")

    # ── Group sequences by gene (first '|'-token of FASTA header) ─────────────
    print(f"Parsing FASTA     : {args.fasta}")
    records_by_gene = defaultdict(list)
    total_sequences = 0
    for record in SeqIO.parse(args.fasta, "fasta"):
        gene_name = record.id.split("|")[0]
        records_by_gene[gene_name].append(record)
        total_sequences += 1
    print(f"  {total_sequences} sequences across {len(records_by_gene)} genes")

    # ── Extract CDS and evaluate completeness ──────────────────────────────────
    valid_count  = defaultdict(int)   # sequences meeting --min-completeness
    total_count  = defaultdict(int)   # sequences with a valid CDS region
    cds_records  = defaultdict(list)  # CDS SeqRecords to write if gene is retained
    skipped_cds  = 0                  # sequences where CDS extends past sequence end

    for gene_id, records in records_by_gene.items():
        if gene_id not in cds_pos:
            continue   # no CDS annotation for this gene

        start_1, end_1 = cds_pos[gene_id]
        py_start = start_1 - 1   # convert to 0-based (Python slice start)
        py_end   = end_1          # Python slice end is exclusive, so end_1 works

        for record in records:
            if py_end > len(record.seq):
                skipped_cds += 1
                continue   # CDS extends beyond this sequence — skip

            cds_seq      = record.seq[py_start:py_end]
            completeness = sequence_completeness(cds_seq)
            total_count[gene_id] += 1

            if completeness >= args.min_completeness:
                valid_count[gene_id] += 1
                cds_rec = SeqRecord(
                    cds_seq,
                    id          = record.id + "|CDS",
                    description = "",
                )
                cds_records[gene_id].append(cds_rec)

    # ── Select genes meeting --min-fraction ────────────────────────────────────
    retained = []
    for gene_id, total in total_count.items():
        if total == 0:
            continue
        ratio = valid_count.get(gene_id, 0) / total
        if ratio >= args.min_fraction:
            retained.append(gene_id)
    retained.sort()

    # ── Print summary ──────────────────────────────────────────────────────────
    print()
    print("Results:")
    print(f"  Genes with CDS annotations : {len(cds_pos)}")
    print(f"  Genes found in FASTA       : {len(records_by_gene)}")
    print(f"  Genes evaluated            : {len(total_count)}")
    print(f"  Genes retained             : {len(retained)}")
    if skipped_cds > 0:
        print(f"  Sequences skipped          : {skipped_cds}  (CDS exceeded sequence length)")
    print()
    print("Filters applied:")
    print(f"  --min-completeness {args.min_completeness:.0%}  "
          f"(each sequence must have >= {args.min_completeness:.0%} non-N bases in CDS)")
    print(f"  --min-fraction     {args.min_fraction:.0%}  "
          f"(at least {args.min_fraction:.0%} of sequences per gene must be complete)")

    # ── Write gene list ────────────────────────────────────────────────────────
    with open(args.gene_list, "w") as fh:
        for gene in retained:
            fh.write(f"{gene}\n")
    print()
    print(f"Gene list written : {args.gene_list}")

    # ── Write one FASTA per retained gene ─────────────────────────────────────
    written = 0
    for gene_id in retained:
        if gene_id in cds_records:
            out_path = os.path.join(args.cds_dir, f"{gene_id}_CDS.fasta")
            with open(out_path, "w") as fh:
                SeqIO.write(cds_records[gene_id], fh, "fasta")
            written += 1
    print(f"CDS FASTA files   : {written} files written to {args.cds_dir}/")


if __name__ == "__main__":
    main()
