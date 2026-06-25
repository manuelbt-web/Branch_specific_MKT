#!/usr/bin/env python3
"""
Compute per-gene polymorphism statistics from codon-aligned FASTA files
(outgroup already removed by remove_outgroup_sequence.py).

Statistics computed:
  lseff         — effective number of analysed nucleotide sites
  S             — number of segregating (polymorphic) sites
  Pi            — total pairwise differences (unnormalised)
  Pi_per_site   — nucleotide diversity (Pi / lseff)
  thetaW        — Watterson's theta (per site)
  Tajimas_D     — Tajima's D neutrality test statistic
  num_pol_S     — number of synonymous polymorphic codon sites
  num_pol_NS    — number of non-synonymous polymorphic codon sites
  pi_S          — mean pairwise differences at synonymous sites
  pi_NS         — mean pairwise differences at non-synonymous sites

EggLib 3.5.2 is used for nucleotide-level statistics (lseff, S, Pi, D, thetaW)
when available (install via pixi — see README.md).  Codon-level statistics
(num_pol_S, num_pol_NS, pi_S, pi_NS) and the folded SFS are always computed
with the built-in pure-Python implementation.

Usage:
  python polymorphisms_stats.py \
      --input-dir   step6_no_outgroup/ \
      --output      step6_stats/polymorphism_stats.tsv \
      [--sfs-output step6_stats/folded_sfs.tsv] \
      [--max-missing 0.0] \
      [--threads 4]

Output:
  --output      TSV table with one row per gene (nucleotide + codon stats)
  --sfs-output  TSV with the folded site frequency spectrum (syn / non-syn)

No mandatory external dependencies.
EggLib 3.5.2 is used automatically when installed (pip/pixi, see README.md).
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import re
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# ── Optional EggLib import ────────────────────────────────────────────────────
import importlib.util
HAS_EGGLIB = importlib.util.find_spec("egglib") is not None

# ── Genetic code ──────────────────────────────────────────────────────────────
_CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_VALID_BASES = frozenset("ACGT")
_ALLELE_SUFFIX = re.compile(r"[._-]?[12aAbB]$")


# ── Haploid downsampling ──────────────────────────────────────────────────────

def _individual_key(header: str) -> str:
    """Extract individual base-ID from a pipeline header 'GENE|species|individual',
    stripping common diploid allele suffixes (_1/_2, _A/_B, .1/.2)."""
    parts = header.split("|")
    ind = parts[2] if len(parts) > 2 else header
    return _ALLELE_SUFFIX.sub("", ind)


def _downsample_to_haploid(
    records: List[Tuple[str, str]], rng: random.Random
) -> Tuple[List[Tuple[str, str]], bool]:
    """Randomly keep one sequence per individual (diploid → n/2 haploid alleles).

    Individuals are identified by the third field of the FASTA header
    (GENE|species|individual), stripping trailing allele suffixes.
    Returns the selected records and a flag indicating whether the sequence
    count was actually reduced (False = no diploid pairing detected).
    """
    groups: Dict[str, List[Tuple[str, str]]] = {}
    for header, seq in records:
        key = _individual_key(header)
        groups.setdefault(key, []).append((header, seq))

    if len(groups) == len(records):
        return list(records), False   # already haploid or no pairing detected

    selected = [rng.choice(alleles) for alleles in groups.values()]
    return selected, True


# ── Pure-Python FASTA reader ──────────────────────────────────────────────────

def _read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header = ""
    parts: List[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(parts)))
                header = line[1:]
                parts = []
            else:
                parts.append(line.strip().upper().replace("U", "T"))
    if header:
        records.append((header, "".join(parts)))
    return records


# ── Temporary FASTA helper ────────────────────────────────────────────────────

def _write_temp_fasta(records: List[Tuple[str, str]]) -> str:
    """Write in-memory records to a temp FASTA file; caller must os.unlink it."""
    fd, path = tempfile.mkstemp(suffix=".fasta")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            for header, seq in records:
                fh.write(f">{header}\n{seq}\n")
    except Exception:
        os.unlink(path)
        raise
    return path


# ── Nucleotide statistics (pure Python) ──────────────────────────────────────

def _pairwise_diffs(chars: Sequence[str]) -> int:
    n = len(chars)
    return sum(
        1 for i in range(n) for j in range(i + 1, n) if chars[i] != chars[j]
    )


def _tajima_d(n: int, s: float, pi: float) -> Optional[float]:
    if n < 2 or s <= 0:
        return None
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))
    b1 = (n + 1) / (3.0 * (n - 1))
    b2 = 2.0 * (n * n + n + 3) / (9.0 * n * (n - 1))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 ** 2)
    var = c1 / a1 * s + c2 / (a1 ** 2 + a2) * s * (s - 1)
    return (pi - s / a1) / (var ** 0.5) if var > 0 else None


def _compute_nt_stats(
    sequences: List[str], max_missing: float
) -> Tuple[int, int, float, Optional[float], float]:
    """Return (lseff, S, Pi_per_site, Tajimas_D, thetaW_per_site)."""
    nseq = len(sequences)
    min_len = min(len(s) for s in sequences)
    lseff = 0
    S = 0
    pi_total = 0

    for i in range(min_len):
        bases = [s[i] for s in sequences if i < len(s) and s[i] in _VALID_BASES]
        if len(bases) / nseq < (1.0 - max_missing):
            continue
        lseff += 1
        if len(set(bases)) > 1:
            S += 1
        pi_total += _pairwise_diffs(bases)

    if lseff == 0:
        return 0, 0, float("nan"), None, float("nan")

    n_pairs = nseq * (nseq - 1) / 2
    pi_per_site = (pi_total / n_pairs) / lseff if n_pairs else float("nan")
    a1 = sum(1.0 / i for i in range(1, nseq)) if nseq > 1 else 1.0
    thetaW = (S / a1) / lseff if a1 and lseff else float("nan")
    d = _tajima_d(nseq, float(S), float(pi_total / n_pairs) if n_pairs else 0.0)
    return lseff, S, pi_per_site, d, thetaW


# ── Codon-level statistics (pure Python) ─────────────────────────────────────

def _codon_type(codons: Sequence[str]) -> Optional[str]:
    """'S' if all variants are synonymous, 'NS' if any amino-acid change."""
    unique = set(codons)
    if len(unique) < 2:
        return None   # monomorphic
    aas = set()
    for c in unique:
        aa = _CODON_TABLE.get(c)
        if aa is None or aa == "*":
            return "NS"
        aas.add(aa)
    return "S" if len(aas) == 1 else "NS"


def _compute_codon_stats(
    sequences: List[str], max_missing: float
) -> Tuple[int, int, float, float, List[Dict]]:
    """Return (pol_S, pol_NS, pi_S, pi_NS, sfs_rows)."""
    nseq = len(sequences)
    codon_limit = min(len(s) for s in sequences)
    codon_limit -= codon_limit % 3

    pol_S = pol_NS = 0
    pi_S_total = pi_NS_total = 0
    sfs_syn: Dict[int, int] = defaultdict(int)
    sfs_ns: Dict[int, int] = defaultdict(int)

    for start in range(0, codon_limit, 3):
        codons = [
            s[start : start + 3]
            for s in sequences
            if start + 3 <= len(s)
            and all(b in _VALID_BASES for b in s[start : start + 3])
        ]
        if len(codons) / nseq < (1.0 - max_missing):
            continue
        ctype = _codon_type(codons)
        if ctype is None:
            continue
        pairs = _pairwise_diffs(codons)
        counts = Counter(codons)
        minor = min(counts.values())
        sample = len(codons)
        folded = min(minor, sample - minor)
        if folded <= 0:
            continue
        if ctype == "S":
            pol_S += 1
            pi_S_total += pairs
            sfs_syn[folded] += 1
        else:
            pol_NS += 1
            pi_NS_total += pairs
            sfs_ns[folded] += 1

    n_pairs = nseq * (nseq - 1) / 2
    pi_S  = (pi_S_total  / n_pairs / pol_S)  if pol_S  and n_pairs else float("nan")
    pi_NS = (pi_NS_total / n_pairs / pol_NS) if pol_NS and n_pairs else float("nan")

    max_bin = nseq // 2
    sfs_rows = [
        {
            "bin": b,
            "frequency": f"{b / nseq:.6f}",
            "syn_count": sfs_syn.get(b, 0),
            "nonsyn_count": sfs_ns.get(b, 0),
        }
        for b in range(1, max_bin + 1)
    ]
    return pol_S, pol_NS, pi_S, pi_NS, sfs_rows


# ── EggLib wrapper (nucleotide stats only) ────────────────────────────────────

def _compute_egglib_nt_stats(
    fasta_path: Path, max_missing: float
) -> Optional[Dict[str, object]]:
    """Use EggLib 3.5.2 for nucleotide-level stats; return None on any failure."""
    if not HAS_EGGLIB:
        return None
    try:
        import egglib
        aln = egglib.io.from_fasta(str(fasta_path), alphabet=egglib.alphabets.DNA)
        cs = egglib.stats.ComputeStats()
        cs.add_stats("lseff", "S", "Pi", "D", "thetaW")
        result = cs.process_align(aln, max_missing=max_missing)
        return {
            "lseff":      result.get("lseff"),
            "S":          result.get("S"),
            "Pi_per_site": (result.get("Pi") / result.get("lseff"))
                            if result.get("lseff") else float("nan"),
            "Tajimas_D":  result.get("D"),
            "thetaW":     result.get("thetaW"),
        }
    except Exception:
        return None


def _compute_egglib_coding_stats(
    fasta_path: Path, max_missing: float, nseq: int
) -> Optional[Dict[str, object]]:
    """Use EggLib CodingDiversity for per-gene codon stats and folded SFS."""
    if not HAS_EGGLIB:
        return None
    try:
        import egglib

        aln  = egglib.io.from_fasta(str(fasta_path), alphabet=egglib.alphabets.DNA)
        cds  = egglib.tools.to_codons(aln)
        cdiv = egglib.stats.CodingDiversity()
        cdiv.process(cds, max_missing=max_missing)

        def _safe_sfs(sites_iter):
            try:
                lst = list(sites_iter)
                if not lst:
                    return []
                result = egglib.stats.SFS(
                    iter(lst), struct=None,
                    max_missing=max_missing, mode="folded",
                    skip_fixed=False,
                )
                return list(result) if result is not None else []
            except Exception:
                return []

        sfs_syn = _safe_sfs(cdiv.sites_S)
        sfs_ns  = _safe_sfs(cdiv.sites_NS)

        maxlen = max(len(sfs_syn), len(sfs_ns), 1)
        sfs_syn += [0] * (maxlen - len(sfs_syn))
        sfs_ns  += [0] * (maxlen - len(sfs_ns))

        sfs_rows = [
            {
                "sample_size":  nseq,
                "bin":          b,
                "frequency":    f"{b / nseq:.6f}" if nseq > 0 else "NA",
                "syn_count":    int(sfs_syn[b]),
                "nonsyn_count": int(sfs_ns[b]),
            }
            for b in range(maxlen)
        ]
        num_pol_S  = sum(int(sfs_syn[b]) for b in range(1, maxlen))
        num_pol_NS = sum(int(sfs_ns[b]) for b in range(1, maxlen))

        return {
            "pi_S":       getattr(cdiv, "Pi_S",    float("nan")),
            "pi_NS":      getattr(cdiv, "Pi_NS",   float("nan")),
            "nseff_S":    getattr(cdiv, "nseff_S", float("nan")),  # eff. synonymous sites
            "nseff_NS":   getattr(cdiv, "nseff",   float("nan")),  # eff. non-syn sites
            "num_pol_S":  num_pol_S,
            "num_pol_NS": num_pol_NS,
            "sfs_rows":   sfs_rows,
        }
    except Exception:
        return None


# ── SDM estimation from folded SFS ───────────────────────────────────────────

def _compute_sdm(sfs_rows: List[Dict], cutoff: float) -> Dict[str, object]:
    """Estimate SDMs from in-memory folded SFS rows.

    Formula: P_SDM = P_N^(<T) - (P_N^(>T) / P_S^(>T)) × P_S^(<T)
    where T = cutoff (default 0.15).
    """
    PnMinus = PnGreater = PsMinus = PsGreater = 0

    for row in sfs_rows:
        b = row["bin"]
        if b == 0:
            continue
        n = row["sample_size"]
        freq = b / n if n > 0 else 0.0
        syn  = row["syn_count"]
        ns   = row["nonsyn_count"]
        if freq <= cutoff:
            PsMinus += syn
            PnMinus += ns
        else:
            PsGreater += syn
            PnGreater += ns

    ratioPs = float(PnGreater) / float(PsGreater) if PsGreater > 0 else float("nan")

    if math.isnan(ratioPs):
        deleterious = pn_neutral = pn_neutral_int = float("nan")
    else:
        deleterious   = float(PnMinus) - ratioPs * float(PsMinus)
        Pn_total      = float(PnMinus + PnGreater)
        pn_neutral    = max(0.0, min(Pn_total - deleterious, Pn_total))
        pn_neutral_int = round(pn_neutral)

    return {
        "PnMinus":              PnMinus,
        "PnGreater":            PnGreater,
        "PsMinus":              PsMinus,
        "PsGreater":            PsGreater,
        "ratioPs":              ratioPs,
        "deleterious":          deleterious,
        "PnNeutral_estimation": pn_neutral,
        "PnNeutral":            pn_neutral_int,
    }


# ── Per-gene computation (called by worker processes) ─────────────────────────

def _fmt(v: object) -> str:
    if v is None or (isinstance(v, float) and v != v):
        return "NA"
    if isinstance(v, float):
        return f"{v:.6f}"
    return str(v)


def process_gene(
    args: Tuple[Path, float, bool, int, float]
) -> Tuple[Dict[str, str], List[Dict]]:
    fasta_path, max_missing, haploid, seed, freq_cutoff = args
    gene = fasta_path.stem
    records = _read_fasta(fasta_path)
    if not records:
        raise ValueError("empty FASTA")

    downsampled = False
    if haploid:
        gene_seed = seed + sum(ord(c) for c in gene)
        records, downsampled = _downsample_to_haploid(records, random.Random(gene_seed))

    sequences = [seq for _, seq in records]
    nseq = len(sequences)

    # When EggLib is available use it for both nucleotide and codon stats.
    # When haploid downsampling was applied, write a temp FASTA so EggLib
    # reads the downsampled sequences rather than the original file.
    nt_result: Optional[Dict] = None
    coding_result: Optional[Dict] = None
    tmp_path: Optional[str] = None

    if HAS_EGGLIB:
        eff_path = fasta_path
        if downsampled:
            tmp_path = _write_temp_fasta(records)
            eff_path = Path(tmp_path)
        try:
            nt_result     = _compute_egglib_nt_stats(eff_path, max_missing)
            coding_result = _compute_egglib_coding_stats(eff_path, max_missing, nseq)
        finally:
            if tmp_path:
                os.unlink(tmp_path)

    # Nucleotide-level stats
    if nt_result:
        lseff, S, pi_site, tajd, thetaW = (
            nt_result["lseff"], nt_result["S"], nt_result["Pi_per_site"],
            nt_result["Tajimas_D"], nt_result["thetaW"],
        )
        nt_engine = "egglib"
    else:
        lseff, S, pi_site, tajd, thetaW = _compute_nt_stats(sequences, max_missing)
        nt_engine = "builtin"

    # Codon / SFS stats
    if coding_result:
        pol_S    = coding_result["num_pol_S"]
        pol_NS   = coding_result["num_pol_NS"]
        pi_S     = coding_result["pi_S"]
        pi_NS    = coding_result["pi_NS"]
        nseff_S  = coding_result["nseff_S"]
        nseff_NS = coding_result["nseff_NS"]
        sfs_rows = coding_result["sfs_rows"]
    else:
        pol_S, pol_NS, pi_S, pi_NS, raw_rows = _compute_codon_stats(sequences, max_missing)
        nseff_S = nseff_NS = float("nan")
        sfs_rows: List[Dict] = [{"sample_size": nseq, "bin": 0, "frequency": "0.000000",
                                  "syn_count": 0, "nonsyn_count": 0}]
        for row in raw_rows:
            row["sample_size"] = nseq
            sfs_rows.append(row)

    # num_codons_eff: prefer EggLib effective-sites sum ÷ 3, else lseff ÷ 3
    if not math.isnan(nseff_S) and not math.isnan(nseff_NS):
        num_codons_eff: object = (nseff_S + nseff_NS) / 3.0
    elif lseff is not None and not math.isnan(float(lseff)):
        num_codons_eff = float(lseff) / 3.0
    else:
        num_codons_eff = float("nan")

    # SDM estimation from the folded SFS (using frequency cutoff T)
    sdm = _compute_sdm(sfs_rows, freq_cutoff)

    stats_row = {
        "gene_id":              gene,
        "number_sites":         _fmt(lseff),
        "num_codons_eff":       _fmt(num_codons_eff),
        "Pi_per_site":          _fmt(pi_site),
        "pi_S":                 _fmt(pi_S),
        "pi_NS":                _fmt(pi_NS),
        "thetaW_per_site":      _fmt(thetaW),
        "Tajimas_D":            _fmt(tajd),
        "num_sites_S_egglib":   _fmt(nseff_S),
        "num_sites_NS_egglib":  _fmt(nseff_NS),
        "polymorphism_count":   str(pol_S + pol_NS),
        "num_pol_S":            str(pol_S),
        "num_pol_NS":           str(pol_NS),
        "sample_size":          str(nseq),
        "PnMinus":              _fmt(sdm["PnMinus"]),
        "PnGreater":            _fmt(sdm["PnGreater"]),
        "PsMinus":              _fmt(sdm["PsMinus"]),
        "PsGreater":            _fmt(sdm["PsGreater"]),
        "ratioPs":              _fmt(sdm["ratioPs"]),
        "PnNeutral_estimation": _fmt(sdm["PnNeutral_estimation"]),
        "PnNeutral":            _fmt(sdm["PnNeutral"]),
        "deleterious":          _fmt(sdm["deleterious"]),
        # Diagnostic columns (engine choice, downsampling flag)
        "downsampled":          "yes" if downsampled else "no",
        "nt_engine":            nt_engine,
        "S":                    _fmt(S),
    }

    for row in sfs_rows:
        row["gene"] = gene

    return stats_row, sfs_rows


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--input-dir", required=True, metavar="DIR",
        help="Directory with codon-aligned FASTA files "
             "(outgroup removed by remove_outgroup_sequence.py)",
    )
    p.add_argument(
        "--output", required=True, metavar="FILE",
        help="Output TSV: one row per gene with all polymorphism statistics",
    )
    p.add_argument(
        "--sfs-output", default=None, metavar="FILE",
        help="Output TSV for the folded site frequency spectrum (optional)",
    )
    p.add_argument(
        "--max-missing", type=float, default=0.30, metavar="FLOAT",
        help="Maximum fraction of missing/gap bases per site to still include "
             "the site in analysis (default: 0.30)",
    )
    p.add_argument(
        "--haploid", action="store_true",
        help="For selfing/diploid species (e.g. Ae. tauschii, T. urartu): "
             "randomly draw one sequence per individual (n/2 alleles) before "
             "computing statistics. Reproducible with --seed.",
    )
    p.add_argument(
        "--seed", type=int, default=42, metavar="INT",
        help="Random seed for haploid downsampling (default: 42). "
             "Each gene uses a deterministic per-gene seed derived from this value.",
    )
    p.add_argument(
        "--freq-cutoff", type=float, default=0.15, metavar="FLOAT",
        help="Folded allele-frequency cutoff T for SDM estimation. Variants with "
             "freq ≤ T are low-frequency (candidate SDMs); variants with freq > T "
             "are high-frequency (effectively neutral). Default: 0.15",
    )
    p.add_argument(
        "--threads", type=int, default=1, metavar="N",
        help="Number of parallel worker processes (default: 1)",
    )
    return p.parse_args()


STATS_FIELDS = [
    # Primary output columns (match target table)
    "gene_id",
    "number_sites", "num_codons_eff",
    "Pi_per_site", "pi_S", "pi_NS",
    "thetaW_per_site", "Tajimas_D",
    "num_sites_S_egglib", "num_sites_NS_egglib",
    "polymorphism_count", "num_pol_S", "num_pol_NS",
    "sample_size",
    "PnMinus", "PnGreater", "PsMinus", "PsGreater",
    "ratioPs", "PnNeutral_estimation", "PnNeutral", "deleterious",
    # Diagnostic columns
    "downsampled", "nt_engine", "S",
]
SFS_FIELDS = ["gene", "sample_size", "bin", "frequency", "syn_count", "nonsyn_count"]


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)

    if not input_dir.is_dir():
        sys.exit(f"ERROR: Input directory not found: {input_dir}")

    fasta_files = sorted(input_dir.glob("*.fasta"))
    if not fasta_files:
        sys.exit(f"ERROR: No .fasta files found in {input_dir}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    engine_note = "EggLib 3.5.2 (nucleotide + codon stats)" if HAS_EGGLIB else "built-in pure Python"
    print("=" * 60)
    print("  Polymorphism Statistics + SDM estimation")
    print("=" * 60)
    print(f"  Input dir    : {input_dir}  ({len(fasta_files)} genes)")
    print(f"  Output       : {args.output}")
    print(f"  Max missing  : {args.max_missing}")
    print(f"  Freq cutoff T: {args.freq_cutoff}  (for SDM estimation)")
    print(f"  Haploid mode : {'yes (seed=' + str(args.seed) + ')' if args.haploid else 'no'}")
    print(f"  Threads      : {args.threads}")
    print(f"  Engine       : {engine_note}")
    if args.haploid and HAS_EGGLIB:
        print("  Note         : haploid sequences written to temp FASTA for EggLib")
    print("=" * 60)

    tasks = [(f, args.max_missing, args.haploid, args.seed, args.freq_cutoff)
             for f in fasta_files]
    stats_rows: List[Dict[str, str]] = []
    all_sfs: List[Dict] = []
    ok = errors = 0

    if args.threads > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as pool:
            futures = {pool.submit(process_gene, t): t[0] for t in tasks}
            for fut in as_completed(futures):
                fpath = futures[fut]
                try:
                    row, sfs = fut.result()
                    stats_rows.append(row)
                    all_sfs.extend(sfs)
                    print(f"  [OK]   {fpath.stem}")
                    ok += 1
                except Exception as exc:
                    print(f"  [ERROR] {fpath.stem}: {exc}")
                    errors += 1
    else:
        for t in tasks:
            fpath = t[0]
            try:
                row, sfs = process_gene(t)
                stats_rows.append(row)
                all_sfs.extend(sfs)
                print(f"  [OK]   {fpath.stem}")
                ok += 1
            except Exception as exc:
                print(f"  [ERROR] {fpath.stem}: {exc}")
                errors += 1

    stats_rows.sort(key=lambda r: r["gene_id"])

    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=STATS_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(stats_rows)

    if args.sfs_output and all_sfs:
        Path(args.sfs_output).parent.mkdir(parents=True, exist_ok=True)
        all_sfs.sort(key=lambda r: (r["gene"], int(r["bin"])))
        with open(args.sfs_output, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=SFS_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(all_sfs)
        print(f"  SFS written  : {args.sfs_output}")

    print()
    print("=" * 60)
    print(f"  Done  —  {ok} genes  |  {errors} errors")
    print(f"  Stats written: {args.output}")
    print("=" * 60)

    if errors > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
