#!/usr/bin/env python3
"""
macse_enrichment.py

Automate MACSE enrichAlignment for orthogroup + polymorphism integration.

Features:
- Reads orthogroup table (TSV) mapping orthogroup -> species gene ids
- Locates orthogroup FASTA and species polymorphism FASTA
- Runs MACSE enrichAlignment per pair, validates *_stats.csv for added;yes
- Logs summary, failed alignments, and macse commands; produces rerun script
- Supports --threads, --dry-run, --resume and --overwrite options

"""
import argparse
import csv
import os
import sys
import subprocess
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Tuple


def read_ortho_table(path: str):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().strip().split('\t')
        cols = hdr
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            rec = dict(zip(cols, parts))
            rows.append(rec)
    return cols, rows


def find_orthofasta(ortho_dir: str, orthogroup: str):
    patterns = [f"{orthogroup}_aln_hmm_codon.fasta", f"{orthogroup}*.fasta"]
    for p in patterns:
        path = os.path.join(ortho_dir, p)
        matches = glob.glob(path)
        if matches:
            return matches[0]
    return None


def find_polyfasta(poly_dir: str, gene_id: str):
    # common naming: {geneid}_CDS_aligned_NT.fasta or {geneid}*.fasta
    patterns = [f"{gene_id}_CDS_aligned_NT.fasta", f"{gene_id}*.fasta"]
    for p in patterns:
        matches = glob.glob(os.path.join(poly_dir, p))
        if matches:
            # if multiple matches, pick first
            return matches[0]
    return None


def parse_stats_csv(stats_path: str) -> Tuple[bool, list]:
    # Return (all_added_bool, list_of_rows_with_added_field)
    if not os.path.exists(stats_path):
        return False, []
    with open(stats_path, 'r', encoding='utf-8', errors='ignore') as fh:
        lines = [l.strip() for l in fh if l.strip()]
    if not lines:
        return False, []
    # detect delimiter
    delim = ',' if ',' in lines[0] and ';' not in lines[0] else ';'
    reader = csv.DictReader(lines, delimiter=delim)
    rows = list(reader)
    if not rows:
        return False, []
    all_added = True
    failed = []
    for r in rows:
        val = r.get('added') or r.get('Added') or r.get('added?')
        if val is None:
            # try to find by header that contains 'added'
            for k in r.keys():
                if 'add' in k.lower():
                    val = r[k]
                    break
        if val is None or val.strip().lower() not in ('yes', 'y', 'true'):
            all_added = False
            failed.append(r)
    return all_added, failed


def build_macse_cmd(java: str, macse_jar: str, polym_fasta: str, ortho_fasta: str, out_prefix: str):
    cmd = [java, '-jar', macse_jar, '-prog', 'enrichAlignment',
           '-align', polym_fasta,
           '-seq', ortho_fasta,
           '-fixed_alignment_ON',
           '-maxINS_inSeq', '1',
           '-maxSTOP_inSeq', '0',
           '-maxFS_inSeq', '1',
           '-out_NT', f'{out_prefix}_NT.fasta',
           '-out_AA', f'{out_prefix}_AA.fasta']
    return cmd


def process_row(rec, args, cols) -> dict:
    orthogroup = rec.get('orthogroup') or rec.get(cols[0])
    species_gene = rec.get(args.species_col)
    result = dict(orthogroup=orthogroup, speltoides_gene=species_gene,
                  ortho_file='-', ortho_found='NO', poly_file='-', poly_found='NO',
                  macse_exit='-', all_added='-', output_files=0)

    ortho_path = find_orthofasta(args.ortho_dir, orthogroup)
    if ortho_path:
        result['ortho_file'] = os.path.basename(ortho_path)
        result['ortho_found'] = 'YES'
    else:
        return result

    poly_path = find_polyfasta(args.poly_dir, species_gene)
    if poly_path:
        result['poly_file'] = os.path.basename(poly_path)
        result['poly_found'] = 'YES'
    else:
        return result

    out_base = os.path.join(args.output_dir, 'enriched', f"{species_gene.replace('.','_')}_{orthogroup}")
    os.makedirs(os.path.dirname(out_base), exist_ok=True)
    nt_out = f"{out_base}_NT.fasta"
    aa_out = f"{out_base}_AA.fasta"
    stats_path = f"{out_base}_stats.csv"

    # resume/skip
    if args.resume and os.path.exists(nt_out) and os.path.exists(aa_out) and os.path.exists(stats_path) and not args.overwrite:
        result['macse_exit'] = 0
        result['output_files'] = 3
        ok, failed = parse_stats_csv(stats_path)
        result['all_added'] = 'YES' if ok else 'NO'
        return result

    cmd = build_macse_cmd(args.java_bin, args.macse_jar, poly_path, ortho_path, out_base)
    # log command
    with open(os.path.join(args.output_dir, 'logs', 'macse_commands.log'), 'a') as mf:
        mf.write(' '.join(cmd) + '\n')

    if args.dry_run:
        result['macse_exit'] = 'DRY'
        return result

    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        result['macse_exit'] = p.returncode
        # write per-job logs
        joblog = os.path.join(args.output_dir, 'logs', f"{species_gene.replace('.','_')}_{orthogroup}.log")
        with open(joblog, 'w') as lf:
            lf.write('STDOUT:\n')
            lf.write(p.stdout + '\n')
            lf.write('STDERR:\n')
            lf.write(p.stderr + '\n')
    except Exception as e:
        result['macse_exit'] = -1
        with open(os.path.join(args.output_dir, 'logs', 'macse_commands.log'), 'a') as mf:
            mf.write(f'ERROR running {orthogroup} {species_gene}: {e}\n')
        return result

    # check outputs
    outs = [nt_out, aa_out]
    out_count = sum(1 for o in outs if os.path.exists(o) and os.path.getsize(o) > 0)
    # attempt to find stats csv
    stats_candidates = [f for f in glob.glob(f"{out_base}*stats*.csv")] + [stats_path]
    stats_file = stats_candidates[0] if stats_candidates else None
    if stats_file and os.path.exists(stats_file):
        ok, failed_rows = parse_stats_csv(stats_file)
        result['all_added'] = 'YES' if ok else 'NO'
    else:
        result['all_added'] = 'NO'
    result['output_files'] = out_count + (1 if stats_file and os.path.exists(stats_file) else 0)
    return result


def main():
    p = argparse.ArgumentParser(description='MACSE enrichment pipeline')
    p.add_argument('--ortho-table', required=True)
    p.add_argument('--ortho-dir', required=True)
    p.add_argument('--poly-dir', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--macse-jar', required=True)
    p.add_argument('--species-col', required=True)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--overwrite', action='store_true')
    p.add_argument('--java-bin', default='java')
    args = p.parse_args()

    os.makedirs(os.path.join(args.output_dir, 'enriched'), exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, 'logs'), exist_ok=True)

    cols, rows = read_ortho_table(args.ortho_table)

    summary_tsv = os.path.join(args.output_dir, 'logs', 'enrichment_summary.tsv')
    failed_txt = os.path.join(args.output_dir, 'logs', 'failed_alignments.txt')
    macse_cmd_log = os.path.join(args.output_dir, 'logs', 'macse_commands.log')

    # Initialize logs early so files exist while the script runs
    with open(summary_tsv, 'w') as sf:
        sf.write('\t'.join(['orthogroup','speltoides_gene','ortho_file','ortho_found','poly_file','poly_found','macse_exit','all_added','output_files']) + '\n')
    # ensure failed file exists (empty)
    open(failed_txt, 'w').close()
    # ensure macse command log exists
    open(macse_cmd_log, 'a').close()

    # run in parallel
    results = []
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = {ex.submit(process_row, r, args, cols): r for r in rows}
        for fut in as_completed(futs):
            r = futs[fut]
            try:
                res = fut.result()
            except Exception as e:
                res = dict(orthogroup=r.get('orthogroup','?'), speltoides_gene=r.get(args.species_col,'?'), ortho_file='-', ortho_found='NO', poly_file='-', poly_found='NO', macse_exit='ERR', all_added='NO', output_files=0)
            results.append(res)

    # write summary and failed lists (append results to existing summary header)
    with open(summary_tsv, 'a') as sf, open(failed_txt, 'a') as ff:
        for r in results:
            sf.write('\t'.join([str(r.get(k,'')) for k in ['orthogroup','speltoides_gene','ortho_file','ortho_found','poly_file','poly_found','macse_exit','all_added','output_files']]) + '\n')
            if r.get('ortho_found') != 'YES' or r.get('poly_found') != 'YES' or str(r.get('macse_exit')) != '0' or r.get('all_added') != 'YES':
                ff.write(f"{r.get('orthogroup')}\t{r.get('speltoides_gene')}\n")
    os.chmod(rerun_sh, 0o755)
    print('Done. Summary written to', summary_tsv)


if __name__ == '__main__':
    main()


#### sbatch 

#!/bin/bash
#SBATCH --job-name=macse_enrich
#SBATCH --output=macse_%j.out
#SBATCH --error=macse_%j.err
#SBATCH --time=142:00:00
#SBATCH --mem=16G
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-long  # À vérifier sur votre cluster

# Charger les modules (vérifier les noms exacts sur votre cluster)
module purge
module load bioinfo-ifb python/3.9 java-jdk  # ou "module load bioinfo-ifb python/3.8"

# === CHEMINS COMPLETS (À ADAPTER À VOTRE CAS) ===
# 1. Table des orthogroupes
ORTHO_TABLE="/home/barrientosm/scratch_barrientosm/non_homologous_sequences/speltoides_mutica_tauschii_orthofinder/orthologs_hmm_cleaner/cleaned_alignments_nucleotides/orthogroup_ids.txt"

# 2. Répertoire des orthologues (fichiers HOG)
ORTHO_DIR="/home/barrientosm/scratch_barrientosm/non_homologous_sequences/speltoides_mutica_tauschii_orthofinder/orthologs_hmm_cleaner/cleaned_alignments_nucleotides"

# 3. Répertoire des polymorphismes (fichiers speltoides)
POLY_DIR="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/02_results/MACSE_alignments/speltoides/cds_sequences_polymorphism"

# 4. Fichier MACSE
MACSE_JAR="/home/barrientosm/scratch_barrientosm/non_homologous_sequences/speltoides_mutica_tauschii_orthofinder/macse_enrichment/macse_v2.07.jar"  # À vérifier où il est

# 5. Répertoire de sortie
OUTPUT_DIR="/home/barrientosm/scratch_barrientosm/non_homologous_sequences/speltoides_mutica_tauschii_orthofinder/macse_enrichment/enrichment_speltoides/enriched"

# === VÉRIFICATIONS ===
echo "=== Vérification des chemins ==="
echo "Ortho table: $ORTHO_TABLE"
[ -f "$ORTHO_TABLE" ] && echo "✓ Existe" || echo "✗ MANQUANT !"
echo "Ortho dir: $ORTHO_DIR"
[ -d "$ORTHO_DIR" ] && echo "✓ Existe ($(ls "$ORTHO_DIR"/*.fasta 2>/dev/null | wc -l) fichiers FASTA)" || echo "✗ MANQUANT !"
echo "Poly dir: $POLY_DIR"
[ -d "$POLY_DIR" ] && echo "✓ Existe ($(ls "$POLY_DIR"/*.fasta 2>/dev/null | wc -l) fichiers FASTA)" || echo "✗ MANQUANT !"
echo "MACSE jar: $MACSE_JAR"
[ -f "$MACSE_JAR" ] && echo "✓ Existe" || echo "✗ MANQUANT !"

# Attendre 5 secondes pour vérifier
sleep 5

# === EXÉCUTION ===
echo "=== Lancement MACSE enrichment ==="
echo "Date: $(date)"
echo "Output: $OUTPUT_DIR"

python3 macse_enrichment.py \
  --ortho-table "$ORTHO_TABLE" \
  --ortho-dir "$ORTHO_DIR" \
  --poly-dir "$POLY_DIR" \
  --output-dir "$OUTPUT_DIR" \
  --macse-jar "$MACSE_JAR" \
  --species-col "Aegilops_speltoides"

echo "=== Terminé ==="
echo "Date: $(date)"
echo "Résultats dans: $OUTPUT_DIR"
