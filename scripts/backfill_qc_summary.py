#!/usr/bin/env python3
"""
backfill_qc_summary.py
======================
Reconstruit un qc_summary.tsv pour les anciens runs traités avant
l'introduction de metrics/ (avant RUN41 environ).

Le script cherche les fichiers aux deux emplacements possibles :
  - Ancien pipeline : pipeline_v0/dup/, pipeline_v0/coverage/
  - Nouveau pipeline : pipeline_v0/metrics/dup/, pipeline_v0/metrics/coverage/

Usage :
  # Un seul run
  python backfill_qc_summary.py --run_path /path/to/20240523_RUN25_...

  # Tous les runs d'un dossier (ignore ceux qui ont déjà un qc_summary.tsv)
  python backfill_qc_summary.py --prod_root /path/to/runs/

  # Forcer la reconstruction même si qc_summary.tsv existe déjà
  python backfill_qc_summary.py --prod_root /path/to/runs/ --force

  # Fournir un fichier TPM pour inclure HBA_total (optionnel)
  python backfill_qc_summary.py --run_path /path/to/run/ --tpm_file /path/to/matrice_gene_tpm_gene.tsv
"""

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import pandas as pd


# ── Helpers ───────────────────────────────────────────────────────────────────

def _find(run_path, *relative_paths):
    """
    Cherche un fichier aux deux emplacements possibles (ancien et nouveau pipeline).
    Retourne le premier trouvé, ou None.
    """
    for rel in relative_paths:
        p = Path(run_path) / "pipeline_v0" / rel
        if p.exists():
            return p
    return None


def _glob_first(pattern):
    hits = sorted(glob.glob(pattern))
    return hits[0] if hits else None


def count_events(filepath):
    """Nombre de lignes de données dans un .tab (header exclu). None si absent."""
    if not filepath or not Path(filepath).exists():
        return None
    try:
        lines = Path(filepath).read_text(errors="replace").splitlines()
        return max(0, len(lines) - 1)
    except Exception:
        return None


def parse_dup_txt(dup_file):
    """
    Parse un fichier Picard MarkDuplicates.
    Retourne (nb_read, mapped, dup_pct) ou (None, None, None).
    Supporte les deux formats : avec et sans ligne d'en-tête 'rnaseq-capture'.
    """
    try:
        for line in Path(dup_file).read_text(errors="replace").splitlines():
            if line.startswith("rnaseq-capture"):
                parts = line.strip().split()
                nb_read = int(parts[1]) + 2 * int(parts[2])
                mapped  = 2 * int(parts[2])
                dup_pct = round(float(parts[8]) * 100, 2)
                return nb_read, mapped, dup_pct
    except Exception as e:
        print(f"[WARN] Impossible de parser {dup_file}: {e}", file=sys.stderr)
    return None, None, None


# ── Reconstruction d'un run ───────────────────────────────────────────────────

def build_qc_summary(run_path, tpm_file=None, target_genes=None):
    """
    Reconstruit le qc_summary.tsv pour un run donné.
    Retourne un DataFrame ou None si aucune donnée trouvée.
    """
    run_path = Path(run_path)
    run_name = run_path.name
    date     = run_name.split("_")[0]
    m        = re.search(r"RUN\d+", run_name, re.IGNORECASE)
    run_id   = m.group(0) if m else None

    if target_genes is None:
        target_genes = ["HBA1", "HBA2", "HBB"]

    # Charger la matrice TPM si fournie
    tpm_df = None
    if tpm_file and Path(tpm_file).exists():
        try:
            tpm_df = pd.read_csv(tpm_file, sep="\t", index_col=0)
        except Exception as e:
            print(f"[WARN] TPM file illisible : {e}", file=sys.stderr)

    # Chercher les fichiers dup — ancien chemin : dup/, nouveau : metrics/dup/
    dup_pattern_old = str(run_path / "pipeline_v0" / "dup" / "*" / "*_dup.txt")
    dup_pattern_new = str(run_path / "pipeline_v0" / "metrics" / "dup" / "*" / "*_dup.txt")
    dup_files = sorted(glob.glob(dup_pattern_new) or glob.glob(dup_pattern_old))

    if not dup_files:
        print(f"[SKIP] Aucun fichier dup trouvé pour {run_name}", file=sys.stderr)
        return None

    rows = []
    for dup_file in dup_files:
        sample_id = Path(dup_file).stem.replace("_dup", "")
        short_id  = re.split(r"[-_]", sample_id)[0]
        row       = {
            "sample id": sample_id,
            "run_name":  run_name,
            "date":      date,
            "run_id":    run_id,
        }

        # Duplication
        nb_read, mapped, dup_pct = parse_dup_txt(dup_file)
        row["nb read"] = nb_read
        row["mapped"]  = mapped
        row["dup"]     = dup_pct

        # On-target — ancien : coverage/, nouveau : metrics/coverage/
        on_target_file = _find(
            run_path,
            f"coverage/{sample_id}/{sample_id}_on_target.txt",
            f"metrics/coverage/{sample_id}/{sample_id}_on_target.txt",
        )
        try:
            row["on target"] = int(on_target_file.read_text().strip()) if on_target_file else None
        except Exception:
            row["on target"] = None

        # OUTRIDER / FRASER normaux
        row["nb event out"]   = count_events(_find(
            run_path,
            f"outrider/filesbysample/{short_id}.outrider.tab",
        ))
        row["nb event fraser"] = count_events(_find(
            run_path,
            f"fraser/filesbysample/{short_id}.fraser.tab",
        ))

        # OUTRIDER / FRASER hyper (absents pour les anciens runs — None)
        row["nb event out hyper"]    = count_events(_find(
            run_path,
            f"outrider_hyper/filesbysample/{short_id}.outrider.tab",
        ))
        row["nb event fraser hyper"] = count_events(_find(
            run_path,
            f"fraser_hyper/filesbysample/{short_id}.fraser.tab",
        ))

        # TPM hémoglobine
        if tpm_df is not None:
            tpm_sample_id = sample_id.replace("-", ".").replace("_", ".")
            for gene in target_genes:
                try:
                    row[gene] = tpm_df.at[gene, tpm_sample_id]
                except KeyError:
                    row[gene] = None
            row["HBA_total"] = sum(row.get(g, 0) or 0 for g in target_genes)
        else:
            for gene in target_genes:
                row[gene] = None
            row["HBA_total"] = None

        rows.append(row)

    if not rows:
        return None

    # Ordre des colonnes cohérent avec recup_metrics.py
    col_order = [
        "sample id", "run_name", "date", "run_id",
        "nb read", "mapped", "dup", "on target",
        "nb event out", "nb event out hyper",
        "nb event fraser", "nb event fraser hyper",
        *target_genes, "HBA_total",
    ]
    df = pd.DataFrame(rows)
    col_order = [c for c in col_order if c in df.columns]
    return df[col_order]


# ── Main ─────────────────────────────────────────────────────────────────────

def process_run(run_path, tpm_file=None, force=False):
    run_path = Path(run_path)
    run_name = run_path.name

    # Chercher aux deux emplacements possibles
    out_new = run_path / "pipeline_v0" / "metrics" / "qc_summary.tsv"
    out_old = run_path / "pipeline_v0" / "qc_summary.tsv"

    if out_new.exists() and not force:
        print(f"[SKIP] {run_name} — qc_summary.tsv déjà présent (--force pour reconstruire)")
        return
    if out_old.exists() and not force:
        print(f"[SKIP] {run_name} — qc_summary.tsv déjà présent (ancien emplacement)")
        return

    # Chercher la matrice TPM si non fournie
    if tpm_file is None:
        tpm_candidates = [
            run_path / "pipeline_v0" / "metrics" / "matrice_gene_tpm_gene.tsv",
            run_path / "pipeline_v0" / "kallisto_bed" / "matrice_gene_tpm_gene.tsv",
        ]
        tpm_file = next((str(p) for p in tpm_candidates if p.exists()), None)

    df = build_qc_summary(run_path, tpm_file=tpm_file)
    if df is None:
        return

    # Écrire dans metrics/ (nouveau chemin) en créant le dossier si nécessaire
    out_new.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_new, sep="\t", index=False)
    print(f"[OK] {run_name} — {len(df)} échantillon(s) → {out_new}")


def main():
    p = argparse.ArgumentParser(
        description="Reconstruit qc_summary.tsv pour les anciens runs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--run_path",  help="Chemin vers un run unique")
    p.add_argument("--prod_root", help="Répertoire racine (traite tous les runs)")
    p.add_argument("--tpm_file",  default=None,
                   help="Matrice TPM gène (optionnel — pour HBA_total)")
    p.add_argument("--force",     action="store_true",
                   help="Reconstruire même si qc_summary.tsv existe déjà")
    args = p.parse_args()

    if args.prod_root:
        runs = sorted(glob.glob(os.path.join(args.prod_root, "20*_RUN*")))
        print(f"[INFO] {len(runs)} run(s) trouvé(s) dans {args.prod_root}")
        for run in runs:
            process_run(run, tpm_file=args.tpm_file, force=args.force)
    elif args.run_path:
        process_run(args.run_path, tpm_file=args.tpm_file, force=args.force)
    else:
        print("[ERROR] Fournir --run_path ou --prod_root", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
