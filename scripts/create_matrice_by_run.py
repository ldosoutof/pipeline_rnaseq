#!/usr/bin/env python3
"""
create_matrice_by_run.py
========================
Assemble per-sample HTSeq count files into a single gene × sample matrix TSV,
ready for OUTRIDER. Scans all historical runs under base_dir.

Usage:
    python create_matrice_by_run.py <base_dir> <output_dir> <blacklist_file> <current_run_tag>
"""
import os
import sys
import re
from datetime import date
from collections import defaultdict

import csv
import pandas as pd
pd.set_option('future.no_silent_downcasting', True)

if len(sys.argv) != 5:
    print("Usage: python create_matrice_by_run.py <base_dir> <output_dir> <blacklist_file> <current_run_tag>")
    sys.exit(1)

base_dir        = sys.argv[1]
out_dir         = sys.argv[2]
blacklist_file  = sys.argv[3]
current_run_tag = sys.argv[4]   # ex: 20260610_RUN50_NextSeq_High_16RNASEQ

# ── Blacklist ─────────────────────────────────────────────────────────────────
blacklist = set()
if os.path.exists(blacklist_file):
    # Format TSV unifié : sample_id <TAB> tool <TAB> reason
    # tool : outrider | fraser | pca | all
    with open(blacklist_file) as _bl:
        reader = csv.DictReader(_bl, delimiter="\t")
        if reader.fieldnames and "sample_id" in reader.fieldnames:
            # Nouveau format unifié
            for row in reader:
                t   = row.get("tool", "all").strip().lower()
                sid = row.get("sample_id", "").strip()
                if sid and not sid.startswith("#") and t in ("outrider", "all"):
                    blacklist.add(sid)
        else:
            # Fallback : ancien format (sample_id reason sur chaque ligne)
            _bl.seek(0)
            for line in _bl:
                line = line.strip()
                if line and not line.startswith("#"):
                    blacklist.add(line.split()[0])
    print(f"[blacklist] {len(blacklist)} samples exclus")
else:
    print("⚠️  Pas de blacklist trouvée, on continue sans filtrage.")

# ── Collect count files ───────────────────────────────────────────────────────
all_counts = None

runs = [d for d in os.listdir(base_dir) if d.startswith("20") and "RUN17" not in d]

# Suivre les échantillons du run courant trouvés
current_run_samples_found = set()

for run in sorted(runs):
    run_number_str = [s for s in run.split("_") if s.startswith("RUN")]
    if not run_number_str:
        continue
    try:
        run_num = int(run_number_str[0].replace("RUN", ""))
    except ValueError:
        continue

    htseq_subfolder = "htseqStrand" if run_num <= 24 else "htseq"
    htseq_dir = os.path.join(base_dir, run, "pipeline_v0", htseq_subfolder)
    if not os.path.isdir(htseq_dir):
        continue

    for root, _, files in os.walk(htseq_dir):
        for file in files:
            if not file.endswith("_gene_counts.txt"):
                continue
            # Filtre sang (MOINS ou PUROMOINS dans le nom)
            if not any(kw in file.upper() for kw in ("MOINS", "PUROMOINS")):
                continue

            sample_name = file.replace("_gene_counts.txt", "")
            short_id = re.split(r"[-_]", sample_name)[0]

            if any(sample_name.startswith(prefix) for prefix in blacklist):
                print(f"[skip] blacklisté : {sample_name}")
                continue

            if run == current_run_tag:
                current_run_samples_found.add(short_id)

            file_path = os.path.join(root, file)
            df = pd.read_csv(file_path, sep="\t", header=None,
                             names=["ENSG", sample_name])
            df = df[df["ENSG"].str.startswith("ENSG")]

            if all_counts is not None:
                all_counts = pd.merge(all_counts, df, on="ENSG", how="outer")
            else:
                all_counts = df

# ── Validation run courant ────────────────────────────────────────────────────
# Reconstruire la liste attendue depuis le dossier htseq du run courant
current_htseq_dir = None
for run in runs:
    if run == current_run_tag:
        run_number_str = [s for s in run.split("_") if s.startswith("RUN")]
        run_num = int(run_number_str[0].replace("RUN", "")) if run_number_str else 0
        subfolder = "htseqStrand" if run_num <= 24 else "htseq"
        current_htseq_dir = os.path.join(base_dir, run, "pipeline_v0", subfolder)
        break

if current_htseq_dir and os.path.isdir(current_htseq_dir):
    # Échantillons attendus = tous les *_gene_counts.txt sang du run courant
    expected = set()
    for root, _, files in os.walk(current_htseq_dir):
        for f in files:
            if f.endswith("_gene_counts.txt") and any(
                kw in f.upper() for kw in ("MOINS", "PUROMOINS")
            ):
                sample_name = f.replace("_gene_counts.txt", "")
                short = re.split(r'[-_]', sample_name)[0]
                if short not in blacklist:
                    expected.add(short)

    missing_current = expected - current_run_samples_found
    if missing_current:
        print(
            f"[ERROR] {len(missing_current)} échantillon(s) du run courant ({current_run_tag}) "
            f"absents de la matrice HTSeq :\n"
            f"  {sorted(missing_current)}\n"
            f"  → Vérifiez que htseq_gene a tourné pour ces échantillons "
            f"(logs sous log/htseq/<sample>/)",
            file=sys.stderr
        )
        sys.exit(2)
    else:
        print(f"[OK] Tous les échantillons du run courant sont présents dans la matrice HTSeq "
              f"({len(expected)} échantillons)")
elif current_run_tag:
    print(f"[WARN] Dossier HTSeq introuvable pour le run courant ({current_run_tag}) "
          f"— validation ignorée", file=sys.stderr)

# ── Final cleanup ─────────────────────────────────────────────────────────────
if all_counts is None:
    print("[ERROR] Aucun fichier de comptage trouvé.", file=sys.stderr)
    sys.exit(1)

all_counts.fillna(0, inplace=True)
all_counts.set_index("ENSG", inplace=True)

# Garder uniquement les échantillons commençant par "1" ou "2", exclure POLYA
all_counts = all_counts.loc[:, all_counts.columns.str.startswith(("1", "2"))]
all_counts = all_counts.loc[:, ~all_counts.columns.str.contains("POLYA")]

# ── Dédupliquer : garder le run le plus récent par sample ID ──────────────────
simplified_names = all_counts.columns.to_series().apply(lambda x: x.split("-")[0])
duplicate_ids = simplified_names[simplified_names.duplicated()].unique().tolist()
name_map = defaultdict(list)
for full_name in all_counts.columns:
    if "-" in full_name:
        short_name = full_name.split("-")[0]
    elif "_" in full_name:
        short_name = full_name.split("_")[0]
    else:
        short_name = full_name
    name_map[short_name].append(full_name)

def get_run_number(colname):
    match = re.search(r"RUN(\d+)", colname)
    return int(match.group(1)) if match else -1

to_drop = []
for dup_id in duplicate_ids:
    cols = name_map[dup_id]
    best_col = max(cols, key=get_run_number)
    drop_cols = [c for c in cols if c != best_col]
    to_drop.extend(drop_cols)
    print(f"[dedup] Gardé {best_col} pour {dup_id}, supprimé {drop_cols}")

all_counts = all_counts.drop(columns=to_drop)

# ── Write matrix ──────────────────────────────────────────────────────────────
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "matrice.txt")
all_counts.to_csv(out_path, sep="\t")
print(f"✅ Matrice HTSeq sauvegardée : {out_path}  {all_counts.shape[0]} gènes × {all_counts.shape[1]} échantillons")
