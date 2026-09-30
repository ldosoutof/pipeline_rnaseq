#!/usr/bin/env python3
"""
backfill_db.py
==============
Peuple la base SQLite depuis tous les qc_summary.tsv existants.
Tout se fait dans un seul processus avec une seule connexion SQLite.

Usage :
    python backfill_db.py \
        --prod_root /path/to/runs/ \
        --db        /path/to/rnaseq_qc.db
"""
import argparse
import glob
import sqlite3
import pandas as pd
from pathlib import Path
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("--prod_root", required=True)
parser.add_argument("--db",        required=True)
args = parser.parse_args()

db_path = Path(args.db)
db_path.parent.mkdir(parents=True, exist_ok=True)

# ── Connexion unique ──────────────────────────────────────────────────────────
con = sqlite3.connect(str(db_path), timeout=60)
cur = con.cursor()
cur.execute("PRAGMA journal_mode=WAL")
cur.execute("PRAGMA synchronous=NORMAL")

cur.execute("""CREATE TABLE IF NOT EXISTS qc_metrics (
    id                    INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id             TEXT NOT NULL,
    run_name              TEXT NOT NULL,
    date                  TEXT,
    run_id                TEXT,
    nb_read               REAL,
    mapped                REAL,
    dup                   REAL,
    on_target             REAL,
    nb_event_out          REAL,
    nb_event_out_hyper    REAL,
    nb_event_fraser       REAL,
    nb_event_fraser_hyper REAL,
    HBA1                  REAL,
    HBA2                  REAL,
    HBB                   REAL,
    HBA_total             REAL,
    nb_DI_green_expressed REAL,
    nb_DI_green_total     REAL,
    pct_DI_green          REAL,
    inserted_at           TEXT,
    UNIQUE(sample_id, run_name)
)""")

cur.execute("""CREATE TABLE IF NOT EXISTS qc_warnings (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id   TEXT NOT NULL,
    run_name    TEXT NOT NULL,
    metric      TEXT NOT NULL,
    value       REAL,
    threshold   REAL,
    message     TEXT,
    inserted_at TEXT
)""")

for idx in [
    "CREATE INDEX IF NOT EXISTS idx_qc_run    ON qc_metrics(run_id)",
    "CREATE INDEX IF NOT EXISTS idx_qc_sample ON qc_metrics(sample_id)",
    "CREATE INDEX IF NOT EXISTS idx_qc_date   ON qc_metrics(date)",
    "CREATE INDEX IF NOT EXISTS idx_warn_run  ON qc_warnings(run_name)",
]:
    cur.execute(idx)

con.commit()

# ── Parcourir les runs ────────────────────────────────────────────────────────
runs = sorted(glob.glob(f"{args.prod_root}/20*_RUN*/"))
print(f"[INFO] {len(runs)} run(s) trouvé(s)")

now = datetime.utcnow().isoformat()
ok = skipped = failed = 0

for run in runs:
    run_name = Path(run).name
    qc_file  = Path(run) / "pipeline_v0" / "metrics" / "qc_summary.tsv"
    warn_file = Path(run) / "pipeline_v0" / "metrics" / "qc_warnings.tsv"

    if not qc_file.exists():
        skipped += 1
        continue

    try:
        df = pd.read_csv(qc_file, sep="\t")
        df = df.rename(columns={
            "sample id":             "sample_id",
            "nb read":               "nb_read",
            "on target":             "on_target",
            "nb event out":          "nb_event_out",
            "nb event out hyper":    "nb_event_out_hyper",
            "nb event fraser":       "nb_event_fraser",
            "nb event fraser hyper": "nb_event_fraser_hyper",
            "pct_DI_green_TPM>10":   "pct_DI_green",
        })

        for _, row in df.iterrows():
            cur.execute("""
                INSERT INTO qc_metrics
                    (sample_id, run_name, date, run_id, nb_read, mapped, dup, on_target,
                     nb_event_out, nb_event_out_hyper, nb_event_fraser, nb_event_fraser_hyper,
                     HBA1, HBA2, HBB, HBA_total,
                     nb_DI_green_expressed, nb_DI_green_total, pct_DI_green, inserted_at)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                ON CONFLICT(sample_id, run_name) DO UPDATE SET
                    nb_read=excluded.nb_read, mapped=excluded.mapped,
                    dup=excluded.dup, on_target=excluded.on_target,
                    nb_event_out=excluded.nb_event_out,
                    nb_event_out_hyper=excluded.nb_event_out_hyper,
                    nb_event_fraser=excluded.nb_event_fraser,
                    nb_event_fraser_hyper=excluded.nb_event_fraser_hyper,
                    HBA1=excluded.HBA1, HBA2=excluded.HBA2, HBB=excluded.HBB,
                    HBA_total=excluded.HBA_total,
                    nb_DI_green_expressed=excluded.nb_DI_green_expressed,
                    nb_DI_green_total=excluded.nb_DI_green_total,
                    pct_DI_green=excluded.pct_DI_green,
                    inserted_at=excluded.inserted_at
            """, (
                row.get("sample_id"), run_name, row.get("date"), row.get("run_id"),
                row.get("nb_read"), row.get("mapped"), row.get("dup"), row.get("on_target"),
                row.get("nb_event_out"), row.get("nb_event_out_hyper"),
                row.get("nb_event_fraser"), row.get("nb_event_fraser_hyper"),
                row.get("HBA1"), row.get("HBA2"), row.get("HBB"), row.get("HBA_total"),
                row.get("nb_DI_green_expressed"), row.get("nb_DI_green_total"),
                row.get("pct_DI_green"), now
            ))

        # Warnings
        if warn_file.exists():
            wdf = pd.read_csv(warn_file, sep="\t")
            wdf = wdf.rename(columns={"sample id": "sample_id"})
            cur.execute("DELETE FROM qc_warnings WHERE run_name = ?", (run_name,))
            for _, row in wdf.iterrows():
                cur.execute("""
                    INSERT INTO qc_warnings
                        (sample_id, run_name, metric, value, threshold, message, inserted_at)
                    VALUES (?,?,?,?,?,?,?)
                """, (
                    row.get("sample_id", row.get("sample id")),
                    run_name, row.get("metric"), row.get("value"),
                    row.get("threshold"), row.get("message"), now
                ))

        con.commit()
        print(f"[OK] {run_name} ({len(df)} samples)")
        ok += 1

    except Exception as e:
        print(f"[ERR] {run_name}: {e}")
        failed += 1

con.close()
print(f"\n✅ {ok} runs insérés, {failed} erreurs, {skipped} sans qc_summary.tsv")
print(f"   Base : {db_path}")
