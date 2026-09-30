#!/usr/bin/env python3
"""
update_db.py
============
Insère ou met à jour les métriques QC d'un run dans la base SQLite.
Appelé automatiquement par la règle generate_metrics à la fin de chaque run.

Usage :
    python update_db.py --qc_summary /path/to/qc_summary.tsv \
                        --warnings   /path/to/qc_warnings.tsv \
                        --db         /path/to/rnaseq_qc.db
"""

import argparse
import sqlite3
import pandas as pd
from pathlib import Path
from datetime import datetime

# ── Arguments ─────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Update QC SQLite database")
parser.add_argument("--qc_summary", required=True)
parser.add_argument("--warnings",   default="")
parser.add_argument("--db",         required=True)
args = parser.parse_args()

db_path = Path(args.db)
db_path.parent.mkdir(parents=True, exist_ok=True)

# ── Connexion ─────────────────────────────────────────────────────────────────
con = sqlite3.connect(db_path, timeout=30)
cur = con.cursor()
# WAL mode — permet les lectures concurrentes sans bloquer les écritures
cur.execute("PRAGMA journal_mode=WAL")
cur.execute("PRAGMA synchronous=NORMAL")

# ── Schéma ────────────────────────────────────────────────────────────────────
for stmt in [
    """CREATE TABLE IF NOT EXISTS qc_metrics (
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
    )""",
    """CREATE TABLE IF NOT EXISTS qc_warnings (
        id         INTEGER PRIMARY KEY AUTOINCREMENT,
        sample_id  TEXT NOT NULL,
        run_name   TEXT NOT NULL,
        metric     TEXT NOT NULL,
        value      REAL,
        threshold  REAL,
        message    TEXT,
        inserted_at TEXT
    )""",
    "CREATE INDEX IF NOT EXISTS idx_qc_run     ON qc_metrics(run_id)",
    "CREATE INDEX IF NOT EXISTS idx_qc_sample  ON qc_metrics(sample_id)",
    "CREATE INDEX IF NOT EXISTS idx_qc_date    ON qc_metrics(date)",
    "CREATE INDEX IF NOT EXISTS idx_warn_run   ON qc_warnings(run_name)",
]:
    cur.execute(stmt)

# ── Insérer les métriques ─────────────────────────────────────────────────────
now = datetime.utcnow().isoformat()

df = pd.read_csv(args.qc_summary, sep="\t")
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

inserted = updated = 0
for _, row in df.iterrows():
    cur.execute("""
        INSERT INTO qc_metrics
            (sample_id, run_name, date, run_id, nb_read, mapped, dup, on_target,
             nb_event_out, nb_event_out_hyper, nb_event_fraser, nb_event_fraser_hyper,
             HBA1, HBA2, HBB, HBA_total,
             nb_DI_green_expressed, nb_DI_green_total, pct_DI_green, inserted_at)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(sample_id, run_name) DO UPDATE SET
            nb_read               = excluded.nb_read,
            mapped                = excluded.mapped,
            dup                   = excluded.dup,
            on_target             = excluded.on_target,
            nb_event_out          = excluded.nb_event_out,
            nb_event_out_hyper    = excluded.nb_event_out_hyper,
            nb_event_fraser       = excluded.nb_event_fraser,
            nb_event_fraser_hyper = excluded.nb_event_fraser_hyper,
            HBA1                  = excluded.HBA1,
            HBA2                  = excluded.HBA2,
            HBB                   = excluded.HBB,
            HBA_total             = excluded.HBA_total,
            nb_DI_green_expressed = excluded.nb_DI_green_expressed,
            nb_DI_green_total     = excluded.nb_DI_green_total,
            pct_DI_green          = excluded.pct_DI_green,
            inserted_at           = excluded.inserted_at
    """, (
        row.get("sample_id"),   row.get("run_name"),  row.get("date"),
        row.get("run_id"),      row.get("nb_read"),   row.get("mapped"),
        row.get("dup"),         row.get("on_target"),
        row.get("nb_event_out"),          row.get("nb_event_out_hyper"),
        row.get("nb_event_fraser"),       row.get("nb_event_fraser_hyper"),
        row.get("HBA1"),  row.get("HBA2"),  row.get("HBB"),  row.get("HBA_total"),
        row.get("nb_DI_green_expressed"), row.get("nb_DI_green_total"),
        row.get("pct_DI_green"),          now
    ))
    if cur.rowcount == 1:
        inserted += 1
    else:
        updated += 1

# ── Insérer les warnings ─────────────────────────────────────────────────────
if args.warnings and Path(args.warnings).exists():
    wdf = pd.read_csv(args.warnings, sep="\t")
    wdf = wdf.rename(columns={"sample id": "sample_id"})
    run_name = df["run_name"].iloc[0] if len(df) > 0 else ""
    # Supprimer les anciens warnings de ce run
    cur.execute("DELETE FROM qc_warnings WHERE run_name = ?", (run_name,))
    for _, row in wdf.iterrows():
        cur.execute("""
            INSERT INTO qc_warnings
                (sample_id, run_name, metric, value, threshold, message, inserted_at)
            VALUES (?,?,?,?,?,?,?)
        """, (
            row.get("sample_id", row.get("sample id")),
            run_name,
            row.get("metric"), row.get("value"), row.get("threshold"),
            row.get("message"), now
        ))

con.commit()
con.close()

print(f"✅ DB mise à jour : {db_path}")
print(f"   {inserted} insérés, {updated} mis à jour")
