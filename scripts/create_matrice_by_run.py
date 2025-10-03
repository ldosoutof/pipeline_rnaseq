#!/usr/bin/env python3
import os
import pandas as pd
pd.set_option('future.no_silent_downcasting', True)
import sys
from datetime import date
from collections import defaultdict
import re

"""
Usage:
    python create_matrice.py <base_dir> <output_dir> <blacklist_file>
"""

if len(sys.argv) != 4:
    print("Usage: python create_matrice.py <base_dir> <output_dir> <blacklist_file>")
    sys.exit(1)

base_dir = sys.argv[1]
out_dir = sys.argv[2]
blacklist_file = sys.argv[3]

# ---------------- Load blacklist ----------------
blacklist = set()
if os.path.exists(blacklist_file):
    bl = pd.read_csv(blacklist_file, sep=r"\s+", header=None, names=["sample", "reason"])
    blacklist = set(bl["sample"])
    print(f"Loaded {len(blacklist)} blacklisted samples")
else:
    print("⚠️ No blacklist file found, continuing without filtering.")

# ---------------- Collect counts ----------------
all_counts = None

runs = [d for d in os.listdir(base_dir) if d.startswith("20") and "RUN17" not in d]
for run in runs:
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
            if file.endswith("_gene_counts.txt") and "MOINS" in file:
                sample_name = file.replace("_gene_counts.txt", "")

                if any(sample_name.startswith(prefix) for prefix in blacklist):
                    print(f"Skipping blacklisted sample: {sample_name}")
                    continue

                file_path = os.path.join(root, file)

                df = pd.read_csv(file_path, sep="\t", header=None,
                                 names=["ENSG", sample_name])
                df = df[df["ENSG"].str.startswith("ENSG")]

                if all_counts is not None:
                    all_counts = pd.merge(all_counts, df, on="ENSG", how="outer")
                else:
                    all_counts = df

# ---------------- Final cleanup ----------------
if all_counts is None:
    print("No matching count files found.")
    sys.exit(0)

all_counts.fillna(0, inplace=True)
all_counts.set_index("ENSG", inplace=True)

# Keep only samples starting with "2" and remove POLYA
all_counts = all_counts.loc[:, all_counts.columns.str.startswith(("1","2"))]
all_counts = all_counts.loc[:, ~all_counts.columns.str.contains("POLYA")]

# ---------------- Resolve duplicate sample IDs ----------------
# Simplified ID = prefix before first "-"
simplified_names = all_counts.columns.to_series().apply(lambda x: x.split("-")[0])
duplicate_ids = simplified_names[simplified_names.duplicated()].unique().tolist()
name_map = defaultdict(list)
for full_name in all_counts.columns:
    #short_name = full_name.split("-")[0]
    if "-" in full_name:
        short_name = full_name.split("-")[0]
    elif "_" in full_name:
        short_name = full_name.split("_")[0]
    else:
        short_name = full_name  # no dash or underscore
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
    print(f"Keeping {best_col} for {dup_id}, dropping {drop_cols}")

all_counts = all_counts.drop(columns=to_drop)

# ---------------- Write dated matrix ----------------
today = date.today().strftime("%Y%m%d")
#file_name = f"{today}_matrice.txt"
file_name = f"matrice.txt"
os.makedirs(out_dir, exist_ok=True)

out_path = os.path.join(out_dir, file_name)
all_counts.to_csv(out_path, sep="\t")

print(f"✅ Matrix saved: {out_path} with shape {all_counts.shape}")

