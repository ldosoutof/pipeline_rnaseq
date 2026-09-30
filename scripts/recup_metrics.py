import pandas as pd
import glob
import os
import argparse
import re

# === PARSE ARGUMENTS ===
parser = argparse.ArgumentParser(description="Generate QC summary for a given run")
parser.add_argument("--run_path",          required=True,  help="Full path to the run folder (e.g. RUN29)")
parser.add_argument("--tpm_file",          required=True,  help="Path to the TPM matrix file")
parser.add_argument("--di_file",           required=True,  help="Path to the panelapp di file")
parser.add_argument("--output",            default="qc_summary.tsv", help="Output filename")
parser.add_argument("--warnings_output",   default="",     help="Output path for the warnings TSV (empty = same dir as --output)")
# Seuils d'alerte (modifiables via le config Snakemake)
parser.add_argument("--max_dup",           type=float, default=40.0,  help="Seuil taux de duplication (%%) — défaut 40")
parser.add_argument("--max_out_hyper",     type=int,   default=20,    help="Seuil nb événements OUTRIDER hyper — défaut 20")
parser.add_argument("--max_fraser_hyper",  type=int,   default=20,    help="Seuil nb événements FRASER hyper — défaut 20")
parser.add_argument("--max_hba_total",     type=float, default=5000.0,help="Seuil TPM HBA_total (contamination érythrocytaire) — défaut 5000")
args = parser.parse_args()

# === CONFIGURATION ===
run_path      = args.run_path
tpm_file      = args.tpm_file
output_file   = args.output
di_green_file = args.di_file
target_genes  = ["HBA1", "HBA2", "HBB"]

# Seuils d'alerte
THRESH_DUP          = args.max_dup
THRESH_OUT_HYPER    = args.max_out_hyper
THRESH_FRASER_HYPER = args.max_fraser_hyper
THRESH_HBA_TOTAL    = args.max_hba_total

# Chemin du fichier warnings (par défaut à côté du qc_summary)
warnings_file = (args.warnings_output if args.warnings_output
                 else os.path.join(os.path.dirname(output_file), "qc_warnings.tsv"))

# === Extract run info ===
run_name = os.path.basename(run_path)   # e.g. '20250526_RUN36_NextSeq_High_16RNASEQ'
print(run_path)
print(run_name)
date   = run_name.split('_')[0]         # e.g. '20250526'
match  = re.search(r'RUN\d+', run_name, re.IGNORECASE)
run_id = match.group(0) if match else None

# === Load DI_green genes ===
try:
    di_df    = pd.read_csv(di_green_file, sep="\t")
    di_genes = set(di_df.loc[di_df["Entity type"] == "gene", "Gene Symbol"].dropna())
    print(f"✅ Loaded {len(di_genes)} DI_green genes")
except Exception as e:
    print(f"⚠️ Could not load DI_green file: {e}")
    di_genes = set()

# === LOAD TPM MATRIX ===
try:
    tpm_df = pd.read_csv(tpm_file, sep="\t", index_col=0)
except Exception as e:
    print(f"❌ Failed to load TPM matrix: {e}")
    exit(1)

# ── Helper ────────────────────────────────────────────────────────────────────
def count_events(filepath, sample_label):
    """Return the number of data lines in a tab file (header excluded).
    Returns None if the file cannot be read, 0 if it is empty."""
    try:
        with open(filepath) as f:
            lines = f.readlines()
        return max(0, len(lines) - 1)
    except Exception as e:
        print(f"⚠️ Could not count events in {filepath} for {sample_label}: {e}")
        return None

# === INIT FINAL TABLE ===
summary_rows = []

# === LOOP THROUGH SAMPLES BASED ON DUP FILES ===
# Chemin mis à jour : metrics/dup/ (remaniement pipeline)
dup_files = glob.glob(f"{run_path}/pipeline_v0/metrics/dup/*/*_dup.txt")
for dup_file in dup_files:
    sample_id  = os.path.basename(dup_file).replace("_dup.txt", "")
    sample_row = {"sample id": sample_id}

    # Run info
    sample_row["run_name"] = run_name
    sample_row["date"]     = date
    sample_row["run_id"]   = run_id

    # === DUPLICATION METRICS ===
    try:
        with open(dup_file) as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith("rnaseq-capture"):
                parts = line.strip().split()
                sample_row["nb read"] = int(parts[1]) + 2 * int(parts[2])
                sample_row["mapped"]  = 2 * int(parts[2])
                sample_row["dup"]     = round(float(parts[8]) * 100, 2)
                break
    except Exception as e:
        print(f"⚠️ Could not parse dup file for {sample_id}: {e}")
        sample_row["nb read"] = sample_row["mapped"] = sample_row["dup"] = None

    # === COVERAGE METRIC ===
    # Chemin mis à jour : metrics/coverage/
    on_target_file = f"{run_path}/pipeline_v0/metrics/coverage/{sample_id}/{sample_id}_on_target.txt"
    try:
        with open(on_target_file) as f:
            sample_row["on target"] = int(f.readline().strip())
    except Exception as e:
        print(f"⚠️ Could not read on-target file for {sample_id}: {e}")
        sample_row["on target"] = None

    # === SHORT ID (7 premiers caractères) ===
    short_id = re.split('[-_]', sample_id)[0]

    # === OUTRIDER EVENT COUNT (pipeline normal) ===
    sample_row["nb event out"] = count_events(
        f"{run_path}/pipeline_v0/outrider/filesbysample/{short_id}.outrider.tab",
        short_id
    )

    # === FRASER EVENT COUNT (pipeline normal) ===
    sample_row["nb event fraser"] = count_events(
        f"{run_path}/pipeline_v0/fraser/filesbysample/{short_id}.fraser.tab",
        short_id
    )

    # === OUTRIDER HYPER EVENT COUNT ===
    sample_row["nb event out hyper"] = count_events(
        f"{run_path}/pipeline_v0/outrider_hyper/filesbysample/{short_id}.outrider.tab",
        short_id
    )

    # === FRASER HYPER EVENT COUNT ===
    sample_row["nb event fraser hyper"] = count_events(
        f"{run_path}/pipeline_v0/fraser_hyper/filesbysample/{short_id}.fraser.tab",
        short_id
    )

    # === TPM VALUES FOR TARGET GENES ===
    tpm_sample_id = sample_id.replace("-", ".").replace("_", ".")
    for gene in target_genes:
        try:
            sample_row[gene] = tpm_df.at[gene, tpm_sample_id]
        except KeyError:
            sample_row[gene] = None

    # === SUM OF HBA1 + HBA2 + HBB ===
    sample_row["HBA_total"] = sum(sample_row.get(g, 0) or 0 for g in ["HBA1", "HBA2", "HBB"])

    # === Percentage of DI_green genes with TPM > 10 ===
    try:
        if di_genes:
            if tpm_sample_id in tpm_df.columns:
                sample_tpms   = tpm_df[tpm_sample_id]
                intersect_genes = di_genes.intersection(sample_tpms.index)
                expressed     = sample_tpms.loc[list(intersect_genes)] > 10
                nb_expressed  = expressed.sum()
                nb_total      = len(intersect_genes)
                pct_expressed = (nb_expressed / nb_total * 100) if nb_total > 0 else 0
                sample_row["nb_DI_green_expressed"]  = nb_expressed
                sample_row["nb_DI_green_total"]      = nb_total
                sample_row["pct_DI_green_TPM>10"]    = round(pct_expressed, 2)
            else:
                sample_row["nb_DI_green_expressed"] = sample_row["nb_DI_green_total"] = sample_row["pct_DI_green_TPM>10"] = None
        else:
            sample_row["nb_DI_green_expressed"] = sample_row["nb_DI_green_total"] = sample_row["pct_DI_green_TPM>10"] = None
    except Exception as e:
        print(f"⚠️ Error calculating DI_green TPM > 10 for {sample_id}: {e}")
        sample_row["nb_DI_green_expressed"] = sample_row["nb_DI_green_total"] = sample_row["pct_DI_green_TPM>10"] = None

    summary_rows.append(sample_row)

# === CREATE AND SAVE FINAL TABLE ===
summary_df = pd.DataFrame(summary_rows)

# Ordre des colonnes
col_order = [
    "sample id", "run_name", "date", "run_id",
    "nb read", "mapped", "dup", "on target",
    "nb event out", "nb event out hyper",
    "nb event fraser", "nb event fraser hyper",
    *target_genes, "HBA_total",
    "nb_DI_green_expressed", "nb_DI_green_total", "pct_DI_green_TPM>10",
]
# Garder uniquement les colonnes présentes (sécurité)
col_order = [c for c in col_order if c in summary_df.columns]
summary_df = summary_df[col_order]

summary_df.to_csv(output_file, sep="\t", index=False)
print(f"✅ QC summary saved to: {output_file}")

# === GENERATE WARNINGS FILE ===
warning_rows = []

for _, row in summary_df.iterrows():
    sample = row.get("sample id", "unknown")

    # Taux de duplication trop élevé
    dup = row.get("dup")
    if dup is not None and dup > THRESH_DUP:
        warning_rows.append({
            "sample id": sample,
            "metric":    "dup",
            "value":     dup,
            "threshold": THRESH_DUP,
            "message":   f"Taux de duplication élevé : {dup:.1f}% (seuil {THRESH_DUP}%)",
        })

    # Nb événements OUTRIDER hyper trop élevé
    out_hyper = row.get("nb event out hyper")
    if out_hyper is not None and out_hyper > THRESH_OUT_HYPER:
        warning_rows.append({
            "sample id": sample,
            "metric":    "nb event out hyper",
            "value":     out_hyper,
            "threshold": THRESH_OUT_HYPER,
            "message":   f"Nombre d'événements OUTRIDER hyper élevé : {out_hyper} (seuil {THRESH_OUT_HYPER})",
        })

    # Nb événements FRASER hyper trop élevé
    fraser_hyper = row.get("nb event fraser hyper")
    if fraser_hyper is not None and fraser_hyper > THRESH_FRASER_HYPER:
        warning_rows.append({
            "sample id": sample,
            "metric":    "nb event fraser hyper",
            "value":     fraser_hyper,
            "threshold": THRESH_FRASER_HYPER,
            "message":   f"Nombre d'événements FRASER hyper élevé : {fraser_hyper} (seuil {THRESH_FRASER_HYPER})",
        })

    # HBA_total trop élevé (contamination érythrocytaire)
    hba = row.get("HBA_total")
    if hba is not None and hba > THRESH_HBA_TOTAL:
        warning_rows.append({
            "sample id": sample,
            "metric":    "HBA_total",
            "value":     round(hba, 1),
            "threshold": THRESH_HBA_TOTAL,
            "message":   f"HBA_total élevé : {hba:.0f} TPM (seuil {THRESH_HBA_TOTAL}) — contamination érythrocytaire ?",
        })

warnings_df = pd.DataFrame(warning_rows, columns=["sample id", "metric", "value", "threshold", "message"])
warnings_df.to_csv(warnings_file, sep="\t", index=False)

if len(warnings_df) == 0:
    print(f"✅ Aucun warning — fichier vide écrit : {warnings_file}")
else:
    print(f"⚠️  {len(warnings_df)} warning(s) écrits dans : {warnings_file}")
    for _, w in warnings_df.iterrows():
        print(f"   [{w['sample id']}] {w['message']}")
