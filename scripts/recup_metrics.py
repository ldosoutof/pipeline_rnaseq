import pandas as pd
import glob
import os
import argparse
import re

# === PARSE ARGUMENTS ===
parser = argparse.ArgumentParser(description="Generate QC summary for a given run")
parser.add_argument("--run_path", required=True, help="Full path to the run folder (e.g. RUN29)")
parser.add_argument("--tpm_file", required=True, help="Path to the TPM matrix file")
parser.add_argument("--di_file", required=True, help="Path to the panelapp di file")
parser.add_argument("--output", default="qc_summary.tsv", help="Output filename")
args = parser.parse_args()

# === CONFIGURATION ===
run_path = args.run_path
tpm_file = args.tpm_file
output_file = args.output
di_green_file = args.di_file
target_genes = ["HBA1", "HBA2", "HBB"]
#target_genes = ["HBA1", "HBA2", "HBB", "SRSF6", "SRSF3"]

# === Extract run info ===
run_name = os.path.basename(run_path)  # e.g. '20250526_RUN36_NextSeq_High_16RNASEQ'
print(run_path)
print(run_name)
date = run_name.split('_')[0]          # e.g. '20250526'

match = re.search(r'RUN\d+', run_name, re.IGNORECASE)
run_id = match.group(0) if match else None

# === Load DI_green genes ===
#di_green_file = "/dataref/bank/human/annotation/GRCh38/current/DI_green.tsv"
try:
    di_df = pd.read_csv(di_green_file, sep="\t")
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

# === INIT FINAL TABLE ===
summary_rows = []

# === LOOP THROUGH SAMPLES BASED ON DUP FILES ===
dup_files = glob.glob(f"{run_path}/pipeline_v0/dup/*/*_dup.txt")
for dup_file in dup_files:
    sample_id = os.path.basename(dup_file).replace("_dup.txt", "")
    sample_row = {"sample id": sample_id}

    # Add run info columns
    sample_row["run_name"] = run_name
    sample_row["date"] = date
    sample_row["run_id"] = run_id

    # === DUPLICATION METRICS ===
    try:
        with open(dup_file) as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith("rnaseq-capture"):
                parts = line.strip().split()
                sample_row["nb read"] = int(parts[1]) + 2 * int(parts[2])
                sample_row["mapped"] = 2 * int(parts[2])
                sample_row["dup"] = round(float(parts[8]) * 100, 2)
                break
    except Exception as e:
        print(f"⚠️ Could not parse dup file for {sample_id}: {e}")
        sample_row["nb read"] = sample_row["mapped"] = sample_row["dup"] = None

    # === COVERAGE METRIC ===
    on_target_file = f"{run_path}/pipeline_v0/coverage/{sample_id}/{sample_id}_on_target.txt"
    try:
        with open(on_target_file) as f:
            sample_row["on target"] = int(f.readline().strip())
    except Exception as e:
        print(f"⚠️ Could not read on-target file for {sample_id}: {e}")
        sample_row["on target"] = None

    # === SHORT ID ===
    #short_id = sample_id.split("-")[0]
    short_id = re.split('[-_]', sample_id)[0]
    # === OUTRIDER EVENT COUNT ===
    outrider_file = f"{run_path}/pipeline_v0/outrider/filesbysample/{short_id}.outrider.tab"
    #outrider_file = f"/datawork/genetique/RNASeq/diag/prod/20250918_RUN40_NextSeq_High_15RNASEQ/pipeline_v0/outrider/filesbysample/{short_id}.outrider.tab"
    try:
        with open(outrider_file) as f:
            lines = f.readlines()
        sample_row["nb event out"] = len(lines) - 1
    except Exception as e:
        print(f"⚠️ Could not count OUTRIDER events for {short_id}: {e}")
        sample_row["nb event out"] = None

    # === FRASER EVENT COUNT ===
    fraser_file = f"{run_path}/pipeline_v0/fraser/filesbysample/{short_id}.fraser.tab"
#    fraser_file = f"/datawork/genetique/RNASeq/diag/prod/20250918_RUN40_NextSeq_High_15RNASEQ/pipeline_v0/fraser/filesbysample/{short_id}.fraser.tab"
    try:
        with open(fraser_file) as f:
            lines = f.readlines()
        sample_row["nb event fraser"] = len(lines) - 1
    except Exception as e:
        print(f"⚠️ Could not count FRASER events for {short_id}: {e}")
        sample_row["nb event fraser"] = None

    # === TPM VALUES FOR TARGET GENES ===
    tpm_sample_id = sample_id.replace("-", ".").replace("_", ".")
    #print(tpm_df)
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
            tpm_sample_id = sample_id.replace("-", ".").replace("_", ".")
            #tpm_sample_id = sample_id
            if tpm_sample_id in tpm_df.columns:
                sample_tpms = tpm_df[tpm_sample_id]
                intersect_genes = di_genes.intersection(sample_tpms.index)
                expressed = sample_tpms.loc[list(intersect_genes)] > 10
                nb_expressed = expressed.sum()
                nb_total = len(intersect_genes)
                pct_expressed = (nb_expressed / nb_total * 100) if nb_total > 0 else 0
                sample_row["nb_DI_green_expressed"] = nb_expressed
                sample_row["nb_DI_green_total"] = nb_total
                sample_row["pct_DI_green_TPM>10"] = round(pct_expressed, 2)
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
summary_df.to_csv(output_file, sep="\t", index=False)
print(f"✅ QC summary saved to: {output_file}")

