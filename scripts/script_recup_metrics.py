import pandas as pd
import glob
import os
import argparse

# === PARSE ARGUMENTS ===
parser = argparse.ArgumentParser(description="Generate QC summary for a given run")
parser.add_argument("--run_path", required=True, help="Full path to the run folder (e.g. RUN29)")
parser.add_argument("--tpm_file", required=True, help="Path to the TPM matrix file")
parser.add_argument("--output", default="qc_summary.tsv", help="Output filename")
args = parser.parse_args()

# === CONFIGURATION ===
run_path = args.run_path
tpm_file = args.tpm_file
output_file = args.output
target_genes = ["HBA1", "HBA2", "HBB"]
#target_genes = ["HBA1", "HBA2", "HBB", "SRSF6", "SRSF3"]

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
    short_id = sample_id.split("-")[0] if "-" in sample_id else sample_id.split("_")[0]


    # === OUTRIDER EVENT COUNT ===
    outrider_file = f"{run_path}/pipeline_v0/outrider/filesbysample/{short_id}.outrider.tab"
    try:
        with open(outrider_file) as f:
            lines = f.readlines()
        sample_row["nb evet out"] = len(lines) - 1
    except Exception as e:
        print(f"⚠️ Could not count OUTRIDER events for {short_id}: {e}")
        sample_row["nb evet out"] = None

    # === FRASER EVENT COUNT ===
    fraser_file = f"{run_path}/pipeline_v0/fraser/filesbysample/{short_id}.fraser.tab"
    try:
        with open(fraser_file) as f:
            lines = f.readlines()
        sample_row["nb evet fraser"] = len(lines) - 1
    except Exception as e:
        print(f"⚠️ Could not count FRASER events for {short_id}: {e}")
        sample_row["nb evet fraser"] = None

    # === TPM VALUES FOR TARGET GENES ===
    #tpm_sample_id = sample_id.replace("-", ".").replace("_", ".")
    for gene in target_genes:
        try:
            sample_row[gene] = tpm_df.at[gene, sample_id]
        except KeyError:
            sample_row[gene] = None

    # === SUM OF HBA1 + HBA2 + HBB ===
    sample_row["HBA_total"] = sum(sample_row.get(g, 0) or 0 for g in ["HBA1", "HBA2", "HBB"])

    summary_rows.append(sample_row)

# === CREATE AND SAVE FINAL TABLE ===
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(output_file, sep="\t", index=False)
print(f"✅ QC summary saved to: {output_file}")

