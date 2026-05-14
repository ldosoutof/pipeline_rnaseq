#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# -----------------------------
# Arguments
# -----------------------------
if len(sys.argv) < 4:
    sys.exit(
        "Usage:\n"
        "  plot_vaf_violin_by_sample.py <output.png> <chrY_mean_expression.tsv> <vaf_table1.tsv> [vaf_table2.tsv ...]"
    )

output_png = sys.argv[1]
sex_table   = sys.argv[2]
vaf_tables  = sys.argv[3:]

# -----------------------------
# Load sex inference table
# -----------------------------
sex_df = pd.read_csv(sex_table, sep="\t")

# Expected columns: sample, mean_chrY_TPM
required_sex_cols = {"sample", "mean_chrY_TPM"}
if not required_sex_cols.issubset(sex_df.columns):
    sys.exit(
        f"Sex table must contain columns {required_sex_cols}, "
        f"found {list(sex_df.columns)}"
    )

# Females = low chrY expression
female_samples = sex_df.loc[
    sex_df["mean_chrY_TPM"] < 1, "sample"
].tolist()

if len(female_samples) == 0:
    sys.exit("No female samples detected (mean_chrY_TPM < 1)")

# -----------------------------
# Load VAF tables
# -----------------------------
dfs = []

for vaf_file in vaf_tables:
    sample = Path(vaf_file).parent.name.replace("-", ".")
    if sample not in female_samples:
        continue

    df = pd.read_csv(vaf_file, sep="\t", header=None,
                 names=["CHROM", "POS", "REF", "ALT", "DP", "VAF", "folded"])

    required_cols = {"CHROM", "POS", "VAF", "folded"}
    if not required_cols.issubset(df.columns):
        sys.exit(
            f"{vaf_file} is missing required columns "
            f"{required_cols}, found {list(df.columns)}"
        )

    df = df[df["CHROM"].isin(["X", "chrX"])]
    df["sample"] = sample
    dfs.append(df)

if not dfs:
    sys.exit("No VAF data available for female samples")

vaf_df = pd.concat(dfs, ignore_index=True)

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(max(10, len(vaf_df["sample"].unique()) * 0.6), 6))

sns.violinplot(
    data=vaf_df,
    x="sample",
    y="folded",
    inner=None,
    cut=0,
    linewidth=1
)

sns.stripplot(
    data=vaf_df,
    x="sample",
    y="folded",
    color="black",
    size=2,
    alpha=0.4
)

# Median per sample (key metric for X-inactivation bias)
medians = vaf_df.groupby("sample")["folded"].median()

for i, (sample, median) in enumerate(medians.items()):
    plt.plot(i, median, "r_", markersize=12, markeredgewidth=2)
    plt.text(i, median + 0.01, f"{median:.2f}",
             ha="center", va="bottom", fontsize=8, color="red")

# Reference lines
plt.axhline(0.5, color="green",  linestyle="--", linewidth=0.8,
            label="Balanced (median=0.5)")
plt.axhline(0.3, color="orange", linestyle="--", linewidth=0.8,
            label="Moderate skew (0.3)")
plt.axhline(0.1, color="red",    linestyle="--", linewidth=0.8,
            label="Extreme skew (0.1)")

plt.ylim(0.0, 0.5)
plt.ylabel("Folded VAF  [ min(VAF, 1-VAF) ]")
plt.xlabel("Sample (females only)")
plt.title("chrX X-inactivation bias\n"
          "Median 0.5 = balanced  |  Median → 0.0 = complete skew")
plt.legend(fontsize=8)

plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(output_png, dpi=300)
plt.close()

print(f"[OK] Violin plot written to {output_png}")

