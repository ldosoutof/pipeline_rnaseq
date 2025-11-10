#!/usr/bin/env python3
# scripts/generate_notebook_pca.py

import nbformat as nbf
import argparse

# =========================
# 0) Parse arguments
# =========================
parser = argparse.ArgumentParser(description="Generate PCA analysis notebook from Snakemake parameters")
parser.add_argument("--base", required=True, help="Base folder with TPM files")
parser.add_argument("--mapping_file", required=True, help="GTF mapping file")
parser.add_argument("--blacklist_file", required=True, help="Blacklist file")
parser.add_argument("--samples_run", required=True, help="Run to include")
parser.add_argument("--excluded_run", required=True, help="Run to exclude")
parser.add_argument("--keyword", required=True, help="Keyword to filter samples")
parser.add_argument("--output_html", required=True, help="Path to save interactive HTML plot")
parser.add_argument("--notebook_path", required=True, help="Path to save the generated notebook")
args = parser.parse_args()

# =========================
# 1) Create notebook object
# =========================
nb = nbf.v4.new_notebook()

# =========================
# 2) Title / description
# =========================
nb.cells.append(nbf.v4.new_markdown_cell("""
# PCA Analysis from Gene-level TPM
**Description**: Build gene-level TPM matrix from Kallisto outputs, perform PCA, and generate an interactive Plotly HTML plot.
Fully Snakemake-ready with blacklist file.
"""))

# =========================
# 3) Imports
# =========================
nb.cells.append(nbf.v4.new_code_cell("""
import pandas as pd
import numpy as np
import glob
import os
import re
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import plotly.io as pio
"""))

# =========================
# 4) Parameters
# =========================
nb.cells.append(nbf.v4.new_code_cell(f"""
# Snakemake-style parameters from script arguments
base = r"{args.base}"
mapping_file = r"{args.mapping_file}"
blacklist_file = r"{args.blacklist_file}"
samples_run = "{args.samples_run}"
excluded_run = "{args.excluded_run}"
keyword = "{args.keyword}"
output_html = r"{args.output_html}"

print("Base folder:", base)
print("Mapping file:", mapping_file)
print("Blacklist file:", blacklist_file)
print("Samples run:", samples_run)
print("Excluded run:", excluded_run)
print("Keyword:", keyword)
print("Output HTML:", output_html)
"""))

# =========================
# 5) Collect TPM files
# =========================
nb.cells.append(nbf.v4.new_code_cell("""
pattern = f"{base}/*/pipeline_v0/kallisto_bed/*/abundance.tsv"
files = [
    f for f in glob.glob(pattern)
    if keyword in os.path.basename(os.path.dirname(f)) and excluded_run not in f
]
print(f"✅ Found {len(files)} files")
"""))

# =========================
# 6) Build TPM matrix
# =========================
nb.cells.append(nbf.v4.new_code_cell("""
tpm_tables = []
standard_index = None
for fp in files:
    sample = os.path.basename(os.path.dirname(fp))
    try:
        df = pd.read_csv(fp, sep="\\t", usecols=["target_id","tpm"], dtype={"target_id": str})
        df["target_id"] = df["target_id"].str.strip()
        df["tpm"] = pd.to_numeric(df["tpm"], errors="coerce")
        df.set_index("target_id", inplace=True)
        df.columns = [sample]
        if standard_index is None:
            standard_index = df.index
        else:
            df = df.reindex(standard_index)
        tpm_tables.append(df)
    except Exception as e:
        print(f"❌ Error reading {fp}: {e}")

tpm_matrix = pd.concat(tpm_tables, axis=1)
print(f"✅ TPM matrix shape: {tpm_matrix.shape}")
"""))

# =========================
# 7) Map transcripts to genes
# =========================
nb.cells.append(nbf.v4.new_code_cell("""
gtf_df = pd.read_csv(
    mapping_file,
    sep="\\t",
    comment="#",
    header=None,
    names=["seqname","source","feature","start","end","score","strand","frame","attribute"]
)

def parse_attributes(attr_str):
    attrs = {}
    for pair in attr_str.strip().split(";"):
        if pair:
            key_value = pair.strip().split(" ", 1)
            if len(key_value) == 2:
                key, value = key_value
                attrs[key] = value.strip('"')
    return attrs

gtf_df["gene_id"] = gtf_df["attribute"].apply(lambda x: parse_attributes(x).get("gene_id"))
gtf_df["transcript_id"] = gtf_df["attribute"].apply(lambda x: parse_attributes(x).get("transcript_id"))

mapping_df = gtf_df[["gene_id", "transcript_id"]].drop_duplicates()
tpm_matrix["transcript_id"] = tpm_matrix.index
merged = pd.merge(tpm_matrix, mapping_df, how="left", on="transcript_id")
merged = merged[~merged["gene_id"].isna()]
gene_tpm = merged.drop(columns=["transcript_id"]).groupby("gene_id").sum(numeric_only=True)
print(f"✅ Aggregated gene-level TPM matrix shape: {gene_tpm.shape}")
"""))

# =========================
# 8) Clean columns and apply blacklist
# =========================
nb.cells.append(nbf.v4.new_code_cell("""
short_names = gene_tpm.columns.to_series().apply(lambda x: re.split(r"[-_\\.]", str(x))[0])
last_occurrence = {short: i for i, short in enumerate(short_names)}
final_indices = sorted(last_occurrence.values())
gene_tpm = gene_tpm.iloc[:, final_indices]
final_short_names = [short_names[i] for i in final_indices]
gene_tpm.columns = final_short_names
print(f"✅ Columns cleaned, final shape: {gene_tpm.shape}")

# Remove blacklisted samples
blacklist_df = pd.read_csv(blacklist_file, sep="\\t", header=None)
blacklist_samples = blacklist_df.iloc[:,0].astype(str).tolist()
df_expr = gene_tpm.drop(columns=[c for c in gene_tpm.columns if c in blacklist_samples], errors="ignore")
df_expr = df_expr.loc[:, df_expr.columns.str.startswith('2')]  # keep only sequencing samples
"""))

# =========================
# 9) PCA computation
# =========================
nb.cells.append(nbf.v4.new_code_cell("""
def min_max_scale_row(row):
    min_val, max_val = row.min(), row.max()
    return (row-min_val)/(max_val-min_val) if max_val != min_val else pd.Series([0]*len(row), index=row.index)

df_norm = df_expr.apply(min_max_scale_row, axis=1)
df_expr_t = df_norm.T.copy()

from sklearn.decomposition import PCA
pca = PCA(n_components=8, random_state=42)
components = pca.fit_transform(df_expr_t)
explained_var = pca.explained_variance_ratio_ * 100
df_pca = pd.DataFrame(components, columns=[f'PC{i+1}' for i in range(8)], index=df_expr_t.index)
df_pca['Sample'] = df_pca.index
df_pca['SampleShort'] = df_pca['Sample'].str.split(r'[-_]').str[0].str[:7]
"""))

# =========================
# 9b) Identify top contributing genes to PC1
# =========================
nb.cells.append(nbf.v4.new_code_cell("""
# Loadings: how much each gene contributes to each PC
loadings = pd.DataFrame(
    pca.components_.T,
    index=df_expr_t.columns,
    columns=[f'PC{i+1}' for i in range(pca.n_components)]
)

# Sort by absolute loading for PC1
top_genes_pc1 = loadings['PC1'].abs().sort_values(ascending=False).head(10)
top_genes_pc1_df = loadings.loc[top_genes_pc1.index, ['PC1']].copy()
top_genes_pc1_df['abs_loading'] = top_genes_pc1
print("Top 10 genes influencing PC1:")
display(top_genes_pc1_df)
"""))

# =========================
# 10) Assign colors
# =========================
nb.cells.append(nbf.v4.new_code_cell("""
name_folders = glob.glob(f"{base}/{samples_run}/pipeline_v0/htseq/*")
#samples_run_list = [os.path.basename(f.strip('/')).split(r'[-_]')[0] for f in name_folders]
samples_run_list = [re.split(r'[-_.]', os.path.basename(f.strip('/')))[0][:7] for f in name_folders]
print(f"samples: {samples_run_list}")
test = df_pca['SampleShort'][0]
print(f"samples in pca: {test}")
df_pca['Color'] = ["red" if s in samples_run_list else "lightgray" for s in df_pca['SampleShort']]
red_samples = df_pca.loc[df_pca['Color'] == "red", 'SampleShort'].tolist()
gray_samples = df_pca.loc[df_pca['Color'] == "lightgray", 'SampleShort'].tolist()
print(f"samples in pca RED: {red_samples}")
"""))
#print(samples_run_list)
# =========================
# 11) Plot PCA
# =========================
nb.cells.append(nbf.v4.new_code_cell("""
import plotly.graph_objects as go
import plotly.io as pio

pc_pairs=[(0,1),(2,3),(4,5),(6,7)]
def pc_label(i): return f"PC{i+1} ({explained_var[i]:.2f}%)"
def title_pair(a,b): return f"PCA - {pc_label(a)} vs {pc_label(b)}"

fig = go.Figure()
for i,(a,b) in enumerate(pc_pairs):
    fig.add_trace(go.Scatter(
        x=df_pca[f'PC{a+1}'], y=df_pca[f'PC{b+1}'],
        mode='markers',
        marker=dict(size=9,color=df_pca['Color'],opacity=0.9),
        hovertext=[f"Sample: {s}<br>{pc_label(a)}: {x:.2f}<br>{pc_label(b)}: {y:.2f}"
                   for s,x,y in zip(df_pca['Sample'], df_pca[f'PC{a+1}'], df_pca[f'PC{b+1}'])],
        hoverinfo='text',
        name=f"{pc_label(a)} vs {pc_label(b)}",
        visible=(i==0)
    ))

dropdown_buttons=[]
for i,(a,b) in enumerate(pc_pairs):
    visible=[j==i for j in range(len(pc_pairs))]
    dropdown_buttons.append(dict(
        label=f"{pc_label(a)} vs {pc_label(b)}",
        method='update',
        args=[{'visible':visible},
              {'title':title_pair(a,b),'xaxis':{'title':pc_label(a)},'yaxis':{'title':pc_label(b)}}]
    ))

fig.update_layout(
    title=dict(text=f"RUN: {samples_run} - {title_pair(0,1)}", x=0.5, xanchor="center"),
    xaxis_title=pc_label(0),
    yaxis_title=pc_label(1),
    plot_bgcolor="white",
    paper_bgcolor="white",
    updatemenus=[dict(type="dropdown",buttons=dropdown_buttons,x=0,y=1.1)]
)

pio.write_html(fig, file=output_html, full_html=True, include_plotlyjs='cdn', auto_open=False)
print(f"✅ PCA interactive plot created: {output_html}")
"""))

# =========================
# 12) Save notebook
# =========================
with open(args.notebook_path, 'w') as f:
    nbf.write(nb, f)

print(f"✅ Notebook created: {args.notebook_path}")

