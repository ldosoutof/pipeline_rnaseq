#!/usr/bin/env python3
import os
import re
import time
import pandas as pd
import numpy as np
import argparse
# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def normalize_sample_id(s):
    """Extract core sample ID (23D1192, X24D1028 → 24D1028)"""
    m = re.search(r'(\d{2,}[A-Z]\d+)', str(s))
    return m.group(1) if m else str(s)

# ------------------------------------------------------------
# Arguments
# ------------------------------------------------------------
parser = argparse.ArgumentParser("Annotate OUTRIDER + FRASER hits")
parser.add_argument('-i', '--outrider', required=True)
parser.add_argument('-r', '--fraser', required=True)
parser.add_argument('-g', '--gtf', required=True)
parser.add_argument('-d', '--di', required=True)
parser.add_argument('-o', '--omim', required=True)
parser.add_argument('-p', '--pli', required=True)
parser.add_argument('-m', '--hpo', required=True)
parser.add_argument('-f', '--output_file', required=True)
args = parser.parse_args()

t0 = time.time()
print("Starting annotation...")

# ------------------------------------------------------------
# Load reference data
# ------------------------------------------------------------
print("Loading annotation databases...")
di = pd.read_csv(args.di, sep='\t')[['EnsemblId(GRch38)','Model_Of_Inheritance','Phenotypes']]
omim = pd.read_csv(args.omim, sep='\t')
omim['omim_id'] = omim['disease_id'].str.extract(r'(OMIM:[0-9]+)')
omim = omim[['gene_symbol','hpo_id','hpo_name','omim_id']]
pli = pd.read_csv(args.pli, sep='\t')[['gene_id','exac_pLI','oe_lof','mis_z']]
hpo = pd.read_csv(args.hpo, sep='\t', comment='#')[['database_id','disease_name']]

# ------------------------------------------------------------
# GTF → ENSG → gene info mapping
# ------------------------------------------------------------
print("Parsing GTF...")
gtf_map = {}
with open(args.gtf) as f:
    for l in f:
        if l.startswith('#'): continue
        c = l.rstrip().split('\t')
        if len(c) < 9 or c[2] not in ('gene','transcript','exon','CDS'): continue
        
        gene_id = re.search(r'gene_id "([^"]+)"', c[8])
        gene_name = re.search(r'gene_name "([^"]+)"', c[8])
        if gene_id:
            gtf_map.setdefault(gene_id.group(1), []).append({
                'chromosome': c[0],
                'start': c[3],
                'gene_name': gene_name.group(1) if gene_name else ''
            })

# ------------------------------------------------------------
# FRASER hits: sample × HGNC → count
# ------------------------------------------------------------
print("Processing FRASER hits...")
fraser_hits = {}
if os.path.exists(args.fraser):
    fr = pd.read_csv(args.fraser, sep='\t', low_memory=False)
    fr.columns = fr.columns.str.replace('"', '')
    fr = fr[fr['pValue'] < 0.0001].dropna(subset=['hgncSymbol'])
    fr['sampleID'] = fr['sampleID'].map(normalize_sample_id)
    fraser_hits = fr.groupby(['sampleID', 'hgncSymbol']).size().to_dict()
print(f"FRASER: {len(fraser_hits)} (sample,gene) pairs")

# ------------------------------------------------------------
# Load & annotate OUTRIDER
# ------------------------------------------------------------
print("Annotating OUTRIDER...")
out = pd.read_csv(args.outrider, sep='\t')
out.columns = out.columns.str.strip('"').str.replace('"', '')

# Drop R row-index column if present (unnamed first column of integers)
if out.columns[0] in ('', 'Unnamed: 0') or out.columns[0].startswith('Unnamed'):
    out = out.iloc[:, 1:]

# Strip quotes from string values (old R write.table format)
for col in out.select_dtypes(include=['object', 'str']).columns:
    out[col] = out[col].astype(str).str.strip('"')

# Normalise column names — outrider_newVersion.R writes padjValue/rawCounts
if 'padjValue' in out.columns and 'padjust' not in out.columns:
    out = out.rename(columns={'padjValue': 'padjust'})
if 'rawCounts' in out.columns and 'rawcounts' not in out.columns:
    out = out.rename(columns={'rawCounts': 'rawcounts'})

rows = []
header = out.columns.tolist() + ['chromosome', 'start', 'gene_name', 'fraser_hits']

for _, r in out.iterrows():
    ensg = str(r['geneID'])
    pval = float(r['pValue'])
    sample = normalize_sample_id(r['sampleID'])
    
    # Filter: significant + GTF mapping
    if ensg not in gtf_map or pval >= 0.05: 
        continue
    
    # Extract gene info from GTF
    gene_names = set()
    chrom, start = '', ''
    for g in gtf_map[ensg]:
        if g['gene_name']: 
            gene_names.add(g['gene_name'])
        if not chrom:
            chrom, start = g['chromosome'], g['start']
    
    # Count FRASER hits for all gene names
    fraser_count = sum(fraser_hits.get((sample, g), 0) for g in gene_names)
    
    rows.append(r.tolist() + [chrom, start, ','.join(sorted(gene_names)), fraser_count])

df = pd.DataFrame(rows, columns=header)

# ------------------------------------------------------------
# Add functional annotations
# ------------------------------------------------------------
print("Adding functional annotations...")

# DI panel inheritance/phenotypes
df['Model_Of_Inheritance'] = df['geneID'].apply(
    lambda g: ';'.join(set(di[di['EnsemblId(GRch38)']==g]['Model_Of_Inheritance'].dropna().astype(str)))
)
df['Phenotypes'] = df['geneID'].apply(
    lambda g: ';'.join(set(di[di['EnsemblId(GRch38)']==g]['Phenotypes'].dropna().astype(str)))
)

# OMIM / HPO
df['omim_id'] = df['gene_name'].apply(
    lambda gns: ';'.join(set(omim[omim['gene_symbol'].isin(str(gns).split(','))]['omim_id'].dropna().astype(str)))
)
df['hpo_name'] = df['gene_name'].apply(
    lambda gns: ';'.join(set(omim[omim['gene_symbol'].isin(str(gns).split(','))]['hpo_name'].dropna().astype(str)))
)

# gnomAD constraint scores
df = df.merge(pli, how='left', left_on='geneID', right_on='gene_id', suffixes=('', '_pli'))
df.drop('gene_id_pli', axis=1, inplace=True, errors='ignore')

# ------------------------------------------------------------
# French Excel formatting
# ------------------------------------------------------------
print("Applying French number formatting...")
numeric_cols = ['pValue', 'padjust', 'zScore', 'normcounts', 'meanCorrected', 'theta', 
                'exac_pLI', 'oe_lof', 'mis_z']
for col in numeric_cols:
    if col in df.columns:
        df[col] = df[col].astype(str).str.replace('.', ',')

# ------------------------------------------------------------
# Save & report
# ------------------------------------------------------------
df.to_csv(args.output_file, sep='\t', index=False)
print(f"[DONE] {len(df)} significant events in {time.time()-t0:.1f}s")
print(f"Output: {args.output_file}")
print("\nPreview:")
print(df[['sampleID', 'geneID', 'gene_name', 'fraser_hits', 'Model_Of_Inheritance']].head())

