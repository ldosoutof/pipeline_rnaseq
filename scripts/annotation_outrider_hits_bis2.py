#!/usr/bin/env python3
"""
annotation_outrider_hits_bis2.py — Annotate OUTRIDER + FRASER hits.

Version optimisée (remplacement direct, mêmes arguments, mêmes colonnes de sortie).

Pourquoi l'ancienne version prenait des heures :
  - annotations DI / OMIM / HPO : pour CHAQUE événement, un filtre pandas sur la
    totalité des tables de référence (O(événements × taille des tables)) ;
  - GTF : un dictionnaire par ligne gene/transcript/exon/CDS (~3 M d'entrées,
    dizaines de Go) alors que seules les infos gène sont utiles ;
  - iterrows() sur toute la table OUTRIDER.

Ici : index précalculés une fois (dict gène -> annotations), GTF réduit à une
entrée par gène, filtrage vectorisé. Résultat identique (l'ordre des valeurs
jointes par ';' dans une cellule est trié, là où l'ancienne version dépendait de
l'ordre arbitraire d'un set).
"""
import os
import re
import time
import argparse
from collections import defaultdict
import pandas as pd

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def normalize_sample_id(s):
    """Extract core sample ID (23D1192, X24D1028 -> 24D1028)"""
    m = re.search(r'(\d{2,}[A-Z]\d+)', str(s))
    return m.group(1) if m else str(s)


def sets_by_key(df, key_col, val_col):
    """{key: set(values non nulles, en str)} en un seul passage (groupby)."""
    sub = df[[key_col, val_col]].dropna(subset=[val_col])
    sub = sub.assign(**{val_col: sub[val_col].astype(str)})
    return sub.groupby(key_col)[val_col].agg(set).to_dict()


def join_sets(index, keys):
    """Union des sets de plusieurs clés, jointe par ';' (triée -> déterministe)."""
    acc = set()
    for k in keys:
        acc |= index.get(k, set())
    return ';'.join(sorted(acc))


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
di = pd.read_csv(args.di, sep='\t')[['EnsemblId(GRch38)', 'Model_Of_Inheritance', 'Phenotypes']]
omim = pd.read_csv(args.omim, sep='\t')
omim['omim_id'] = omim['disease_id'].str.extract(r'(OMIM:[0-9]+)')
omim = omim[['gene_symbol', 'hpo_id', 'hpo_name', 'omim_id']]
pli = pd.read_csv(args.pli, sep='\t')[['gene_id', 'exac_pLI', 'oe_lof', 'mis_z']]
hpo = pd.read_csv(args.hpo, sep='\t', comment='#')[['database_id', 'disease_name']]  # chargé comme avant (non utilisé en aval)

# Index précalculés (une seule passe chacun) — remplace les filtres par ligne
di_moi_idx  = sets_by_key(di, 'EnsemblId(GRch38)', 'Model_Of_Inheritance')
di_phen_idx = sets_by_key(di, 'EnsemblId(GRch38)', 'Phenotypes')
omim_id_idx = sets_by_key(omim, 'gene_symbol', 'omim_id')
hpo_name_idx = sets_by_key(omim, 'gene_symbol', 'hpo_name')
print(f"  index: DI={len(di_moi_idx)} genes, OMIM/HPO={len(omim_id_idx)} symbols  ({time.time()-t0:.1f}s)")

# ------------------------------------------------------------
# GTF -> ENSG -> gene info (une entrée par gène : chrom, start, gene_name)
# ------------------------------------------------------------
# L'ancienne version stockait toutes les lignes gene/transcript/exon/CDS puis
# prenait, par gène, le chrom/start de la 1re ligne rencontrée et l'ensemble
# des gene_name. On conserve cette sémantique : 1re occurrence -> chrom/start ;
# gene_name collectés (identiques pour un même gene_id dans un GTF Ensembl).
print("Parsing GTF...")
gtf_first = {}                     # ensg -> (chrom, start)
gtf_names = defaultdict(set)       # ensg -> {gene_name}
gid_re  = re.compile(r'gene_id "([^"]+)"')
name_re = re.compile(r'gene_name "([^"]+)"')
wanted = ('gene', 'transcript', 'exon', 'CDS')
with open(args.gtf) as f:
    for l in f:
        if l.startswith('#'):
            continue
        c = l.split('\t', 9)
        if len(c) < 9 or c[2] not in wanted:
            continue
        m = gid_re.search(c[8])
        if not m:
            continue
        gid = m.group(1)
        if gid not in gtf_first:
            gtf_first[gid] = (c[0], c[3])
        mn = name_re.search(c[8])
        if mn:
            gtf_names[gid].add(mn.group(1))
print(f"  GTF: {len(gtf_first)} genes  ({time.time()-t0:.1f}s)")

# ------------------------------------------------------------
# FRASER hits: sample x HGNC -> count
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
# Load & annotate OUTRIDER (vectorisé)
# ------------------------------------------------------------
print("Annotating OUTRIDER...")
out = pd.read_csv(args.outrider, sep='\t')
out.columns = out.columns.str.strip('"').str.replace('"', '')

# Drop R row-index column if present
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

orig_cols = out.columns.tolist()
n_total = len(out)

# Filtre : significatif + présent dans le GTF (même critère qu'avant)
ensg = out['geneID'].astype(str)
pval = pd.to_numeric(out['pValue'], errors='coerce')
keep = (pval < 0.05) & ensg.isin(gtf_first.keys())
df = out.loc[keep].copy()
df['geneID'] = df['geneID'].astype(str)

# Infos gène (lookups O(1))
df['chromosome'] = df['geneID'].map(lambda g: gtf_first[g][0])
df['start']      = df['geneID'].map(lambda g: gtf_first[g][1])
df['gene_name']  = df['geneID'].map(lambda g: ','.join(sorted(gtf_names.get(g, set()))))

# FRASER hits par (sample, chacun des gene_name)
samples = df['sampleID'].map(normalize_sample_id)
df['fraser_hits'] = [
    sum(fraser_hits.get((s, g), 0) for g in gn.split(',') if g)
    for s, gn in zip(samples, df['gene_name'])
]
df = df[orig_cols + ['chromosome', 'start', 'gene_name', 'fraser_hits']]
print(f"  kept {len(df)} / {n_total} events  ({time.time()-t0:.1f}s)")

# ------------------------------------------------------------
# Add functional annotations (lookups dans les index précalculés)
# ------------------------------------------------------------
print("Adding functional annotations...")
df['Model_Of_Inheritance'] = df['geneID'].map(lambda g: join_sets(di_moi_idx, [g]))
df['Phenotypes']           = df['geneID'].map(lambda g: join_sets(di_phen_idx, [g]))
df['omim_id']  = df['gene_name'].map(lambda gns: join_sets(omim_id_idx,  str(gns).split(',')))
df['hpo_name'] = df['gene_name'].map(lambda gns: join_sets(hpo_name_idx, str(gns).split(',')))

# gnomAD constraint scores (jointure, comme avant)
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
        df[col] = df[col].astype(str).str.replace('.', ',', regex=False)

# ------------------------------------------------------------
# Save & report
# ------------------------------------------------------------
df.to_csv(args.output_file, sep='\t', index=False)
print(f"[DONE] {len(df)} significant events in {time.time()-t0:.1f}s")
print(f"Output: {args.output_file}")
print("\nPreview:")
print(df[['sampleID', 'geneID', 'gene_name', 'fraser_hits', 'Model_Of_Inheritance']].head())
