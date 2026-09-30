import os
import re
import pandas as pd
import argparse
import numpy as np

from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from collections import defaultdict

# =============================
# Arguments
# =============================

parser = argparse.ArgumentParser(description="Annotate Fraser file")
parser.add_argument('-f', '--fraser', required=True)
parser.add_argument('-g', '--gtf', required=True)
parser.add_argument('-d', '--di', required=True)
parser.add_argument('-o', '--omim', required=True)
parser.add_argument('-p', '--pli', required=True)
parser.add_argument('-m', '--hpo', required=True)
parser.add_argument('-oR', '--outrider', required=True)
parser.add_argument('--output',  required=True,  help='Path to the output file')
parser.add_argument('--magnis',  default='',      help='Chemin vers le fichier Excel Magnis (optionnel, pour enrichissement ID_SPICE)')

args = parser.parse_args()

# =============================
# Load annotation tables
# =============================

di = pd.read_csv(args.di, sep='\t', na_values=[''])
di_subset = di[['EnsemblId(GRch38)', 'Model_Of_Inheritance', 'Phenotypes']]

omim = pd.read_csv(args.omim, sep='\t', na_values=[''])
omim['omim_id'] = omim['disease_id'].str.extract(r'(OMIM:[0-9]+)')
omim_subset = omim[['gene_symbol', 'hpo_id', 'hpo_name', 'omim_id']]

pli = pd.read_csv(args.pli, sep='\t', na_values=[''])
pli_subset = pli[['gene_id', 'exac_pLI', 'oe_lof', 'mis_z']]

hpo = pd.read_csv(args.hpo, sep='\t', na_values=[''], comment='#')
hpo_subset = hpo[['database_id', 'disease_name']]

# =============================
# Load OUTRIDER table
# =============================

outrider = pd.read_csv(args.outrider, sep='\t', na_values=[''])
outrider.columns = outrider.columns.str.strip('"')

# Drop R row-index column if present
if outrider.columns[0] in ('', 'Unnamed: 0') or outrider.columns[0].startswith('Unnamed'):
    outrider = outrider.iloc[:, 1:]

# Strip quotes from string values (old R write.table format)
for col in outrider.select_dtypes(include=['object', 'str']).columns:
    outrider[col] = outrider[col].astype(str).str.strip('"')

# Normalise column names — outrider_newVersion.R writes padjValue/rawCounts
if 'padjValue' in outrider.columns and 'padjust' not in outrider.columns:
    outrider = outrider.rename(columns={'padjValue': 'padjust'})
if 'rawCounts' in outrider.columns and 'rawcounts' not in outrider.columns:
    outrider = outrider.rename(columns={'rawCounts': 'rawcounts'})

def normalize_outrider_sample_id(s):
    """Extract core sample ID from any known format:
      25D2693.STEMC.PUROMOINS.AVITI  (R output, dots)
      25D2693-STEMC-PUROMOINS-AVITI  (raw name, dashes)
      X25D2693                        (older R output with X prefix)
    """
    s = str(s)
    if s.startswith('X'):
        s = s[1:]
    return re.split(r'[.\-]', s)[0]

outrider['sampleID_norm'] = outrider['sampleID'].apply(normalize_outrider_sample_id)

# significance cutoff for OUTRIDER
outrider_sig = outrider[outrider['padjust'] < 0.05].copy()

outrider_hits_by_sample = defaultdict(set)
for _, r in outrider_sig.iterrows():
    ou_samp = r['sampleID_norm']
    ou_gene = str(r['geneID'])
    outrider_hits_by_sample[ou_samp].add(ou_gene)

# =============================
# Load Excel for ID_SPICE
# =============================

_magnis_path = args.magnis
if _magnis_path and os.path.isfile(_magnis_path):
    df_magnis = pd.read_excel(_magnis_path, engine='openpyxl')
    df_with_id_spice = df_magnis[df_magnis['ID_SPICE'].notna()]
else:
    if _magnis_path:
        print(f"[WARN] Fichier Magnis introuvable : {_magnis_path} — ID_SPICE non renseigné")
    df_with_id_spice = pd.DataFrame(columns=['N° Genno', 'ID_SPICE'])

def get_id_spice(sample_id):
    result = df_with_id_spice[df_with_id_spice['N° Genno'] == sample_id]['ID_SPICE']
    return result.values[0] if not result.empty else None

# =============================
# Parse GTF and build interval trees (fast O(n log m) lookup)
# =============================

import numpy as np
try:
    from ncls import NCLS
    USE_NCLS = True
except ImportError:
    USE_NCLS = False
    print("[WARN] ncls not available — falling back to linear GTF scan (slow for large files)")

gtf_entries  = {}   # chrom -> list of (start, end, gene_id, gene_name)
hgnc_to_ensg = {}

with open(args.gtf, 'r') as gtf_file:
    for line in gtf_file:
        if line.startswith('#'):
            continue
        data = line.strip().split('\t')
        if len(data) < 9 or data[2] not in ['gene', 'transcript', 'exon', 'CDS']:
            continue
        chromosome = data[0]
        start_pos  = int(data[3])
        end_pos    = int(data[4])
        attributes = data[8]

        gene_id = gene_name = None
        for attribute in attributes.split(';'):
            attribute = attribute.strip()
            if attribute.startswith('gene_id') and '"' in attribute:
                gene_id = attribute.split('"')[1]
            elif attribute.startswith('gene_name') and '"' in attribute:
                gene_name = attribute.split('"')[1]

        if chromosome not in gtf_entries:
            gtf_entries[chromosome] = []
        gtf_entries[chromosome].append((start_pos, end_pos, gene_id, gene_name))

        if gene_name and gene_id:
            hgnc_to_ensg.setdefault(gene_name, set()).add(gene_id)

# Build interval trees per chromosome
if USE_NCLS:
    gtf_trees = {}
    for chrom, entries in gtf_entries.items():
        starts = np.array([e[0] for e in entries], dtype=np.int64)
        ends   = np.array([e[1] for e in entries], dtype=np.int64)
        ids    = np.arange(len(entries), dtype=np.int64)
        gtf_trees[chrom] = (NCLS(starts, ends, ids), entries)
    print(f"[GTF] Built interval trees for {len(gtf_trees)} chromosomes", flush=True)
else:
    gtf_trees = None

def query_gtf(chrom, start, end):
    """Return list of (gene_id, gene_name) overlapping [start, end]."""
    if USE_NCLS and gtf_trees and chrom in gtf_trees:
        tree, entries = gtf_trees[chrom]
        hits = list(tree.find_overlap(start, end))
        results = []
        for _, _, idx in hits:
            e = entries[idx]
            if e[2] or e[3]:
                results.append((e[2], e[3]))
        return results
    elif chrom in gtf_entries:
        return [
            (e[2], e[3]) for e in gtf_entries[chrom]
            if e[0] <= end and e[1] >= start
        ]
    return []

# =============================
# Read FRASER file and annotate
# =============================

print(f"[INFO] Loading FRASER file: {args.fraser}", flush=True)

# Auto-detect format:
# - New format (fraser_newVer.R fwrite): first column = 'sampleID', no row index
# - Old format (fraser.R write.table): first column = '"seqnames"' with R row index prepended
with open(args.fraser) as _fh:
    _first_col = _fh.readline().split('\t')[0].strip().strip('"')
_has_row_index = _first_col != 'sampleID'
_index_col = 0 if _has_row_index else None
print(f"[INFO] Format: {'old (row index)' if _has_row_index else 'new (no row index)'} — first col='{_first_col}'", flush=True)

# Pre-filter during load using chunked reading to avoid OOM on large files (e.g. fraser_results_all.tsv)
# Only keep rows with pValue < 0.001 — the annotation loop skips everything else anyway
PVAL_THRESHOLD = 0.001
chunk_size = 500_000
chunks = []
total_rows = 0
kept_rows  = 0

for chunk in pd.read_csv(args.fraser, sep='\t', na_values=[''], low_memory=False,
                          index_col=_index_col, chunksize=chunk_size):
    # Strip quotes from column names on first chunk
    chunk.columns = chunk.columns.str.strip('"')
    if total_rows == 0:
        print(f"[INFO] Columns: {list(chunk.columns)}", flush=True)
    # Strip quotes from string values
    for col in chunk.select_dtypes(include=['object', 'str']).columns:
        chunk[col] = chunk[col].astype(str).str.strip('"')
    # Normalise pValue column name
    if 'pValue' not in chunk.columns and 'pvalue' in chunk.columns:
        chunk = chunk.rename(columns={'pvalue': 'pValue'})
    # Filter early
    if 'pValue' in chunk.columns:
        chunk['pValue'] = pd.to_numeric(chunk['pValue'], errors='coerce')
        chunk = chunk[chunk['pValue'] < PVAL_THRESHOLD]
    total_rows += chunk_size
    kept_rows  += len(chunk)
    chunks.append(chunk)
    if total_rows % 5_000_000 == 0:
        print(f"[INFO] Read {total_rows:,} rows, kept {kept_rows:,} so far...", flush=True)

fraser_df = pd.concat(chunks, ignore_index=False) if chunks else pd.DataFrame()

# If sampleID was used as index (old format detection), restore it as a column
if _has_row_index and fraser_df.index.name == 'sampleID':
    fraser_df = fraser_df.reset_index()
elif 'sampleID' not in fraser_df.columns and fraser_df.index.name:
    fraser_df = fraser_df.reset_index()

print(f"[INFO] Loaded {kept_rows:,} rows (pValue < {PVAL_THRESHOLD}) from {total_rows:,} total", flush=True)
print(f"[INFO] Columns after load: {list(fraser_df.columns)}", flush=True)

# Normalise column names — fraser_newVer.R uses padjValue, older scripts use padjust
if 'padjValue' in fraser_df.columns and 'padjust' not in fraser_df.columns:
    fraser_df = fraser_df.rename(columns={'padjValue': 'padjust'})

# Normalise chromosome column
if 'seqnames' not in fraser_df.columns and 'CHROM' in fraser_df.columns:
    fraser_df = fraser_df.rename(columns={'CHROM': 'seqnames'})

# Normalise hgncSymbol
if 'hgncSymbol' not in fraser_df.columns and 'gene_name' in fraser_df.columns:
    fraser_df = fraser_df.rename(columns={'gene_name': 'hgncSymbol'})

# Convert numeric columns
for col in ['pValue', 'padjust', 'deltaPsi', 'psiValue']:
    if col in fraser_df.columns:
        fraser_df[col] = pd.to_numeric(fraser_df[col], errors='coerce')

rows = []
extra_cols = ['ENSG', 'gene_name', 'ID_SPICE', 'ucsc_link', 'outrider_hit']
out_cols = list(fraser_df.columns) + extra_cols

for _, row in fraser_df.iterrows():
    chromosome  = str(row.get('seqnames', '')).strip('"')
    start_pos   = int(row.get('start', 0))
    end_pos     = int(row.get('end', 0))
    key_id      = f"{chromosome}:{start_pos}-{end_pos}"
    type_fraser = str(row.get('type', ''))
    pValue      = float(row.get('pValue', 1.0))
    sample_id   = str(row.get('sampleID', '')).strip('"')
    hgnc_val    = row.get('hgncSymbol', '')
    hgnc_symbol = str(hgnc_val).strip('"') if pd.notna(hgnc_val) else ''

    matching_genes = []
    matching_names = []
    link = ''
    outrider_hit_flag = 'no'

    if len(chromosome) <= 2:
        hits = query_gtf(chromosome, start_pos, end_pos)
        for gene_id_hit, gene_name_hit in hits:
            if type_fraser == 'theta':
                s = start_pos - 100
                e = end_pos + 100
                link = (
                    f'"https://genome-euro.ucsc.edu/cgi-bin/hgTracks?pix=1243&textSize=12&textFont=Helvetica'
                    f'&hgt.labelWidth=20&hgS_otherUserName=View&hgS_otherUserSessionName=Clinical_SNVs_hg38'
                    f'&hgS_otherUserSessionLabel=Clinical%20SNVs&hgS_otherUserSessionDesc='
                    f'Assess%20potential%20disease%20contributions%20of%20single%20nucleotide%20variants%20in%20coding%20regions'
                    f'&hgS_doOtherUser=submit&position=chr{chromosome}:{s}-{e}&highlight=chr{key_id}"'
                )
            else:
                length = max(end_pos - start_pos, 1)
                s = start_pos - length
                e = end_pos + length
                link = (
                    f'"https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default'
                    f'&lastVirtModeExtraState=&virtModeType=default&virtMode=0'
                    f'&position=chr{chromosome}:{s}-{e}&highlight=chr{key_id}";"{key_id}"'
                )

            if gene_id_hit and gene_id_hit not in matching_genes:
                matching_genes.append(gene_id_hit)
            if gene_name_hit and gene_name_hit not in matching_names:
                matching_names.append(gene_name_hit)

        # OUTRIDER integration via HGNC → ENSG
        ensg_from_hgnc = hgnc_to_ensg.get(hgnc_symbol, set()) if hgnc_symbol else set()
        ou_genes = outrider_hits_by_sample.get(sample_id, set())
        if any(g in ou_genes for g in ensg_from_hgnc):
            outrider_hit_flag = 'yes'

        id_spice = get_id_spice(sample_id)
        data_row = list(row) + [
            ','.join(matching_genes),
            ','.join(filter(None, matching_names)),
            id_spice or '',
            link,
            outrider_hit_flag
        ]
        rows.append(data_row)
df = pd.DataFrame(rows, columns=out_cols)

# =============================
# Reorder and annotate DataFrame
# =============================

# Move ucsc_link after "end" column
ucsc_link_col = df.pop('ucsc_link')
end_col = 'end' if 'end' in df.columns else '"end"'
if end_col in df.columns:
    df.insert(df.columns.get_loc(end_col) + 1, 'ucsc_link', ucsc_link_col)
else:
    df['ucsc_link'] = ucsc_link_col

def merge_model(row):
    ensg_values = row['ENSG'].split(',')
    model_names = di_subset[di_subset['EnsemblId(GRch38)'].isin(ensg_values)]['Model_Of_Inheritance'].tolist()
    model_names = list(set([str(name) for name in model_names if not pd.isnull(name)]))
    return ';'.join(model_names)

df['Model_Of_Inheritance'] = df.apply(merge_model, axis=1)

def pheno_model(row):
    ensg_values = row['ENSG'].split(',')
    pheno_names = di_subset[di_subset['EnsemblId(GRch38)'].isin(ensg_values)]['Phenotypes'].tolist()
    pheno_names = list(set([str(name) for name in pheno_names if not pd.isnull(name)]))
    return ';'.join(pheno_names)

df['Phenotypes'] = df.apply(pheno_model, axis=1)

def omim_id_merge(row):
    gene_values = row['gene_name'].split(',')
    omim_names = omim_subset[omim_subset['gene_symbol'].isin(gene_values)]['omim_id'].tolist()
    omim_names = list(set([str(name) for name in omim_names if not pd.isnull(name)]))
    return ';'.join(omim_names)

df['omim_id'] = df.apply(omim_id_merge, axis=1)

def create_first_omim_link(omim_ids):
    if pd.isna(omim_ids) or not isinstance(omim_ids, str) or omim_ids.strip() == '':
        return ''
    omim_id_list = [x.strip() for x in omim_ids.split(',') if x.strip()]
    if not omim_id_list:
        return ''
    first_omim_id = omim_id_list[0]
    parts = first_omim_id.split(':')
    if len(parts) == 2:
        return f'=LIEN_HYPERTEXTE("<https://www.omim.org/entry/{parts[1]}"; "{first_omim_id}")'
    return ''

df['first_omim_link'] = df['omim_id'].apply(create_first_omim_link)

def merge_hpo_names(row):
    gene_values = row['gene_name'].split(',')
    hpo_names = omim_subset[omim_subset['gene_symbol'].isin(gene_values)]['hpo_name'].tolist()
    hpo_names = list(set([str(name) for name in hpo_names if not pd.isnull(name)]))
    return ';'.join(hpo_names)

def merge_hpo_id(row):
    gene_values = row['gene_name'].split(',')
    hpo_id = omim_subset[omim_subset['gene_symbol'].isin(gene_values)]['hpo_id'].tolist()
    hpo_id = list(set([str(name) for name in hpo_id if not pd.isnull(name)]))
    return ';'.join(hpo_id)

df['hpo_id'] = df.apply(merge_hpo_id, axis=1)
df['hpo_name'] = df.apply(merge_hpo_names, axis=1)

def exac_id_merge(row):
    ensg_values = row['ENSG'].split(',')
    exac_names = pli_subset[pli_subset['gene_id'].isin(ensg_values)]['exac_pLI'].tolist()
    exac_names = list(set([str(name) for name in exac_names if not pd.isnull(name)]))
    return ';'.join(exac_names)

df['gnomAD_pLI'] = df.apply(exac_id_merge, axis=1)

def lof_id_merge(row):
    ensg_values = row['ENSG'].split(',')
    lof_names = pli_subset[pli_subset['gene_id'].isin(ensg_values)]['oe_lof'].tolist()
    lof_names = list(set([str(name) for name in lof_names if not pd.isnull(name)]))
    return ';'.join(lof_names)

df['oe_lof'] = df.apply(lof_id_merge, axis=1)

def mis_id_merge(row):
    ensg_values = row['ENSG'].split(',')
    mis_names = pli_subset[pli_subset['gene_id'].isin(ensg_values)]['mis_z'].tolist()
    mis_names = list(set([str(name) for name in mis_names if not pd.isnull(name)]))
    return ';'.join(mis_names)

df['mis_z'] = df.apply(mis_id_merge, axis=1)

def hpo_merge(row):
    omim_values = row['omim_id'].split(',')
    hpo_names = hpo_subset[hpo_subset['database_id'].isin(omim_values)]['disease_name'].tolist()
    hpo_names = list(set([str(name) for name in hpo_names if not pd.isnull(name)]))
    return ';'.join(hpo_names)

df['disease_name'] = df.apply(hpo_merge, axis=1)

inheritance_mapping = {
    'HP:0000006': 'Autosomal Dominant',
    'HP:0012275': 'Autosomal Dominant with Imprinting (Maternal)',
    'HP:0012274': 'Autosomal Dominant with Imprinting (Paternal)',
    'HP:0000007': 'Autosomal Recessive',
    'HP:0001417': 'X-Linked',
    'HP:0001419': 'X-Linked Recessive',
    'HP:0001423': 'X-Linked Dominant'
}

def extract_inheritance(row):
    hpo_ids = row['hpo_id'].split(',')
    inheritance = [inheritance_mapping[hpo] for hpo in hpo_ids if hpo in inheritance_mapping]
    return ';'.join(inheritance)

df['OMIM_Inheritance'] = df.apply(extract_inheritance, axis=1)

# =============================
# Final formatting
# =============================

df = df.drop_duplicates().reset_index(drop=True)

columns_to_exclude = df.columns[df.columns != 'first_omim_link']
df[columns_to_exclude] = df[columns_to_exclude].replace('"', '', regex=True)
df.columns = df.columns.str.replace('"', '')

def categorize_row(row):
    try:
        if float(row['deltaPsi']) >= 0.2:
            inheritance_values = str(row['Model_Of_Inheritance']).lower()
            if 'monoallelic' in inheritance_values or 'x-linked' in inheritance_values:
                return 'high'
            gnomAD_pLI_values = str(row['gnomAD_pLI']).split(';')
            gnomAD_pLI_values = [v for v in gnomAD_pLI_values if v not in ['', 'nan']]
            if len(gnomAD_pLI_values) == 0:
                return None
            max_pLI = float(max(gnomAD_pLI_values))
            if (float(row['pValue']) <= 0.0001 and
                pd.notna(row['gnomAD_pLI']) and row['gnomAD_pLI'] != '' and
                max_pLI >= 0.9):
                if pd.notna(row['omim_id']) and row['omim_id'] != '':
                    return 'mid'
                else:
                    return 'research'
    except Exception:
        return None
    return None

df['outrider'] = df.apply(categorize_row, axis=1)

for col in ['pValue', 'padjust', 'psiValue', 'deltaPsi',
            'meanCounts', 'meanTotalCounts', 'gnomAD_pLI',
            'oe_lof', 'mis_z']:
    if col in df.columns:
        df[col] = df[col].astype(str).str.replace('.', ',', regex=False)

print(df)
df.to_csv(args.output, sep="\t", index=False)

