import os
import pandas as pd
import argparse
import numpy as np


import pathlib
import duckdb # version 0.10.1
import polars # version 0.20.7
import glob

parser = argparse.ArgumentParser(description="Annotate output file from outrider")
parser.add_argument('-i', '--outrider', required=True)
parser.add_argument('-g', '--gtf', required=True)
parser.add_argument('-d', '--di', required=True)
parser.add_argument('-o', '--omim', required=True)
parser.add_argument('-p', '--pli', required=True)
parser.add_argument('-m', '--hpo', required=True)
parser.add_argument('-f', '--output_file', required=True, help="Output file name")

args = parser.parse_args()

outrider = pd.read_csv(args.outrider, sep='\t')
di = pd.read_csv(args.di, sep='\t', na_values=[''])
di_subset = di[['EnsemblId(GRch38)', 'Model_Of_Inheritance', 'Phenotypes']]
omim = pd.read_csv(args.omim, sep='\t', na_values=[''])
omim['omim_id'] = omim['disease_id'].str.extract(r'(OMIM:[0-9]+)')
omim_subset = omim[['gene_symbol', 'hpo_id', 'hpo_name', 'omim_id']]
pli = pd.read_csv(args.pli, sep='\t', na_values=[''])
pli_subset = pli[['gene_id', 'exac_pLI', 'oe_lof', 'mis_z']]
hpo = pd.read_csv(args.hpo, sep='\t', na_values=[''], comment='#')
hpo_subset = hpo[['database_id', 'disease_name']]


# Load the Excel file (replace 'your_file.xlsx' with the actual file path)
file_path = '/datawork/genetique/RNASeq/diag/prod/Results_LymphoRNA_Magnis_NS.xlsx'

# Read the Excel file into a pandas DataFrame
df_magnis = pd.read_excel(file_path, engine='openpyxl')

# Check if required columns exist
#if 'ID_SPICE' not in df_magnis.columns or 'N° Genno' not in df_magnis.columns:
#    raise ValueError("Required columns 'ID_SPICE' and 'N° Genno' not found in the Excel file.")

df_with_id_spice = df_magnis[df_magnis['ID_SPICE'].notna()]
#print(df_with_id_spice['ID_SPICE'])
# Function to retrieve ID_SPICE based on sample ID
def get_id_spice(sample_id):
    result = df_with_id_spice[df_with_id_spice['N° Genno'] == sample_id]['ID_SPICE']
    return result.values[0] if not result.empty else None


# Step 1: Parse the GTF file to create a mapping of genomic positions to gene IDs and transcript IDs
gtf_mapping = {}
with open(args.gtf, 'r') as gtf_file:
    for line in gtf_file:
        if line.startswith('#'):
            continue
        data = line.strip().split('\t')
        if data[2] in ['gene', 'transcript', 'exon', 'CDS']:
            chromosome = data[0]
            start_pos = int(data[3])
            end_pos = int(data[4])
            attributes = data[8]
            gene_id = None
            transcript_id = None
            gene_name = None
            for attribute in attributes.split(';'):
                if 'gene_id' in attribute:
                    gene_id = attribute.split('"')[1]
                elif 'transcript_id' in attribute:
                    transcript_id = attribute.split('"')[1]
                elif 'gene_name' in attribute:
                    gene_name = attribute.split('"')[1]
            if gene_id not in gtf_mapping:
                gtf_mapping[gene_id] = []
            gtf_mapping[gene_id].append({
                'chromosome': chromosome,
                'start': start_pos,
                'end': end_pos,
                'gene_name': gene_name,
                'gene_id': gene_id,
                'transcript_id': transcript_id
            })

# Create empty lists to store data
rows = []

# Step 2: Iterate through the Fraser file and find matching positions
with open(args.outrider, 'r') as fraser_file:
    lines = fraser_file.readlines()
    header = lines[0].strip().split('\t')
    header.extend(['chromosome', 'start', 'gene_name', 'ID_SPICE'])
    header = header[1:]    
    #print(header)    
#    header = [h.strip('"') for h in header]  # Remove quotes from header
#   rows = [line.strip().split('\t') for line in lines[1:]]
    for line in lines[1:]:
        data = line.strip().split('\t')
        ensg = data[1].split('"')[1]
        pValue = float(data[3])
        matching_transcripts = []
        matching_genes = []
        matching_names = []
        chromosome = ''
        start_pos = ''
        sample_id = data[2].split('"')[1]
        #print(sample_id)        
        if ensg in gtf_mapping and pValue < 0.05:
            for entry in gtf_mapping[ensg]:
                if chromosome == '' and start_pos == '':
                    chromosome = entry['chromosome']
                    start_pos = str(entry['start'])
                if not entry['gene_name'] in matching_names:
                    matching_names.append(entry['gene_name'])
            # Retrieve the corresponding ID_SPICE based on the sample ID
            id_spice = get_id_spice(sample_id)

            data.extend([chromosome, start_pos, ','.join(filter(None, matching_names)), id_spice or ''])
            rows.append(data[1:])

df = pd.DataFrame(rows, columns=header)
df = df.replace('"', '', regex=True)
df.columns = df.columns.str.replace('"', '')

def merge_model(row):
    ensg_values = row['geneID'].split(',')
    model_names = di_subset[di_subset['EnsemblId(GRch38)'].isin(ensg_values)]['Model_Of_Inheritance'].tolist()
    model_names = list(set([str(name) for name in model_names if not pd.isnull(name)]))
    return ';'.join(model_names)

df['Model_Of_Inheritance'] = df.apply(merge_model, axis=1)

def pheno_model(row):
    ensg_values = row['geneID'].split(',')
    pheno_names = di_subset[di_subset['EnsemblId(GRch38)'].isin(ensg_values)]['Phenotypes'].tolist()
    pheno_names = list(set([str(name) for name in pheno_names if not pd.isnull(name)]))
    return ';'.join(pheno_names)

df['Phenotypes'] = df.apply(pheno_model, axis=1)

def omim_id_merge(row):
    ensg_values = row['gene_name'].split(',')
    omim_names = omim_subset[omim_subset['gene_symbol'].isin(ensg_values)]['omim_id'].tolist()
    omim_names = list(set([str(name) for name in omim_names if not pd.isnull(name)]))
    return ';'.join(omim_names)

df['omim_id'] = df.apply(omim_id_merge, axis=1)

def create_first_omim_link(omim_ids):
    omim_id_list = omim_ids.split(',')
    if omim_id_list:
        first_omim_id = omim_id_list[0]
        parts = first_omim_id.strip().split(':')
        if len(parts) == 2:
            return f'=LIEN_HYPERTEXTE(\"https://www.omim.org/entry/{parts[1]}\"; \"{first_omim_id}\")'
        else:
            return ''
    else:
        return ''

df['first_omim_link'] = df['omim_id'].apply(create_first_omim_link)

def merge_hpo_names(row):
    ensg_values = row['gene_name'].split(',')
    hpo_names = omim_subset[omim_subset['gene_symbol'].isin(ensg_values)]['hpo_name'].tolist()
    hpo_names = list(set([str(name) for name in hpo_names if not pd.isnull(name)]))
    return ';'.join(hpo_names)

def merge_hpo_id(row):
    ensg_values = row['gene_name'].split(',')
    hpo_id = omim_subset[omim_subset['gene_symbol'].isin(ensg_values)]['hpo_id'].tolist()
    hpo_id = list(set([str(name) for name in hpo_id if not pd.isnull(name)]))
    return ';'.join(hpo_id)

df['hpo_id'] = df.apply(merge_hpo_id, axis=1)
df['hpo_name'] = df.apply(merge_hpo_names, axis=1)

def exac_id_merge(row):
    ensg_values = row['geneID'].split(',')
    exac_names = pli_subset[pli_subset['gene_id'].isin(ensg_values)]['exac_pLI'].tolist()
    exac_names = list(set([str(name) for name in exac_names if not pd.isnull(name)]))
    return ';'.join(exac_names)

df['gnomAD_pLI'] = df.apply(exac_id_merge, axis=1)

def lof_id_merge(row):
    ensg_values = row['geneID'].split(',')
    lof_names = pli_subset[pli_subset['gene_id'].isin(ensg_values)]['oe_lof'].tolist()
    lof_names = list(set([str(name) for name in lof_names if not pd.isnull(name)]))
    return ';'.join(lof_names)

df['oe_lof'] = df.apply(lof_id_merge, axis=1)

def mis_id_merge(row):
    ensg_values = row['geneID'].split(',')
    mis_names = pli_subset[pli_subset['gene_id'].isin(ensg_values)]['mis_z'].tolist()
    mis_names = list(set([str(name) for name in mis_names if not pd.isnull(name)]))
    return ';'.join(mis_names)

df['mis_z'] = df.apply(mis_id_merge, axis=1)

def hpo_merge(row):
    ensg_values = row['omim_id'].split(',')
    hpo_names = hpo_subset[hpo_subset['database_id'].isin(ensg_values)]['disease_name'].tolist()
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
    return '; '.join(inheritance)

df['OMIM_Inheritance'] = df.apply(extract_inheritance, axis=1)

# Drop duplicates from df to keep unique rows
df = df.drop_duplicates().reset_index(drop=True)

columns_to_exclude = df.columns[df.columns != 'first_omim_link']

df[columns_to_exclude] = df[columns_to_exclude].replace('"', '', regex=True)
df.columns = df.columns.str.replace('"', '')

# Function to categorize row based on the conditions
def categorize_row(row):
    if float(row['l2fc']) < -0.3:
        inheritance_values = str(row['Model_Of_Inheritance']).lower()
        if 'monoallelic' in inheritance_values or 'x-linked' in inheritance_values:
            return 'high'
        if (float(row['pValue']) <= 0.0001) and (pd.notna(row['gnomAD_pLI']) and row['gnomAD_pLI'] != '' and (float(row['gnomAD_pLI']) >= 0.9)):
            if pd.notna(row['omim_id']):
                return 'mid'
            else:
                return 'research'
# Apply the conditional filters to each row in the DataFrame
df['category'] = df.apply(categorize_row, axis=1)



# Replace ',' with '.' in 'pValue' column
df['pValue'] = df['pValue'].astype(str).str.replace('.', ',')
df['padjust'] = df['padjust'].astype(str).str.replace('.', ',')
df['zScore'] = df['zScore'].astype(str).str.replace('.', ',')
df['normcounts'] = df['normcounts'].astype(str).str.replace('.', ',')
df['meanCorrected'] = df['meanCorrected'].astype(str).str.replace('.', ',')
df['theta'] = df['theta'].astype(str).str.replace('.', ',')
df['gnomAD_pLI'] = df['gnomAD_pLI'].astype(str).str.replace('.', ',')
df['oe_lof'] = df['oe_lof'].astype(str).str.replace('.', ',')
df['mis_z'] = df['mis_z'].astype(str).str.replace('.', ',')

print(df)

df.to_csv(args.output_file, sep = "\t",index=False)
