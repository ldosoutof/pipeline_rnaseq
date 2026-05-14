import os
import pandas as pd
import argparse
import numpy as np

from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

# Create a workbook and select the active worksheet
wb = Workbook()
ws = wb.active


parser = argparse.ArgumentParser(description="Annotate Fraser file")
parser.add_argument('-f', '--fraser', required=True)
parser.add_argument('-g', '--gtf', required=True)
parser.add_argument('-d', '--di', required=True)
parser.add_argument('-o', '--omim', required=True)
parser.add_argument('-p', '--pli', required=True)
parser.add_argument('-m', '--hpo', required=True)

parser.add_argument('--output', required=True, help='Path to the output file')

args = parser.parse_args()


di = pd.read_csv(args.di, sep='\t', na_values=[''])
di_subset = di[['EnsemblId(GRch38)', 'Model_Of_Inheritance', 'Phenotypes']]
#print(di_subset)
omim= pd.read_csv(args.omim, sep='\t', na_values=[''])
omim['omim_id']=omim['disease_id'].str.extract(r'(OMIM:[0-9]+)')
omim_subset = omim[['gene_symbol', 'hpo_id', 'hpo_name', 'omim_id']]
#print(omim)
pli= pd.read_csv(args.pli, sep='\t', na_values=[''])
pli_subset = pli[['gene_id', 'exac_pLI', 'oe_lof','mis_z']]
#print(pli)
hpo= pd.read_csv(args.hpo, sep='\t', na_values=[''], comment='#')
hpo_subset=hpo[['database_id','disease_name']]
#print(fraser)


# Load the Excel file (replace 'your_file.xlsx' with the actual file path)
file_path = '/datawork2/genetique/RNASeq/diag/prod/Results_LymphoRNA_Magnis_NS.xlsx'

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
            # Extract gene and transcript IDs
            gene_id = None
            transcript_id = None
            gene_name = None    
            for attribute in attributes.split(';'):
               if 'gene_id' in attribute:
                    gene_id = attribute.split('"')[1]
               elif 'transcript_id' in attribute:
                    #print(attribute.split('"')[1])
                    transcript_id = attribute.split('"')[1]
               elif 'gene_name' in attribute:
                    #print(attribute.split('"')[1])
                    gene_name = attribute.split('"')[1]       
            # Store mapping
            if chromosome not in gtf_mapping:
                gtf_mapping[chromosome] = []
            gtf_mapping[chromosome].append({
                    'start': start_pos,
                    'end': end_pos,
                    'gene_name': gene_name,    
                    'gene_id': gene_id,
                    'transcript_id': transcript_id
                })
            #print(gtf_mapping)


# Create empty lists to store data
rows = []

# Step 2: Iterate through the Fraser file and find matching positions
with open(args.fraser, 'r') as fraser_file:
    lines = fraser_file.readlines()
    header = lines[0].strip().split('\t')
    header.extend(['ENSG', 'gene_name', 'ID_SPICE', 'ucsc_link'])
    link=''
    for line in lines[1:]:
        data = line.strip().split('\t')
        #print(data)	
        chromosome = data[1].split('"')[1]
        start_pos = int(data[2])
        end_pos = int(data[3])
        key_id = chromosome+":"+str(start_pos)+"-"+str(end_pos)
        #type_fraser = data[9].split('"')[1]   #for FRASER
        type_fraser = data[8]	# for FRASER2
        #pValue = float(data[10])
        pValue = float(data[9])	# for FRASER2
        #print(pValue)
        sample_id = data[6].split('"')[1]
        matching_transcripts =[] 
        matching_genes = []
        matching_names = []
        if pValue < 0.001:
            if chromosome in gtf_mapping and len(chromosome) <=2:
                for entry in gtf_mapping[chromosome]:
                    #print(str(pValue)+"entry")			
                    if entry['start'] <= start_pos <= entry['end'] or entry['start'] <= end_pos <= entry['end']:
                        #print(str(start_pos)+"start")
                        #print(str(end_pos)+"end")			
                        #print(str(pValue)+"in")
                        if type_fraser=='theta':
                            start=start_pos-100
                            end=end_pos+100
                            link_base = "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?pix=1243&textSize=12&textFont=Helvetica&hgt.labelWidth=20"
                            link_params = "&hgS_otherUserName=View&hgS_otherUserSessionName=Clinical_SNVs_hg38&hgS_otherUserSessionLabel=Clinical%20SNVs&hgS_otherUserSessionDesc="
                            position = f"chr{chromosome}:{start}-{end}"
                            full_link = f"{link_base}&{link_params}{position}"

			    # Create the Excel formula using CONCATENATE
                        
                            link=f'"https://genome-euro.ucsc.edu/cgi-bin/hgTracks?pix=1243&textSize=12&textFont=Helvetica&hgt.labelWidth=20&hgS_otherUserName=View&hgS_otherUserSessionName=Clinical_SNVs_hg38&hgS_otherUserSessionLabel=Clinical%20SNVs&hgS_otherUserSessionDesc=Assess%20potential%20disease%20contributions%20of%20single%20nucleotide%20variants%20in%20coding%20regions&hgS_doOtherUser=submit&position=chr{chromosome}:{start}-{end}&highlight=chr{key_id}"'
                            #print(link)			    
                        else:
                            length=end_pos-start_pos
                            start=start_pos-length
                            end=end_pos+length
                            #link=f'"https://genome-euro.ucsc.edu/cgi-bin/hgTracks?pix=1243&textSize=12&textFont=Helvetica&hgt.labelWidth=20&hgS_otherUserName=View&hgS_otherUserSessionName=Clinical_SNVs_hg38&hgS_otherUserSessionLabel=Clinical%20SNVs&hgS_otherUserSessionDesc=Assess%20potential%20disease%20contributions%20of%20single%20nucleotide%20variants%20in%20coding%20regions&hgS_doOtherUser=submit&position=chr{chromosome}:{start}-{end}&highlight=chr{key_id}"'
                            link=f'"https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&position=chr{chromosome}:{start}-{end}&highlight=chr{key_id}";"{key_id}"'
                        if not entry['gene_id'] in matching_genes:
                            matching_genes.append(entry['gene_id'])
                        if not entry['gene_name'] in matching_names:                 
                            matching_names.append(entry['gene_name'])
                data.extend([','.join(matching_genes)])
                #data.extend([','.join(filter(None,matching_names))])
                id_spice = get_id_spice(sample_id)
                data.extend([','.join(filter(None, matching_names)), id_spice or ''])                
                #print('\t'.join(data))
                data.extend([link])
                # Retrieve the corresponding ID_SPICE based on the sample ID
                rows.append(data[1:])  
                #print(header)
    df = pd.DataFrame(rows, columns=header)
    #print(df)
    #print(header)
# Insert the 'ucsc_link' column after the 'end' column
ucsc_link_col = df.pop('ucsc_link')  # Remove the column to reinsert it
df.insert(df.columns.get_loc('"end"') + 1, 'ucsc_link', ucsc_link_col)  # Insert at 'end' column position + 1

def merge_model(row):
    ensg_values = row['ENSG'].split(',')  # Split multiple ENSG values in a row
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

#def create_omim_links(omim_ids):
#    links = []
#    for omim_id in omim_ids.split(','):
#        parts = omim_id.strip().split(':')
#        if len(parts) == 2:
#            links.append(f'https://www.omim.org/entry/{parts[1]}')
#        else:
#            links.append(f"Invalid OMIM ID: {omim_id}")
#    return links
#
#def omim_id_merge_with_links(row):
#    ensg_values = row['gene_name'].split(',')
#    omim_ids = omim_subset[omim_subset['gene_symbol'].isin(ensg_values)]['omim_id'].tolist()
#    omim_ids = list(set([str(name) for name in omim_ids if not pd.isnull(name)]))
#    omim_links = create_omim_links(','.join(omim_ids))
#    return ','.join(omim_ids), ','.join(omim_links)
#
#df['omim_id'], df['omim_links'] = zip(*df.apply(omim_id_merge_with_links, axis=1))

## Assuming 'omim_id' column contains OMIM IDs in your DataFrame 'df'
#df['OMIM_Links'] = df['omim_id'].apply(create_omim_links)
#
#
#df['OMIM_Links_String'] = df['OMIM_Links'].apply(lambda x: ','.join(x))
#
## Drop duplicates based on the new string column
#df_unique = df.drop_duplicates(subset='OMIM_Links_String')
#
## Drop the intermediate string column if not needed
#df_unique = df_unique.drop(columns=['OMIM_Links_String'])
#
def merge_hpo_names(row):
    ensg_values = row['gene_name'].split(',')  # Split multiple ENSG values in a row
    hpo_names = omim_subset[omim_subset['gene_symbol'].isin(ensg_values)]['hpo_name'].tolist()
    hpo_names = list(set([str(name) for name in hpo_names if not pd.isnull(name)]))
    return ';'.join(hpo_names)

def merge_hpo_id(row):
    ensg_values = row['gene_name'].split(',')  # Split multiple ENSG values in a row
    hpo_id = omim_subset[omim_subset['gene_symbol'].isin(ensg_values)]['hpo_id'].tolist()
    hpo_id = list(set([str(name) for name in hpo_id if not pd.isnull(name)]))
    return ';'.join(hpo_id)
df['hpo_id'] = df.apply(merge_hpo_id, axis=1)

# Apply the custom function to each row to merge hpo_names based on ENSG
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


#merged_df = pd.merge(merged_df, hpo_subset, how='left', left_on='omim_id', right_on='database_id')

def hpo_merge(row):
    ensg_values = row['omim_id'].split(',')
    hpo_names = hpo_subset[hpo_subset['database_id'].isin(ensg_values)]['disease_name'].tolist()
    hpo_names = list(set([str(name) for name in hpo_names if not pd.isnull(name)]))
    return ';'.join(hpo_names)

df['disease_name'] = df.apply(hpo_merge, axis=1)

# Your mapping
inheritance_mapping = {
    'HP:0000006': 'Autosomal Dominant',
    'HP:0012275': 'Autosomal Dominant with Imprinting (Maternal)',
    'HP:0012274': 'Autosomal Dominant with Imprinting (Paternal)',
    'HP:0000007': 'Autosomal Recessive',
    'HP:0001417': 'X-Linked',
    'HP:0001419': 'X-Linked Recessive',
    'HP:0001423': 'X-Linked Dominant'
    }

# Function to extract inheritance patterns

def extract_inheritance(row):
    hpo_ids = row['hpo_id'].split(',')
    inheritance = [inheritance_mapping[hpo] for hpo in hpo_ids if hpo in inheritance_mapping]
    return ';'.join(inheritance)
   
   
   
   
   
   

df['OMIM_Inheritance'] = df.apply(extract_inheritance, axis=1)

# Drop duplicates from df to keep unique rows
df = df.drop_duplicates().reset_index(drop=True)

columns_to_exclude = df.columns[df.columns != 'first_omim_link']

# Perform df.replace on all columns except 'first_omim_link'
df[columns_to_exclude] = df[columns_to_exclude].replace('"', '', regex=True)
#df = df.replace('"', '', regex=True)
df.columns = df.columns.str.replace('"', '')
print(df)


def categorize_row(row):
    if float(row['deltaPsi']) >= 0.2:
        inheritance_values = str(row['Model_Of_Inheritance']).lower()
        if 'monoallelic' in inheritance_values or 'x-linked' in inheritance_values:
            return 'high'
        gnomAD_pLI_values = str(row['gnomAD_pLI']).split(';')
        gnomAD_pLI = float(gnomAD_pLI_values[0]) if gnomAD_pLI_values[0] else None
        #print(gnomAD_pLI, gnomAD_pLI_values)
        max_pLI=max(gnomAD_pLI_values)
        if (float(row['pValue']) <= 0.0001) and (pd.notna(row['gnomAD_pLI']) and row['gnomAD_pLI'] != '' and (float(max_pLI) >= 0.9)):
            if pd.notna(row['omim_id']):
                return 'mid'
            else:
                return 'research'
	    
# Apply the conditional filters to each row in the DataFrame
df['category'] = df.apply(categorize_row, axis=1)


# Replace ',' with '.' in 'pValue' column
df['pValue'] = df['pValue'].astype(str).str.replace('.', ',')

df['padjust'] = df['padjust'].astype(str).str.replace('.', ',')
#df['zScore'] = df['zScore'].astype(str).str.replace('.', ',') # for FRASER
df['psiValue'] = df['psiValue'].astype(str).str.replace('.', ',')
df['deltaPsi'] = df['deltaPsi'].astype(str).str.replace('.', ',')
df['meanCounts'] = df['meanCounts'].astype(str).str.replace('.', ',')
df['meanTotalCounts'] = df['meanTotalCounts'].astype(str).str.replace('.', ',')
df['gnomAD_pLI'] = df['gnomAD_pLI'].astype(str).str.replace('.', ',')
df['oe_lof'] = df['oe_lof'].astype(str).str.replace('.', ',')
df['mis_z'] = df['mis_z'].astype(str).str.replace('.', ',')

#def categorize_row(row):
#    if row['delta_psi'] >= 0.2:
#        if df['Model_Of_Inheritance'].str.contains('monoallelic|x-linked', case=False):
#            return 'high'
#        elif (df['pValue'] <= 0.0001) & (df['gnomAD_pLI'] >= 0.9) & (df['omim_id'].notna()):
#            return 'mid'
#        else:
#            return 'research'

# Apply the conditional filters to each row in the DataFrame
#df['new_column'] = df.apply(categorize_row, axis=1)
print(df)
df.to_csv(args.output,sep = "\t",index=False)



#for r_idx, row in enumerate(dataframe_to_rows(df, index=False, header=True), 1):
#    for c_idx, value in enumerate(row, 1):
#        cell = ws.cell(row=r_idx, column=c_idx, value=value)
#
#wb.save('output.xlsx')







##Splittinggthe rows with multiple 'ENSG' values into separate rows
#df['ENSG'] = df['ENSG'].str.split(',')
#df = df.explode('ENSG')
#
## Initialize an empty dictionary to store 'ENSG' and 'hpo_name' associations
#ensg_hpo_dict = {}
#for index, row in df.iterrows():
#    ensg_values = row['gene_name']
#    hpo_names = []
#    for ensg in ensg_values:
#        matching_hpo = omim_subset.loc[omim_subset['gene_symbol'] == ensg, 'hpo_name'].tolist()
#        hpo_names.extend(matching_hpo)
#        ensg_hpo_dict[row['ENSG']] = hpo_names
#
## Update the DataFrame with aggregated 'hpo_name' values for each 'ENSG'
#df['hpo_name'] = df['ENSG'].map(lambda x: ','.join(set(ensg_hpo_dict[x])) if x in ensg_hpo_dict else '')
#
#print(df)
#
#
#
## Update the DataFrame with aggregated 'hpo_name' values for each 'ENSG'
#df['hpo_name'] = df['gene_name'].map(lambda x: ','.join(set(ensg_hpo_dict[x])) if x in ensg_hpo_dict else '')
#
## Group by the initial index to join 'ENSG' values back into single rows
#result_df = df.groupby(df.index).agg({
#    'seqnames': 'first',
#    'start': 'first',
#    'end': 'first',
#    'width': 'first',
#    'counts': 'first',
#    'totalCounts': 'first',
#    'ENSG': ','.join,
#    'gene_name': ','.join,
#    'hpo_name': 'first'  # Take the first value of hpo_name, as it has been aggregated
#    }).reset_index(drop=True)
#
#print(result_df)
#
#
#
#
#
#
#
#
#
#
#df['ENSG'] = df['ENSG'].str.split(',')
#df = df.explode('ENSG')
#merged_df = pd.merge(df, di_subset, how='left', left_on='ENSG', right_on='EnsemblId(GRch38)')
#
#merged_df['Model_Of_Inheritance'] = merged_df['Model_Of_Inheritance'].astype(str)
#merged_df['Model_Of_Inheritance'] = merged_df.groupby(merged_df.index)['Model_Of_Inheritance'].agg(','.join)
#
#
## Function to merge hpo_names for multiple ENSG values in a row
#def merge_hpo_names(row):
#    gene_values = row['gene_name'].split(',')  # Split ENSG values in the row
#    hpo_names = omim_subset[omim_subset['gene_symbol'].isin(gene_values)]['hpo_name'].tolist()
#    return ','.join(set(hpo_names))  # Join unique hpo_names using a comma
#
#    # Applying the function to each row
#df['hpo_name'] = df.apply(merge_hpo_names, axis=1)
#
#    # Dropping duplicates to keep unique rows
#df = df.drop_duplicates()
#
#print(df)
#
#
#

#merged_df['Phenotypes'] = merged_df['Phenotypes'].astype(str)
#merged_df['Phenotypes'] = merged_df.groupby(merged_df.index)['Phenotypes'].agg(','.join)
#
#
#merged_df = pd.merge(merged_df, omim_subset, how='left', left_on='gene_name', right_on='gene_symbol')
#
#merged_df['omim_id'] = merged_df['omim_id'].astype(str)
#merged_df['omim_id'] = merged_df.groupby(merged_df.index)['omim_id'].agg(','.join)
#
#merged_df['hpo_name'] = merged_df['hpo_name'].astype(str)
#merged_df['hpo_name'] = merged_df.groupby(merged_df.index)['hpo_name'].agg(','.join)
#
##'gene_id', 'exac_pLI', 'oe_lof','mis_z'
#merged_df = pd.merge(merged_df, pli_subset, how='left', left_on='ENSG', right_on='gene_id')
#
#merged_df['exac_pLI'] = merged_df['exac_pLI'].astype(str)
#merged_df['exac_pLI'] = merged_df.groupby(merged_df.index)['exac_pLI'].agg(','.join)
#
#merged_df['oe_lof'] = merged_df['oe_lof'].astype(str)
#merged_df['oe_lof'] = merged_df.groupby(merged_df.index)['oe_lof'].agg(','.join)
#
#merged_df['mis_z'] = merged_df['mis_z'].astype(str)
#merged_df['mis_z'] = merged_df.groupby(merged_df.index)['mis_z'].agg(','.join)
#
##'database_id','disease_namei'
#
#merged_df = pd.merge(merged_df, hpo_subset, how='left', left_on='omim_id', right_on='database_id')
#
#merged_df['disease_name'] = merged_df['disease_name'].astype(str)
#merged_df['disease_name'] = merged_df.groupby(merged_df.index)['disease_name'].agg(','.join)
#print(merged_df)
#merged_df.to_csv("/data2/laura_home/pipelines/fraser_tab_test_annot.tsv",sep = "\t",index=False)
### Combining additional_info for rows with the same index
#merged_df = merged_df.groupby(merged_df.index).agg({
#    '"seqnames"': 'first',
#    '"start"': 'first', 
#    '"end"': 'first', 
#    '"width"': 'first', 
#    '"strand"': 'first', 
#    '"sampleID"': 'first',
#    '"hgncSymbol"': 'first', 
#    '"addHgncSymbols"': 'first',
#    '"type"': 'first', 
#    '"pValue"': 'first', 
#    '"padjust"': 'first', 
#    '"zScore"': 'first', 
#    '"psiValue"': 'first', 
#    '"deltaPsi"': 'first', 
#    '"meanCounts"': 'first', 
#    '"meanTotalCounts"': 'first', 
#    '"counts"': 'first', 
#    '"totalCounts"': 'first',
#    'ENSG': lambda x: ','.join(x.unique()),
#    'gene_name': lambda x: ','.join(x.unique()),
#    'Model_Of_Inheritance': lambda x: ','.join(x.unique())
#   # 'hpo_name':lambda x: ','.join(x.unique()),
#   # 'omim_id':lambda x: ','.join(x.unique()),
#   # 'exac_pLI':lambda x: ','.join(x.unique()),
#   # 'oe_lof':lambda x: ','.join(x.unique()),
#   # 'mis_z':lambda x: ','.join(x.unique()),
#   # 'disease_name':lambda x: ','.join(x.unique())
#    }).reset_index(drop=True)
#merged_df = merged_df.drop_duplicates()
#print(merged_df)
### Grouping the additional_info based on the original index
###merged_df['Model_Of_Inheritance'] = merged_df['Model_Of_Inheritance'].astype(str)
###merged_df['Model_Of_Inheritance'] = merged_df.groupby(merged_df.index)['Model_Of_Inheritance'].agg(','.join)
### Displaying the merged DataFrame
##G#print(merged_df['Model_Of_Inheritancenitial index to join 'ENSG' values back into single rows
#result_df = df.groupby(df.index).agg({
#    'seqnames': 'first',
#    'start': 'first',
#    'end': 'first',
#    'width': 'first',
#    'counts': 'first',
#    'totalCounts': 'first',
#    'ENSG': ','.join,
#    'gene_name': ','.join,
#    'hpo_name': 'first'  # Take the first value of hpo_name, as it has been aggregated
#    }).reset_index(drop=True)
#
#print(result_df)
