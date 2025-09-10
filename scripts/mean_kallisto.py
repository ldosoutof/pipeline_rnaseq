import os
import pandas as pd
import re
from statistics import mean, stdev
import subprocess


##  Calculate mean for sample of different type (puroplus, puromoins, paxgene)

def parse_gtf(gtf_file):
    gene_mapping = {}
    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            
            data = line.strip().split('\t')
            if len(data)<9:
                continue		    
	
            if data[2] == 'gene':
                attributes = data[8].split('; ')
                gene_id = attributes[0].split(' ')[1].strip('"')
                gene_name = attributes[4].split(' ')[1].strip('"')
                gene_mapping[gene_id] = gene_name
        return gene_mapping

def ensembl_to_symbol_from_gtf(gtf_file, ensembl_ids):
    ensembl_to_symbol = {}

    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            data = line.strip().split('\t')

            if len(data)<9:
                continue

            if data[2] == 'gene':
                attributes = data[8].split(';')
                gene_id = None
                gene_name = None
                for attr in attributes:
                    if 'gene_id' in attr:
                        gene_id = attr.split('"')[1]
                    elif 'gene_name' in attr:
                        gene_name = attr.split('"')[1]
                if gene_id and gene_name:
                    ensembl_to_symbol[gene_id] = gene_name

    ensembl_ids = set(ensembl_ids)
    ensembl_to_symbol_filtered = {ensembl_id: ensembl_to_symbol.get(ensembl_id) for ensembl_id in ensembl_ids}
    return ensembl_to_symbol_filtered



# Example usage:
gtf_path = '/data2/laura_home/test_data/Homo_sapiens.GRCh38.106.gtf'
#ensembl_to_hugo_mapping = parse_gtf(gtf_file_path)

def run_tximport(abundance_file, sample_path, sample_name):
    r_script = f'''
    library(GenomicFeatures)
    library(tximport)
    library(ggplot2)
    library(rhdf5)

    txdb <- makeTxDbFromGFF(file="/data2/Exome_analysis_BC/Exome_pipeline_datas/GRCh38/Homo_sapiens.GRCh38.106.gtf")
    k <- keys(txdb, keytype = "TXNAME")
    tx2gene <- select(txdb, k, "GENEID", "TXNAME")
    #print("{abundance_file}")
   
    txi <- "{abundance_file}"

    # Run tximport for gene-level counts
    txi_gene <- tximport(txi, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
    
    # Return the gene-level counts or save them to a file
    write.table(txi_gene$counts, file="{sample_path}/abundance_gene_level_counts.tsv", sep="\t", quote=FALSE)
    '''

    with open('tximport_script.R', 'w') as file:
        file.write(r_script)

    try:
        subprocess.run(['/home/ldosouto/miniconda3/envs/test/bin/Rscript', 'tximport_script.R'], check=True)
    except subprocess.CalledProcessError as e:
        print("R script execution failed!")
        print(e.output.decode("utf-8"))


# Function to extract sample names with "puromoins"
def aggregate_tpm(root_directory, type_cult):
    gene_tpm = {}
    with open('/data2/laura_home/pipelines/bam_config_xths1.txt', 'r') as sample_file:
        next(sample_file)  # Skip header
        sample_names_in_file = {line.split()[1] for line in sample_file}


    for root, dirs, files in os.walk(root_directory):
        if 'pipeline_v0' in dirs:  # Checking for 'pipeline_v0' folders
            pipeline_path = os.path.join(root, 'pipeline_v0')
            #print(pipeline_path)    
            for subroot, subdirs, subfiles in os.walk(pipeline_path):
                if 'kallisto' in subdirs:  # Checking for 'kallisto' folders inside 'pipeline_v0'
                    kallisto_path = os.path.join(subroot, 'kallisto')#dirs.remove('kallisto')  # Prevent descending into 'kallisto' folders
                    #print(kallisto_path)
                    for subroot, subdirs, subfiles in os.walk(kallisto_path):
                        #print(subdirs)
                        for dir_sample in subdirs:
                            #print(dir_sample)
                            if not 'ONCO' in subroot:			    
                                if type_cult in dir_sample or type_cult.lower() in dir_sample:
                                    #print('puromoins')
                                    sample_path =  os.path.join(subroot,dir_sample)
                                    sample_name = sample_path.split('/')[-1]
                                    print(sample_path)    
                                    if sample_name in sample_names_in_file:
                                        abundance_file = os.path.join(sample_path, 'abundance_gene_level_counts.tsv')
                                        if not os.path.isfile(abundance_file):
                                            abund_file = os.path.join(sample_path, 'abundance.h5')
                                            print(abund_file)
                                            sample_name = sample_path.split('/')[-1]
                                            run_tximport(abund_file, sample_path, sample_name)
                                            #sample_name = sample_path.split('/')[-1]
        #print(abundance_file)
                                        sample_name = sample_path.split('/')[-1]
                                        if os.path.exists(abundance_file):
                                            run_id = subroot.split('/')[-3]  # Extracting run_id from path
                                            print(run_id)
                                            sample_name = sample_path.split('/')[-1]  # Extracting sample name from path
                                            print(sample_name)
                                            df = pd.read_csv(abundance_file, sep='\t')
                                            print(df)
                                            for gene_id, tpm_value in df.iterrows():
                                                #print(gene_id, tpm_value['V1'])
                                                tpm = tpm_value['V1']
                                                if gene_id not in gene_tpm:
                                                    gene_tpm[gene_id] = [tpm]
                                                else:
                                                    gene_tpm[gene_id].append(tpm)
    gene_stats = {}
    #print(gene_tpm)
    for gene, tpms in gene_tpm.items():
        mean_tpm = mean(tpms)
        sd_tpm = stdev(tpms) if len(tpms) > 1 else 0
        gene_stats[gene] = {'Mean TPM': mean_tpm, 'SD TPM': sd_tpm, 'TPM':tpms}

    return pd.DataFrame.from_dict(gene_stats, orient='index')

# Directory where to search for 'kallisto' folders
root_directory = '/data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq/'

#df_paxgene = aggregate_tpm(root_directory, 'PAX')


## Assuming df_paxgene is your DataFrame with Ensembl IDs in the index
#ensembl_ids_list = df_paxgene.index.tolist()
#ensembl_to_symbol_mapping = ensembl_to_symbol_from_gtf(gtf_path, ensembl_ids_list)
#
## Update your DataFrame with gene symbols
#df_paxgene['HUGO_symbol'] = df_paxgene.index.map(ensembl_to_symbol_mapping)
##print(df_paxgene)

##print(df_paxgene)
#df_paxgene.to_csv("/data2/laura_home/pipelines/mean_sd_kallisto_paxgene.tsv",sep = "\t",index=True)


df_puromoins = aggregate_tpm(root_directory, 'PUROMOINS')
ensembl_ids_list = df_puromoins.index.tolist()
ensembl_to_symbol_mapping = ensembl_to_symbol_from_gtf(gtf_path, ensembl_ids_list)

df_puromoins['HUGO_symbol'] = df_puromoins.index.map(ensembl_to_symbol_mapping)
df_puromoins.to_csv("/data2/laura_home/pipelines/mean_sd_kallisto_puromoins.tsv",sep = "\t",index=True)

df_puroplus = aggregate_tpm(root_directory, 'PUROPLUS')
ensembl_ids_list = df_puroplus.index.tolist()
ensembl_to_symbol_mapping = ensembl_to_symbol_from_gtf(gtf_path, ensembl_ids_list)

df_puroplus['HUGO_symbol'] = df_puroplus.index.map(ensembl_to_symbol_mapping)
#print(df_puroplus)
df_puroplus.to_csv("/data2/laura_home/pipelines/mean_sd_kallisto_puroplus.tsv",sep = "\t",index=True)

