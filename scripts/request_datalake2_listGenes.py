import pathlib
import duckdb  # version 0.10.1
import polars as pl  # version 0.20.7
import pandas as pd
import argparse



# Set up argument parsing for input and output file
parser = argparse.ArgumentParser(description="Process genetic data")
parser.add_argument("--input", required=True, help="Input outrider file path")
parser.add_argument("--fraser", required=True, help="Input fraser file path")
parser.add_argument("--output", required=True, help="Output file path for outrider with custom IDs")
args = parser.parse_args()


# Define constants
data_type = "germline"  # Options: "germline", "somatic"

# Paths
dl_path = pathlib.Path("/data2/seqoia/seqoia_lake/")  # Update path as needed
aggregations_path = dl_path / "aggregations"
annotations_path = dl_path / "annotations"
data_path = dl_path / data_type
partitions_path = data_path / "genotypes" / "partitions"
transmission_path = data_path / "genotypes" / "transmissions"

# Read the Outrider data
#outrider_file = '/data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq/20240718_RUN28_NextSeq_High_16RNASEQ/fastq/../pipeline_v0/outrider/sample_files2/23D2655.outrider.tab'
outrider_file = args.input
outrider = pd.read_csv(outrider_file, sep='\t', dtype={'ID_SPICE': str})
outrider.sort_values('pValue')

outrider = outrider[outrider['pValue'].str.replace(',', '.').astype(float) <= 0.01]
#outrider = outrider.head(10)
#print(outrider)
sample_id_of_kindex = outrider['ID_SPICE'].astype(str).unique().tolist()


# New: Read the Outrider file for comparison
fraser_file = args.fraser
fraser = pd.read_csv(fraser_file, sep='\t', dtype={'gene_name': str})

# Compare Fraser 'gene_name' with Outrider 'gene_name' and add a new 'comparison' column
fraser_genes = set(fraser['gene_name'])  # Genes in Fraser
outrider_genes = set(outrider['gene_name'])  # Genes in Outrider

outrider['comparison'] = outrider['gene_name'].apply(lambda gene: 2 if gene in fraser_genes else 1)




# Check if 'ID_SPICE' is empty
if outrider['ID_SPICE'].isna().all() or outrider['ID_SPICE'].str.strip().eq("").all():
    print("ID_SPICE is empty. Writing the file as it is without custom processing.")
    output_file = args.output
    outrider.to_csv(output_file, sep='\t', index=False)
    exit()



sample_filename = sample_id_of_kindex[0][:-1] + "0"
print(sample_id_of_kindex,sample_filename)
sample_list = [sample_filename]


# Sample file path (you can modify this as needed)
sample_file_path = pathlib.Path(data_path / "genotypes" / "samples" / f"{sample_filename}.parquet")

# Check if the file exists
if not sample_file_path.exists():
    print(f"No files found that match the pattern '{sample_file_path}'. Exiting the script.")
    output_file = args.output
    outrider.to_csv(output_file, sep='\t', index=False)
    exit()  # Exit 


gnomad_ac_threshold = 5

gene_names = outrider['gene_name'].unique().tolist()

gene_names_str = ', '.join(f"'{gene}'" for gene in gene_names)


chromosomes = outrider['chromosome'].unique().tolist()

all_gnomaded_results = []


# Initialize new columns for custom IDs
outrider['custom_ids_high'] = None
outrider['custom_ids_moderate'] = None
outrider['custom_ids_low'] = None
outrider['custom_ids_modifier'] = None

for chrom in chromosomes:
    
    # Query to get variant annotations for all genes in gene_names_str
    query = f"""
        select
            snpeff.id,
            snpeff.impact as snpeff_impact,
            snpeff.hgvs_c as snpeff_hgvs_c,
            snpeff.hgvs_p as snpeff_hgvs_p,
            snpeff.gene as snpeff_gene,
            snpeff.geneid as snpeff_geneid,
            *
        from
            read_parquet('{annotations_path / "snpeff" / data_type / "v4.3t" / f"{chrom}.parquet"}') as snpeff
        where
            snpeff.gene IN ({gene_names_str})
        """
    # Run the query
    snpeffed = duckdb.query(query).pl()
    print(snpeffed)
        # GnomAD query for each chromosome file
    query = f"""
        select
            snpeffed.*, gnomad.AC as gnomad_ac
        from
            snpeffed
        left join
            read_parquet('{annotations_path / "gnomad" / "3.1.2" / f"{chrom}.parquet"}') as gnomad
        on
            snpeffed.id = gnomad.id
    """
    
    # Execute the GnomAD query for the current chromosome and process the results
    gnomaded = (
        duckdb.query(query)
        .pl()
        .with_columns(
            gnomad_ac=pl.col("gnomad_ac").fill_null(0).list.get(0)  
        )
        .filter(pl.col("gnomad_ac") < gnomad_ac_threshold)  
    )
    
    # Append the results to the list
    all_gnomaded_results.append(gnomaded)

# Concatenate the results from all chromosomes
gnomaded_df = pl.concat(all_gnomaded_results)
print(gnomaded_df)
# Register the concatenated gnomaded_df as a DuckDB view so we can use it in the next query
duckdb.register('gnomaded_df', gnomaded_df)

# Variant data query to retrieve the variant info (assuming variants.parquet is a single file)
query = f"""
    select
        v.chr, v.pos, v.ref, v.alt, gnomaded_df.*
    from
        gnomaded_df
    left join
        read_parquet('{data_path / "variants.parquet"}') as v
    on
        gnomaded_df.id = v.id
"""

# Run the final query to get variants
variants = duckdb.query(query).pl()
#print(variants)

# Function to get sample-specific data
def worker(data: pl.DataFrame, target_sample_id: str) -> pl.DataFrame:
    query = f"""
    select
        data.*, sample.*
    from
        data
    join
        read_parquet('{data_path / "genotypes" / "samples" / f"{target_sample_id}.parquet"}') as sample
    on
        data.id == sample.id
    """
    return duckdb.query(query).pl()

results_list = [worker(variants, target_sample_id) for target_sample_id in sample_list]

results = pl.concat(results_list).drop("id_1")  # Remove duplicate columns

print(results)
# Group by 'id' (variant ID) to count unique samples
sample_count_per_variant = (
    results.group_by("id")
    .agg(pl.col("sample").n_unique().alias("sample_count"))  # Count unique samples per variant
)
print(sample_count_per_variant)
results_with_count = results.join(sample_count_per_variant, on="id")

results_with_custom_id = results_with_count.with_columns(
    (
        (pl.col("chr").cast(str) + "_" +
         pl.col("pos").cast(str) + "_" +
         pl.col("ref") + "_" +
         pl.col("alt") + ":" +
         pl.col("snpeff_hgvs_c") +
         pl.when(pl.col("sample_count") == 1).then(pl.lit("(1)")).otherwise(pl.lit("(2)"))).alias("custom_id")
    )
)

filtered_results = results_with_custom_id.filter(
    pl.col("sample").is_in(sample_id_of_kindex)
)

#print(sample_id_of_kindex)
#print(filtered_results)

for index, row in outrider.iterrows():
    
    gene_name = row['gene_name']


# Filter the results for the current gene_name (matching snpeff_gene)
    gene_filtered_results = filtered_results.filter(pl.col("snpeff_gene").cast(pl.Utf8) == str(gene_name))
    #print(gene_filtered_results)
    # Filter for high-impact variants
    high_impact_results = gene_filtered_results.filter(pl.col("snpeff_impact") == "HIGH")
    id_list_high = high_impact_results["custom_id"].to_list()

    # Filter for moderate impact variants
    moderate_impact_results = gene_filtered_results.filter(pl.col("snpeff_impact") == "MODERATE")
    id_list_moderate = moderate_impact_results["custom_id"].to_list()

    # Filter for low-impact variants
    low_impact_results = gene_filtered_results.filter(pl.col("snpeff_impact") == "LOW")
    id_list_low = low_impact_results["custom_id"].to_list()
     
    # Filter for modifier-impact variants
    modifier_impact_results = gene_filtered_results.filter(pl.col("snpeff_impact") == "MODIFIER")
    id_list_modifier = modifier_impact_results["custom_id"].to_list()


    outrider.at[index, 'custom_ids_high'] = id_list_high if id_list_high else None
    outrider.at[index, 'custom_ids_moderate'] = id_list_moderate if id_list_moderate else None
    outrider.at[index, 'custom_ids_low'] = id_list_low if id_list_low else None
    outrider.at[index, 'custom_ids_modifier'] = id_list_modifier if id_list_modifier else None

#output_file = '23D2655_outrider_with_custom_ids.tsv'
output_file = args.output
outrider.to_csv(output_file, sep='\t', index=False)

