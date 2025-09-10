import pathlib
import duckdb  # version 0.10.1
import polars as pl  # version 0.20.7
import pandas as pd
import altair
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import argparse



THREADS = 10

altair.data_transformers.disable_max_rows()

# Define constants
data_type = "germline"  # Options: "germline", "somatic"

# Paths
dl_path = pathlib.Path("/data2/seqoia/seqoia_lake/")  # Update path as needed
aggregations_path = dl_path / "aggregations"
annotations_path = dl_path / "annotations"
data_path = dl_path / data_type
partitions_path = data_path / "genotypes" / "partitions"
transmission_path = data_path / "genotypes" / "transmissions"



# Set up argument parsing for input and output file
parser = argparse.ArgumentParser(description="Process genetic data")
parser.add_argument("--input", required=True, help="Input Fraser file path")
parser.add_argument("--outrider", required=True, help="Input Outrider file path")
parser.add_argument("--output", required=True, help="Output file path for Fraser with custom IDs")
args = parser.parse_args()

# Read the fraser data
#fraser_file = '/data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq/20240718_RUN28_NextSeq_High_16RNASEQ/fastq/../pipeline_v0/fraser/sample_files/24D0157.fraser.tab'
fraser_file = args.input
fraser = pd.read_csv(fraser_file, sep='\t', dtype={'ID_SPICE': str})
fraser.sort_values('pValue')

# Filter genes based on pValue
fraser = fraser[fraser['pValue'].str.replace(',', '.').astype(float) <= 0.01]
sample_id_of_kindex = fraser['ID_SPICE'].astype(str).unique().tolist()
print("Sample ID(s):", sample_id_of_kindex)





# New: Read the Outrider file for comparison
outrider_file = args.outrider
outrider = pd.read_csv(outrider_file, sep='\t', dtype={'gene_name': str})

# Compare Fraser 'gene_name' with Outrider 'gene_name' and add a new 'comparison' column
fraser_genes = set(fraser['gene_name'])  # Genes in Fraser
outrider_genes = set(outrider['gene_name'])  # Genes in Outrider

fraser['comparison'] = fraser['gene_name'].apply(lambda gene: 2 if gene in outrider_genes else 1)





# Check if 'ID_SPICE' is empty
if fraser['ID_SPICE'].isna().all() or fraser['ID_SPICE'].str.strip().eq("").all():
    print("ID_SPICE is empty. Writing the file as it is without custom processing.")
    output_file = args.output
    fraser.to_csv(output_file, sep='\t', index=False)
    exit()





#print(fraser)
gnomad_ac_threshold = 5

# Initialize new columns for custom IDs
fraser['custom_ids_high'] = None
fraser['custom_ids_moderate'] = None
fraser['custom_ids_low'] = None
fraser['custom_ids_modifier'] = None

# Store results for all rows
all_gnomaded_results = []

# Prepare the sample filename

## Query all variants for the sample
#all_variants_query = f"""
#select
#    v.*, g.sample, g.gt
#from
#    read_parquet('{ data_path / "genotypes" / "samples" / f"{sample_filename}.parquet" }') as g
#join
#    read_parquet('{ data_path / "variants.parquet" }') as v
#on
#    v.id == g.id
#where
#    g.sample = '{sample_id_of_kindex[0]}' OR g.sample = '{sample_filename1}' OR g.sample = '{sample_filename2}'
#"""
##g.sample = '{sample_filename1}' OR g.sample = '{sample_filename2}'
#
#all_variants=duckdb.query(all_variants_query).pl()

# Prepare the sample filename
sample_filename = sample_id_of_kindex[0][:-1] + "0"
sample_filename1 = sample_id_of_kindex[0][:-1] + "1"
sample_filename2 = sample_id_of_kindex[0][:-1] + "2"

# Sample file path (you can modify this as needed)
sample_file_path = pathlib.Path(data_path / "genotypes" / "samples" / f"{sample_filename}.parquet")


# Check if the file exists
if not sample_file_path.exists():
    print(f"No files found that match the pattern '{sample_file_path}'. Exiting the script.")
    output_file = args.output
    fraser.to_csv(output_file, sep='\t', index=False)
    exit()  # Exit


# Query for sample_id_of_kindex[0] first (main sample)
main_query = f"""
select
    v.*, g.sample, g.gt
from
    read_parquet('{ data_path / "genotypes" / "samples" / f"{sample_filename}.parquet" }') as g
join
    read_parquet('{ data_path / "variants.parquet" }') as v
on
    v.id == g.id
where
    g.sample = '{sample_id_of_kindex[0]}'
"""

main_variants = duckdb.query(main_query).pl()
#print(main_variants)
# Now query for variants in sample_filename1 and sample_filename2, and filter to keep only those in the main sample
combined_query = f"""
select
    v.*, g.sample, g.gt
from
    read_parquet('{ data_path / "genotypes" / "samples" / f"{sample_filename}.parquet" }') as g
join
    read_parquet('{ data_path / "variants.parquet" }') as v
on
    v.id == g.id
where
    g.sample IN ('{sample_filename1}', '{sample_filename2}')
and
    v.id in (select id from main_variants)  -- Keep only variants that are present in the main sample
"""

combined_variants = duckdb.query(combined_query).pl()

# Sort the combined_variants DataFrame by 'id'
sorted_combined_variants = combined_variants.sort("id")

# Print the sorted DataFrame
print(sorted_combined_variants)

# You can now merge these results or perform any further filtering as needed

print(combined_variants)

# Concatenate the two DataFrames
all_variants = pl.concat([main_variants, combined_variants], how="vertical")

#results_with_count_columns = all_variants.select(["pos", "alt", "ref", "chr", "id","sample"])
#print(results_with_count_columns)

intervals_df = fraser[['seqnames', 'start', 'end']].rename(columns={'seqnames': 'chromosome'})

# Helper function to process each target interval
def process_interval(target):
    print(f"Processing interval: {target['chromosome']}:{target['start']}-{target['end']}")
    # Query the variants within the interval
    query = f"""
    select
        a.*
    from
        all_variants as a
    where
        a.chr = '{target['chromosome']}'
    and
        a.pos > {target['start']}-1000
    and
        a.pos < {target['end']}+1000
    """
    variants = duckdb.query(query).pl()

    # Return early if no variants are found
    if variants.is_empty():
        return None

    # Build and run the gnomad query
    query = f"""
    select
      variants.*, gnomad.AC as gnomad_ac
    from
      variants
    left join
      read_parquet('{annotations_path / "gnomad" / "3.1.2" / f"{target['chromosome']}.parquet"}') as gnomad
    on
      variants.id = gnomad.id
    """
    gnomaded = (
        duckdb.query(query).pl()
        .with_columns(
            gnomad_ac=pl.col("gnomad_ac").fill_null([0]).list.get(0)
        )
        .filter(pl.col("gnomad_ac") < gnomad_ac_threshold)
    )

    # If gnomaded is empty, return None
    if gnomaded.is_empty():
        return None

    # Build and run the snpeff query
    query = f"""
    select
      gnomaded.*,
      snpeff.impact as snpeff_impact, snpeff.hgvs_c as snpeff_hgvs_c
    from
      gnomaded
    left join
      read_parquet('{annotations_path / "snpeff" / data_type / "v4.3t" / f"{target['chromosome']}.parquet"}') as snpeff
    on
      gnomaded.id = snpeff.id
    """
    annotated = duckdb.query(query).pl()

    return annotated
# Use sequential processing for debugging
#results = []
#for target in tqdm(intervals_df.to_dict('records')):
#    result = process_interval(target)
#    if result is not None:
#        results.append(result)
#

# Initialize an empty Polars DataFrame to store results
results = pl.DataFrame()

# Iterate through intervals and process them
for target in tqdm(intervals_df.to_dict('records')):
    result = process_interval(target)
    
    # Only append non-empty results
    if result is not None and not result.is_empty():
        # Concatenate the current result with the accumulated results
        results = pl.concat([results, result], how="vertical")


print(results)

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
        data.id = sample.id
    """
    return duckdb.query(query).pl()

# Retrieve sample-specific data
results_list = [worker(results, target_sample_id) for target_sample_id in sample_id_of_kindex]

# Combine results from all samples
results = pl.concat(results_list).drop("id_1").unique()  # Remove duplicate columns

# Group by 'id' (variant ID) to count unique samples
sample_count_per_variant = (
    results.group_by("id")
    .agg(pl.col("sample").n_unique().alias("sample_count"))
)
print(sample_count_per_variant)
# Merge this sample count back to the results DataFrame
results_with_count = results.join(sample_count_per_variant, on="id")

# Determine custom IDs based on sample count
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
# Filter to keep only lines where sample matches sample_id_of_kindex
filtered_results = results_with_custom_id.filter(
    pl.col("sample").is_in(sample_id_of_kindex)
).unique()

#print(type(pl.col('chr')))
#print(type(pl.col('pos')))        
#results_with_count_columns = results.select(["pos", "alt", "ref", "chr", "id","sample"])
#print("Results ", results)
#print("Results count",results_with_count_columns)

selected_columns = filtered_results.select(["pos", "alt", "ref", "chr", "id","sample"])

# Print the selected columns
print("Results filt ", selected_columns)

# Iterate over each row in the Fraser DataFrame to update custom_ids based on the intervals
for index, row in fraser.iterrows():
    # Get the relevant chromosome, start, and end from the Fraser data
    chromosome = str(row['seqnames'])
    start = row['start']
    end = row['end']
    #print(type(chromosome))
    #print(type(start))

    # Filter the variants that fall within this interval
    matching_variants = filtered_results.filter(
        (pl.col('chr') == chromosome) &
        (pl.col('pos').cast(int) >= start-1000) &
        (pl.col('pos').cast(int) <= end+1000)
    )
    #print("matching",matching_variants.select(["pos", "alt", "ref", "chr", "id","sample"])) 
    # For each impact level, create the custom IDs
    if not matching_variants.is_empty():
        # High-impact variants
        high_impact_variants = matching_variants.filter(pl.col("snpeff_impact") == "HIGH")
        id_list_high = list(set(high_impact_variants["custom_id"].to_list()))

        # Moderate or modifier-impact variants
        moderate_impact_variants = matching_variants.filter(
            (pl.col("snpeff_impact") == "MODERATE") | (pl.col("snpeff_impact") == "MODIFIER")
        )
        id_list_moderate = list(set(moderate_impact_variants["custom_id"].to_list()))

        # Low-impact variants
        low_impact_variants = matching_variants.filter(pl.col("snpeff_impact") == "LOW")
        id_list_low = list(set(low_impact_variants["custom_id"].to_list()))

        # Update the Fraser DataFrame with the custom IDs, only if they exist
        fraser.at[index, 'custom_ids_high'] = id_list_high if id_list_high else None
        fraser.at[index, 'custom_ids_moderate'] = id_list_moderate if id_list_moderate else None
        fraser.at[index, 'custom_ids_low'] = id_list_low if id_list_low else None

    else:
        fraser.at[index, 'custom_ids_high'] = None
        fraser.at[index, 'custom_ids_moderate'] = None
        fraser.at[index, 'custom_ids_low'] = None

print(fraser)


# Write the updated DataFrame to a new file
#output_file = '24D0157_fraser_with_custom_ids.tsv'
# Write the modified Fraser DataFrame to the output file provided as a parameter
output_file = args.output
fraser.to_csv(output_file, sep='\t', index=False)




