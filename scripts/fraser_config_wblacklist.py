import os
import string
import sys
import pandas as pd

def generate_config_fraser(output_filename, pattern, blacklist_file=None):
    # Define the root directory
    root_dir = '/datawork2/genetique/RNASeq/diag/prod/'
    
    # Initialize group increment
    group_increment = 1

    # Initialize gene increment
    gene_sequence = (
        f"{letter}{number}"
        for number in range(1, sys.maxsize)  # Cycle through numbers indefinitely
        for letter in string.ascii_uppercase
    )

    # Initialize an empty list to collect rows
    rows = []

    # Load blacklist if provided
    blacklist = set()
    if blacklist_file:
        with open(blacklist_file, 'r') as bf:
            blacklist = set(line.strip() for line in bf)
    
    # Walk through each directory in the root directory
    for dirpath, dirnames, filenames in os.walk(root_dir):
        if 'star' in dirpath:
            for filename in filenames:
                # Check if the filename contains 'Aligned.sortedByCoord.out.bam'
                if filename.endswith('Aligned.sortedByCoord.out.bam') and pattern in filename:
                    # Extract sampleID and other information
                    sample_id = filename.split('-')[0]
                    
                    # Check if the sample_id is in the blacklist
                    if sample_id in blacklist or 'RUN17' in dirpath:
                        continue
                    
                    group = group_increment
                    try:
                        gene = next(gene_sequence)
                    except StopIteration:
                        # If the sequence reaches the end, recreate the sequence
                        gene_sequence = (
                            f"{letter}{number}"
                            for number in range(1, sys.maxsize)
                            for letter in string.ascii_uppercase
                        )
                        gene = next(gene_sequence)
                    
                    paired_end = 'TRUE'
                    
                    # Append the information to the rows list
                    rows.append([sample_id, os.path.join(dirpath, filename), group, gene, paired_end])
                    
                    group_increment += 1

    # Convert the list to a DataFrame
    df = pd.DataFrame(rows, columns=["sampleID", "bamFile", "group", "gene", "pairedEnd"])
    
    # Write the DataFrame to the output file
    df.to_csv(output_filename, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python script_name.py output_filename pattern [blacklist_file]")
    else:
        output_filename = sys.argv[1]
        pattern = sys.argv[2]
        blacklist_file = sys.argv[3] if len(sys.argv) == 4 else None
        generate_config_fraser(output_filename, pattern, blacklist_file)

