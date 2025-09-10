import os
import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Split input file by sampleID and save to separate files")
parser.add_argument('-i', '--input_file', required=True, help="Path to the input file")
parser.add_argument('-b', '--basename', required=True, help="Base name for the output files")
args = parser.parse_args()

# Read the input file
df = pd.read_csv(args.input_file, sep='\t', index_col=False, header=0)
print(df)


# Ensure the folder from the basename exists
output_dir = os.path.dirname(args.basename)  # Extract the directory part from the basename
if output_dir and not os.path.exists(output_dir):  # Check if directory exists, and create if it doesn't
    os.makedirs(output_dir, exist_ok=True)


# Group by sampleID and write to separate files
for sample_id, group in df.groupby('sampleID'):
    # Create a filename based on the sampleID and the provided basename
    output_file = os.path.join(output_dir, f'{output_dir}/{sample_id}.fraser.tab')

    # Write the group to the file
    group.to_csv(output_file, sep='\t', index=False)

    print(f'Created file: {output_file}')
