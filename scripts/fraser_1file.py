import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Split input file by sampleID and save to separate files")
parser.add_argument('-i', '--input_file', required=True, help="Path to the input file")
parser.add_argument('-b', '--basename',   required=True, help="Base name for the output files")
args = parser.parse_args()

df = pd.read_csv(args.input_file, sep='\t', index_col=False, header=0)

# Strip quotes from column names (old R format)
df.columns = df.columns.str.strip('"')

# Find sampleID column — handle variants
if 'sampleID' not in df.columns:
    candidates = [c for c in df.columns if 'sample' in c.lower()]
    if candidates:
        df = df.rename(columns={candidates[0]: 'sampleID'})
        print(f"[INFO] Renamed '{candidates[0]}' to 'sampleID'")
    else:
        print(f"[ERROR] No sampleID column found. Columns: {list(df.columns)}")
        raise KeyError("sampleID")

print(f"[INFO] {len(df)} rows, {df['sampleID'].nunique()} samples")

output_dir = os.path.dirname(args.basename)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

for sample_id, group in df.groupby('sampleID'):
    output_file = os.path.join(output_dir, f'{sample_id}.fraser.tab')
    group.to_csv(output_file, sep='\t', index=False)
    print(f'Created file: {output_file}')
