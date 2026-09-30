import os
import re
import pandas as pd
import argparse

def extract_short_id(s):
    """
    Extrait l'ID court depuis n'importe quel format connu :
      26D0643                              -> 26D0643
      26D0643-MOINS                        -> 26D0643
      25D2693-STEMC-PUROMOINS-AVITI        -> 25D2693
      X26D0643                             -> 26D0643  (ancien format R avec préfixe X)
      25D2693.STEMC.PUROMOINS.AVITI        -> 25D2693  (format R avec points)
    """
    s = str(s).lstrip('X')
    return re.split(r'[.\-]', s)[0]

parser = argparse.ArgumentParser(description="Split OUTRIDER output by sampleID into per-sample files")
parser.add_argument('-i', '--input_file', required=True, help="Path to the input file")
parser.add_argument('-b', '--basename',   required=True, help="Base name / output directory")
args = parser.parse_args()

df = pd.read_csv(args.input_file, sep='\t', index_col=False, header=0)

# Strip quotes from column names (old R write.table format)
df.columns = df.columns.str.strip('"')

# Normalise sampleID column name
if 'sampleID' not in df.columns:
    candidates = [c for c in df.columns if 'sample' in c.lower()]
    if candidates:
        df = df.rename(columns={candidates[0]: 'sampleID'})
        print(f"[INFO] Renamed '{candidates[0]}' to 'sampleID'")
    else:
        raise KeyError(f"No sampleID column found. Columns: {list(df.columns)}")

print(f"[INFO] {len(df)} rows, {df['sampleID'].nunique()} unique sampleIDs")

output_dir = os.path.dirname(args.basename)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

for sample_id, group in df.groupby('sampleID'):
    short_id = extract_short_id(sample_id)
    output_file = os.path.join(output_dir, f'{short_id}.outrider.tab')
    group.to_csv(output_file, sep='\t', index=False)
    print(f'[OK] {sample_id} -> {output_file}')
