import os
import sys

# Define a function to extract PERCENT_DUPLICATION from a file
def extract_duplication_percentage(file_path):
    with open(file_path, 'r') as file:
        #print(file_path)    
        sample_name = os.path.basename(file_path).split('_')[0]  # Extract sample name from filename
        for line in file:
            if line.startswith('rnaseq-capture'):
                data = line.strip().split('\t')
                percent_duplication = float(data[8])
                return sample_name, percent_duplication

    # Directory where the folders containing samples are located

# Check if the correct number of arguments is provided
if len(sys.argv) != 2:
    print("Usage: python script_name.py parent_folder_path")
    sys.exit(1)

parent_folder = sys.argv[1]


#parent_folder = '/data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq/20231012-RUN15-NextSeq_Mid_07RNASEQ/pipeline_v0/dup/'

    # Process each folder and extract duplication percentage
duplication_data = []
for folder_name in os.listdir(parent_folder):
    folder_path = os.path.join(parent_folder, folder_name)
    if os.path.isdir(folder_path):
            # List all files ending with 'dup.txt' in the current folder
        file_list = [f for f in os.listdir(folder_path) if f.endswith('_dup.txt')]
        for file_name in file_list:
            file_path = os.path.join(folder_path, file_name)
            result = extract_duplication_percentage(file_path)
            if result:
                duplication_data.append(result)

# Displaying the results
print("Sample Name\tPercent Duplication")
for sample, percent_duplication in duplication_data:
    print(f"{sample}\t{percent_duplication}")
