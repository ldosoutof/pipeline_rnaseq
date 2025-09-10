import os
import pandas as pd
import re
from statistics import mean, stdev
import subprocess


# Function to extract sample names with "puromoins"
def aggregate_matrice(root_directory):
    matrix = {'SampleID': [], 'Percent_On_Target': [], 'SamplePath': []}
    data = []
    for root, dirs, files in os.walk(root_directory):
        if 'pipeline_v0' in dirs:  # Checking for 'pipeline_v0' folders
            pipeline_path = os.path.join(root, 'pipeline_v0')
            #print(pipeline_path)    
            for subroot, subdirs, subfiles in os.walk(pipeline_path):
                if 'coverage' in subdirs:  # Checking for 'kallisto' folders inside 'pipeline_v0'
                    coverage_path = os.path.join(subroot, 'coverage')
                    #print(kallisto_path)
                    for subroot, subdirs, subfiles in os.walk(coverage_path):
                        #print(subdirs)
                        for dir_sample in subdirs:
                            #print(dir_sample)
                            if not 'ONCO' in subroot:    
                                #print('puromoins')
                                sample_path =  os.path.join(subroot,dir_sample)
                                sample_name = sample_path.split('/')[-1]
                                #print(sample_path)    
                                on_target_file = os.path.join(sample_path, f"{sample_name}_on_target.txt")
                                with open(on_target_file, 'r') as on_target:
                                   on_target_count = int(on_target.read().strip())
                                   #print(on_target_count)				   
                                stats_file = os.path.join(sample_path, f"{sample_name}_stats.txt")
                                with open(stats_file, 'r') as stats:
                                    for line in stats:
                                        #print(line)					    
                                        if line.startswith('SN\tsequences:'):
                                            #print(line)
                                            sequences_count = int(line.split('\t')[2])
                                            percent_on_target = (on_target_count / sequences_count) * 100
                                            matrix['SampleID'].append(sample_name)
                                            matrix['Percent_On_Target'].append(percent_on_target)
                                            matrix['SamplePath'].append(sample_path)
                                            data.append({'SampleID': sample_name, 'Percent_On_Target': percent_on_target, 'SamplePath': sample_path})					   
                                            break

                            # Calculate percent of on-target
                                percent_on_target = (on_target_count / sequences_count) * 100

                            # Add to matrix
                                #matrix[sample_name] = percent_on_target

    #df = pd.DataFrame(list(matrix.items()), columns=['SampleID', 'Percent_On_Target'])
    df = pd.DataFrame(data)
    return df

# Example usage:
root_directory = '/data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq/'

result_df = aggregate_matrice(root_directory)
print(result_df)


output_csv_path = 'percent_on_target.csv'
result_df.to_csv(output_csv_path, index=False)

print(f"Matrix saved to {output_csv_path}")
