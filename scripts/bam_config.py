import os
import string
import sys


### Create file with info of path sample and condition and run and technique (XTHS1 ou XTHS2) ###

def generate_config(output_filename):

    # Define the root directory
    root_dir = '/data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq/'
    # Open a file to write the output
    with open(output_filename, 'w') as output_file:
        output_file.write("path\tsample\tcondition\trun\ttechnique\n")
    
        # Initialize group increment
        group_increment = 1
     
        # Initialize gene increment
        gene_sequence = iter(string.ascii_uppercase + string.digits)
    
        # Walk through each directory in the root directory
        for root, dirnames, filenames in os.walk(root_dir):
            if 'pipeline_v0' in dirnames:  # Checking for 'pipeline_v0' folders
                pipeline_path = os.path.join(root, 'pipeline_v0')
                #print(pipeline_path)
                for subroot, subdirs, subfiles in os.walk(pipeline_path):
                    #print(subdirs)		    
                    if 'kallisto_bed' in subdirs:
                        kallisto_path = os.path.join(subroot, 'kallisto_bed')
                        for subroot, subdirs, subfiles in os.walk(kallisto_path):
                            if not 'ONCO' in subroot:
                                for dir_sample in subdirs:
                                    if 'PAX' in dir_sample or 'PAX'.lower() in dir_sample or 'PX' in dir_sample:
                                        sample_path =  os.path.join(subroot,dir_sample)
                                        print(sample_path)
                                        abund_file = os.path.join(sample_path, 'abundance.tsv')
                                        if os.path.exists(abund_file):
                                            run_id = subroot.split('/')[-3]
                                            # Extract sampleID and other information
                                            sample = sample_path.split('/')[-1]
                                            sample_id = sample.split('-')[0]                                            
                                            # Write the information to the output file
                                            output_file.write(
                                                f"{abund_file}\t{sample_id}\tPAXGENE\t{run_id}\tXTHS1\n"
                                            )
                        
                                    if 'PUROMOINS' in dir_sample or 'PUROMOINS'.lower() in dir_sample:
                                        sample_path =  os.path.join(subroot,dir_sample)
                                        print(sample_path)
                                        abund_file = os.path.join(sample_path, 'abundance.tsv')
                                        if os.path.exists(abund_file):
                                            run_id = subroot.split('/')[-3] 
                                            # Extract sampleID and other information
                                            sample = sample_path.split('/')[-1]       
                                            sample_id = sample.split('-')[0]					    
                                            # Write the information to the output file
                                            output_file.write(
                                                f"{abund_file}\t{sample_id}\tPUROMOINS\t{run_id}\tXTHS1\n"
                                            )
         
                                    if 'PUROPLUS' in dir_sample or 'PUROPLUS'.lower() in dir_sample:
                                        sample_path =  os.path.join(subroot,dir_sample)
                                        print(sample_path)
                                        abund_file = os.path.join(sample_path, 'abundance.tsv')
                                        if os.path.exists(abund_file):
                                            run_id = subroot.split('/')[-3]
                                            # Extract sampleID and other information
                                            sample = sample_path.split('/')[-1]
                                            sample_id = sample.split('-')[0]					    
         
                                            # Write the information to the output file
                                            output_file.write(
                                                f"{abund_file}\t{sample_id}\tPUROPLUS\t{run_id}\tXTHS1\n"
                                        )




if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py output_filename")
    else:
        output_filename = sys.argv[1]
        generate_config(output_filename)

