import os
import string
import sys


### Create config file for fraser ###

def generate_config_fraser(output_filename, pattern):
	
    # Define the root directory
    root_dir = '/data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq/'
    # Open a file to write the output
    with open(output_filename, 'w') as output_file:
        output_file.write("sampleID\tbamFile\tgroup\tgene\tpairedEnd\n")
    
        # Initialize group increment
        group_increment = 1
     
        # Initialize gene increment
        #gene_sequence = iter(string.ascii_uppercase + string.digits)
        gene_sequence = (
           f"{letter}{number}"
           for number in range(1, sys.maxsize)  # Cycle through numbers indefinitely
           for letter in string.ascii_uppercase
        )
    
        # Walk through each directory in the root directory
        for dirpath, dirnames, filenames in os.walk(root_dir):
            if 'star' in dirpath:
                for filename in filenames:
                        # Check if the filename contains 'Aligned.sortedByCoord.out.bam'
                    if filename.endswith('Aligned.sortedByCoord.out.bam') and pattern in filename:
                            # Extract sampleID and other information
                            sample_id = filename.split('-')[0]
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
                            
                       #  gene = next(gene_sequence)
		       #  # If the sequence reaches the end, reset to the start with an incremented number
                       #     if gene == 'Z':
                       #         gene_sequence = iter([f'{letter}1' for letter in string.ascii_uppercase])
        
                            paired_end = 'TRUE' 
                            # Write the information to the output file
                            output_file.write(
                                f"{sample_id}\t{os.path.join(dirpath, filename)}\t{group}\t{gene}\t{paired_end}\n"
                        )
        
                            group_increment += 1 

if __name__ == "__main__":
    if len(sys.argv) != 3:
         print("Usage: python script_name.py output_filename pattern")
    else:
        output_filename = sys.argv[1]
        pattern = sys.argv[2]
        generate_config_fraser(output_filename,pattern)
