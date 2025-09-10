#!/bin/bash


# Function to display usage message
usage() {
    echo "Usage: $0 -r <run_name> -i <input_dir>"
        exit 1
}

# Parse command-line options
while getopts "r:i:" opt; do
    case ${opt} in
	    r )
                run_name=$OPTARG
		;;
		i) 
                input_dir=$OPTARG
		;;
		*) 
                usage
		;;
        esac
done

# Check if both parameters are set
if [ -z "${run_name}" ] || [ -z "${input_dir}" ]; then
    usage
fi

#run_name="20231122_RUN17_NextSeq_Mid_8RNASEQ"
#input_dir="/data2/Server_S_nextseq/231122_NB551465_0295_AH3W2MAFX5/Alignment_1/20231123_034200/Fastq"


# Define the output directory
output_dir="/data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq"

mkdir ${output_dir}/${run_name}
mkdir ${output_dir}/${run_name}/fastq

# Loop through each R1 FASTQ file in the specified directory
for i in ${input_dir}/2*_L001_R1_001.fastq.gz; do
    # Extract the base filename from the path
    filename=$(basename "${i}")
    
    # Extract the 'end' part after the first underscore
    end="${filename#*_}"
    
    # Extract the 'index' part before the first underscore in 'end'
    index="${end%%_*}"
    
    # Extract the 'id' part before the first underscore in 'filename'
    id="${filename%%_*}"
    
    # Extract the directory path
    dir="${i%/*}"
    
    
    # Construct and launch the zcat commands for R1 and R2 files in the background
    zcat ${input_dir}/${id}_${index}_L001_R1_001.fastq.gz ${input_dir}/${id}_${index}_L002_R1_001.fastq.gz ${input_dir}/${id}_${index}_L003_R1_001.fastq.gz ${input_dir}/${id}_${index}_L004_R1_001.fastq.gz | gzip -c > ${output_dir}/${run_name}/fastq/${id}_R1.fastq.gz &
    zcat ${input_dir}/${id}_${index}_L001_R2_001.fastq.gz ${input_dir}/${id}_${index}_L002_R2_001.fastq.gz ${input_dir}/${id}_${index}_L003_R2_001.fastq.gz ${input_dir}/${id}_${index}_L004_R2_001.fastq.gz | gzip -c > ${output_dir}/${run_name}/fastq/${id}_R2.fastq.gz &
    done

    # Wait for all background jobs to finish
    wait
