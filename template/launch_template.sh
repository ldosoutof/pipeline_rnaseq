#!/bin/bash

###Script launch snakemake creating a rule graph and a file graph for analysis to see input output and path of directory

export CONDA_PKGS_DIRS=/data1/rnaseq/try/conda_pkgs
export CONDA_ENVS_DIRS=/data1/rnaseq/try/conda_envs


source /app/miniconda/2024.10-1/bin/activate 
conda acti vate /app/envs_laura/snakemake

snakefile="ROOT_PIPELINE/snakemake/pipeline.smk"
config_file="PATH_TO_CONFIG"
path_to_launch=`echo ${config_file} |xargs dirname`


echo "Creating filegraph..."
snakemake -s $snakefile --configfile $config_file --filegraph | dot -Tpng > ${path_to_launch}"/filegraph.png"

echo "Creating rulegraph..."
snakemake -s $snakefile --configfile $config_file --rulegraph | dot -Tpng > ${path_to_launch}"/rulegraph.png"

echo "Launching pipeline"
snakemake \
            -s  $snakefile \
	                    --configfile $config_file \
			                        --latency-wait 60 \
						                        -j 100 \
									                        --resources mem_gb=60 \
												                            --restart-times 2 \
															                                --rerun-incomplete \
																			                            --use-conda \
																						                                    --jobname {rulename}_{wildcards}.{jobid}.sh


