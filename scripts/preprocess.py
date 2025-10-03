#!/usr/bin/python3

"""
Script creating symbolik link for data and renaming them / create also from template a launching and config file for each analysis
"""
import argparse
import os
from os import makedirs
from pathlib import Path
import glob
import sys
import re
import fnmatch
import codecs
import csv
import yaml

"""
Options required to launch script 
- path : where you launch the script but mostly where are going to go your launcher directory with launch & config file 
- info : tabulated file with run and condition of runs
- workdir : where you want the analysis and results to go
- dataDir : directory where are the runs 
"""

parser = argparse.ArgumentParser(prog="preprocess.py", description='')
parser.add_argument('-p', '--path', required=True)
parser.add_argument('-w', '--workDir', required=True)
parser.add_argument('-d', '--dataDir', required=True)

args = parser.parse_args()

path_pipeline = Path(args.path).absolute()
path_data = Path(args.dataDir).absolute()
path_workdir = Path(args.workDir).absolute()
workdir = os.path.basename(path_workdir)


folder=path_data
subFolders = [ f.path for f in os.scandir(folder) if f.is_dir() ]

"""
    -  
    - 
"""
makedirs(path_pipeline / workdir , exist_ok=True)

print(path_pipeline)
with codecs.open(path_pipeline / "template/launch_template.sh",  encoding='utf-8', errors='ignore') as fp:
    lines = "".join(fp.readlines())

launch_folder = path_pipeline / workdir/ 'launch_folder'
#launch_folder = path_pipeline / (workdir / 'launch_folder')

makedirs(launch_folder, exist_ok=True)
path_to_config = Path(launch_folder) / "config.yml"
lines = lines.replace("PATH_TO_CONFIG", str(path_to_config.absolute()))
lines = lines.replace("ROOT_PIPELINE", str(path_pipeline.absolute()))

with open(launch_folder / "launch.sh", "w+") as fp:
   fp.write(lines)


#filename = "/data2/Exome_analysis_BC/Server_S_exome/MAGNIS/Runs_nextseq/20230602_RUN11_NextSeq_MidOutput_8RNASEQ/FDR-NS-20230601-RNA12-PAT.csv"

print(os.listdir(os.path.join(path_data,'fastq')))

sample_ids =[]

for fastq_files in os.listdir(os.path.join(path_data,'fastq')):
	if fastq_files.endswith('R1.fastq.gz'):
	    if not fastq_files.startswith('Undetermined'):
                #print(fastq_files.rsplit('_',1)[0])
                sample_ids.append(fastq_files.rsplit('_',1)[0])




"""
Retrieve path for config file and other info needed to complete it
Load yaml (config file template)
Write new config file with info into launch folder for analysis
"""

new_config = launch_folder / "config.yml"
config = yaml.safe_load(open(path_pipeline / "template/config_template.yml"))

config["configuration"]["samples"] = sample_ids
config["configuration"]["outputDir"] = str(path_workdir)
#config["configuration"]["inputDir"] = "/tmp/Data"
config["configuration"]["fastq_dir"] = str(path_data/ 'fastq')
config["configuration"]["pipeline_dir"] = str(path_pipeline)

with open(new_config, "w+") as fp:
    yaml.dump(config, fp)


