shell.executable('bash')
from os.path import join
import os
import re
import subprocess
from time import strftime, localtime

#configfile: "config.yml"

##--------------------------------------------------------------------------------------##
## Auteur : Laura DO SOUTO FERREIRA (dosoutoferreira.laura@chu-nantes.fr)
## Affiliation : CHU Nantes
## But : Fichier Snakemake pour le pipeline de rnaseq
## 
## Latest modification : 27/01/2019
## Latest modification : 11/10/2023  
## 
##--------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------##
## Declaration des constantes
##--------------------------------------------------------------------------------------##

RUNS_DIR = config['runs_dir']
FASTQ_DIR = config['configuration']['fastq_dir']
OUTPUT_REP = config['configuration']['outputDir']
PIPELINE_DIR = config['configuration']['pipeline_dir']
SAMPLES = config['configuration']["samples"] 
GENOME = config["genome"]
DI_BED = config["di_bed"]
PANELAPP = config["panelapp"]
HPO = config["hpo"]
PHENO = config["pheno"]
PLI = config["pli"]
BED = config["bed"]
fraser_count = config["fraser_count"]
PCA_blacklist = config["pca_blacklist"]
outrider_blacklist = config["outrider_blacklist"]
fraser_blacklist = config["fraser_blacklist"]
RSEQ_BED = config["rseq_bed"]
STAR_GENOME = config["star_genome"]
RSEM_GENOME = config["rsem_genome"]
KALLISTO_IDX = config["kallisto_idx"]
PADDED = config["padded"]
GTF = config["gtfFile"]
MATRICES = config["matrices"]
TPM = config["matrice_tpm"]
run_annot_sake = config["run_annot_sake"]

# Set environment variables
#os.environ["CONDARC"] = OUTPUT_REP + ".condarc"

include: '../rules/01_trim_fastqc.smk'
include: '../rules/02_alignement.smk'
include: '../rules/03_comptage.smk'
include: '../rules/04_outrider_fraser.smk'
include: '../rules/05_metrics.smk'

workdir: OUTPUT_REP


SAMPLES_ID = [s[:7] for s in SAMPLES]

rule all:
    input:
        expand(FASTQ_DIR+"/{sample}_R1.fastq.gz",sample=SAMPLES),
        expand(FASTQ_DIR+"/{sample}_R2.fastq.gz",sample=SAMPLES),
        expand(rules.fastqc_report.output, sample=SAMPLES),
        expand(rules.fastp.output.html, sample=SAMPLES),
        expand(rules.fastp.output.R1, sample=SAMPLES),
        expand(rules.fastp.output.R2, sample=SAMPLES),
        expand(rules.fastqc_trim_report.output, sample=SAMPLES),
        expand(rules.alignment_star.output.bamg, sample=SAMPLES),
        expand(rules.index_bam.output.baig, sample=SAMPLES),
        expand(rules.htseq_gene.output.gene, sample=SAMPLES),
        expand(rules.matrix.output, sample=SAMPLES),
        expand(rules.matrix_tpm.output.gene, sample=SAMPLES),
        expand(rules.kallistoBed.output.h5, sample=SAMPLES),
        expand(rules.kallisto2gene.output, sample=SAMPLES),
        expand(rules.bam_stats.output.on_target, sample=SAMPLES),
        expand(rules.rseqc.output, sample=SAMPLES),
        expand(rules.multiqc.output, sample=SAMPLES), 
        expand(rules.outrider.output.out_file, sample=SAMPLES),
        expand(rules.fraser_config.output, sample=SAMPLES),
        expand(rules.fraser.output.fraser, sample=SAMPLES),
        expand(rules.fraser_annot_rare.output.fraser_rare, samples_id=SAMPLES_ID),
        expand(rules.outrider_annot_rare.output.outrider_rare, samples_id=SAMPLES_ID),
        expand(rules.volcano.output, samples_id=SAMPLES_ID),
        expand(rules.boxplot.output.filt, samples_id=SAMPLES_ID),
        expand(rules.fraser_boxplot.output.filt, samples_id=SAMPLES_ID),
        rules.generate_and_run_param_notebook.output.executed_nb,
        rules.generate_metrics.output.metrics



#addresses = ["laura.dosoutoferreira@chu-nantes.fr"]
#t = "ERREUR ROUTINE !"

#onerror:
 #   for mail in addresses :
  #      shell('mail -s "an error occurred" {mail} ')
