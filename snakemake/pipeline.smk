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


def get_version_from_env(env_yml, cmd):
    """
    Detect a tool version by searching for its binary inside the local conda_env folder.
    """
    import subprocess
    from pathlib import Path

    print(f"[DEBUG] get_version_from_env called with: {env_yml}, {cmd}")
    bin_name = cmd.split()[0]

    # 🔧 Ensure PIPELINE_DIR is a Path, even if it's a string globally
    base_dir = Path(PIPELINE_DIR) / "conda_env"
    print(f"[DEBUG] Searching envs in {base_dir}")

    found_binary = None

    # Recursively look for binary under all conda_env subdirectories
    for env_dir in sorted(base_dir.glob("*")):
        if not env_dir.is_dir():
            continue
        for binary in env_dir.glob(f"**/bin/{bin_name}"):
            if binary.exists():
                print(f"[DEBUG] ✅ Found binary {binary}")
                found_binary = binary
                break
        if found_binary:
            break

    if not found_binary:
        print(f"[WARN] ❌ Could not find {bin_name} in any conda_env directory")
        return f"{bin_name}: not found"

    # Run the version command, capturing both stdout and stderr
    try:
        result = subprocess.run(
            f"{found_binary} {cmd[len(bin_name):]}",
            shell=True,
            check=True,
            capture_output=True,
            text=True,
            executable="/bin/bash"
        )
        output = (result.stdout + result.stderr).strip()
        version = output.split("\n")[0]
        print(f"[DEBUG] ✅ {bin_name} version detected: {version}")
        return version if version else f"{bin_name}: version not found"
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to run {cmd} in {found_binary}: {e}")
        return f"/bin/sh: 1: {bin_name}: not found"


# Set environment variables
#os.environ["CONDARC"] = OUTPUT_REP + ".condarc"

include: '../rules/01_trim_fastqc.smk'
include: '../rules/02_alignement.smk'
include: '../rules/03_comptage.smk'
include: '../rules/04_outrider_fraser.smk'
include: '../rules/05_metrics.smk'
include: '../rules/06_versions.smk'

workdir: OUTPUT_REP


SAMPLES_ID = [s[:7] for s in SAMPLES]

def active_samples(blacklist_file, samples=SAMPLES_ID):
    """
    Return list of sample IDs excluding any in the blacklist file.
    If the blacklist file does not exist, return all samples.
    """
    if os.path.exists(blacklist_file):
        with open(blacklist_file) as f:
            excluded = set(f.read().split())
    else:
        excluded = set()
    return [s for s in samples if s not in excluded]
ACTIVE_FRASER = active_samples(fraser_blacklist)
ACTIVE_OUTRIDER = active_samples(outrider_blacklist)
ACTIVE_PCA = active_samples(PCA_blacklist)

rule all:
    input:
        expand(FASTQ_DIR+"/{sample}_R1.fastq.gz",sample=SAMPLES),
        expand(FASTQ_DIR+"/{sample}_R2.fastq.gz",sample=SAMPLES),
        expand(rules.fastqc_report.output, sample=SAMPLES),
        expand(rules.fastp.output.html, sample=SAMPLES),
        expand(rules.fastp.output.R1, sample=SAMPLES),
        expand(rules.fastp.output.R2, sample=SAMPLES),
        expand(rules.fastqc_trim_report.output, sample=SAMPLES),
        expand(rules.alignment_star.output.bam, sample=SAMPLES),
        expand(rules.index_bam.output.bai, sample=SAMPLES),
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
        rules.fraser.output.fraser,
        expand(rules.fraser.output.fraser),
        expand(rules.fraser_boxplot.output.filt, samples_id=ACTIVE_FRASER),
        expand(rules.volcano.output, samples_id=ACTIVE_PCA),
        expand(rules.boxplot.output.filt, samples_id=ACTIVE_PCA),
        rules.generate_and_run_param_notebook.output.executed_nb,
        rules.generate_metrics.output.metrics,
        "benchmarks/versions/pipeline_versions.tsv"



addresses = ["laura.dosoutoferreira@chu-nantes.fr"]
#t = "ERREUR ROUTINE !"
#
#onerror:
#    for mail in addresses :
#        shell('mail -s "an error occurred" {mail} ')
#
#onsuccess:
#    for mail in addresses:
#        shell(f'mail -s "Pipeline completed successfully" {mail} <<< "The Snakemake pipeline finished without errors."')

addresses = ["laura.dosoutoferreira@chu-nantes.fr"]

import smtplib
from email.mime.text import MIMEText
from datetime import datetime

# -----------------------------
# CONFIGURATION
# -----------------------------
SMTP_HOST =   "172.27.162.183"
SMTP_PORT = 25  

FROM = "laura.dosoutoferreira@chu-nantes.fr"  
TO = ["laura.dosoutoferreira@chu-nantes.fr"]   

# -----------------------------
# EMAIL FUNCTION
# -----------------------------
def send_email(subject, body):
    try:
        msg = MIMEText(body)
        msg["Subject"] = subject
        msg["From"] = FROM
        msg["To"] = ", ".join(TO)

        with smtplib.SMTP(SMTP_HOST, SMTP_PORT, timeout=10) as server:
            server.send_message(msg)

        print(f"[INFO] Email sent to {TO} from {FROM}")
    except Exception as e:
        print(f"[WARN] Could not send email: {e}")

# -----------------------------
# envoie mail
# -----------------------------
onsuccess:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    send_email(
        subject="✅ Snakemake pipeline completed successfully",
        body=f"The pipeline finished successfully at {now}.\nAll rules completed without errors."
    )

onerror:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    send_email(
        subject="❌ Snakemake pipeline failed",
        body=f"An error occurred during the pipeline at {now}.\nCheck the logs for details."
    )

