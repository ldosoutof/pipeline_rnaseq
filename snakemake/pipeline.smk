shell.executable('bash')
from os.path import join
import os
import re
import subprocess
import sys
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

RUNS_DIR  = config['runs_dir'].rstrip("/")
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
fraser_count       = config["fraser_count"]
fraser_count_hyper = config.get("fraser_count_hyper",
                                config["fraser_count"].rstrip("/") + "_hyper/")
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

# --- featureCounts / make_zip globals ----------------------------------------
# PROD_ROOT  : root directory of all run folders (same as runs_dir)
PROD_ROOT       = config["runs_dir"].rstrip("/")
GTF_REFSEQ      = config.get("gtf_refseq", "")
GTF_ENSEMBL     = config.get("gtfFile", GTF)     # reuse Ensembl GTF already loaded
GENEID_MAP      = config.get("geneid_map", PROD_ROOT + "/geneid_to_ensg.tsv")
GNOMAD          = config.get("gnomad", "")
MENDELIOME      = config.get("mendeliome", "")

# RUN_SAMPLE is built after SAMPLES_ID is defined — see below
import glob as _glob

def _build_run_sample(prod_root, samples_id):
    """
    Scan prod_root for aligned BAMs and return
    [(run_tag, sample_id, bam_path), ...] for each active sample.
    """
    triples = []
    pattern = os.path.join(prod_root, "20*", "pipeline_v0", "star",
                           "*_Aligned.sortedByCoord.out.bam")
    for bam in sorted(_glob.glob(pattern)):
        parts = bam.split(os.sep)
        run_tag = next((p for p in parts if p.startswith("20")), None)
        if run_tag is None:
            continue
        basename = os.path.basename(bam)
        sample_id = basename.split("-")[0] if "-" in basename else basename.split("_")[0]
        if sample_id in samples_id:
            triples.append((run_tag, sample_id, bam))
    return triples



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

SAMPLES_ID = [s[:7] for s in SAMPLES]

# FRASER and OUTRIDER are only meaningful for blood RNA-seq
FRASER_KEYWORDS = tuple(
    kw.upper() for kw in
    config.get("fraser_keywords", "MOINS PUROMINS").split()
)

# Try to detect blood samples from full sample names first.
# If SAMPLES only contains short IDs (e.g. "26D0198"), fall back to
# scanning FASTQ filenames to recover the full name with the suffix.
def _get_full_sample_names(samples, fastq_dir):
    """
    If sample names already contain a keyword, use them as-is.
    Otherwise scan fastq_dir for R1 files to recover full names.
    """
    import glob as _g
    full = []
    for s in samples:
        if any(kw in s.upper() for kw in FRASER_KEYWORDS):
            full.append(s)
        else:
            # look for <s>*_R1.fastq.gz in fastq_dir
            pattern = os.path.join(fastq_dir, f"{s}*_R1.fastq.gz")
            hits = _g.glob(pattern)
            if hits:
                # derive full name from filename: strip _R1.fastq.gz
                fname = os.path.basename(hits[0])
                full_name = fname.replace("_R1.fastq.gz", "")
                full.append(full_name)
            else:
                full.append(s)  # keep as-is if not found
    return full

SAMPLES_FULL     = _get_full_sample_names(SAMPLES, FASTQ_DIR)
SAMPLES_BLOOD    = [s for s in SAMPLES_FULL if any(kw in s.upper() for kw in FRASER_KEYWORDS)]
SAMPLES_ID_BLOOD = [s[:7] for s in SAMPLES_BLOOD]

print(f"[INFO] FRASER_KEYWORDS: {FRASER_KEYWORDS}", file=sys.stderr)
print(f"[INFO] SAMPLES ({len(SAMPLES)}): {SAMPLES[:3]}...", file=sys.stderr)
print(f"[INFO] SAMPLES_BLOOD ({len(SAMPLES_BLOOD)}): {SAMPLES_BLOOD[:3]}...", file=sys.stderr)
print(f"[INFO] SAMPLES_ID_BLOOD ({len(SAMPLES_ID_BLOOD)}): {SAMPLES_ID_BLOOD[:3]}...", file=sys.stderr)

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

ACTIVE_FRASER   = active_samples(fraser_blacklist,   samples=SAMPLES_ID_BLOOD)
ACTIVE_OUTRIDER = active_samples(outrider_blacklist,  samples=SAMPLES_ID_BLOOD)
ACTIVE_PCA      = active_samples(PCA_blacklist)

include: '../rules/logging.smk'
include: '../rules/01_trim_fastqc.smk'
include: '../rules/02_alignement.smk'
include: '../rules/03_comptage.smk'
include: '../rules/03bis_featureCount.smk'
include: '../rules/04_outrider_fraser.smk'
include: '../rules/04_outrider_fraser_hyper.smk'
include: '../rules/04_make_zip.smk'
include: '../rules/05_metrics.smk'
include: '../rules/06_versions.smk'
include: '../rules/07_variantCalling_bis.smk'

workdir: OUTPUT_REP

# RUN_SAMPLE: built here, after SAMPLES_ID is available
RUN_SAMPLE = _build_run_sample(PROD_ROOT, set(SAMPLES_ID))

rule all:
    input:
        expand(FASTQ_DIR+"/{sample}_R1.fastq.gz",sample=SAMPLES),
        expand(FASTQ_DIR+"/{sample}_R2.fastq.gz",sample=SAMPLES),
        #expand(rules.fastqc_report.output, sample=SAMPLES),
        #expand(rules.fastp.output.html, sample=SAMPLES),
        #expand(rules.fastp.output.R1, sample=SAMPLES),
        #expand(rules.fastp.output.R2, sample=SAMPLES),
        #expand(rules.fastqc_trim_report.output, sample=SAMPLES),
        #expand(rules.alignment_star.output.bam, sample=SAMPLES),
        #expand(rules.index_bam.output.bai, sample=SAMPLES),
        #expand(rules.htseq_gene.output.gene, sample=SAMPLES),
        #expand(rules.matrix.output, sample=SAMPLES),
        #expand(rules.matrix_tpm.output.gene, sample=SAMPLES),
        #expand(rules.kallistoBed.output.h5, sample=SAMPLES),
        #expand(rules.kallisto2gene.output, sample=SAMPLES),
        #expand(rules.bam_stats.output.on_target, sample=SAMPLES),
        #expand(rules.rseqc.output, sample=SAMPLES),
        #expand(rules.multiqc.output, sample=SAMPLES), 
        expand(rules.outrider.output.out_file, sample=SAMPLES),
        rules.annotation_outrider.output.annot,
        expand(rules.fraser_config.output, sample=SAMPLES),
        rules.fraser.output.fraser,
        expand(rules.fraser.output.fraser),
        expand(rules.fraser_boxplot.output.filt, samples_id=ACTIVE_FRASER),
        expand(rules.volcano.output, samples_id=ACTIVE_OUTRIDER),
        expand(rules.boxplot.output.filt, samples_id=ACTIVE_OUTRIDER),
        # featureCounts (03bis) — one output per (run, sample) pair
        [rules.featurecounts_gene.output.gene.format(run=r, sample=s)
         for r, s, _ in RUN_SAMPLE],
        [rules.map_refseq_to_ensembl.output.ensembl_counts.format(run=r, sample=s)
         for r, s, _ in RUN_SAMPLE],
        rules.matrix_featurecounts.output.matrix,
        # hyper pipeline (outrider + fraser on featureCounts matrix)
        rules.annotation_outrider_hyper.output.annot,
        rules.annotation_fraser_hyper.output.annot_fraser,
        rules.rnaseq_per_sample_hyper.output.zip_file,
        # ZIP bundling + per-sample analysis (04_make_zip)
        [rules.make_analysis_zip.output.zip.format(run=r)
         for r in sorted({r for r, s, _ in RUN_SAMPLE})],
        [rules.run_rnaseq_analysis.output.result_zip.format(run=r)
         for r in sorted({r for r, s, _ in RUN_SAMPLE})],
        rules.generate_and_run_param_notebook.output.executed_nb,
        rules.generate_metrics.output.metrics,
        "benchmarks/versions/pipeline_versions.tsv",
        rules.mean_chrY_expression.output.tsv,
        rules.vaf_violin_plot_run_females.output.plot


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
    failures = collect_failed_logs(log_dir="log")
    if failures:
        details = []
        for f in failures:
            details.append(
                f"Rule      : {f['rule']}\n"
                f"Exit code : {f['exit_code']}\n"
                f"Timestamp : {f['timestamp']}\n"
                f"Log file  : {f['log_path']}\n"
                f"--- last 30 lines ---\n{f['tail']}\n"
            )
        body = (
            f"An error occurred during the pipeline at {now}.\n\n"
            + "\n" + "="*60 + "\n"
            + ("\n" + "="*60 + "\n").join(details)
        )
    else:
        body = (
            f"An error occurred during the pipeline at {now}.\n"
            "No structured log entries found — check Snakemake's own output.\n"
            "Log directory: log/"
        )
    send_email(
        subject="❌ Snakemake pipeline failed",
        body=body
    )

