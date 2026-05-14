import os
import subprocess

# ACTIVE_FRASER, ACTIVE_OUTRIDER, ACTIVE_PCA, SAMPLES_ID_BLOOD
# are defined in pipeline.smk with blood-sample filtering.
# Do NOT redefine them here.

# ----------------------
# TOOL VERSIONS
# ----------------------
TOOL_VERSIONS = {
    "OUTRIDER": subprocess.getoutput('Rscript -e "packageVersion(\'OUTRIDER\')"'),
    "FRASER": subprocess.getoutput('Rscript -e "packageVersion(\'FRASER\')"'),
    "python": subprocess.getoutput("python --version | cut -d' ' -f2"),
}


# --------------------------
# OUTRIDER main rule
# --------------------------
rule outrider_hyper:
    version: TOOL_VERSIONS["OUTRIDER"]
    input:
        matrice = rules.matrix_featurecounts.output.matrix
    output:
        out_file     = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq.tab",
        out_file_all = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq_all.tsv"
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    params:
        scripts         = PIPELINE_DIR,
        yaml_config     = config.get("outrider_yaml", ""),
        current_samples = lambda wc: ",".join(SAMPLES_ID_BLOOD)
    benchmark:
        "benchmarks/outrider_hyper/benchmark_outrider.tsv"
    log:
        run_info = "log/outrider_hyper/log.txt",
        time     = "log/outrider_hyper/time.txt"
    threads: 1
    resources:
        single_job = 8
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time}
        OUTRIDER_CURRENT_SAMPLES="{params.current_samples}" \
        Rscript {params.scripts}/scripts/outrider_newVersion.R \
            {input.matrice} \
            {output.out_file} \
            {params.yaml_config}
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''

# --------------------------
# OUTRIDER annotation rule
# --------------------------
rule annotation_outrider_hyper:
    version: TOOL_VERSIONS["OUTRIDER"]
    input:
        out_file = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq_all.tsv",
        fraser   = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_results_all.tsv"
    output:
        annot = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq_annot.tsv",
        output_files_dir = directory(FASTQ_DIR + "/../pipeline_v0/outrider_hyper/filesbysample/"),
        files = touch(expand(FASTQ_DIR + "/../pipeline_v0/outrider_hyper/filesbysample/{samples_id}.outrider.tab", samples_id=ACTIVE_OUTRIDER))
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    params:
        gtf = GTF,
        panel = PANELAPP,
        pli = PLI,
        hpo = HPO,
        pheno = PHENO,
        scripts = PIPELINE_DIR
    log:
        run_info = "log/outrider_hyper/log_annotation.txt",
        time = "log/outrider_hyper/time_annotation.txt"
    threads: 1
    resources:
        single_job = 4
    shell:
        """
        set -euo pipefail
        SNAKEMAKE_LOG="{log.run_info}"
        echo "start : $(date +"%d-%m-%y %T")" > {log.time}
        python {params.scripts}/scripts/annotation_outrider_hits_bis2.py \
            -i {input.out_file} -g {params.gtf} -d {params.panel} \
            -p {params.pli} -o {params.pheno} -m {params.hpo} \
            -f {output.annot} --fraser {input.fraser} >> "$SNAKEMAKE_LOG" 2>&1
        python {params.scripts}/scripts/outrider_1file.py \
            -i {output.annot} -b {output.output_files_dir}/ >> "$SNAKEMAKE_LOG" 2>&1
        for f in {output.files}; do [ -f "$f" ] || touch "$f"; done
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        """

# --------------------------
# FRASER config generation
# --------------------------
rule fraser_config_hyper:
    input:
        bam = expand(rules.alignment_star.output.bam, sample=SAMPLES),
        bai = expand(FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        dir = FASTQ_DIR
    output:
        fraser = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_config.txt"
    conda:
        PIPELINE_DIR + "/envs/fraser_env.yml"
    params:
        blacklist = fraser_blacklist,
        fraser = fraser_count_hyper,
        scripts = PIPELINE_DIR
    benchmark:
        "benchmarks/fraser_config/benchmark.tsv"
    log:
        run_info = "log/fraser_config/log.txt",
        time = "log/fraser_config/time.txt"
    threads: 1
    resources:
        single_job = 8
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time}
        python {params.scripts}/scripts/fraser_config_create.py {output.fraser} MOINS {input.dir}/../.. {params.fraser} {params.blacklist} > {log.run_info} 2>&1
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''

# --------------------------
# FRASER main rule (core)
# --------------------------
rule fraser_hyper:
    version: TOOL_VERSIONS["FRASER"]
    input:
        config = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_config.txt"
    output:
        fraser     = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_results_aberrant.tsv",
        fraser_all = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_results_all.tsv"
    conda:
        PIPELINE_DIR + "/envs/fraser_env.yml"
    params:
        dir         = fraser_count_hyper,
        scripts     = PIPELINE_DIR,
        yaml_config = config.get("fraser_yaml", "")
    benchmark:
        "benchmarks/fraser_hyper/benchmark_fraser.tsv"
    log:
        run_info = "log/fraser_hyper/log_fraser.txt",
        time     = "log/fraser_hyper/time_fraser.txt"
    threads: 10
    resources:
        single_job = 8
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time}
        python {params.scripts}/scripts/backup_fraser_counts.py {params.dir}
        FRASER_THREADS={threads} Rscript {params.scripts}/scripts/fraser_newVer.R \
            {params.dir} \
            {input.config} \
            {params.yaml_config}
        mkdir -p $(dirname {output.fraser})
        cp {params.dir}/fraser_results_aberrant.tsv {output.fraser}
        cp {params.dir}/fraser_results_all.tsv      {output.fraser_all}
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''

# --------------------------
# FRASER annotation rule
# --------------------------
rule annotation_fraser_hyper:
    version: TOOL_VERSIONS["FRASER"]
    input:
        fraser   = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_results_all.tsv",
        outrider = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq_all.tsv"
    output:
        annot_fraser = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_annot.tsv",
        files_dir = directory(FASTQ_DIR + "/../pipeline_v0/fraser_hyper/filesbysample/"),
        files = touch(expand(FASTQ_DIR + "/../pipeline_v0/fraser_hyper/filesbysample/{samples_id}.fraser.tab", samples_id=ACTIVE_FRASER))
    conda:
        PIPELINE_DIR + "/envs/fraser_env.yml"
    params:
        blacklist = fraser_blacklist,
        gtf = GTF,
        panel = PANELAPP,
        pli = PLI,
        hpo = HPO,
        pheno = PHENO,
        scripts = PIPELINE_DIR
    log:
        run_info = "log/fraser_hyper/log_annotation_fraser.txt",
        time = "log/fraser_hyper/time_annotation_fraser.txt"
    threads: 1
    resources:
        single_job = 4
    shell:
        """
        set -euo pipefail
        SNAKEMAKE_LOG="{log.run_info}"
        echo "start : $(date +"%d-%m-%y %T")" > {log.time}
        python {params.scripts}/scripts/annotation_fraser_test4.py \
            -f {input.fraser} -g {params.gtf} -d {params.panel} \
            -p {params.pli} -m {params.hpo} -o {params.pheno} \
            --outrider {input.outrider} --output {output.annot_fraser} >> "$SNAKEMAKE_LOG" 2>&1
        python {params.scripts}/scripts/fraser_1file.py \
            -i {output.annot_fraser} -b {output.files_dir}/ >> "$SNAKEMAKE_LOG" 2>&1
        for f in {output.files}; do [ -f "$f" ] || touch "$f"; done
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        """

# --------------------------
# Per-sample output rule
# --------------------------
rule rnaseq_per_sample_hyper:
    input:
        fraser   = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_results_all.tsv",
        outrider = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq_all.tsv"
    output:
        files_dir = directory(FASTQ_DIR + "/../pipeline_v0/per_sample_hyper/"),
        zip_file  = touch(FASTQ_DIR + "/../pipeline_v0/per_sample_hyper/results.zip")
    conda:
        PIPELINE_DIR + "/envs/python_plot_env.yml"
    params:
        gtf        = GTF,
        gnomad     = config.get("gnomad", ""),
        mendeliome = config.get("mendeliome", ""),
        samples    = config.get("samples_file", ""),
        mode       = config.get("per_sample_mode", "all"),
        pvalue     = config.get("per_sample_pvalue", ""),
        workers    = config.get("per_sample_workers", 4),
        scripts    = PIPELINE_DIR
    log:
        run_info = "log/per_sample_hyper/log.txt",
        time     = "log/per_sample_hyper/time.txt"
    threads: 4
    resources:
        single_job = 4
    shell:
        """
        set -euo pipefail
        SNAKEMAKE_LOG="{log.run_info}"
        echo "start : $(date +"%d-%m-%y %T")" > {log.time}
        mkdir -p {output.files_dir}
        python {params.scripts}/scripts/rnaseq_analysis_per_sample.py \
            --fraser   {input.fraser} \
            --outrider {input.outrider} \
            --gtf      {params.gtf} \
            --output   {output.files_dir} \
            --mode     {params.mode} \
            --workers  {params.workers} \
            $([ -n "{params.samples}"    ] && echo "--samples {params.samples}")       \
            $([ -n "{params.gnomad}"     ] && echo "--gnomad {params.gnomad}")         \
            $([ -n "{params.mendeliome}" ] && echo "--mendeliome {params.mendeliome}") \
            $([ -n "{params.pvalue}"     ] && echo "--pvalue {params.pvalue}")         \
            >> "$SNAKEMAKE_LOG" 2>&1
        [ -f "{output.zip_file}" ] || \
            mv $(ls -t {output.files_dir}/*.zip 2>/dev/null | head -1) \
               {output.zip_file} 2>/dev/null || \
            touch "{output.zip_file}"
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        """

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#






















