import os
import subprocess

SAMPLES_ID = [s[:7] for s in SAMPLES]

# Fallback for active sample lists
ACTIVE_OUTRIDER = ACTIVE_OUTRIDER if 'ACTIVE_OUTRIDER' in globals() else SAMPLES_ID
ACTIVE_FRASER   = ACTIVE_FRASER   if 'ACTIVE_FRASER' in globals() else SAMPLES_ID
ACTIVE_PCA      = ACTIVE_PCA      if 'ACTIVE_PCA' in globals() else SAMPLES_ID

# ----------------------
# TOOL VERSIONS
# ----------------------
TOOL_VERSIONS = {
    "picard": subprocess.getoutput("picard MarkDuplicates --version 2>&1 | head -1"),
    "samtools": subprocess.getoutput("samtools --version | head -1 | cut -d' ' -f2"),
    "mosdepth": subprocess.getoutput("mosdepth --version 2>&1 | head -1 | cut -d' ' -f3"),
    "multiqc": subprocess.getoutput("multiqc --version | head -1 | cut -d' ' -f2"),
    "python": subprocess.getoutput("python --version | cut -d' ' -f2"),
}

# --------------------------
# markDuplicate
# --------------------------
rule markDuplicate:
    input:
        aln = rules.alignment_star.output.bam,
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        bam = temp(FASTQ_DIR + "/../pipeline_v0/dup/{sample}/{sample}_dup.bam"),
        txt = FASTQ_DIR + "/../pipeline_v0/dup/{sample}/{sample}_dup.txt",
    conda:
        PIPELINE_DIR + "/envs/mark_env.yml"
    version: TOOL_VERSIONS["picard"]
    threads: 8
    resources:
        single_job = 2,
        tmpdir = OUTPUT_REP + "/dup/tmp"
    benchmark:
        "benchmarks/dup/{sample}.tsv"
    log:
        run_info = "log/dup/{sample}/log.txt",
        time = "log/dup/{sample}/time.txt"
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        picard MarkDuplicates \
            INPUT={input.aln} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.txt} \
            VALIDATION_STRINGENCY=LENIENT \
            REMOVE_DUPLICATES=false \
            TMP_DIR={resources.tmpdir} > {log.run_info} 2>&1 &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''

# --------------------------
# volcano
# --------------------------
rule volcano:
    input:
        outrider = FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab"
    output:
        png = FASTQ_DIR + "/../pipeline_v0/outrider/plot/{samples_id}.volcano.png",
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    params:
        scripts = PIPELINE_DIR
    benchmark:
        "benchmarks/volcano/{samples_id}.tsv"
    log:
        run_info = "log/volcano/{samples_id}/log.txt",
        time = "log/volcano/{samples_id}/time.txt"
    threads: 2
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        python {params.scripts}/scripts/create_volcano.py {input.outrider} {output.png} > {log.run_info} 2>&1 &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''

# --------------------------
# boxplot
# --------------------------
rule boxplot:
    input:
        outrider = expand(
            FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab",
            samples_id = ACTIVE_OUTRIDER
        )
    output:
        box = FASTQ_DIR + "/../pipeline_v0/outrider/plot/boxplot.png",
        filt = FASTQ_DIR + "/../pipeline_v0/outrider/plot/boxplot_filt.png",
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    params:
        dir = FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/",
        samples_id_str = ",".join(ACTIVE_OUTRIDER),
        scripts = PIPELINE_DIR
    threads: 2
    benchmark:
        "benchmarks/boxplot/boxplot.tsv"
    log:
        run_info = "log/boxplot/log.txt",
        time = "log/boxplot/time.txt"
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        python {params.scripts}/scripts/create_outrider_boxplots2.py {params.dir} {output.box} {output.filt} {params.samples_id_str} > {log.run_info} 2>&1 &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''

# --------------------------
# fraser_boxplot
# --------------------------
rule fraser_boxplot:
    input:
        fraser = expand(
            FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser.tab",
            samples_id = ACTIVE_FRASER
        )
    output:
        box = FASTQ_DIR + "/../pipeline_v0/fraser/plot/boxplot.png",
        filt = FASTQ_DIR + "/../pipeline_v0/fraser/plot/boxplot_filt.png",
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    params:
        dir = FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/",
        samples_id_str = ",".join(ACTIVE_FRASER),
        scripts = PIPELINE_DIR
    threads: 2
    benchmark:
        "benchmarks/boxplot_fraser/boxplot.tsv"
    log:
        run_info = "log/boxplot_fraser/log.txt",
        time = "log/boxplot_fraser/time.txt"
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        python {params.scripts}/scripts/create_fraser_boxplots5.py {params.dir} {output.box} {output.filt} {params.samples_id_str} > {log.run_info} 2>&1 &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''

# --------------------------
# rseqc
# --------------------------
rule rseqc:
    input:
        bam = rules.alignment_star.output.bam,
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        FASTQ_DIR + "/../pipeline_v0/rseqc/{sample}/{sample}.rseqc.results",
    params:
        bed = RSEQ_BED
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    threads: 8
    benchmark:
        "benchmarks/rseqc/{sample}.tsv"
    log:
        run_info = "log/rseqc/{sample}/log.txt",
        time = "log/rseqc/{sample}/time.txt"
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        read_distribution.py -i {input.bam} -r {params.bed} > {output} 2> {log.run_info} &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''

# --------------------------
# bam_stats
# --------------------------
rule bam_stats:
    input:
        bam = rules.alignment_star.output.bam,
        bai = FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        on_target = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_on_target.txt",
        padded = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_padded.txt",
        stats = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_stats.txt",
        hist = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_hist.txt",
        #mos = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}.mosdepth.global.dist.txt",
    params:
        bed = BED,
        padded_bed = PADDED,
        DI_bed = DI_BED
    threads: 12
    conda:
        PIPELINE_DIR + "/envs/comptage_env.yml"
    resources:
        single_job = 16,
        tmpdir = OUTPUT_REP + "/coverage/tmp",
        mem_gb = 500
    version: TOOL_VERSIONS["samtools"]
    log:
        on_target = "log/{sample}/stats_on_target.log",
        padded = "log/{sample}/stats_padded.log",
        hist = "log/{sample}/stats_hist.log",
        insert_size = "log/{sample}/stats_insert_size.log",
        time = "log/coverage/{sample}/time.txt"
    benchmark:
        "benchmarks/sam_stats/{sample}.tsv"
    shell:
        """
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        samtools view {input.bam} -b -@ {threads} -F 260 | \
        tee >(bedtools intersect -bed -u -abam stdin -b {params.bed} | wc -l > {output.on_target} 2>>{log.on_target}) \
            >(bedtools intersect -bed -u -abam stdin -b {params.padded_bed} | wc -l > {output.padded} 2>>{log.padded}) \
            >(bedtools coverage -hist -abam stdin -b {params.bed} | grep all > {output.hist} 2>>{log.hist}) \
        1>/dev/null &&
        samtools stats -F 4 -@ {threads} {input.bam} > {output.stats} 2>>{log.insert_size} &&
        #mosdepth --by {params.DI_bed} --threads {threads} --thresholds 1,10,20,30 {input.dir}/../pipeline_v0/coverage/{wildcards.sample}/{wildcards.sample} {input.bam} &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        """
# --------------------------
# multiqc
# --------------------------
rule multiqc:
    input:
        bam = expand(rules.alignment_star.output.bam, sample=SAMPLES),
        mark = expand(rules.markDuplicate.output.bam, sample=SAMPLES),
        dir = FASTQ_DIR
    output:
        html = FASTQ_DIR + "/../pipeline_v0/multiqc/multiqc_report.html",
        data = directory(FASTQ_DIR + "/../pipeline_v0/multiqc/multiqc_data"),
    params:
        bed = RSEQ_BED
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    threads: 8
    resources:
        single_job = 4
    benchmark:
        "benchmarks/multiqc/multiqc.tsv"
    log:
        run_info = "log/multiqc/log.txt",
        time = "log/multiqc/time.txt"
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        export TMPDIR={input.dir}/tmp &&
        multiqc --force {input.dir}/../pipeline_v0 -o {FASTQ_DIR}/../pipeline_v0/multiqc > {log.run_info} 2>&1 &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''

# --------------------------
# generate_and_run_param_notebook (PCA)
# --------------------------
RUN_NAME = os.path.basename(OUTPUT_REP)
rule generate_and_run_param_notebook:
    input:
        runs_folder = RUNS_DIR,
        mapping_file = GTF,
        blacklist  = PCA_blacklist,
        htseq = expand(rules.htseq_gene.output.gene, sample=SAMPLES), 
    output:
        executed_nb = FASTQ_DIR + f"/../pipeline_v0/notebooks/notebook_pca_{RUN_NAME}.ipynb",
    params:
        samples_run  = RUN_NAME,
        excluded_run = "20231122_RUN17_NextSeq_Mid_8RNASEQ",
        keyword      = "MOINS",
        output_html  = f"PCA-{RUN_NAME}.html",
        scripts = PIPELINE_DIR
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    log:
        run_info = f"log/notebook/log_{RUN_NAME}.txt",
        time = f"log/notebook/time_{RUN_NAME}.txt"
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        python {params.scripts}/scripts/generate_pca_notebook.py \
            --base {input.runs_folder} \
            --mapping_file {input.mapping_file} \
            --blacklist_file {input.blacklist} \
            --samples_run {params.samples_run} \
            --excluded_run {params.excluded_run} \
            --keyword {params.keyword} \
            --output_html {params.output_html} \
            --notebook_path {output.executed_nb} > {log.run_info} 2>&1 &&
        jupyter nbconvert --to notebook --execute --inplace {output.executed_nb} >> {log.run_info} 2>&1 &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''

# --------------------------
# generate_metrics
# --------------------------
rule generate_metrics:
    input:
        run_dir = os.path.dirname(FASTQ_DIR),
        tpm = rules.matrix_tpm.output.gene_gene,
        fraser = expand(FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser.tab", samples_id = ACTIVE_FRASER),
        outrider = expand(FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab", samples_id = ACTIVE_OUTRIDER)
    output:
        metrics = FASTQ_DIR + "/../pipeline_v0/metrics/qc_summary.tsv",
    params:
        scripts = PIPELINE_DIR,
        di_file = PANELAPP
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    log:
        run_info = "log/metrics/log.txt",
        time = "log/metrics/time.txt"
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        python {params.scripts}/scripts/recup_metrics.py \
            --run_path {input.run_dir} \
            --tpm_file {input.tpm} \
            --di {params.di_file} \
            --output {output.metrics} > {log.run_info} 2>&1 &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''

