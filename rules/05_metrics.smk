# ../rules/05_metrics.smk
import os
import subprocess

SAMPLES_ID = [s[:7] for s in SAMPLES]

# Expect ACTIVE_OUTRIDER, ACTIVE_FRASER, ACTIVE_PCA to be defined earlier (pipeline.smk or included file).
# If not defined, fall back to SAMPLES_ID
try:
    ACTIVE_OUTRIDER
except NameError:
    ACTIVE_OUTRIDER = SAMPLES_ID
try:
    ACTIVE_FRASER
except NameError:
    ACTIVE_FRASER = SAMPLES_ID
try:
    ACTIVE_PCA
except NameError:
    ACTIVE_PCA = SAMPLES_ID


rule markDuplicate:
    """
    Marquage des duplicats
    """
    input:
        aln = FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam",
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        bam = temp(FASTQ_DIR + "/../pipeline_v0/dup/{sample}/{sample}_dup.bam"),
        txt = FASTQ_DIR + "/../pipeline_v0/dup/{sample}/{sample}_dup.txt"
    conda:
        PIPELINE_DIR + "/envs/mark_env.yml"
    version:
        # Use picard to get a version string if available
        subprocess.getoutput("picard MarkDuplicates --version 2>&1 | head -1")
    resources:
        single_job = 2,
        tmpdir = OUTPUT_REP + "/dup/tmp"
    benchmark:
        "benchmarks/dup/{sample}.tsv"
    log:
        run_info = "log/dup/{sample}/log.txt",
        time = "log/dup/{sample}/time.txt"
    threads: 8
    shell:
        '''
        echo "start : $(date +"%d-%m-%y   %T")" > {log.time} && \
        mkdir -p {input.dir}/../pipeline_v0/dup/{wildcards.sample} && \
        mkdir -p {input.out}/dup/{wildcards.sample} && \
        picard MarkDuplicates \
            INPUT={input.aln} \
            OUTPUT={input.dir}/../pipeline_v0/dup/{wildcards.sample}/{wildcards.sample}_dup.bam \
            METRICS_FILE={input.dir}/../pipeline_v0/dup/{wildcards.sample}/{wildcards.sample}_dup.txt \
            VALIDATION_STRINGENCY=LENIENT \
            REMOVE_DUPLICATES=false \
            TMP_DIR={resources.tmpdir} > {log.run_info} 2>&1 && \
        echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
        '''


rule volcano:
    """
    Per-sample volcano plot from OUTRIDER per-sample results.
    This rule is instantiated per {samples_id}.
    """
    input:
        outrider = FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab"
    output:
        png = FASTQ_DIR + "/../pipeline_v0/outrider/plot/{samples_id}.volcano.png"
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
        mkdir -p $(dirname {output.png}) && \
        python {params.scripts}/scripts/create_volcano.py {input.outrider} {output.png}
        '''


rule boxplot:
    """
    Generate OUTRIDER boxplots comparing sample counts for the run and other samples.
    This is a single aggregated job using ACTIVE_PCA samples as reference.
    """
    input:
        # list of outrider per-sample files for ACTIVE_PCA
        outrider = expand(
            FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab",
            samples_id = ACTIVE_PCA
        )
    output:
        box = FASTQ_DIR + "/../pipeline_v0/outrider/plot/boxplot.png",
        filt = FASTQ_DIR + "/../pipeline_v0/outrider/plot/boxplot_filt.png"
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    benchmark:
        "benchmarks/boxplot/boxplot.tsv"
    log:
        run_info = "log/boxplot/log.txt",
        time = "log/boxplot/time.txt"
    params:
        outrider_files = ",".join(expand(FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab", samples_id = ACTIVE_PCA)),
        dir = FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/",
        samples_id_str = ",".join(ACTIVE_PCA),
        scripts = PIPELINE_DIR
    threads: 2
    shell:
        '''
        mkdir -p $(dirname {output.box}) && \
        python {params.scripts}/scripts/create_outrider_boxplots2.py {params.dir} {output.box} {output.filt} {params.samples_id_str}
        '''


rule fraser_boxplot:
    """
    Generate FRASER boxplots comparing sample counts for the run and other samples.
    Aggregated single job using ACTIVE_FRASER list.
    """
    input:
        fraser = expand(
            FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser.tab",
            samples_id = ACTIVE_FRASER
        )
    output:
        box = FASTQ_DIR + "/../pipeline_v0/fraser/plot/boxplot.png",
        filt = FASTQ_DIR + "/../pipeline_v0/fraser/plot/boxplot_filt.png"
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    benchmark:
        "benchmarks/boxplot_fraser/boxplot.tsv"
    log:
        run_info = "log/boxplot_fraser/log.txt",
        time = "log/boxplot_fraser/time.txt"
    params:
        fraser_files = ",".join(expand(FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser.tab", samples_id = ACTIVE_FRASER)),
        dir = FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/",
        samples_id_str = ",".join(ACTIVE_FRASER),
        scripts = PIPELINE_DIR
    threads: 2
    shell:
        '''
        mkdir -p $(dirname {output.box}) && \
        python {params.scripts}/scripts/create_fraser_boxplots5.py {params.dir} {output.box} {output.filt} {params.samples_id_str}
        '''


rule rseqc:
    """
    distribution type
    """
    input:
        bam = FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam",
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        FASTQ_DIR + "/../pipeline_v0/rseqc/{sample}/{sample}.rseqc.results"
    params:
        bed = RSEQ_BED
    resources:
        single_job = 2
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    benchmark:
        "benchmarks/rseqc/{sample}.tsv"
    log:
        run_info = "log/rseqc/{sample}/log.txt",
        time = "log/rseqc/{sample}/time.txt"
    threads: 8
    shell:
        'read_distribution.py -i {input.bam} -r {params.bed} > {output}'


rule bam_stats:
    """
    Generate coverage/stats files per sample.
    """
    input:
        bam = FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam",
        bai = FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        on_target = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_on_target.txt",
        padded = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_padded.txt",
        stats = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_stats.txt",
        hist = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_hist.txt",
        mos = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}.mosdepth.global.dist.txt"
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
    version:
        subprocess.getoutput("samtools --version | head -1 | cut -d' ' -f2")
    log:
        on_target = "log/{sample}/stats_on_target.log",
        padded = "log/{sample}/stats_padded.log",
        hist = "log/{sample}/stats_hist.log",
        insert_size = "log/{sample}/stats_insert_size.log"
    benchmark:
        "benchmarks/sam_stats/{sample}.tsv"
    message:
        "--------- stats : {wildcards.sample} ---------"
    shell:
        r"""
        mkdir -p {input.dir}/../pipeline_v0/coverage/{wildcards.sample} && \
        samtools view {input.bam} -b -@ {threads} -F 260 | \
        tee >(bedtools intersect -bed -u -abam stdin -b {params.bed} | wc -l > {input.dir}/../pipeline_v0/coverage/{wildcards.sample}/{wildcards.sample}_on_target.txt 2>>{log.on_target}) \
            >(bedtools intersect -bed -u -abam stdin -b {params.padded_bed} | wc -l > {input.dir}/../pipeline_v0/coverage/{wildcards.sample}/{wildcards.sample}_padded.txt 2>>{log.padded}) \
            >(bedtools coverage -hist -abam stdin -b {params.bed} | grep all > {input.dir}/../pipeline_v0/coverage/{wildcards.sample}/{wildcards.sample}_hist.txt 2>>{log.hist}) \
        1>/dev/null && \
        samtools stats -F 4 -@ {threads} {input.bam} > {input.dir}/../pipeline_v0/coverage/{wildcards.sample}/{wildcards.sample}_stats.txt 2>>{log.insert_size} && \
        mosdepth --by {params.DI_bed} --threads 8 --thresholds 1,10,20,30 {input.dir}/../pipeline_v0/coverage/{wildcards.sample}/{wildcards.sample} {input.bam}
        """


rule multiqc:
    """
    multiqc
    """
    input:
        bam = expand(FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        mark = expand(rules.markDuplicate.output.bam, sample=SAMPLES),
        dir = FASTQ_DIR
    output:
        html = FASTQ_DIR + "/../pipeline_v0/multiqc/multiqc_report.html",
        data = directory(FASTQ_DIR + "/../pipeline_v0/multiqc/multiqc_data")
    params:
        bed = RSEQ_BED
    resources:
        single_job = 4
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    benchmark:
        "benchmarks/multiqc/multiqc.tsv"
    log:
        run_info = "log/multiqc/log.txt",
        time = "log/multiqc/time.txt"
    threads: 8
    shell:
        """
        export TMPDIR={input.dir}/tmp && \
        multiqc --force {input.dir}/../pipeline_v0 -o {FASTQ_DIR}/../pipeline_v0/multiqc
        """


RUN_NAME = os.path.basename(OUTPUT_REP)
rule generate_and_run_param_notebook:
    input:
        runs_folder = RUNS_DIR,
        mapping_file = GTF,
        blacklist  = PCA_blacklist
    output:
        executed_nb = FASTQ_DIR + f"/../pipeline_v0/notebooks/notebook_pca_{RUN_NAME}.ipynb"
    params:
        samples_run  = RUN_NAME,
        excluded_run = "20231122_RUN17_NextSeq_Mid_8RNASEQ",
        keyword      = "MOINS",
        output_html  = f"PCA-{RUN_NAME}.html",
        scripts = PIPELINE_DIR
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    shell:
        """
        python {params.scripts}/scripts/generate_pca_notebook.py \
            --base {input.runs_folder} \
            --mapping_file {input.mapping_file} \
            --blacklist_file {input.blacklist} \
            --samples_run {params.samples_run} \
            --excluded_run {params.excluded_run} \
            --keyword {params.keyword} \
            --output_html {params.output_html} \
            --notebook_path {output.executed_nb}

        jupyter nbconvert --to notebook --execute --inplace {output.executed_nb}
        """


rule generate_metrics:
    """
    Aggregate QC metrics for the run.
    """
    input:
        run_dir = os.path.dirname(FASTQ_DIR),
        tpm = FASTQ_DIR + "/../pipeline_v0/kallisto_bed/marice_gene_tpm_gene.tsv",
        fraser = expand(FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser.tab", samples_id = ACTIVE_FRASER),
        outrider = expand(FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab", samples_id = ACTIVE_OUTRIDER)
    output:
        metrics = FASTQ_DIR + "/../pipeline_v0/metrics/qc_summary.tsv"
    params:
        scripts = PIPELINE_DIR
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    shell:
        """
        python {params.scripts}/scripts/script_recup_metrics.py --run_path {input.run_dir} --tpm_file {input.tpm} --output {output.metrics}
        """

