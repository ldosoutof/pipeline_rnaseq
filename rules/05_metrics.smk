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
    params:
        log_start = lambda wc, input, threads: log_start("markDuplicate", wc, threads),
        log_end   = LOG_END
    shell:
        """
        set -euo pipefail
        {params.log_start}
        picard MarkDuplicates \
            INPUT={input.aln} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.txt} \
            VALIDATION_STRINGENCY=LENIENT \
            REMOVE_DUPLICATES=false \
            TMP_DIR={resources.tmpdir} >> {log.run_info} 2>&1
        {params.log_end}
        """

# --------------------------
# volcano
# --------------------------
rule volcano:
    input:
        outrider = FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab"
    output:
        png = touch(FASTQ_DIR + "/../pipeline_v0/outrider/plot/{samples_id}.volcano.png"),
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    params:
        scripts = PIPELINE_DIR,
        log_start = lambda wc, input, threads: log_start("volcano", wc, threads),
        log_end   = LOG_END
    benchmark:
        "benchmarks/volcano/{samples_id}.tsv"
    log:
        run_info = "log/volcano/{samples_id}/log.txt",
        time = "log/volcano/{samples_id}/time.txt"
    threads: 2
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p $(dirname {output.png})
        if [ ! -s {input.outrider} ]; then
            echo "[INFO] No OUTRIDER events for {wildcards.samples_id} — creating empty placeholder" >> {log.run_info}
            touch {output.png}
        else
            python {params.scripts}/scripts/create_volcano.py \
                {input.outrider} {output.png} >> {log.run_info} 2>&1
        fi
        {params.log_end}
        """

# --------------------------
# boxplot
# --------------------------
rule boxplot:
    input:
        outrider = lambda wc: expand(
            FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab",
            samples_id=ACTIVE_OUTRIDER
        )
    output:
        box  = touch(FASTQ_DIR + "/../pipeline_v0/outrider/plot/boxplot.png"),
        filt = touch(FASTQ_DIR + "/../pipeline_v0/outrider/plot/boxplot_filt.png"),
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    params:
        dir            = FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/",
        samples_id_str = lambda wc: ",".join(ACTIVE_OUTRIDER),
        scripts        = PIPELINE_DIR,
        log_start      = lambda wc, input, threads: log_start("boxplot", wc, threads),
        log_end   = LOG_END
    threads: 2
    benchmark:
        "benchmarks/boxplot/boxplot.tsv"
    log:
        run_info = "log/boxplot/log.txt",
        time = "log/boxplot/time.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p $(dirname {output.box})
        python {params.scripts}/scripts/create_outrider_boxplots2.py \
            {params.dir} {output.box} {output.filt} \
            {params.samples_id_str} >> {log.run_info} 2>&1 || \
            (echo "[WARN] boxplot failed — creating empty placeholders" >> {log.run_info} && \
             touch {output.box} {output.filt})
        {params.log_end}
        """

# --------------------------
# fraser_boxplot
# --------------------------
rule fraser_boxplot:
    input:
        fraser = lambda wc: expand(
            FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser.tab",
            samples_id=ACTIVE_FRASER
        )
    output:
        box  = touch(FASTQ_DIR + "/../pipeline_v0/fraser/plot/boxplot.png"),
        filt = touch(FASTQ_DIR + "/../pipeline_v0/fraser/plot/boxplot_filt.png"),
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    params:
        dir            = FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/",
        samples_id_str = lambda wc: ",".join(ACTIVE_FRASER),
        scripts        = PIPELINE_DIR,
        log_start      = lambda wc, input, threads: log_start("fraser_boxplot", wc, threads),
        log_end        = LOG_END
    threads: 2
    benchmark:
        "benchmarks/boxplot_fraser/boxplot.tsv"
    log:
        run_info = "log/boxplot_fraser/log.txt",
        time = "log/boxplot_fraser/time.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p $(dirname {output.box})
        python {params.scripts}/scripts/create_fraser_boxplots5.py \
            {params.dir} {output.box} {output.filt} \
            {params.samples_id_str} >> {log.run_info} 2>&1 || \
            (echo "[WARN] fraser_boxplot failed — creating empty placeholders" >> {log.run_info} && \
             touch {output.box} {output.filt})
        {params.log_end}
        """

# --------------------------
# rseqc
# --------------------------
rule rseqc:
    input:
        bam = rules.alignment_star.output.bam,
    output:
        FASTQ_DIR + "/../pipeline_v0/rseqc/{sample}/{sample}.rseqc.results",
    params:
        bed = RSEQ_BED,
        log_start = lambda wc, input, threads: log_start("rseqc", wc, threads),
        log_end   = LOG_END
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    threads: 8
    benchmark:
        "benchmarks/rseqc/{sample}.tsv"
    log:
        run_info = "log/rseqc/{sample}/log.txt",
        time = "log/rseqc/{sample}/time.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        read_distribution.py -i {input.bam} -r {params.bed} \
            > {output} 2>> {log.run_info}
        {params.log_end}
        """

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
        on_target = touch(FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_on_target.txt"),
        padded    = touch(FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_padded.txt"),
        stats     = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_stats.txt",
        hist      = touch(FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_hist.txt"),
        #mos = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}.mosdepth.global.dist.txt",
    params:
        bed = BED,
        padded_bed = PADDED,
        DI_bed = DI_BED,
        log_start = lambda wc, input, threads: log_start("bam_stats", wc, threads),
        log_end   = LOG_END
    threads: 12
    conda:
        PIPELINE_DIR + "/envs/comptage_env.yml"
    resources:
        single_job = 16,
        tmpdir = OUTPUT_REP + "/coverage/tmp",
        mem_gb = 500
    version: TOOL_VERSIONS["samtools"]
    log:
        on_target  = "log/{sample}/stats_on_target.log",
        padded     = "log/{sample}/stats_padded.log",
        hist       = "log/{sample}/stats_hist.log",
        insert_size = "log/{sample}/stats_insert_size.log",
        run_info   = "log/coverage/{sample}/log.txt",
        time       = "log/coverage/{sample}/time.txt"
    benchmark:
        "benchmarks/sam_stats/{sample}.tsv"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p $(dirname {output.on_target})

        # Run tee with process substitutions; capture PIDs for explicit wait
        samtools view {input.bam} -b -@ {threads} -F 260 | \
        tee \
            >(bedtools intersect -bed -u -abam stdin -b {params.bed} \
                | wc -l > {output.on_target} 2>>{log.on_target} ; ) \
            >(bedtools intersect -bed -u -abam stdin -b {params.padded_bed} \
                | wc -l > {output.padded} 2>>{log.padded} ; ) \
            >(bedtools coverage -hist -abam stdin -b {params.bed} \
                | grep all > {output.hist} 2>>{log.hist} ; ) \
        1>/dev/null
        # Wait for all background subshells to flush their output files
        wait

        # Ensure output files exist even if no reads intersected
        touch {output.on_target} {output.padded} {output.hist}

        samtools stats -F 4 -@ {threads} {input.bam} \
            > {output.stats} 2>>{log.insert_size}
        {params.log_end}
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
        bed = RSEQ_BED,
        fastq_dir = FASTQ_DIR,
        log_start = lambda wc, input, threads: log_start("multiqc", wc, threads),
        log_end   = LOG_END
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
        """
        set -euo pipefail
        {params.log_start}
        export TMPDIR={params.fastq_dir}/tmp
        multiqc --force {params.fastq_dir}/../pipeline_v0 \
            -o {FASTQ_DIR}/../pipeline_v0/multiqc >> {log.run_info} 2>&1
        {params.log_end}
        """

# --------------------------
# generate_and_run_param_notebook (PCA)
# --------------------------
RUN_NAME = os.path.basename(OUTPUT_REP)
rule generate_and_run_param_notebook:
    input:
        runs_folder  = RUNS_DIR,
        mapping_file = GTF,
#        htseq = expand(rules.htseq_gene.output.gene, sample=SAMPLES),
        htseq = expand(FASTQ_DIR + "/../pipeline_v0/htseq/{sample}/{sample}_gene_counts.txt", sample=SAMPLES)
    output:
        executed_nb = FASTQ_DIR + f"/../pipeline_v0/notebooks/notebook_pca_{RUN_NAME}.ipynb",
    params:
        samples_run  = RUN_NAME,
        excluded_run = "20231122_RUN17_NextSeq_Mid_8RNASEQ",
        keyword      = "MOINS",
        output_html  = f"PCA-{RUN_NAME}.html",
        blacklist    = PCA_blacklist,
        scripts      = PIPELINE_DIR,
        log_start    = lambda wc, input, threads: log_start("generate_and_run_param_notebook", wc, threads),
        log_end      = LOG_END
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    log:
        run_info = f"log/notebook/log_{RUN_NAME}.txt",
        time = f"log/notebook/time_{RUN_NAME}.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        BL="{params.blacklist}"
        if [ ! -f "$BL" ]; then
            mkdir -p "$(dirname "$BL")"
            touch "$BL"
        fi
        python {params.scripts}/scripts/generate_pca_notebook.py \
            --base {input.runs_folder} \
            --mapping_file {input.mapping_file} \
            --blacklist_file "$BL" \
            --samples_run {params.samples_run} \
            --excluded_run {params.excluded_run} \
            --keyword {params.keyword} \
            --output_html {params.output_html} \
            --notebook_path {output.executed_nb} >> {log.run_info} 2>&1
        jupyter nbconvert --to notebook --execute --inplace \
            {output.executed_nb} >> {log.run_info} 2>&1
        {params.log_end}
        """

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
        scripts     = PIPELINE_DIR,
        di_file     = PANELAPP,
        gateway     = config.get("prometheus_gateway", "http://localhost:9091"),
        log_start   = lambda wc, input, threads: log_start("generate_metrics", wc, threads),
        log_end     = LOG_END
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    log:
        run_info = "log/metrics/log.txt",
        time = "log/metrics/time.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        python {params.scripts}/scripts/recup_metrics.py \
            --run_path {input.run_dir} \
            --tpm_file {input.tpm} \
            --di {params.di_file} \
            --output {output.metrics} >> {log.run_info} 2>&1
        python {params.scripts}/monitoring/push_metrics.py \
            --run_path {input.run_dir} \
            --gateway  {params.gateway} >> {log.run_info} 2>&1 || \
            echo "[WARN] Pushgateway unavailable — metrics not pushed" >> {log.run_info}
        {params.log_end}
        """

