# --- fetch tool versions at parse time ---
TOOL_VERSIONS = {
    "fastqc": get_version_from_env("fastqc --version | head -1 | cut -d' ' -f2"),
    "fastp": get_version_from_env("fastp --version"),
}

# ----------------------
rule fastqc_report:
    """
    Run FastQC on raw FASTQ data
    """
    version: TOOL_VERSIONS["fastqc"]
    input:
        R1 = FASTQ_DIR + "/{sample}_R1.fastq.gz",
        R2 = FASTQ_DIR + "/{sample}_R2.fastq.gz"
    output:
        report_dir = directory(FASTQ_DIR + "/../pipeline_v0/metrics/qc/fastqc/{sample}"),
        R1 = FASTQ_DIR + "/../pipeline_v0/metrics/qc/fastqc/{sample}/{sample}_R1_fastqc.html",
        R2 = FASTQ_DIR + "/../pipeline_v0/metrics/qc/fastqc/{sample}/{sample}_R2_fastqc.html"
    conda:
        PIPELINE_DIR + "/envs/fastqc_env.yml"
    log:
        run_info = "log/fastqc/{sample}/log.txt",
        time     = "log/fastqc/{sample}/time.txt"
    threads: 8
    params:
        log_start = lambda wc, input, threads: log_start("fastqc_report", wc, threads),
        log_end   = LOG_END
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p {output.report_dir}
        fastqc -o {output.report_dir} -t {threads} {input.R1} {input.R2} >> {log.run_info} 2>&1
        {params.log_end}
        """

# ----------------------
rule fastp:
    """
    Trim and filter reads with quality control using fastp
    """
    version: TOOL_VERSIONS["fastp"]
    input:
        R1 = FASTQ_DIR + "/{sample}_R1.fastq.gz",
        R2 = FASTQ_DIR + "/{sample}_R2.fastq.gz"
    output:
        R1   = temp(FASTQ_DIR + "/../pipeline_v0/trimmed_fastq/{sample}_R1_trimmed.fastq.gz"),
        R2   = temp(FASTQ_DIR + "/../pipeline_v0/trimmed_fastq/{sample}_R2_trimmed.fastq.gz"),
        html = FASTQ_DIR + "/../pipeline_v0/metrics/qc/fastp/{sample}/report.html",
        json = FASTQ_DIR + "/../pipeline_v0/metrics/qc/fastp/{sample}/report.json"
    params:
        cut_right_size          = 10,
        cut_right_mean_quality  = 15,
        qualified_quality_phred = 15,
        minimum_read_length     = 30,
        log_start = lambda wc, input, threads: log_start("fastp", wc, threads),
        log_end   = LOG_END
    threads: 2
    conda:
        PIPELINE_DIR + "/envs/fastp_env.yml"
    log:
        run_info = "log/fastp/{sample}/fastp.log",
        time     = "log/fastp/{sample}/fastp_time.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        fastp \
            -i {input.R1} \
            -I {input.R2} \
            -o {output.R1} \
            -O {output.R2} \
            -w {threads} \
            -l {params.minimum_read_length} \
            --detect_adapter_for_pe \
            --trim_front1 5 \
            --trim_front2 5 \
            --cut_right \
            --cut_right_window_size {params.cut_right_size} \
            --cut_right_mean_quality {params.cut_right_mean_quality} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            -h {output.html} \
            -j {output.json} \
            >> {log.run_info} 2>&1
        {params.log_end}
        """

# ----------------------
rule fastqc_trim_report:
    """
    FastQC report after trimming reads
    """
    version: TOOL_VERSIONS["fastqc"]
    input:
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2
    output:
        report_dir = directory(FASTQ_DIR + "/../pipeline_v0/metrics/qc/fastqc_trim/{sample}"),
        R1 = FASTQ_DIR + "/../pipeline_v0/metrics/qc/fastqc_trim/{sample}/{sample}_R1_trimmed_fastqc.html",
        R2 = FASTQ_DIR + "/../pipeline_v0/metrics/qc/fastqc_trim/{sample}/{sample}_R2_trimmed_fastqc.html"
    conda:
        PIPELINE_DIR + "/envs/fastqc_env.yml"
    log:
        run_info = "log/fastqc_trim/{sample}/log.txt",
        time     = "log/fastqc_trim/{sample}/time.txt"
    threads: 8
    params:
        log_start = lambda wc, input, threads: log_start("fastqc_trim_report", wc, threads),
        log_end   = LOG_END
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p {output.report_dir}
        fastqc -o {output.report_dir} -t {threads} {input.R1} {input.R2} >> {log.run_info} 2>&1
        {params.log_end}
        """
