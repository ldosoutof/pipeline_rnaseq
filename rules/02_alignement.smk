TOOL_VERSIONS = {
    "star": get_version_from_env(
        PIPELINE_DIR + "/envs/star_env.yml",
        "STAR --version | head -1 | cut -d' ' -f2"
    ),
    "samtools": get_version_from_env(
        PIPELINE_DIR + "/envs/samtools_env.yml",
        "samtools --version | head -1 | cut -d' ' -f2"
    ),
}

# ----------------------
rule alignment_star:
    """
    Align reads to the genome using STAR
    """
    version: TOOL_VERSIONS["star"]
    input:
        R1 = rules.fastp.output.R1, 
        R2 = rules.fastp.output.R2,
    output:
        bam = os.path.abspath(FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam"),
        temps= temp(directory(FASTQ_DIR + "/../pipeline_v0/star/{sample}_tmp"))
    conda:
        PIPELINE_DIR + "/envs/comptage_env.yml"
    params:
        star_genome = STAR_GENOME,
        prefix      = lambda wc: os.path.abspath(
            FASTQ_DIR + f"/../pipeline_v0/star/{wc.sample}_"),
        log_start = lambda wc, input, threads: log_start("alignment_star", wc, threads),
        log_end   = LOG_END
    resources:
        tmpdir= temp(OUTPUT_REP + "/star/tmp"),
        single_job=8,
    benchmark:
        "benchmarks/alignement/{sample}.tsv"
    log:
        run_info = "log/star/{sample}/star.log",
        time = "log/star/{sample}/star_time.txt"
    threads: 16
    shell:
        """
        set -euo pipefail
        {params.log_start}
        rm -rf {output.temps}
        STAR --runThreadN {threads} --genomeDir {params.star_genome} \
            --readFilesIn {input.R1} {input.R2} \
            --outSAMtype BAM SortedByCoordinate \
            --chimSegmentMin 20 \
            --twopassMode Basic \
            --outFileNamePrefix {params.prefix} \
            --readFilesCommand zcat \
            --outSAMunmapped Within \
            --outSAMattrRGline ID:4 LB:rnaseq-capture PL:ILLUMINA SM:20 PU:unit1 \
            --outTmpDir {output.temps} \
            --quantMode GeneCounts \
            >> {log.run_info} 2>&1
        {params.log_end}
        """


# ----------------------
rule index_bam:
    """
    Index BAM files using samtools
    """
    version: TOOL_VERSIONS["samtools"]
    input:
        bam = rules.alignment_star.output.bam
    output:
        bai = FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam.bai"
    conda:
        PIPELINE_DIR + "/envs/comptage_env.yml"
    log:
        run_info = "log/samtools/{sample}/index.log",
        time = "log/samtools/{sample}/index_time.txt"
    threads: 4
    params:
        log_start = lambda wc, input, threads: log_start("index_bam", wc, threads),
        log_end   = LOG_END
    shell:
        """
        set -euo pipefail
        {params.log_start}
        samtools index -b {input.bam} >> {log.run_info} 2>&1
        {params.log_end}
        """

