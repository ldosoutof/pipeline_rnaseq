# ----------------------
# 03_comptage_save.smk
# ----------------------

from datetime import datetime

DATE = datetime.now().strftime("%Y-%m-%d")

# ----------------------
# TOOL VERSION DICTIONARY
# ----------------------
TOOL_VERSIONS = {
    "htseq": get_version_from_env(PIPELINE_DIR + "/envs/htseq_env.yml",
                                  "htseq-count --version | head -1 | cut -d' ' -f2"),
    "kallisto": get_version_from_env(PIPELINE_DIR + "/envs/count_env.yml",
                                     "kallisto version"),
    "tx2gene": get_version_from_env(PIPELINE_DIR + "/envs/htseq_env.yml",
                                    "Rscript --version | head -1"),
    "tpm": get_version_from_env(PIPELINE_DIR + "/envs/htseq_env.yml",
                                "Rscript --version | head -1")
}

# ----------------------
# RULES
# ----------------------

rule htseq_gene:
    """
    Count genes using HTSeq (absolute counts)
    """
    version: TOOL_VERSIONS["htseq"]
    input:
        bamg = rules.alignment_star.output.bam,
    output:
        gene = FASTQ_DIR + "/../pipeline_v0/htseq/{sample}/{sample}_gene_counts.txt",
    conda:
        PIPELINE_DIR + "/envs/htseq_env.yml"
    params:
        gtf=GTF,
        log_start = lambda wc, input, threads: log_start("htseq_gene", wc, threads),
        log_end   = LOG_END
    benchmark:
        "benchmarks/htseq/{sample}.tsv"
    log:
        run_info = "log/htseq/{sample}/log.txt",
        time = "log/htseq/{sample}/time.txt"
    threads: 2
    shell:
        """
        set -euo pipefail
        {params.log_start}
        htseq-count -s reverse -r pos -f bam {input.bamg} {params.gtf} -i gene_id \
            > {output.gene}.tmp 2>> {log.run_info}
        grep '^ENSG' {output.gene}.tmp > {output.gene}
        rm {output.gene}.tmp
        {params.log_end}
        """

rule matrix:
    """
    Create gene count matrix aggregating ALL historical runs under PROD_ROOT.
    Includes current run samples + all previous runs for robust OUTRIDER correction.
    """
    input:
        htseq    = expand(rules.htseq_gene.output.gene, sample=SAMPLES),
    output:
        matrix = FASTQ_DIR + "/../pipeline_v0/htseq/matrice.txt"
    params:
        out_dir   = FASTQ_DIR + "/../pipeline_v0/htseq",
        blacklist = outrider_blacklist,
        prod_root = lambda wc: PROD_ROOT,
        scripts   = PIPELINE_DIR,
        log_start = lambda wc, input, threads: log_start("matrix", wc, threads),
        log_end   = LOG_END
    conda:
        PIPELINE_DIR + "/envs/htseq_env.yml"
    log:
        run_info = "log/matrix/log.txt",
        time = "log/matrix/time.txt"
    threads: 8
    shell:
        """
        set -euo pipefail
        {params.log_start}
        BL="{params.blacklist}"
        if [ ! -f "$BL" ]; then
            mkdir -p "$(dirname "$BL")"
            touch "$BL"
        fi
        python {params.scripts}/scripts/create_matrice_by_run.py \
            {params.prod_root} {params.out_dir} "$BL" >> {log.run_info} 2>&1
        {params.log_end}
        """

rule kallistoBed:
    """
    Kallisto quantification (stranded protocol)
    """
    version: TOOL_VERSIONS["kallisto"]
    input:
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
        dir = FASTQ_DIR
    output:
        abund = FASTQ_DIR + "/../pipeline_v0/kallisto_bed/{sample}/abundance.tsv",
        h5 = FASTQ_DIR + "/../pipeline_v0/kallisto_bed/{sample}/abundance.h5"
    conda:
        PIPELINE_DIR + "/envs/count_env.yml"
    params:
        index=KALLISTO_IDX,
        fastq_dir = FASTQ_DIR,
        log_start = lambda wc, input, threads: log_start("kallistoBed", wc, threads),
        log_end   = LOG_END
    benchmark:
        "benchmarks/kallisto_bed/{sample}.tsv"
    log:
        run_info = "log/kallisto_bed/{sample}/log.txt",
        time = "log/kallisto_bed/{sample}/time.txt"
    threads: 8
    shell:
        """
        set -euo pipefail
        {params.log_start}
        kallisto quant -t {threads} -i {params.index} --rf-stranded \
            -o {params.fastq_dir}/../pipeline_v0/kallisto_bed/{wildcards.sample} \
            {input.R1} {input.R2} >> {log.run_info} 2>&1
        {params.log_end}
        """

rule kallisto2gene:
    """
    Convert Ensembl gene IDs to gene names
    """
    version: TOOL_VERSIONS["tx2gene"]
    input:
        h5 = rules.kallistoBed.output.h5
    output:
        gene_level = FASTQ_DIR + "/../pipeline_v0/kallisto_bed/{sample}/abundance_gene_level_counts.tsv"
    conda:
        PIPELINE_DIR + "/envs/htseq_env.yml"
    params:
        scripts = PIPELINE_DIR,
        log_start = lambda wc, input, threads: log_start("kallisto2gene", wc, threads),
        log_end   = LOG_END
    log:
        run_info = "log/tx2gene/{sample}/log.txt",
        time = "log/tx2gene/{sample}/time.txt"
    benchmark:
        "benchmarks/tx2gene/{sample}.tsv"
    threads: 8
    shell:
        """
        set -euo pipefail
        {params.log_start}
        Rscript {params.scripts}/scripts/tx2gene.R {input.h5} >> {log.run_info} 2>&1
        {params.log_end}
        """

rule matrix_tpm:
    """
    Create TPM matrix for Kallisto results
    """
    version: TOOL_VERSIONS["tpm"]
    input:
        counts = expand(rules.kallisto2gene.output.gene_level, sample=SAMPLES),
        dir = FASTQ_DIR
    output:
        gene = FASTQ_DIR + "/../pipeline_v0/kallisto_bed/matrice_gene_tpm.tsv",
        gene_gene = FASTQ_DIR + "/../pipeline_v0/kallisto_bed/matrice_gene_tpm_gene.tsv"
    conda:
        PIPELINE_DIR + "/envs/htseq_env.yml"
    params:
        scripts = PIPELINE_DIR,
        matrix_dir = TPM,
        gtf = GTF,
        log_start = lambda wc, input, threads: log_start("matrix_tpm", wc, threads),
        log_end   = LOG_END
    log:
        run_info = "log/matrixtpm/log.txt",
        time = "log/matrixtpm/time.txt"
    benchmark:
        "benchmarks/matrixtpm/mat.tsv"
    threads: 8
    shell:
        """
        set -euo pipefail
        {params.log_start}
        Rscript {params.scripts}/scripts/create_matrice_tpm_gene_by_run.R \
            FASTQ_DIR/../ {params.matrix_dir}/ {params.gtf} >> {log.run_info} 2>&1
        {params.log_end}
        """
