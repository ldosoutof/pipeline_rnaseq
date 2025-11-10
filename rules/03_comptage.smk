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
        gtf=GTF
    benchmark:
        "benchmarks/htseq/{sample}.tsv"
    log:
        run_info = "log/htseq/{sample}/log.txt",
        time = "log/htseq/{sample}/time.txt"
    threads: 2
    shell:
        """
        echo "start : $(date +"%d-%m-%y   %T")" > {log.time}
        htseq-count -s reverse -r pos -f bam {input.bamg} {params.gtf} -i gene_id \
            > {output.gene}.tmp 2> {log.run_info}
        sed -n '/^ENSG/,$p' {output.gene}.tmp > {output.gene} 2>&1 | tee -a {log.run_info}
        rm {output.gene}.tmp
        echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
        """

rule matrix:
    """
    Create gene count matrix for all samples, apply blacklist filtering, remove duplicates
    """
    input:
        htseq = expand(rules.htseq_gene.output.gene, sample=SAMPLES),
        blacklist = outrider_blacklist,
        base_dir = FASTQ_DIR + "/../.."
    output:
        matrix = FASTQ_DIR + "/../pipeline_v0/htseq/matrice.txt"
    params:
        out_dir = FASTQ_DIR + "/../pipeline_v0/htseq",
        scripts = PIPELINE_DIR
    conda:
        PIPELINE_DIR + "/envs/htseq_env.yml"
    log:
        run_info = "log/matrix/log.txt",
        time = "log/matrix/time.txt"
    threads: 8
    shell:
        """
        echo "start : $(date +"%d-%m-%y   %T")" > {log.time}
        python {params.scripts}/scripts/create_matrice_by_run.py {input.base_dir} {params.out_dir} {input.blacklist} >> {log.run_info} 2>&1
        echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
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
        index=KALLISTO_IDX
    benchmark:
        "benchmarks/kallisto_bed/{sample}.tsv"
    log:
        run_info = "log/kallisto_bed/{sample}/log.txt",
        time = "log/kallisto_bed/{sample}/time.txt"
    threads: 8
    shell:
        """
        echo "start : $(date +"%d-%m-%y   %T")" > {log.time}
        kallisto quant -t {threads} -i {params.index} --rf-stranded -o {input.dir}/../pipeline_v0/kallisto_bed/{wildcards.sample} \
            {input.R1} {input.R2} > {log.run_info} 2>&1
        echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
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
        scripts = PIPELINE_DIR
    log:
        run_info = "log/tx2gene/{sample}/log.txt",
        time = "log/tx2gene/{sample}/time.txt"
    benchmark:
        "benchmarks/tx2gene/{sample}.tsv"
    threads: 8
    shell:
        """
        echo "start : $(date +"%d-%m-%y   %T")" > {log.time}
        Rscript {params.scripts}/scripts/tx2gene.R {input.h5}
        echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
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
        gtf = GTF
    log:
        run_info = "log/matrixtpm/log.txt",
        time = "log/matrixtpm/time.txt"
    benchmark:
        "benchmarks/matrixtpm/mat.tsv"
    threads: 8
    shell:
        """
        echo "start : $(date +"%d-%m-%y   %T")" > {log.time}
        Rscript {params.scripts}/scripts/create_matrice_tpm_gene_by_run.R {input.dir}/../ {params.matrix_dir}/ {params.gtf}
        echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
        """
