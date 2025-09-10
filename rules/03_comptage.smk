rule htseq_gene:
    """
    comptage gene avec htseq, comptage absolu
    """
    input:
        bamg= FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam",
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        gene = FASTQ_DIR + "/../pipeline_v0/htseq/{sample}/{sample}_gene_counts.txt",
    conda:
        "../envs/htseq_env.yml"
    params:
        gtf=GTF
    benchmark:
        "benchmarks/htseq_g/{sample}.tsv"
    version: # pour récupérer la version de l'outil avec une commande shell
        subprocess.getoutput(
            "htseq-count --version | "
            "head -1 | "
            "cut -d' ' -f2"
        )
    log:
        run_info = "logs/htseq_g/{sample}/log.txt",
        time = "logs/htseq_g/{sample}/time.txt"
    threads:2
    resources:
        single_job=16
    shell:
        'echo "start : $(date +"%d-%m-%y   %T")" > {log.time} && '

        "htseq-count -s reverse -r pos -f bam {input.bamg} {params.gtf} -i gene_id > {input.dir}/../pipeline_v0/htseq/{wildcards.sample}/{wildcards.sample}_gene_counts.txt.tmp 2> {log.run_info} && "
        "sed -n '/^ENSG/,$p' {input.dir}/../pipeline_v0/htseq/{wildcards.sample}/{wildcards.sample}_gene_counts.txt.tmp > {input.dir}/../pipeline_v0/htseq/{wildcards.sample}/{wildcards.sample}_gene_counts.txt 2>&1 | tee -a {log.run_info} && "
        "rm {input.dir}/../pipeline_v0/htseq/{wildcards.sample}/{wildcards.sample}_gene_counts.txt.tmp && "
        'echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}'
rule matrix:
    """
    Create gene count matrix for all PUROMOINS samples in a run,
    apply blacklist filtering, remove POLYA samples, and resolve duplicates.
    """
    input:
        base_dir = FASTQ_DIR + "/../..",
        blacklist = "config/blacklist.txt"
    output:
        matrix = FASTQ_DIR + "/../pipeline_v0/htseq/{run}_matrice.txt"
    params:
        out_dir = FASTQ_DIR + "/../pipeline_v0/htseq"
    conda:
        "../envs/python_env.yml"
    threads: 8
    shell:
        """
        python ../scripts/create_by_run.py {input.base_dir} {params.out_dir} {input.blacklist}
        """
#rule matrix:
#    """
#    creer matrice avec tous les échantillons
#    """
#    input:
#        htseq = expand(FASTQ_DIR + "/../pipeline_v0/htseq/{sample}/{sample}_gene_counts.txt", sample=SAMPLES),
#        dir = FASTQ_DIR,
#        out = OUTPUT_REP
#    output:
#        gene = FASTQ_DIR + "/../pipeline_v0/htseq/matrice_gene_counts.tsv"
#    conda:
#        "../envs/htseq_env.yml"
#    params:
#        matrix=MATRICES
#    benchmark:
#        "benchmarks/matrix/mat.tsv"
#    log:
#        run_info = "logs/matrix/log.txt",
#        time = "logs/matrix/time.txt"
#    threads:8
#    resources:
#        single_job=4
#    shell:
#        'echo "start : $(date +"%d-%m-%y   %T")" > {log.time} && '
#        "Rscript /home/ldosoutoferreira/pipeline/RNASEQ/routine/scripts/create_matrice.R {input.dir}/../ {params.matrix} && "
#        'echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}'
rule kallistoBed:
    """
    kallisto : comptage des genes avec l'option strand (protocole stranded)
    """
    input:
        R1=FASTQ_DIR + "/../pipeline_v0/trimmed_fastq/{sample}_R1_trimmed.fastq.gz",
        R2=FASTQ_DIR + "/../pipeline_v0/trimmed_fastq/{sample}_R2_trimmed.fastq.gz",
        #R1 = rules.fastp.output.R1,
        #R2 = rules.fastp.output.R2,
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        abund=FASTQ_DIR + "/../pipeline_v0/kallisto_bed/{sample}/abundance.tsv",
        h5=FASTQ_DIR + "/../pipeline_v0/kallisto_bed/{sample}/abundance.h5"
    conda:
        "../envs/count_env.yml"
    params:
        index=KALLISTO_IDX
    version: # pour récupérer la version de l'outil avec une commande shell
        subprocess.getoutput(
            "kallisto | "
            "head -1 | "
            "cut -d' ' -f2"
        )
    resources:
        single_job=16
    benchmark:
        "benchmarks/kallisto_bed/{sample}.tsv"
    log:
        run_info = "logs/kallisto_bed/{sample}/log.txt",
        time = "logs/kallisto_bed/{sample}/time.txt"
    threads:8
    shell:
        'echo "start : $(date +"%d-%m-%y   %T")" > {log.time} && '
        "kallisto quant -t {threads} -i {params.index} --rf-stranded -o {input.dir}/../pipeline_v0/kallisto_bed/{wildcards.sample} {input.R1} {input.R2} > {log.run_info} 2>&1 && "
        'echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}'
rule kallisto2gene:
    """
    modifie gene ensembl en gene name en passant par le gtf
    """
    input:
         h5=FASTQ_DIR + "/../pipeline_v0/kallisto_bed/{sample}/abundance.h5"
    output:
         FASTQ_DIR + "/../pipeline_v0/kallisto_bed/{sample}/abundance_gene_level_counts.tsv"
    conda:
        "../envs/htseq_env.yml"
    benchmark:
        "benchmarks/tx2gene/{sample}.tsv"
    log:
        run_info = "logs/tx2gene/{sample}/log.txt",
        time = "logs/tx2gene/{sample}/time.txt"
    threads:8
    resources:
        single_job=4
    shell:
        'echo "start : $(date +"%d-%m-%y   %T")" > {log.time} && '
        "Rscript /home/ldosoutoferreira/pipeline/RNASEQ/routine/scripts/tx2gene.R {input.h5} && "
        'echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}'
rule matrix_tpm:
    """
    creer matrice pour kallisto (utilisation pour cibersortx)
    """
    input:
        htseq = expand(FASTQ_DIR + "/../pipeline_v0/kallisto_bed/{sample}/abundance_gene_level_counts.tsv", sample=SAMPLES),
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        gene = FASTQ_DIR + "/../pipeline_v0/kallisto_bed/marice_gene_tpm.tsv",
    conda:
        "../envs/htseq_env.yml"
    params:
        matrix=TPM,
        gtf=GTF
    benchmark:
        "benchmarks/matrixtpm/mat.tsv"
    log:
        run_info = "logs/matrixtpm/log.txt",
        time = "logs/matrixtpm/time.txt"
    threads:8
    resources:
        single_job=4
    shell:
        'echo "start : $(date +"%d-%m-%y   %T")" > {log.time} && '
        "Rscript /home/ldosoutoferreira/pipeline/RNASEQ/routine/scripts/create_matrice_tpm_gene.R {input.dir}/../ {params.matrix}/ {params.gtf} && "
        'echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}'

