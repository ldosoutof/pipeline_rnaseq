rule alignment_star:
    """
    Alignement des fastq avec l'outil star
    2 pass
    bam genome et transcriptome (sort pour le transcriptome)
    """
    input:
        R1=FASTQ_DIR + "/../pipeline_v0/trimmed_fastq/{sample}_R1_trimmed.fastq.gz",
        R2=FASTQ_DIR + "/../pipeline_v0/trimmed_fastq/{sample}_R2_trimmed.fastq.gz",
        dir = FASTQ_DIR,
        out = OUTPUT_REP        
    output:
        bamg= FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam",
    conda:
        "../envs/comptage_env.yml"
    params:
        star_genome=STAR_GENOME
    resources:
        tmpdir= OUTPUT_REP + "/star/tmp",
        single_job=8,
        mem_gb=60  # Set the memory resource limit to 60GB
    benchmark:
        "benchmarks/alignement/{sample}.tsv"
    log:
        run_info = "logs/star/{sample}/log.txt",
        time = "logs/star/{sample}/time.txt"
    version: # pour récupérer la version de l'outil avec une commande shell 
        subprocess.getoutput(
            "STAR --version | "
            "head -1 | "
            "cut -d' ' -f2"
        )
    threads:12
    shell:
        'echo "start : $(date +"%d-%m-%y   %T")" > {log.time} && '
        "if [ ! -d {input.dir}/../pipeline_v0/star ]; then mkdir {input.dir}/../pipeline_v0/star; fi && "
        "if [ ! -d {input.out}/star ]; then mkdir {input.out}/star; fi && "
        "STAR --runThreadN {threads} --genomeDir {params.star_genome} "
        "--readFilesIn {input.R1} {input.R2} "
        "--outSAMtype BAM SortedByCoordinate "
        "--chimSegmentMin 20 "
        "--twopassMode Basic "
        "--outFileNamePrefix {input.dir}/../pipeline_v0/star/{wildcards.sample}_ "
        "--readFilesCommand zcat "
        "--outSAMunmapped Within "
        "--outSAMattrRGline ID:4 LB:rnaseq-capture PL:ILLUMINA SM:20 PU:unit1 "
        "--outTmpDir {input.out}/star/{wildcards.sample}_tmp "
        "--quantMode GeneCounts "
        "> {log.run_info} 2>&1 && "
        'echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}'

rule index_bam:
    """
    indexation des bam
    """
    input:
        bamg= FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam",
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        baig=FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
    conda:
        "../envs/comptage_env.yml"
    benchmark:
        "benchmarks/alignement/{sample}_bai.tsv"
    version: # pour récupérer la version de l'outil avec une commande shell
        subprocess.getoutput(
            "samtools --version | "
            "head -1 | "
            "cut -d' ' -f2"
        )
    log:
        run_info = "logs/index/{sample}/log.txt",
        time = "logs/index/{sample}/time.txt"
    threads:8
    resources:
        single_job=2
    shell:
        'echo "start : $(date +"%d-%m-%y   %T")" > {log.time} && '

        "samtools index -b {input.bamg} && "
        'echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}'





























































