rule fastqc_report:
    """
    rapport fastqc sur les raw data (fastq)
    """
    input:
        R1 = FASTQ_DIR + "/{sample}_R1.fastq.gz",
        R2 = FASTQ_DIR + "/{sample}_R2.fastq.gz",
        dir = FASTQ_DIR,
        out = OUTPUT_REP 
    output:
        directory(FASTQ_DIR + "/../pipeline_v0/fastqc/{sample}")
    conda:
        "../envs/fastqc_env.yml"
    log: 
        run_info = "logs/fastqc/{sample}/log.txt",
        time = "logs/fastqc/{sample}/time.txt"
    benchmark:
        "benchmarks/fastqc/{sample}.tsv"
    version:
        subprocess.getoutput(
            "fastqc --version | "
            "head -1 | "
            "cut -d' ' -f2"
        )
    threads:8 
    shell:
        'echo "start : $(date +"%d-%m-%y   %T")" > {log.time} && '
        "mkdir -p {output} && "
        "fastqc "
        "-o {output} "
        "-t {threads} "
        "{input.R1} {input.R2} "
        "> {log.run_info} 2>&1 "
        'echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}'
rule fastp:
    """
    Règle filtrant les reads et les trimmant, sur la base de leur qualité
    ET RETIRE UMI
    """
    input:
        R1 = FASTQ_DIR + "/{sample}_R1.fastq.gz",
        R2 = FASTQ_DIR + "/{sample}_R2.fastq.gz",
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        R1=FASTQ_DIR + "/../pipeline_v0/trimmed_fastq/{sample}_R1_trimmed.fastq.gz",
        R2=FASTQ_DIR + "/../pipeline_v0/trimmed_fastq/{sample}_R2_trimmed.fastq.gz",
        html=FASTQ_DIR + "/../pipeline_v0/fastp/{sample}/report.html",
        json=FASTQ_DIR + "/../pipeline_v0/fastp/{sample}/report.json"
    params:
        cut_right_size = 10,
        cut_right_mean_quality = 15,
        qualified_quality_phred = 15,
        max_percent_unqualified = 15,
        minimum_read_length = 30
    threads: 2 # fastp n'en utilise pas plus de 16
    version:    
        subprocess.getoutput(
            "fastp --version | "
            "head -1 | "
            "cut -d' ' -f2"
        )

    conda:
        "../envs/fastp_env.yml"
    log:
        run_info = "logs/fastp/{sample}/fastp.log"
    benchmark:
        "benchmarks/fastp/{sample}.tsv"
    message:
        "--------- fastp : {wildcards.sample} ---------"
    shell:
        "fastp "
            "-i {input.R1} "
            "-I {input.R2} "
            "-o {input.dir}/../pipeline_v0/trimmed_fastq/{wildcards.sample}_R1_trimmed.fastq.gz "
            "-O {input.dir}/../pipeline_v0/trimmed_fastq/{wildcards.sample}_R2_trimmed.fastq.gz "
            "-w {threads} "
            "-l {params.minimum_read_length} "                             # minimum read length
            "--detect_adapter_for_pe "                                     # by default, the adapter sequence auto-detection is enabled for SE data only, turn on this option to enable it for PE data
            "--trim_front1 5 "
            "--trim_front2 5 "
            "--cut_right "                                                 # move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
            "--cut_right_window_size {params.cut_right_size} "
            "--cut_right_mean_quality {params.cut_right_mean_quality} "    # the mean quality requirement option for cut_right, default to cut_mean_quality if not specified 
            "--qualified_quality_phred {params.qualified_quality_phred} "  # the quality value that a base is qualified
            #"-u {params.max_percent_unqualified} "                         # how many percents of bases are allowed to be unqualified (0~100)
            "-h {input.dir}/../pipeline_v0/fastp/{wildcards.sample}/report.html "
            "-j {input.dir}/../pipeline_v0/fastp/{wildcards.sample}/report.json "
        # "--dont_overwrite "
        "> {log.run_info} 2>&1  "
rule fastqc_trim_report:
    """
    Rapport fastqc après trimming de reads
    """
    input:
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        directory(FASTQ_DIR + "/../pipeline_v0/fastqc_trim/{sample}")
    conda:
        "../envs/fastqc_env.yml"
    log:
        run_info = "logs/fastqc_trim/{sample}/log.txt",
        time = "logs/fastqc_trim/{sample}/time.txt"
    version:
        subprocess.getoutput(
            "fastqc --version | "
            "head -1 | "
            "cut -d' ' -f2"
        )
    threads:8
    shell:
        'echo "start : $(date +"%d-%m-%y   %T")" > {log.time} && '
        "mkdir -p {output} && "
        "fastqc "
        "-o {output} "
        "-t {threads} "
        "{input.R1} {input.R2}"
        "> {log.run_info} 2>&1 "
        'echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}'
