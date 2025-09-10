rule markDuplicate:
    """
    marquage des duplicats
    """
    input:
        aln=FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam",
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        bam=temp(FASTQ_DIR + "/../pipeline_v0/dup/{sample}/{sample}_dup.bam"),
        txt=FASTQ_DIR + "/../pipeline_v0/dup/{sample}/{sample}_dup.txt"
    conda:
        "../envs/mark_env.yml"
    version: # pour récupérer la version de l'outil avec une commande shell
        subprocess.getoutput(
            "kallisto | "
            "head -1 | "
            "cut -d' ' -f2"
        )
    resources:
        single_job=2,
        tmpdir= OUTPUT_REP + "/dup/tmp"
    benchmark:
        "benchmarks/dup/{sample}.tsv"
    log:
        run_info = "logs/dup/{sample}/log.txt",
        time = "logs/dup/{sample}/time.txt"
    threads:8
    shell:
        'echo "start : $(date +"%d-%m-%y   %T")" > {log.time} && '

        "if [ ! -d {input.dir}/../pipeline_v0/dup ]; then mkdir {input.dir}/../pipeline_v0/dup ;fi && "
        "if [ ! -d {input.dir}/../pipeline_v0/dup/{wildcards.sample} ]; then mkdir {input.dir}/../pipeline_v0/dup/{wildcards.sample};fi && "
        "if [ ! -d {input.out}/dup ];then mkdir {input.out}/dup ;fi && "
        "if [ ! -d {input.out}/dup/{wildcards.sample} ];then mkdir {input.out}/dup/{wildcards.sample} ;fi && "
        "picard MarkDuplicates "
        "INPUT={input.aln} "
        "OUTPUT={input.dir}/../pipeline_v0/dup/{wildcards.sample}/{wildcards.sample}_dup.bam "
        "METRICS_FILE={input.dir}/../pipeline_v0/dup/{wildcards.sample}/{wildcards.sample}_dup.txt "
        "VALIDATION_STRINGENCY=LENIENT "
        "REMOVE_DUPLICATES=false "
        "TMP_DIR=tmp >{log.run_info} 2>&1 && "
        'echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}'
SAMPLES_ID = [s[:7] for s in SAMPLES]
rule volcano:
    """
    graph volcano
    """
    input:
        outrider = FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab",
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        FASTQ_DIR + "/../pipeline_v0/outrider/plot/{samples_id}.volcano.png"
    conda:
        "../envs/metrics_env.yml"
    benchmark:
        "benchmarks/volcano/{samples_id}.tsv"
    log:
        run_info = "logs/volcano/{samples_id}/log.txt",
        time = "logs/volcano/{samples_id}/time.txt"
    threads:2
    shell:
        'Rscript ../scripts/volcano2.R {input.outrider} {output} '
rule boxplot:
    """
    Generate boxplots comparing sample counts for the actual run and other samples.
    """
    input:
        outrider = expand(FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab", samples_id=SAMPLES_ID)
    output:
        box=FASTQ_DIR + "/../pipeline_v0/outrider/plot/boxplot.png",
        filt=FASTQ_DIR + "/../pipeline_v0/outrider/plot/boxplot_filt.png"
    conda:
        "../envs/metrics_env.yml"
    benchmark:
        "benchmarks/boxplot/boxplot.tsv"
    log:
        run_info = "logs/boxplot/log.txt",
        time = "logs/boxplot/time.txt"
    params:
        outrider_files = ",".join(expand(FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab", samples_id=SAMPLES_ID)),
        dir = FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/",
        samples_id_str = ",".join(SAMPLES_ID)
    threads: 2
    shell:
        """
        echo '{input.outrider}' &&
        Rscript ../scripts/boxplots4444.R {params.dir} {output.box} {output.filt} {params.samples_id_str}
        """
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
        bed=RSEQ_BED
    resources:
        single_job=2
    conda:
        "../envs/metrics_env.yml"
    benchmark:
        "benchmarks/rseqc/{sample}.tsv"
    log:
        run_info = "logs/rseqc/{sample}/log.txt",
        time = "logs/rseqc/{sample}/time.txt"
    threads:8
    shell:
        'read_distribution.py -i {input.bam} -r {params.bed} > {output} '
rule bam_stats:
    """
    Règle générant des fichiers de statistiques sur le sam telles que la longueur moyenne des inserts, le nombre de reads "on target" etc..
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
#        coverage = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_coverage.tsv",
        hist = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_hist.txt",
        mos = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}.mosdepth.global.dist.txt",
#        di = FASTQ_DIR + "/../pipeline_v0/coverage/{sample}/{sample}_DI.txt"
    params:
        bed = BED,
        padded_bed = PADDED,
        DI_bed = DI_BED
    threads: 12
    conda:
        "../envs/comptage_env.yml"
    resources:
        single_job=16,
        tmpdir= OUTPUT_REP + "/coverage/tmp",
        mem_gb=500  # Set the memory resource limit to 60GB
    version: # pour récupérer la version de l'outil avec une commande shell
        subprocess.getoutput(
            "samtools --version | "
            "head -1 | "
            "cut -d' ' -f2"
        )
    log:
        on_target="log/{sample}/stats_on_target.log",
        padded="log/{sample}/stats_padded.log",
#        coverage="log/{sample}/stats_coverage.log",
        hist="log/{sample}/stats_hist.log",
        insert_size="log/{sample}/stats_insert_size.log",
 #       di="log/{sample}/di.log"
    benchmark:
        "benchmarks/sam_stats/{sample}.tsv"
    message:
        "--------- stats : {wildcards.sample} ---------"
    shell:
        #"if [ ! -d {input.out}/coverage ]; then mkdir {input.out}/coverage; fi && "
        #"if [ ! -d {input.out}/coverage/{wildcards.sample} ]; then mkdir {input.out}/coverage/{wildcards.sample}; fi && "
        "samtools view {input.bam} -b -@ {threads} -F 260 | "

        # avec tee, redirection de stdout vers le stdin des différents bedtools
        "tee >(bedtools intersect -bed -u -abam stdin -b {params.bed} | wc -l > {input.dir}/../pipeline_v0/coverage/{wildcards.sample}/{wildcards.sample}_on_target.txt 2>>{log.on_target}) " # compte le nombre de reads "on target"
        ">(bedtools intersect -bed -u -abam stdin -b {params.padded_bed} | wc -l > {input.dir}/../pipeline_v0/coverage/{wildcards.sample}/{wildcards.sample}_padded.txt 2>>{log.padded}) " # compte le nombre de reads "padded"
#        ">(bedtools coverage -abam stdin -b {params.bed} -d > {output.coverage} 2>>{log.coverage}) " # crée un fichier tsv contenant la depth toutes les positions couvertes par le panel
        ">(bedtools coverage -hist -abam stdin -b {params.bed} | grep all > {input.dir}/../pipeline_v0/coverage/{wildcards.sample}/{wildcards.sample}_hist.txt 2>>{log.hist}) " # crée un histogramme du nombre de reads en fonction de la profondeur
#        ">(bedtools coverage -abam stdin -b {params.DI_bed} -d > {output.di} 2>>{log.di}) "
        "1>/dev/null &&  "
        # calcule les stats du sam (ici, nous sommes intéressés par la longueur des inserts) --> http://www.htslib.org/doc//ssamtools-stats.html
        "samtools stats -F 4 -@ {threads} {input.bam} > {input.dir}/../pipeline_v0/coverage/{wildcards.sample}/{wildcards.sample}_stats.txt 2>>{log.insert_size} && "
        "mosdepth --by {params.DI_bed} --threads 8 --thresholds 1,10,20,30 {input.dir}/../pipeline_v0/coverage/{wildcards.sample}/{wildcards.sample} {input.bam} && "
        "if [ ! -d {input.dir}/../pipeline_v0/coverage ]; then mkdir {input.dir}/../pipeline_v0/coverage; fi &&"
        "if [ ! -d {input.dir}/../pipeline_v0/coverage/{wildcards.sample} ]; then mkdir {input.dir}/../pipeline_v0/coverage/{wildcards.sample}; fi &&"
        "echo 'Le fichier stats.txt a été généré avec samtools, cf http://www.htslib.org/doc/samtools-stats.html' 1>>{log.insert_size}"
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
        "../envs/metrics_env.yml"
    benchmark:
        "benchmarks/multiqc/multiqc.tsv"
    log:
        run_info = "logs/multiqc/log.txt",
        time = "logs/multiqc/time.txt"
    threads: 8
    shell:
        """
        export TMPDIR={input.dir}/tmp && \
        multiqc --force {input.dir}/../pipeline_v0 -o {FASTQ_DIR}/../pipeline_v0/multiqc
        """

