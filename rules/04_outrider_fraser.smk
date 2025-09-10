SAMPLES_ID = [s[:7] for s in SAMPLES]
rule outrider:
    """
    Run OUTRIDER gene expression analysis.
    """
    input:
        htseq_matrice = FASTQ_DIR + "/../pipeline_v0/htseq/matrice_gene_counts.tsv",
        matrix = MATRICES,
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        out_file = FASTQ_DIR + "/../pipeline_v0/outrider/outrider_htseq.tab",
        annot = FASTQ_DIR + "/../pipeline_v0/outrider/outrider_htseq_annot.tsv",
        output_files_dir = directory(FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/"),
        output_files = expand(FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab", samples_id=SAMPLES_ID)
    conda:
        "../envs/outrider_env.yml"
    params:
        gtf = GTF,
        panel = PANELAPP,
        pli = PLI,
        hpo = HPO,
        pheno = PHENO,
        date = lambda wildcards: strftime("%Y%m%d", localtime())
    benchmark:
        "benchmarks/outrider/benchmark_outrider.tsv"
    log:
        run_info = "logs/outrider/log.txt",
        time = "logs/outrider/time.txt"
    threads: 1
    resources:
        single_job = 8
    shell:
        '''
        echo "start : $(date +"%d-%m-%y   %T")" > {log.time} &&
        LAST_FILE=$(ls -t {input.matrix} | head -n 1) &&
        Rscript ../scripts/outrider_v3.R {input.matrix}/$LAST_FILE {output.out_file} &&
        python ../scripts/annotation_outrider_v2_rare.py -i {output.out_file} -g {params.gtf} -d {params.panel} -p {params.pli} -o {params.pheno} -m {params.hpo} -f {output.annot} &&
        python ../scripts/outrider_1file.py -i {output.annot} -b {output.output_files_dir}/ &&
        echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
        '''
rule fraser_config:
    """
    Generate FRASER configuration file.
    """
    input:
        bam = expand(FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        dir = FASTQ_DIR,
        out = OUTPUT_REP
    output:
        fraser = FASTQ_DIR + "/../pipeline_v0/fraser/fraser_config.txt",
    conda:
        "../envs/fraser_env.yml"
    params:
        blacklist = blacklist, 
    benchmark:
        "benchmarks/fraser_config/benchmark.tsv"
    log:
        run_info = "logs/fraser_config/log.txt",
        time = "logs/fraser_config/time.txt"
    threads: 1
    resources:
        single_job = 8
    shell:
        '''
        echo "start : $(date +"%d-%m-%y   %T")" > {log.time} &&
        python ../scripts/fraser_config_create.py {output.fraser} PUROMOINS {params.blacklist}  > {log.run_info} 2>&1 &&
        echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
        '''
rule fraser:
    """
    Run FRASER splicing analysis.
    """
    input:
        mat = MATRICES,
        dir = FASTQ_DIR,
        out = OUTPUT_REP,
        config = FASTQ_DIR + "/../pipeline_v0/fraser/fraser_config.txt"
    output:
        fraser = FASTQ_DIR + "/../pipeline_v0/fraser/fraser.tab",
        annot_fraser = FASTQ_DIR + "/../pipeline_v0/fraser/fraser_annot.tsv",
        files_dir = directory(FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/"),
        output_files = expand(FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser.tab", samples_id=SAMPLES_ID)
    conda:
        "../envs/fraser_env.yml"
    params:
        gtf = GTF,
        panel = PANELAPP,
        pli = PLI,
        hpo = HPO,
        pheno = PHENO,
        dir = fraser_count,
        date = lambda wildcards: strftime("%Y%m%d", localtime())
    benchmark:
        "benchmarks/benchmark_fraser.tsv"
    log:
        run_info = "logs/fraser/log_fraser.txt",
        time = "logs/fraser/time_fraser.txt"
    threads: 1
    resources:
        single_job = 8
    shell:
        '''
        echo "start : $(date +"%d-%m-%y   %T")" > {log.time} &&
        LAST_FILE=$(ls -t {input.mat} | head -n 1) &&
        Rscript ../scripts/fraser.R {params.dir} {input.config} {output.fraser}
        python ../scripts/annotation_test_fraser.py -f {output.fraser} -g {params.gtf} -d {params.panel} -p {params.pli} -m {params.hpo} -o {params.pheno} --output {output.annot_fraser} &&
        python ../scripts/fraser_1file.py -i {output.annot_fraser} -b {output.files_dir}/ &&
        echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
        '''
if config.get("run_annot_sake", False):
    rule fraser_annot_rare:
        """
        Run FRASER splicing analysis.
        """
        input:
            fraser=FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser.tab",
            outrider=FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab"
        output:
            fraser_rare=FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser_annot_seqoia.tab"
        params:
            date=lambda wildcards: strftime("%Y%m%d", localtime()),
            run=run_annot_sake
        conda:
            "../envs/fraser_env.yml"
        benchmark:
            "benchmarks/benchmark_fraser_rare_{samples_id}.tsv"
        log:
            run_info = "logs/fraser_rare/{samples_id}_log.txt",
            time = "logs/fraser_rare/{samples_id}_time.txt"
        threads: 1
        resources:
            single_job = 8
        shell:
            '''
            echo "start : $(date +"%d-%m-%y   %T")" > {log.time} &&
            python ../scripts/request_fraser_datalake2_listGenes_v4_multithreads.py --input {input.fraser} --outrider {input.outrider} --output {output.fraser_rare} &&
            echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
            '''
    rule outrider_annot_rare:
        """
        Run FRASER splicing analysis.
        """
        input:
            fraser = FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser.tab",
            outrider =  FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab"
        output:
            outrider_rare = FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider_annot_seqoia.tab"
            #annot_fraser = FASTQ_DIR + "/../pipeline_v0/fraser/fraser_annot.tsv",
            #files_dir = directory(FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/")
        conda:
            "../envs/outrider_env.yml"
        params:
            date = lambda wildcards: strftime("%Y%m%d", localtime())
        benchmark:
            "benchmarks/benchmark_outrider_rare_{samples_id}.tsv"
        log:
            run_info = "logs/outrider_rare/{samples_id}_log.txt",
            time = "logs/outrider_rare/{samples_id}_time.txt"
        threads: 1
        resources:
            single_job = 8
        shell:
            '''
            echo "start : $(date +"%d-%m-%y   %T")" > {log.time} &&
            python ../scripts/request_datalake2_listGenes.py --input {input.outrider} --fraser {input.fraser} --output {output.outrider_rare} &&
            echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
            '''

