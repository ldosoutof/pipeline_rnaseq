SAMPLES_ID = [s[:7] for s in SAMPLES]

def active_samples(blacklist_file, samples=SAMPLES_ID):
    """
    Return list of sample IDs excluding any in the blacklist file.
    """
    if os.path.exists(blacklist_file):
        with open(blacklist_file) as f:
            excluded = set(f.read().split())
    else:
        excluded = set()
    return [s for s in samples if s not in excluded]


# Active sample lists
ACTIVE_FRASER   = active_samples(fraser_blacklist)
ACTIVE_OUTRIDER = active_samples(outrider_blacklist)
ACTIVE_PCA      = active_samples(PCA_blacklist)


# --------------------------
# OUTRIDER main rule
# --------------------------
rule outrider:
    """
    Run OUTRIDER gene expression analysis.
    """
    input:
        htseq_matrice = FASTQ_DIR + "/../pipeline_v0/htseq/matrice.txt"
    output:
        out_file = FASTQ_DIR + "/../pipeline_v0/outrider/outrider_htseq.tab",
        annot = FASTQ_DIR + "/../pipeline_v0/outrider/outrider_htseq_annot.tsv",
        output_files_dir = directory(FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/"),
        files = expand(
        FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab",
        samples_id=ACTIVE_OUTRIDER
    )
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    params:
        gtf = GTF,
        panel = PANELAPP,
        pli = PLI,
        hpo = HPO,
        pheno = PHENO,
        scripts = PIPELINE_DIR
    benchmark:
        "benchmarks/outrider/benchmark_outrider.tsv"
    log:
        run_info = "log/outrider/log.txt",
        time = "log/outrider/time.txt"
    threads: 1
    resources:
        single_job = 8
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        Rscript {params.scripts}/scripts/outrider_new.R {input.htseq_matrice} {output.out_file} &&
        python {params.scripts}/scripts/annotation_outrider_v2_rare.py \
            -i {output.out_file} -g {params.gtf} -d {params.panel} \
            -p {params.pli} -o {params.pheno} -m {params.hpo} -f {output.annot} &&
        python {params.scripts}/scripts/outrider_1file.py \
            -i {output.annot} -b {output.output_files_dir}/ &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''


# --------------------------
# FRASER config generation
# --------------------------
rule fraser_config:
    """
    Generate FRASER configuration file.
    """
    input:
        bam = expand(FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        bai = expand(FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        dir = FASTQ_DIR
    output:
        fraser = FASTQ_DIR + "/../pipeline_v0/fraser/fraser_config.txt"
    conda:
        PIPELINE_DIR + "/envs/fraser_env.yml"
    params:
        blacklist = fraser_blacklist,
        fraser = "/datawork/genetique/RNASeq/diag/prod/fraser_2_2",
        scripts = PIPELINE_DIR
    benchmark:
        "benchmarks/fraser_config/benchmark.tsv"
    log:
        run_info = "log/fraser_config/log.txt",
        time = "log/fraser_config/time.txt"
    threads: 1
    resources:
        single_job = 8
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        python {params.scripts}/scripts/fraser_config_create.py {output.fraser} MOINS {input.dir}/../.. {params.fraser} {params.blacklist} > {log.run_info} 2>&1 &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''


# --------------------------
# FRASER main rule
# --------------------------
rule fraser:
    """
    Run FRASER splicing analysis.
    """
    input:
        mat = MATRICES,
        config = FASTQ_DIR + "/../pipeline_v0/fraser/fraser_config.txt"
    output:
        fraser = FASTQ_DIR + "/../pipeline_v0/fraser/fraser.tab",
        annot_fraser = FASTQ_DIR + "/../pipeline_v0/fraser/fraser_annot.tsv",
        files_dir = directory(FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/"),
        files = expand(
        FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser.tab",
        samples_id=ACTIVE_FRASER
    )
    conda:
        PIPELINE_DIR + "/envs/fraser_env.yml"
    params:
        blacklist = fraser_blacklist,
        gtf = GTF,
        panel = PANELAPP,
        pli = PLI,
        hpo = HPO,
        pheno = PHENO,
        dir = fraser_count,
        scripts = PIPELINE_DIR
    benchmark:
        "benchmarks/fraser/benchmark_fraser.tsv"
    log:
        run_info = "log/fraser/log_fraser.txt",
        time = "log/fraser/time_fraser.txt"
    threads: 1
    resources:
        single_job = 8
    shell:
        '''
        echo "start : $(date +"%d-%m-%y %T")" > {log.time} &&
        python {params.scripts}/scripts/backup_fraser_counts.py {params.dir} &&
        Rscript {params.scripts}/scripts/fraser.R {params.dir} {input.config} {output.fraser} &&
        python {params.scripts}/scripts/annotation_test_fraser.py \
            -f {output.fraser} -g {params.gtf} -d {params.panel} \
            -p {params.pli} -m {params.hpo} -o {params.pheno} --output {output.annot_fraser} &&
        python {params.scripts}/scripts/fraser_1file.py -i {output.annot_fraser} -b {output.files_dir}/ &&
        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        '''


