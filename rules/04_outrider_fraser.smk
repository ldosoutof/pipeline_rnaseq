import os
import subprocess

# ACTIVE_FRASER, ACTIVE_OUTRIDER, ACTIVE_PCA, SAMPLES_ID_BLOOD
# are defined in pipeline.smk with blood-sample filtering.
# Do NOT redefine them here.

# ----------------------
# TOOL VERSIONS
# ----------------------
def _r_pkg_version(pkg):
    """Extrait la version propre depuis la sortie de packageVersion() R."""
    raw = subprocess.getoutput(f'Rscript -e "packageVersion(\'{pkg}\')"' )
    for line in raw.splitlines():
        line = line.strip()
        if line.startswith("[1]"):
            return line.replace("[1]", "").strip().strip('"\' ')
    return raw.strip().split("\n")[0].strip('"\' ') if raw.strip() else "unknown"

TOOL_VERSIONS = {
    "OUTRIDER": _r_pkg_version("OUTRIDER"),
    "FRASER":   _r_pkg_version("FRASER"),
    "python":   subprocess.getoutput("python --version | cut -d' ' -f2"),
}


# --------------------------
# OUTRIDER main rule
# --------------------------
rule outrider:
    version: TOOL_VERSIONS["OUTRIDER"]
    input:
        htseq_matrice = FASTQ_DIR + "/../pipeline_v0/htseq/matrice.txt"
    output:
        out_file     = FASTQ_DIR + "/../pipeline_v0/outrider/outrider_htseq.tab",
        out_file_all = FASTQ_DIR + "/../pipeline_v0/outrider/outrider_htseq_all.tsv"
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    params:
        scripts      = PIPELINE_DIR,
        yaml_config  = config.get("outrider_yaml", ""),
        current_samples = lambda wc: ",".join(SAMPLES_ID_BLOOD),
        log_start    = lambda wc, input, threads: log_start("outrider", wc, threads),
        log_end      = LOG_END
    benchmark:
        "benchmarks/outrider/benchmark_outrider.tsv"
    log:
        run_info = "log/outrider/log.txt",
        time     = "log/outrider/time.txt"
    threads: 8
    resources:
        single_job = 8
    shell:
        """
        set -euo pipefail
        {params.log_start}
        OUTRIDER_THREADS={threads} \
        OUTRIDER_CURRENT_SAMPLES="{params.current_samples}" \
        Rscript {params.scripts}/scripts/outrider_newVersion.R \
            {input.htseq_matrice} \
            {output.out_file} \
            {params.yaml_config} >> {log.run_info} 2>&1
        {params.log_end}
        """

# --------------------------
# OUTRIDER annotation rule
# --------------------------
rule annotation_outrider:
    """
    Annotate OUTRIDER hits with GTF, PanelApp DI, PLI, HPO.
    Enrichit chaque hit avec le nombre de jonctions FRASER aberrantes
    co-localisées (colonne fraser_hits) — signal de convergence expression+épissage.

    La dépendance à fraser.tab est intentionnelle mais non bloquante :
    si fraser.tab est absent, le script produit fraser_hits=0 pour tous les hits.
    Le chemin est passé via params (pas input) pour ne pas bloquer le DAG
    si FRASER échoue indépendamment.
    """
    version: TOOL_VERSIONS["OUTRIDER"]
    input:
        out_file = FASTQ_DIR + "/../pipeline_v0/outrider/outrider_htseq_all.tsv",
    output:
        annot = FASTQ_DIR + "/../pipeline_v0/outrider/outrider_htseq_annot.tsv",
        output_files_dir = directory(FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/"),
        files = touch(expand(FASTQ_DIR + "/../pipeline_v0/outrider/filesbysample/{samples_id}.outrider.tab", samples_id=ACTIVE_OUTRIDER))
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    params:
        gtf     = GTF,
        panel   = PANELAPP,
        pli     = PLI,
        hpo     = HPO,
        pheno   = PHENO,
        scripts = PIPELINE_DIR,
        # Chemin vers fraser.tab passé en params (non bloquant dans le DAG).
        # Le script charge ce fichier seulement s'il existe (os.path.exists).
        fraser  = FASTQ_DIR + "/../pipeline_v0/fraser/fraser.tab",
        log_start = lambda wc, input, threads: log_start("annotation_outrider", wc, threads),
        log_end   = LOG_END
    log:
        run_info = "log/outrider/log_annotation.txt",
        time = "log/outrider/time_annotation.txt"
    threads: 1
    resources:
        single_job = 4
    shell:
        """
        set -euo pipefail
        {params.log_start}
        python {params.scripts}/scripts/annotation_outrider_hits_bis2.py \
            -i {input.out_file} -g {params.gtf} -d {params.panel} \
            -p {params.pli} -o {params.pheno} -m {params.hpo} \
            -f {output.annot} --fraser {params.fraser} >> {log.run_info} 2>&1
        python {params.scripts}/scripts/outrider_1file.py \
            -i {output.annot} -b {output.output_files_dir}/ >> {log.run_info} 2>&1
        for f in {output.files}; do [ -f "$f" ] || touch "$f"; done
        {params.log_end}
        """
# --------------------------
rule fraser_config:
    input:
        bam = expand(rules.alignment_star.output.bam, sample=SAMPLES),
        bai = expand(FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        dir = FASTQ_DIR
    output:
        fraser = FASTQ_DIR + "/../pipeline_v0/fraser/fraser_config.txt"
    conda:
        PIPELINE_DIR + "/envs/fraser_env.yml"
    params:
        blacklist      = fraser_blacklist,
        fraser         = fraser_count,
        scripts        = PIPELINE_DIR,
        excluded_runs  = config.get("fraser_excluded_runs", ""),
        log_start = lambda wc, input, threads: log_start("fraser_config", wc, threads),
        log_end   = LOG_END
    benchmark:
        "benchmarks/fraser_config/benchmark.tsv"
    log:
        run_info = "log/fraser_config/log.txt",
        time = "log/fraser_config/time.txt"
    threads: 1
    resources:
        single_job = 8
    shell:
        """
        set -euo pipefail
        {params.log_start}
        python {params.scripts}/scripts/fraser_config_create.py \
            --output        {output.fraser}          \
            --pattern       MOINS                    \
            --root_dir      {input.dir}/../..        \
            --fraser_dir    {params.fraser}          \
            --blacklist     {params.blacklist}       \
            --excluded_runs "{params.excluded_runs}" \
            >> {log.run_info} 2>&1
        {params.log_end}
        """

# --------------------------
# FRASER main rule (core)
# --------------------------
rule fraser:
    version: TOOL_VERSIONS["FRASER"]
    input:
        config = FASTQ_DIR + "/../pipeline_v0/fraser/fraser_config.txt"
    output:
        fraser = FASTQ_DIR + "/../pipeline_v0/fraser/fraser.tab"
    conda:
        PIPELINE_DIR + "/envs/fraser_env.yml"
    params:
        dir         = fraser_count,
        scripts     = PIPELINE_DIR,
        yaml_config = config.get("fraser_yaml", ""),
        log_start = lambda wc, input, threads: log_start("fraser", wc, threads),
        log_end   = LOG_END
    benchmark:
        "benchmarks/fraser/benchmark_fraser.tsv"
    log:
        run_info = "log/fraser/log_fraser.txt",
        time = "log/fraser/time_fraser.txt"
    threads: 10
    resources:
        single_job = 8
    shell:
        """
        set -euo pipefail
        {params.log_start}
        python {params.scripts}/scripts/backup_fraser_counts.py {params.dir} >> {log.run_info} 2>&1
        FRASER_THREADS={threads} \
        Rscript {params.scripts}/scripts/fraser.R \
            {params.dir} {input.config} {output.fraser} \
            {params.yaml_config} >> {log.run_info} 2>&1
        {params.log_end}
        """

# --------------------------
# FRASER annotation rule
# --------------------------
rule annotation_fraser:
    """
    Annotate FRASER hits with GTF, PanelApp DI, PLI, HPO.
    Enrichit chaque hit avec un flag outrider_hit indiquant si le même gène
    est aberrant en expression chez le même patient — signal de convergence.

    La dépendance à outrider_htseq.tab est déclarée en input (couplage DAG) :
    Snakemake garantit qu'OUTRIDER est produit avant cette annotation.
    """
    version: TOOL_VERSIONS["FRASER"]
    input:
        fraser = FASTQ_DIR + "/../pipeline_v0/fraser/fraser.tab",
        outrider = FASTQ_DIR + "/../pipeline_v0/outrider/outrider_htseq.tab",
    output:
        annot_fraser = FASTQ_DIR + "/../pipeline_v0/fraser/fraser_annot.tsv",
        files_dir = directory(FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/"),
        files = touch(expand(FASTQ_DIR + "/../pipeline_v0/fraser/filesbysample/{samples_id}.fraser.tab", samples_id=ACTIVE_FRASER))
    conda:
        PIPELINE_DIR + "/envs/fraser_env.yml"
    params:
        blacklist = fraser_blacklist,
        gtf     = GTF,
        panel   = PANELAPP,
        pli     = PLI,
        hpo     = HPO,
        pheno   = PHENO,
        scripts = PIPELINE_DIR,
        magnis  = config.get("magnis_file", ""),
        # Chemin vers outrider_htseq.tab déclaré en input (dépendance DAG).
        log_start = lambda wc, input, threads: log_start("annotation_fraser", wc, threads),
        log_end   = LOG_END
    log:
        run_info = "log/fraser/log_annotation_fraser.txt",
        time = "log/fraser/time_annotation_fraser.txt"
    threads: 1
    resources:
        single_job = 4
    shell:
        """
        set -euo pipefail
        {params.log_start}
        python {params.scripts}/scripts/annotation_fraser_test4.py \
            -f {input.fraser} -g {params.gtf} -d {params.panel} \
            -p {params.pli} -m {params.hpo} -o {params.pheno} \
            --outrider {input.outrider} --output {output.annot_fraser} \
            --magnis "{params.magnis}" >> {log.run_info} 2>&1
        python {params.scripts}/scripts/fraser_1file.py \
            -i {output.annot_fraser} -b {output.files_dir}/ >> {log.run_info} 2>&1
        for f in {output.files}; do [ -f "$f" ] || touch "$f"; done
        {params.log_end}
        """

# --------------------------
# Per-sample output rule
# --------------------------
rule rnaseq_per_sample:
    input:
        fraser   = FASTQ_DIR + "/../pipeline_v0/fraser/fraser.tab",
        outrider = FASTQ_DIR + "/../pipeline_v0/outrider/outrider_htseq_all.tsv"
    output:
        files_dir = directory(FASTQ_DIR + "/../pipeline_v0/per_sample/"),
        zip_file  = touch(FASTQ_DIR + "/../pipeline_v0/per_sample/results.zip")
    conda:
        PIPELINE_DIR + "/envs/python_plot_env.yml"
    params:
        gtf        = GTF,
        gnomad     = config.get("gnomad", ""),
        mendeliome = config.get("mendeliome", ""),
        samples    = config.get("samples_file", ""),
        mode       = config.get("per_sample_mode", "all"),
        pvalue     = config.get("per_sample_pvalue", ""),
        workers    = config.get("per_sample_workers", 4),
        scripts    = PIPELINE_DIR,
        log_start  = lambda wc, input, threads: log_start("rnaseq_per_sample", wc, threads),
        log_end    = LOG_END
    log:
        run_info = "log/per_sample/log.txt",
        time     = "log/per_sample/time.txt"
    threads: 4
    resources:
        single_job = 4
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p {output.files_dir}
        python {params.scripts}/scripts/rnaseq_analysis_per_sample.py \
            --fraser   {input.fraser} \
            --outrider {input.outrider} \
            --gtf      {params.gtf} \
            --output   {output.files_dir} \
            --mode     {params.mode} \
            --workers  {params.workers} \
            $([ -n "{params.samples}"    ] && echo "--samples {params.samples}")       \
            $([ -n "{params.gnomad}"     ] && echo "--gnomad {params.gnomad}")         \
            $([ -n "{params.mendeliome}" ] && echo "--mendeliome {params.mendeliome}") \
            $([ -n "{params.pvalue}"     ] && echo "--pvalue {params.pvalue}")         \
            >> {log.run_info} 2>&1
        [ -f "{output.zip_file}" ] || \
            mv $(ls -t {output.files_dir}/*.zip 2>/dev/null | head -1) \
               {output.zip_file} 2>/dev/null || \
            touch "{output.zip_file}"
        {params.log_end}
        """

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#






















