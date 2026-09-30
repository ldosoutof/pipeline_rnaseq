import os
import subprocess

# ACTIVE_FRASER, ACTIVE_OUTRIDER, ACTIVE_PCA, SAMPLES_ID_BLOOD
# are defined in pipeline.smk with blood-sample filtering.
# Do NOT redefine them here.

# ── Tool versions (evaluated once at parse time) ──────────────────────────────
def _r_pkg_version(pkg):
    """Extrait la version propre depuis la sortie de packageVersion() R."""
    raw = subprocess.getoutput(f'Rscript -e "packageVersion(\'{pkg}\')"')
    for line in raw.splitlines():
        line = line.strip()
        if line.startswith("[1]"):
            return line.replace("[1]", "").strip().strip('"\'  ')
    return raw.strip().split("\n")[0].strip('"\'  ') if raw.strip() else "unknown"

TOOL_VERSIONS = {
    "OUTRIDER": _r_pkg_version("OUTRIDER"),
    "FRASER":   _r_pkg_version("FRASER"),
    "python":   subprocess.getoutput("python --version | cut -d' ' -f2"),
}


# ─────────────────────────────────────────────────────────────────────────────
# OUTRIDER hyper
# ─────────────────────────────────────────────────────────────────────────────
rule outrider_hyper:
    """
    Run OUTRIDER on the featureCounts CDS matrix (hyper pipeline).
    """
    version: TOOL_VERSIONS["OUTRIDER"]
    input:
        matrice = rules.matrix_featurecounts.output.matrix
    output:
        out_file     = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq.tab",
        out_file_all = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq_all.tsv"
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    params:
        scripts         = PIPELINE_DIR,
        yaml_config     = config.get("outrider_yaml", ""),
        current_samples = lambda wc: ",".join(SAMPLES_ID_BLOOD),
        log_start = lambda wc, input, threads: log_start("outrider_hyper", wc, threads),
        log_end   = LOG_END
    benchmark:
        "benchmarks/outrider_hyper/benchmark_outrider.tsv"
    log:
        run_info = "log/outrider_hyper/log.txt",
        time     = "log/outrider_hyper/time.txt"
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
            {input.matrice} \
            {output.out_file} \
            {params.yaml_config} >> {log.run_info} 2>&1
        {params.log_end}
        """


# ─────────────────────────────────────────────────────────────────────────────
# OUTRIDER annotation hyper
# ─────────────────────────────────────────────────────────────────────────────
rule annotation_outrider_hyper:
    """
    Annotate OUTRIDER hyper hits.
    La dépendance à fraser_results_all.tsv est déclarée en input (colonne fraser_hits,
    couplage DAG : Snakemake attend FRASER avant l'annotation).
    """
    version: TOOL_VERSIONS["OUTRIDER"]
    input:
        out_file = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq_all.tsv",
    output:
        annot = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq_annot.tsv",
        output_files_dir = directory(FASTQ_DIR + "/../pipeline_v0/outrider_hyper/filesbysample/"),
        files = touch(expand(
            FASTQ_DIR + "/../pipeline_v0/outrider_hyper/filesbysample/{samples_id}.outrider.tab",
            samples_id=ACTIVE_OUTRIDER))
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    params:
        gtf     = GTF,
        panel   = PANELAPP,
        pli     = PLI,
        hpo     = HPO,
        pheno   = PHENO,
        scripts = PIPELINE_DIR,
        fraser  = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_results_all.tsv",
        log_start = lambda wc, input, threads: log_start("annotation_outrider_hyper", wc, threads),
        log_end   = LOG_END
    log:
        run_info = "log/outrider_hyper/log_annotation.txt",
        time     = "log/outrider_hyper/time_annotation.txt"
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


# ─────────────────────────────────────────────────────────────────────────────
# FRASER config hyper
# ─────────────────────────────────────────────────────────────────────────────
rule fraser_config_hyper:
    """
    Generate the FRASER2 sample config for the hyper pipeline.
    """
    input:
        bam = expand(rules.alignment_star.output.bam, sample=SAMPLES),
        bai = expand(FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
                     sample=SAMPLES),
        dir = FASTQ_DIR
    output:
        fraser = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_config.txt"
    conda:
        PIPELINE_DIR + "/envs/fraser_env.yml"
    params:
        blacklist     = fraser_blacklist,
        fraser        = fraser_count_hyper,
        scripts       = PIPELINE_DIR,
        excluded_runs = config.get("fraser_excluded_runs", ""),
        log_start = lambda wc, input, threads: log_start("fraser_config_hyper", wc, threads),
        log_end   = LOG_END
    benchmark:
        "benchmarks/fraser_config/benchmark_hyper.tsv"
    log:
        run_info = "log/fraser_config_hyper/log.txt",
        time     = "log/fraser_config_hyper/time.txt"
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


# ─────────────────────────────────────────────────────────────────────────────
# FRASER hyper
# ─────────────────────────────────────────────────────────────────────────────
rule fraser_hyper:
    """
    Run FRASER2 for the hyper pipeline (featureCounts CDS BAMs).

    Le script R écrit maintenant directement dans les chemins déclarés par
    output: via l'argument out_dir explicite — supprime le cp post-run et le
    risque de copie de fichiers partiels si R échoue en cours d'exécution.

    Ordre des opérations :
      1. backup_fraser_counts.py  — sauvegarde les comptes RNA précédents
                                    et purge les fraser_results_*.tsv périmés
      2. fraser_newVer.R          — analyse FRASER2, écrit dans out_dir
      3. Snakemake valide output  — les deux .tsv existent et sont non-vides
    """
    version: TOOL_VERSIONS["FRASER"]
    input:
        config = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_config.txt"
    output:
        fraser     = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_results_aberrant.tsv",
        fraser_all = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_results_all.tsv"
    conda:
        PIPELINE_DIR + "/envs/fraser_env.yml"
    params:
        dir         = fraser_count_hyper,
        scripts     = PIPELINE_DIR,
        out_dir     = FASTQ_DIR + "/../pipeline_v0/fraser_hyper",
        yaml_config = config.get("fraser_yaml", ""),
        log_start = lambda wc, input, threads: log_start("fraser_hyper", wc, threads),
        log_end   = LOG_END
    benchmark:
        "benchmarks/fraser_hyper/benchmark_fraser.tsv"
    log:
        run_info = "log/fraser_hyper/log_fraser.txt",
        time     = "log/fraser_hyper/time_fraser.txt"
    threads: 10
    resources:
        single_job = 8
    shell:
        """
        set -euo pipefail
        {params.log_start}
        python {params.scripts}/scripts/backup_fraser_counts.py \
            {params.dir} >> {log.run_info} 2>&1
        FRASER_THREADS={threads} \
        Rscript {params.scripts}/scripts/fraser_newVer.R \
            {params.dir}        \
            {input.config}      \
            {params.out_dir}    \
            {params.yaml_config} >> {log.run_info} 2>&1
        {params.log_end}
        """


# ─────────────────────────────────────────────────────────────────────────────
# FRASER annotation hyper
# ─────────────────────────────────────────────────────────────────────────────
rule annotation_fraser_hyper:
    """
    Annotate FRASER hyper hits.
    La dépendance à outrider_htseq_all.tsv est déclarée en input (flag outrider_hit,
    couplage DAG : Snakemake attend la table OUTRIDER avant l'annotation).
    """
    version: TOOL_VERSIONS["FRASER"]
    input:
        fraser = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_results_all.tsv",
        outrider = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq_all.tsv",
    output:
        annot_fraser = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_annot.tsv",
        files_dir = directory(FASTQ_DIR + "/../pipeline_v0/fraser_hyper/filesbysample/"),
        files = touch(expand(
            FASTQ_DIR + "/../pipeline_v0/fraser_hyper/filesbysample/{samples_id}.fraser.tab",
            samples_id=ACTIVE_FRASER))
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
        log_start = lambda wc, input, threads: log_start("annotation_fraser_hyper", wc, threads),
        log_end   = LOG_END
    log:
        run_info = "log/fraser_hyper/log_annotation_fraser.txt",
        time     = "log/fraser_hyper/time_annotation_fraser.txt"
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


# ─────────────────────────────────────────────────────────────────────────────
# Per-sample output hyper
# ─────────────────────────────────────────────────────────────────────────────
rule rnaseq_per_sample_hyper:
    input:
        fraser   = FASTQ_DIR + "/../pipeline_v0/fraser_hyper/fraser_results_all.tsv",
        outrider = FASTQ_DIR + "/../pipeline_v0/outrider_hyper/outrider_htseq_all.tsv"
    output:
        files_dir = directory(FASTQ_DIR + "/../pipeline_v0/per_sample_hyper/"),
        zip_file  = touch(FASTQ_DIR + "/../pipeline_v0/per_sample_hyper/results.zip"),
        # ── Marqueur du zip filtré (le zip lui-même vit DANS per_sample_hyper/,
        #    donc couvert par le directory() ci-dessus ; ce marqueur hors-dossier
        #    permet aux règles avales d'en dépendre sans chevaucher directory()). ─
        zip_current = touch(FASTQ_DIR + "/../pipeline_v0/.per_sample_hyper_current.done"),
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
        # Sample IDs du run courant (pour le zip filtré)
        current_ids = " ".join(SAMPLES_ID),
        run_tag     = _CURRENT_RUN_TAG,
        # Dossier de travail de la passe 2 HORS du directory() output (évite conflit Snakemake)
        current_dir = FASTQ_DIR + "/../pipeline_v0/.per_sample_hyper_current_tmp",
        samples_txt = FASTQ_DIR + "/../pipeline_v0/.per_sample_hyper_current_tmp/current_samples.txt",
        log_start = lambda wc, input, threads: log_start("rnaseq_per_sample_hyper", wc, threads),
        log_end   = LOG_END
    log:
        run_info = "log/per_sample_hyper/log.txt",
        time     = "log/per_sample_hyper/time.txt"
    threads: 4
    resources:
        single_job = 4
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p {output.files_dir}

        # ── Passe 1 : analyse complète (mode all ou config) ────────────────
        python {params.scripts}/scripts/rnaseq_analysis_per_sample.py \
            --fraser   {input.fraser}    \
            --outrider {input.outrider}  \
            --gtf      {params.gtf}      \
            --output   {output.files_dir} \
            --mode     {params.mode}     \
            --workers  {params.workers}  \
            $([ -n "{params.samples}"    ] && echo "--samples {params.samples}")       \
            $([ -n "{params.gnomad}"     ] && echo "--gnomad {params.gnomad}")         \
            $([ -n "{params.mendeliome}" ] && echo "--mendeliome {params.mendeliome}") \
            $([ -n "{params.pvalue}"     ] && echo "--pvalue {params.pvalue}")         \
            >> {log.run_info} 2>&1
        [ -f "{output.zip_file}" ] || \
            mv $(ls -t {output.files_dir}/*.zip 2>/dev/null | head -1) \
               {output.zip_file} 2>/dev/null || \
            touch "{output.zip_file}"

        # ── Passe 2 : analyse restreinte aux échantillons du run courant ───
        # Écrire la liste des sample IDs du run (un par ligne)
        mkdir -p {params.current_dir}
        echo "{params.current_ids}" | tr ' ' '\\n' | grep -v '^$' > {params.samples_txt}
        echo "[INFO] Zip hyper current - samples du run :" >> {log.run_info}
        cat {params.samples_txt} >> {log.run_info}

        python {params.scripts}/scripts/rnaseq_analysis_per_sample.py \
            --fraser   {input.fraser}    \
            --outrider {input.outrider}  \
            --gtf      {params.gtf}      \
            --output   {params.current_dir} \
            --mode     samples           \
            --samples  {params.samples_txt} \
            --partial-match              \
            --workers  {params.workers}  \
            $([ -n "{params.gnomad}"     ] && echo "--gnomad {params.gnomad}")         \
            $([ -n "{params.mendeliome}" ] && echo "--mendeliome {params.mendeliome}") \
            $([ -n "{params.pvalue}"     ] && echo "--pvalue {params.pvalue}")         \
            >> {log.run_info} 2>&1

        # Récupérer le zip produit par la passe 2 (nommé run_<timestamp>.zip)
        # -> le copier DANS per_sample_hyper/ sous le nom livrable
        _z=$(ls -t {params.current_dir}/*.zip 2>/dev/null | head -1)
        if [ -n "$_z" ]; then
            cp "$_z" "{output.files_dir}/{params.run_tag}_results_hyper_current.zip"
            echo "[OK] zip filtré : {output.files_dir}/{params.run_tag}_results_hyper_current.zip" >> {log.run_info}
        else
            echo "[WARN] aucun zip produit par la passe 2 (samples du run)" >> {log.run_info}
        fi
        # nettoyer le dossier de travail temporaire de la passe 2
        rm -rf {params.current_dir}
        {params.log_end}
        """
