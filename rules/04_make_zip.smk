# ── Pack FRASER + OUTRIDER results into ZIPs for downstream analysis ──────────
#
# Two independent pairs of rules (can run in parallel):
#   • normal : fraser/        + outrider/        → analysis_output/
#   • hyper  : fraser_hyper/  + outrider_hyper/  → analysis_output_hyper/
#
# Downstream consumers:
#   • rnaseq_analysis_per_sample.py  -- called with --fraser / --outrider flags
#   • analyze_from_zip_per_sample.py -- called with --zip (auto-detects files)
# ─────────────────────────────────────────────────────────────────────────────


# =============================================================================
# NORMAL
# =============================================================================

rule make_analysis_zip:
    """
    Bundle fraser.tab + outrider.tab (normal run) into a single ZIP archive.
    Inputs are resolved via PROD_ROOT/{run}/ so each run gets its own data.
    """
    input:
        # Pour le run courant : référencer directement les outputs des règles
        # (Snakemake chaîne le DAG). Pour les runs historiques (re-bundling),
        # construire le chemin via PROD_ROOT/{run}.
        fraser   = lambda wc: (
            rules.fraser.output.fraser
            if wc.run == _CURRENT_RUN_TAG
            else f"{PROD_ROOT}/{wc.run}/pipeline_v0/fraser/fraser.tab"
        ),
        outrider = lambda wc: (
            rules.outrider.output.out_file
            if wc.run == _CURRENT_RUN_TAG
            else f"{PROD_ROOT}/{wc.run}/pipeline_v0/outrider/outrider_htseq.tab"
        ),
    output:
        zip = f"{PROD_ROOT}/{{run}}/pipeline_v0/analysis_input/{{run}}_fraser_outrider.zip",
    log:
        f"{PROD_ROOT}/{{run}}/log/make_analysis_zip/{{run}}.log",
    run:
        import zipfile, os, datetime
        os.makedirs(os.path.dirname(output.zip), exist_ok=True)
        with open(log[0], "w") as lf:
            lf.write(f"start : {datetime.datetime.now():%d-%m-%y   %T}\n")
            with zipfile.ZipFile(output.zip, "w", zipfile.ZIP_DEFLATED) as zf:
                zf.write(input.fraser,   arcname=os.path.basename(input.fraser))
                lf.write(f"  added : {input.fraser}\n")
                zf.write(input.outrider, arcname=os.path.basename(input.outrider))
                lf.write(f"  added : {input.outrider}\n")
            lf.write(f"  zip : {output.zip}\n")
            lf.write(f"end   : {datetime.datetime.now():%d-%m-%y   %T}\n")


rule run_rnaseq_analysis:
    """
    Run analyze_from_zip_per_sample.py on the normal ZIP (--mode all).

    Equivalent CLI:
      python analyze_from_zip_per_sample.py
          --zip        {run}_fraser_outrider.zip
          --mode       all
          --output     {PROD_ROOT}/{run}/pipeline_v0/analysis_output/
          --gtf        references/gencode.v44.annotation.gtf
          --gnomad     references/gnomad_v4_by_gene.tsv
          --mendeliome references/mendeliome_australia.json
          --workers    48

    Required config keys:
        GTF, GNOMAD, MENDELIOME, ANALYSIS_WORKERS (default 4)
    Optional:
        PVALUE_FILTER (e.g. 0.05)
    """
    input:
        zip        = rules.make_analysis_zip.output.zip,
    output:
        result_zip = f"{PROD_ROOT}/{{run}}/pipeline_v0/analysis_output/{{run}}_results.zip",
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    params:
        output_dir = f"{PROD_ROOT}/{{run}}/pipeline_v0/analysis_output",
        gtf        = GTF,
        gnomad     = GNOMAD,
        mendeliome = MENDELIOME,
        workers    = config.get("ANALYSIS_WORKERS", 4),
        pvalue     = config.get("PVALUE_FILTER", ""),
        scripts    = PIPELINE_DIR,
        log_start  = lambda wc, input, threads: log_start("run_rnaseq_analysis", wc, threads),
        log_end    = LOG_END
    log:
        run_info = f"{PROD_ROOT}/{{run}}/log/rnaseq_analysis/{{run}}.log",
        time     = f"{PROD_ROOT}/{{run}}/log/rnaseq_analysis/{{run}}_time.txt",
    threads: config.get("ANALYSIS_WORKERS", 4)
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p {params.output_dir}

        PVAL_FLAG=""
        [ -n "{params.pvalue}" ] && PVAL_FLAG="--pvalue {params.pvalue}"

        GNOMAD_FLAG=""
        [ -f "{params.gnomad}" ] && GNOMAD_FLAG="--gnomad {params.gnomad}"

        MEND_FLAG=""
        [ -f "{params.mendeliome}" ] && MEND_FLAG="--mendeliome {params.mendeliome}"

        python {params.scripts}/scripts/analyze_from_zip_per_sample.py \
            --zip        {input.zip}         \
            --mode       all                 \
            --output     {params.output_dir} \
            --gtf        {params.gtf}        \
            --workers    {params.workers}    \
            $GNOMAD_FLAG                     \
            $MEND_FLAG                       \
            $PVAL_FLAG                       \
            >> {log.run_info} 2>&1

        RESULT=$(ls -t {params.output_dir}/run_*.zip 2>/dev/null | head -1)
        if [ -z "$RESULT" ]; then
            echo "ERROR: no output ZIP found in {params.output_dir}" >> {log.run_info}
            exit 1
        fi
        mv "$RESULT" {output.result_zip}
        {params.log_end}
        """


# =============================================================================
# NORMAL — échantillons du run courant uniquement
# =============================================================================

rule run_rnaseq_analysis_current_run:
    """
    Analyse per-sample restreinte aux échantillons du run courant.

    Utilise --mode samples + un fichier texte listant les sample IDs
    du run (un par ligne), afin d'exclure les échantillons des runs
    précédents inclus dans le ZIP global (make_analysis_zip).

    Sortie : analysis_output_current/{run}_results_current.zip
    """
    input:
        zip = rules.make_analysis_zip.output.zip,
    output:
        result_zip   = f"{PROD_ROOT}/{{run}}/pipeline_v0/analysis_output_current/{{run}}_results_current.zip",
        samples_file = f"{PROD_ROOT}/{{run}}/pipeline_v0/analysis_input/{{run}}_current_samples.txt",
    params:
        output_dir   = f"{PROD_ROOT}/{{run}}/pipeline_v0/analysis_output_current",
        gtf          = GTF,
        gnomad       = GNOMAD,
        mendeliome   = MENDELIOME,
        workers      = config.get("ANALYSIS_WORKERS", 4),
        pvalue       = config.get("PVALUE_FILTER", ""),
        scripts      = PIPELINE_DIR,
        # Sample IDs du run courant séparés par des espaces (pour le shell)
        sample_ids   = lambda wc: " ".join(_RUN_TO_SAMPLES.get(wc.run, [])),
        log_start    = lambda wc, input, threads: log_start("run_rnaseq_analysis_current_run", wc, threads),
        log_end      = LOG_END
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    log:
        run_info = f"{PROD_ROOT}/{{run}}/log/rnaseq_analysis_current_run/{{run}}.log",
        time     = f"{PROD_ROOT}/{{run}}/log/rnaseq_analysis_current_run/{{run}}_time.txt",
    threads: config.get("ANALYSIS_WORKERS", 4)
    shell:
        """
        set -euo pipefail
        {params.log_start}

        mkdir -p $(dirname {output.samples_file})
        mkdir -p {params.output_dir}

        # Écrire le fichier samples (un ID par ligne)
        echo "{params.sample_ids}" | tr ' ' '\n' | grep -v '^$' > {output.samples_file}

        if [ ! -s {output.samples_file} ]; then
            echo "[ERROR] Aucun sample trouvé pour le run {wildcards.run}" >&2
            exit 1
        fi

        echo "[INFO] Samples du run courant :" >> {log.run_info}
        cat {output.samples_file} >> {log.run_info}

        # Construire les arguments optionnels
        PVALUE_ARG=""
        [ -n "{params.pvalue}" ] && PVALUE_ARG="--pvalue {params.pvalue}"

        GNOMAD_ARG=""
        [ -f "{params.gnomad}" ] && GNOMAD_ARG="--gnomad {params.gnomad}"

        MENDELIOME_ARG=""
        [ -f "{params.mendeliome}" ] && MENDELIOME_ARG="--mendeliome {params.mendeliome}"

        python {params.scripts}/scripts/analyze_from_zip_per_sample.py \
            --zip     {input.zip} \
            --mode    samples \
            --samples {output.samples_file} \
            --output  {params.output_dir} \
            --gtf     {params.gtf} \
            --workers {threads} \
            $PVALUE_ARG $GNOMAD_ARG $MENDELIOME_ARG \
            >> {log.run_info} 2>&1

        # Renommer le ZIP produit vers le nom attendu
        RESULT=$(ls -t {params.output_dir}/run_*.zip 2>/dev/null | head -1)
        if [ -z "$RESULT" ]; then
            echo "[ERROR] Aucun ZIP trouvé dans {params.output_dir}" >&2
            exit 1
        fi
        mv "$RESULT" {output.result_zip}

        {params.log_end}
        """


# =============================================================================
# HYPER
# =============================================================================

rule make_analysis_zip_hyper:
    """
    Bundle fraser_hyper all results + outrider_hyper into a single ZIP archive.
    Includes both aberrant-only and full FRASER results.
    Inputs are resolved via PROD_ROOT/{run}/ so each run gets its own data.
    """
    input:
        fraser       = lambda wc: (
            rules.fraser_hyper.output.fraser
            if wc.run == _CURRENT_RUN_TAG
            else f"{PROD_ROOT}/{wc.run}/pipeline_v0/fraser_hyper/fraser_results_aberrant.tsv"
        ),
        fraser_all   = lambda wc: (
            rules.fraser_hyper.output.fraser_all
            if wc.run == _CURRENT_RUN_TAG
            else f"{PROD_ROOT}/{wc.run}/pipeline_v0/fraser_hyper/fraser_results_all.tsv"
        ),
        outrider     = lambda wc: (
            rules.outrider_hyper.output.out_file
            if wc.run == _CURRENT_RUN_TAG
            else f"{PROD_ROOT}/{wc.run}/pipeline_v0/outrider_hyper/outrider_htseq.tab"
        ),
        outrider_all = lambda wc: (
            rules.outrider_hyper.output.out_file_all
            if wc.run == _CURRENT_RUN_TAG
            else f"{PROD_ROOT}/{wc.run}/pipeline_v0/outrider_hyper/outrider_htseq_all.tsv"
        ),
    output:
        zip = f"{PROD_ROOT}/{{run}}/pipeline_v0/analysis_input/{{run}}_fraser_outrider_hyper.zip",
    log:
        f"{PROD_ROOT}/{{run}}/log/make_analysis_zip_hyper/{{run}}.log",
    run:
        import zipfile, os, datetime
        os.makedirs(os.path.dirname(output.zip), exist_ok=True)
        with open(log[0], "w") as lf:
            lf.write(f"start : {datetime.datetime.now():%d-%m-%y   %T}\n")
            with zipfile.ZipFile(output.zip, "w", zipfile.ZIP_DEFLATED) as zf:
                for src in [input.fraser, input.fraser_all,
                            input.outrider, input.outrider_all]:
                    if os.path.exists(src) and os.path.getsize(src) > 0:
                        zf.write(src, arcname=os.path.basename(src))
                        lf.write(f"  added : {src}\n")
                    else:
                        lf.write(f"  skipped (empty): {src}\n")
            lf.write(f"  zip : {output.zip}\n")
            lf.write(f"end   : {datetime.datetime.now():%d-%m-%y   %T}\n")


rule run_rnaseq_analysis_hyper:
    """
    Run analyze_from_zip_per_sample.py on the hyper ZIP (--mode all).

    Equivalent CLI:
      python analyze_from_zip_per_sample.py
          --zip        {run}_fraser_outrider_hyper.zip
          --mode       all
          --output     {PROD_ROOT}/{run}/pipeline_v0/analysis_output_hyper/
          --gtf        references/gencode.v44.annotation.gtf
          --gnomad     references/gnomad_v4_by_gene.tsv
          --mendeliome references/mendeliome_australia.json
          --workers    48

    Uses the same config keys as run_rnaseq_analysis.
    Optional: PVALUE_FILTER_HYPER overrides PVALUE_FILTER for the hyper run only.
    """
    input:
        zip        = rules.make_analysis_zip_hyper.output.zip,
    output:
        result_zip = f"{PROD_ROOT}/{{run}}/pipeline_v0/analysis_output_hyper/{{run}}_results_hyper.zip",
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    params:
        output_dir = f"{PROD_ROOT}/{{run}}/pipeline_v0/analysis_output_hyper",
        gtf        = GTF,
        gnomad     = GNOMAD,
        mendeliome = MENDELIOME,
        workers    = config.get("ANALYSIS_WORKERS", 4),
        pvalue     = config.get("PVALUE_FILTER_HYPER", config.get("PVALUE_FILTER", "")),
        scripts    = PIPELINE_DIR,
        log_start  = lambda wc, input, threads: log_start("run_rnaseq_analysis_hyper", wc, threads),
        log_end    = LOG_END
    log:
        run_info = f"{PROD_ROOT}/{{run}}/log/rnaseq_analysis_hyper/{{run}}.log",
        time     = f"{PROD_ROOT}/{{run}}/log/rnaseq_analysis_hyper/{{run}}_time.txt",
    threads: config.get("ANALYSIS_WORKERS", 4)
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p {params.output_dir}

        PVAL_FLAG=""
        [ -n "{params.pvalue}" ] && PVAL_FLAG="--pvalue {params.pvalue}"

        GNOMAD_FLAG=""
        [ -f "{params.gnomad}" ] && GNOMAD_FLAG="--gnomad {params.gnomad}"

        MEND_FLAG=""
        [ -f "{params.mendeliome}" ] && MEND_FLAG="--mendeliome {params.mendeliome}"

        python {params.scripts}/scripts/analyze_from_zip_per_sample.py \
            --zip        {input.zip}         \
            --mode       all                 \
            --output     {params.output_dir} \
            --gtf        {params.gtf}        \
            --workers    {params.workers}    \
            $GNOMAD_FLAG                     \
            $MEND_FLAG                       \
            $PVAL_FLAG                       \
            >> {log.run_info} 2>&1

        RESULT=$(ls -t {params.output_dir}/run_*.zip 2>/dev/null | head -1)
        if [ -z "$RESULT" ]; then
            echo "ERROR: no output ZIP found in {params.output_dir}" >> {log.run_info}
            exit 1
        fi
        mv "$RESULT" {output.result_zip}
        {params.log_end}
        """


# =============================================================================
# CLEANUP FINAL
# =============================================================================
# Supprime les sorties intermédiaires/redondantes à la toute fin du pipeline,
# une fois que TOUTES les sorties à conserver sont produites.
#
# Supprimé :
#   - analysis_input/          (zips techniques FRASER+OUTRIDER bruts)
#   - analysis_output/         (analyse normale globale, toute la cohorte)
#   - analysis_output_current/ (analyse normale filtrée au run)
#
# Conservé :
#   - per_sample_hyper/        (dont {run}_results_hyper_current.zip, le livrable)
#   - analysis_output_hyper/   (analyse hyper globale)
#
# Dépend des sorties à conserver pour garantir l'ordre (cleanup en dernier).
# =============================================================================

rule cleanup_final:
    """
    Nettoyage final : retire analysis_input / analysis_output /
    analysis_output_current après que les livrables hyper soient produits.
    Produit un marqueur .cleanup_done pour que Snakemake trace l'exécution.
    """
    input:
        # Dépendances = tout ce qu'on CONSERVE (garantit que le cleanup est le dernier)
        hyper_current = rules.rnaseq_per_sample_hyper.output.zip_current,
        hyper_global  = f"{PROD_ROOT}/{_CURRENT_RUN_TAG}/pipeline_v0/analysis_output_hyper/{_CURRENT_RUN_TAG}_results_hyper.zip",
    output:
        marker = f"{PROD_ROOT}/{_CURRENT_RUN_TAG}/pipeline_v0/.cleanup_done",
    params:
        base = f"{PROD_ROOT}/{_CURRENT_RUN_TAG}/pipeline_v0",
    log:
        run_info = f"{PROD_ROOT}/{_CURRENT_RUN_TAG}/log/cleanup_final/log.txt",
    run:
        import shutil, os, datetime
        os.makedirs(os.path.dirname(log.run_info), exist_ok=True)
        to_remove = ["analysis_input", "analysis_output", "analysis_output_current"]
        with open(log.run_info, "w") as lf:
            lf.write(f"start : {datetime.datetime.now():%d-%m-%y %T}\n")
            for d in to_remove:
                path = os.path.join(params.base, d)
                if os.path.isdir(path):
                    shutil.rmtree(path)
                    lf.write(f"  removed : {path}\n")
                else:
                    lf.write(f"  absent  : {path}\n")
            lf.write(f"end   : {datetime.datetime.now():%d-%m-%y %T}\n")
        # marqueur de fin
        with open(output.marker, "w") as mf:
            mf.write(f"cleanup done {datetime.datetime.now():%d-%m-%y %T}\n")
