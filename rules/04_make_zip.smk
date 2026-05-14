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
    """
    input:
        fraser   = rules.fraser.output.fraser,
        outrider = rules.outrider.output.out_file,
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
    params:
        output_dir = f"{PROD_ROOT}/{{run}}/pipeline_v0/analysis_output",
        gtf        = GTF,
        gnomad     = GNOMAD,
        mendeliome = MENDELIOME,
        workers    = config.get("ANALYSIS_WORKERS", 4),
        pvalue     = config.get("PVALUE_FILTER", ""),
        scripts    = PIPELINE_DIR,
    conda:
        "outrider"
    log:
        f"{PROD_ROOT}/{{run}}/log/rnaseq_analysis/{{run}}.log",
    threads: config.get("ANALYSIS_WORKERS", 4)
    shell:
        """
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
            &> {log}

        RESULT=$(ls -t {params.output_dir}/run_*.zip 2>/dev/null | head -1)
        if [ -z "$RESULT" ]; then
            echo "ERROR: no output ZIP found in {params.output_dir}" >> {log}
            exit 1
        fi
        mv "$RESULT" {output.result_zip}
        """


# =============================================================================
# HYPER
# =============================================================================

rule make_analysis_zip_hyper:
    """
    Bundle fraser_hyper all results + outrider_hyper into a single ZIP archive.
    Includes both aberrant-only and full FRASER results.
    """
    input:
        fraser         = rules.fraser_hyper.output.fraser,
        fraser_all     = rules.fraser_hyper.output.fraser_all,
        outrider       = rules.outrider_hyper.output.out_file,
        outrider_all   = rules.outrider_hyper.output.out_file_all,
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
    params:
        output_dir = f"{PROD_ROOT}/{{run}}/pipeline_v0/analysis_output_hyper",
        gtf        = GTF,
        gnomad     = GNOMAD,
        mendeliome = MENDELIOME,
        workers    = config.get("ANALYSIS_WORKERS", 4),
        pvalue     = config.get("PVALUE_FILTER_HYPER", config.get("PVALUE_FILTER", "")),
        scripts    = PIPELINE_DIR,
    conda:
        "outrider"
    log:
        f"{PROD_ROOT}/{{run}}/log/rnaseq_analysis_hyper/{{run}}.log",
    threads: config.get("ANALYSIS_WORKERS", 4)
    shell:
        """
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
            &> {log}

        RESULT=$(ls -t {params.output_dir}/run_*.zip 2>/dev/null | head -1)
        if [ -z "$RESULT" ]; then
            echo "ERROR: no output ZIP found in {params.output_dir}" >> {log}
            exit 1
        fi
        mv "$RESULT" {output.result_zip}
        """
