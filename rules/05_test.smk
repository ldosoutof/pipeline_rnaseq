RUN_NAME = os.path.basename(OUTPUT_REP)
rule generate_and_run_param_notebook:
    input:
        tpm_folder = "/datawork/genetique/RNASeq/diag/prod/",
        mapping_file = GTF,
        blacklist  = PCA_blacklist
    output:
        executed_nb = f"notebooks/notebook_pca_{RUN_NAME}.ipynb"
    params:
        samples_run  = RUN_NAME,
        excluded_run = "20231122_RUN17_NextSeq_Mid_8RNASEQ",
        keyword      = "MOINS",
        output_html  = f"PCA-{RUN_NAME}.html",
        scripts = PIPELINE_DIR
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    shell:
        """
        python {params.scripts}/scripts/generate_pca_notebook.py \
            --base {input.tpm_folder} \
            --mapping_file {input.mapping_file} \
            --blacklist_file {input.blacklist} \
            --samples_run {params.samples_run} \
            --excluded_run {params.excluded_run} \
            --keyword {params.keyword} \
            --output_html {params.output_html} \
            --notebook_path {output.executed_nb}

        jupyter nbconvert --to notebook --execute --inplace {output.executed_nb}
        """

