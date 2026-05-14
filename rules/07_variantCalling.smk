rule mean_chrY_expression:
    """
    Compute mean TPM of chrY lymphocyte genes (sex inference)
    """
    input:
        matrix = FASTQ_DIR + "/../pipeline_v0/kallisto_bed/matrice_gene_tpm.tsv",
        genes  = FASTQ_DIR + "/../chrY_top10_lymphocyte_genes.txt"
    output:
        tsv = FASTQ_DIR + "/../pipeline_v0/qc/chrY_mean_expression.tsv"
    conda:
        PIPELINE_DIR + "/envs/htseq_env.yml"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript {PIPELINE_DIR}/scripts/mean_chrY_expression.R \
            {input.matrix} {input.genes} {output.tsv}
        """
rule vaf_violin_plot_run_females:
    """
    One violin plot per run, chrX VAF, females only
    """
    input:
        vaf_tables = expand(
            FASTQ_DIR + "/../pipeline_v0/deepvariant/{sample}/vaf_table.tsv",
            sample=SAMPLES
        ),
        sex = FASTQ_DIR + "/../pipeline_v0/qc/chrY_mean_expression.tsv"
    output:
        plot = FASTQ_DIR + "/../pipeline_v0/deepvariant/vaf_violin_chrX_females.png"
    conda:
        PIPELINE_DIR + "/envs/python_plot_env.yml"
    shell:
        """
        python {PIPELINE_DIR}/scripts/plot_vaf_violin_by_sample.py \
            {output.plot} \
            {input.vaf_tables}
        """
rule deepvariant_vaf_plot:
    """
    Run DeepVariant (WES mode) using Docker, filter PASS SNPs with DP ≥ 20,
    extract Variant Allele Frequencies (VAF), and plot chrX VAF scatter.
    """
    input:
        ref = GENOME,
        bam = FASTQ_DIR + "/../pipeline_v0/star/{sample}_Aligned.sortedByCoord.out.bam",
        script = PIPELINE_DIR + "/scripts/plot_vaf.py"
    output:
        vcf = FASTQ_DIR + "/../pipeline_v0/deepvariant_chrX/{sample}.vcf.gz",
        vcf_filtered = FASTQ_DIR + "/../pipeline_v0/deepvariant_chrX/snps_pass_dp20_{sample}.vcf.gz",
        vaf_table = FASTQ_DIR + "/../pipeline_v0/deepvariant_chrX/vaf_table_{sample}.tsv",
        plot = FASTQ_DIR + "/../pipeline_v0/deepvariant_chrX/vaf_chrX_{sample}_scatter.png"
    params:
        logs = FASTQ_DIR + "/../pipeline_v0/deepvariant_chrX/logs",
        intermediate = FASTQ_DIR + "/../pipeline_v0/deepvariant_chrX/intermediate_results_dir_all",
        outdir = FASTQ_DIR + "/../pipeline_v0/deepvariant_chrX",
        docker_image = "google/deepvariant:1.8.0"
    threads: 16
    log:
        run_info = "log/deepvariant/{sample}/log.txt",
        time = "log/deepvariant/{sample}/time.txt"
    version:
        "DeepVariant 1.8.0 | bcftools 1.17 | matplotlib 3.8 | pandas 2.2"
    conda:
        PIPELINE_DIR + "/envs/deepvariant_env.yml"
    shell:
        r"""
        echo "start : $(date +"%d-%m-%y %T")" > {log.time}

        # Step 1: DeepVariant (Docker run)
        sudo docker run \
            -v /datawork:/datawork \
            -v /dataref:/dataref \
            -w /datawork \
            {params.docker_image} \
            run_deepvariant \
                --model_type=WES \
                --ref={input.ref} \
                --reads={input.bam} \
                --output_vcf={output.vcf} \
                --num_shards={threads} \
                --intermediate_results_dir={params.intermediate} \
                --make_examples_extra_args="channels='',split_skip_reads=true" \
                --logging_dir={params.logs} \
            2>&1 | tee {params.outdir}/deepvariant_{wildcards.sample}.log

        # Step 2: Filter SNPs (PASS, DP ≥ 20)
        bcftools view -f PASS -v snps -i 'FORMAT/DP>=20' {output.vcf} -Oz -o {output.vcf_filtered}

        # Step 3: Extract Variant Allele Frequencies (VAF)
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD]\n' {output.vcf_filtered} \
            | awk 'BEGIN{{OFS="\t"}} {{split($6,a,","); if($5>0) print $1,$2,$3,$4,$5,a[2]/$5}}' \
            > {output.vaf_table}

        # Step 4: Plot chrX VAF scatter
        python {input.script} {output.vaf_table} {output.plot} >> {log.run_info} 2>&1

        echo "end : $(date +"%d-%m-%y %T")" >> {log.time}
        """
