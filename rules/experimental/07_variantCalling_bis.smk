# ── MODULE OPTIONNEL — activable via run_variant_calling: true ────────────────
# Analyse du biais d'inactivation du chromosome X.
#
# Règles :
#   mean_chrY_expression      — inférence du sexe par expression chrY (Kallisto TPM)
#   deepvariant_chrX_WES      — appel de variants sur chrX (DeepVariant WES)
#   bgzip_and_index_vcf       — compression + indexation du VCF
#   annotate_chrX_population  — annotation rsID (1000G) et AF (gnomAD)
#   vaf_table_chrX            — extraction des SNPs hétérozygotes PASS
#   vaf_violin_plot_run_females — violin plot du biais d'inactivation X
#
# Activation :
#   Ajouter dans le config.yml du run :
#       run_variant_calling: true
#   Le include: et les cibles dans rule all sont gérés conditionnellement
#   par pipeline.smk.
#
# Prérequis système :
#   Singularity/Apptainer installé sur le nœud d'exécution.
#   Lancer Snakemake avec --use-singularity.
#   Configurer singularity_args dans pipeline.smk si les montages par défaut
#   ne couvrent pas GENOME, BAM et BED (voir config docker_mounts).
# ─────────────────────────────────────────────────────────────────────────────

rule mean_chrY_expression:
    """
    Compute mean TPM of chrY lymphocyte genes (sex inference)
    """
    input:
        matrix = FASTQ_DIR + "/../pipeline_v0/kallisto_bed/matrice_gene_tpm_gene.tsv",
        genes  = FASTQ_DIR + "/../../chrY_top10_lymphocyte_genes.txt"
    output:
        tsv = FASTQ_DIR + "/../pipeline_v0/qc/chrY_mean_expression.tsv"
    conda:
        PIPELINE_DIR + "/envs/htseq_env.yml"
    params:
        scripts   = PIPELINE_DIR,
        log_start = lambda wc, input, threads: log_start("mean_chrY_expression", wc, threads),
        log_end   = LOG_END
    log:
        run_info = "log/mean_chrY_expression/log.txt",
        time     = "log/mean_chrY_expression/time.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p $(dirname {output.tsv})
        Rscript {params.scripts}/scripts/mean_chrY_expression.R \
            {input.matrix} {input.genes} {output.tsv} >> {log.run_info} 2>&1
        {params.log_end}
        """
rule deepvariant_chrX_WES:
    """
    Run DeepVariant on chrX using WES model and capture BED.
    Requires Singularity/Apptainer — launch Snakemake with --use-singularity.
    """
    input:
        bam = rules.alignment_star.output.bam,
        bai = rules.index_bam.output.bai
    output:
        vcf = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/{sample}_chrX.vcf.gz"
    container:
        config.get("deepvariant_image", "docker://google/deepvariant:1.8.0")
    params:
        ref  = GENOME,
        bed  = config.get("chrx_bed", ""),
        tmp  = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/intermediate_results_dir",
        logs = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/logs",
        log_start = lambda wc, input, threads: log_start("deepvariant_chrX_WES", wc, threads),
        log_end   = LOG_END
    threads: 16
    log:
        run_info = "log/deepvariant_WES/{sample}/log.txt",
        time     = "log/deepvariant_WES/{sample}/time.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p $(dirname {output.vcf}) {params.tmp} {params.logs}
        run_deepvariant \
            --model_type=WES \
            --ref={params.ref} \
            --reads={input.bam} \
            --output_vcf={output.vcf} \
            --num_shards={threads} \
            --regions={params.bed} \
            --intermediate_results_dir={params.tmp} \
            --logging_dir={params.logs} \
            --make_examples_extra_args="channels='',split_skip_reads=true" \
            >> {log.run_info} 2>&1
        {params.log_end}
        """
rule bgzip_and_index_vcf:
    """
    Ensure DeepVariant VCF is bgzip-compressed and indexed
    """
    input:
        vcf = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/{sample}_chrX.vcf.gz"
    output:
        vcf = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/{sample}_chrX.bgz.vcf.gz",
        tbi = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/{sample}_chrX.bgz.vcf.gz.tbi"
    conda:
        PIPELINE_DIR + "/envs/bcftools_env.yml"
    params:
        log_start = lambda wc, input, threads: log_start("bgzip_and_index_vcf", wc, threads),
        log_end   = LOG_END
    log:
        run_info = "log/bgzip_vcf/{sample}/log.txt",
        time     = "log/bgzip_vcf/{sample}/time.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        gunzip -c {input.vcf} | bgzip -c > {output.vcf}
        tabix -f -p vcf {output.vcf}
        {params.log_end}
        """
rule annotate_chrX_population:
    """
    Annotate chrX VCF with 1000G rsID and gnomAD AF
    """
    input:
        vcf = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/{sample}_chrX.bgz.vcf.gz"
    output:
        vcf = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/{sample}_chrX_annot.vcf.gz",
        tbi = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/{sample}_chrX_annot.vcf.gz.tbi"
    params:
        tmp    = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/tmp_annot.vcf.gz",
        g1000  = config.get("g1000_chrx", ""),
        gnomad = config.get("gnomad_chrx", ""),
        log_start = lambda wc, input, threads: log_start("annotate_chrX_population", wc, threads),
        log_end   = LOG_END
    conda:
        PIPELINE_DIR + "/envs/bcftools_env.yml"
    log:
        run_info = "log/deepvariant/annotate/{sample}/log.txt",
        time     = "log/deepvariant/annotate/{sample}/time.txt"
    shell:
        r"""
        set -euo pipefail
        {params.log_start}

        echo "[INFO] Index annotation VCFs if needed"
        bcftools index -f -t {params.g1000}
        bcftools index -f -t {params.gnomad}

        echo "[INFO] Step 1: add rsID from 1000G"
        bcftools annotate \
          -a {params.g1000} \
          -c ID \
          {input.vcf} \
        | bgzip -c > {params.tmp}

        bcftools index -t {params.tmp}

        echo "[INFO] Step 2: add gnomAD AF"
        bcftools annotate \
          -a {params.gnomad} \
          -c INFO/AF,INFO/AF_non_v2_XX \
          -h <(cat <<'EOF'
##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD allele frequency">
##INFO=<ID=gnomAD_AF_XX,Number=A,Type=Float,Description="gnomAD allele frequency in XX samples">
EOF
          ) \
          {params.tmp} \
        | bgzip -c > {output.vcf}

        bcftools index -t {output.vcf}

        rm -f {params.tmp} {params.tmp}.tbi
        {params.log_end}
        """


rule vaf_table_chrX:
    """
    Extract PASS het SNPs (DP>=20, 0.2<=VAF<=0.8) for X-inactivation bias.
    VAF folded to upper half: VAF = VAF if VAF>=0.5 else 1-VAF.
    Median sits at 0.5 for balanced inactivation, rises toward 1.0 for skew.
    """
    input:
        vcf = rules.annotate_chrX_population.output.vcf
    output:
        tsv = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/vaf_table.tsv"
    conda:
        PIPELINE_DIR + "/envs/bcftools_env.yml"
    params:
        log_start = lambda wc, input, threads: log_start("vaf_table_chrX", wc, threads),
        log_end   = LOG_END
    log:
        run_info = "log/vaf_table/{sample}/log.txt",
        time     = "log/vaf_table/{sample}/time.txt"
    shell:
        r"""
        set -euo pipefail
        {params.log_start}
        bcftools view -f PASS -v snps \
            -i 'FORMAT/DP>=20 && FORMAT/VAF>=0.2 && FORMAT/VAF<=0.8' \
            {input.vcf} | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%VAF]\n' | \
        awk 'BEGIN{{OFS="\t"}} {{
            vaf=$6+0;
            folded = (vaf <= 0.5) ? vaf : 1-vaf;
            print $1,$2,$3,$4,$5,vaf,folded
        }}' \
        > {output.tsv}
        {params.log_end}
        """

rule vaf_violin_plot_run_females:
    """
    One violin plot per run, chrX VAF, females only, blood samples only.
    Uses minor allele VAF (1-VAF when VAF>0.5) to show X-inactivation bias.

    Contrôle positif (optionnel) :
      Renseignez dans le config :
        x_bias_control_vaf   : chemin vers la vaf_table.tsv d'un run précédent
        x_bias_control_label : identifiant affiché sur le plot (ex: "25D1234-ctrl")
      Si non renseignés, le plot est produit sans contrôle.
    """
    input:
        vaf_tables = lambda wc: expand(
            FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/vaf_table.tsv",
            sample=SAMPLES
        ),
        sex = FASTQ_DIR + "/../pipeline_v0/qc/chrY_mean_expression.tsv"
    output:
        plot = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/vaf_violin_chrX_females.png"
    params:
        scripts       = PIPELINE_DIR,
        control_vaf   = config.get("x_bias_control_vaf",   ""),
        control_label = config.get("x_bias_control_label", ""),
        log_start = lambda wc, input, threads: log_start("vaf_violin_plot_run_females", wc, threads),
        log_end   = LOG_END
    conda:
        PIPELINE_DIR + "/envs/python_plot_env.yml"
    log:
        run_info = "log/vaf_violin/log.txt",
        time     = "log/vaf_violin/time.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        CTRL_FLAG=""
        if [ -n "{params.control_vaf}" ] && [ -f "{params.control_vaf}" ]; then
            CTRL_FLAG="--control {params.control_vaf} {params.control_label}"
        fi
        python {params.scripts}/scripts/plot_vaf_violin_by_sample.py \
            {output.plot} \
            {input.sex} \
            $CTRL_FLAG \
            {input.vaf_tables} >> {log.run_info} 2>&1
        {params.log_end}
        """

