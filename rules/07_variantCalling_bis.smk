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
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript {PIPELINE_DIR}/scripts/mean_chrY_expression.R \
            {input.matrix} {input.genes} {output.tsv}
        """
rule deepvariant_chrX_WES:
    """
    Run DeepVariant on chrX using WES model and capture BED
    """
    input:
        bam = rules.alignment_star.output.bam,
        bai = rules.index_bam.output.bai
    output:
        vcf = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/{sample}_chrX.vcf.gz"
    params:
        ref = GENOME,
        bed = "/dataref/bank/human/bed/SureSelect_V8/S33266340_Padded_wo_chr_onlyX.bed",
        tmp = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/intermediate_results_dir",
        logs = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/logs",
        image = "google/deepvariant:1.8.0"
    threads: 16
    log:
        "log/deepvariant_WES/{sample}/run.log"
    shell:
        r"""
        mkdir -p $(dirname {output.vcf}) {params.tmp} {params.logs}

        echo "Labelleetlabete1991!" | sudo -S docker run \
          -v /datawork2:/datawork2 \
          -v /dataref:/dataref \
          -w /datawork2 \
          {params.image} \
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
          > {log} 2>&1
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
    shell:
        """
        gunzip -c {input.vcf} | bgzip -c > {output.vcf}
        tabix -f -p vcf {output.vcf}
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
        tmp = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/tmp_annot.vcf.gz",
        g1000 = "/dataref/bank/human/annotation/GRCh38/1000G/ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz",
        gnomad = "/dataref/bank/human/annotation/GRCh38/gnomad/gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz"
    conda:
        PIPELINE_DIR + "/envs/bcftools_env.yml"
    log:
        "log/deepvariant/annotate/{sample}.log"
    shell:
        r"""
        set -euo pipefail

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
        """

#rule annotate_chrX_population:
#    """
#    Annotate chrX VCF with 1000G rsID and gnomAD AF
#    """
#    input:
#        vcf = rules.deepvariant_chrX_WES.output.vcf
#    output:
#        vcf = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/{sample}_chrX_annot.vcf.gz"
#    conda:
#        PIPELINE_DIR + "/envs/bcftools_env.yml"
#    shell:
#        """
#        bcftools annotate \
#          -a /dataref/bank/human/annotation/GRCh38/1000G/ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
#          -c ID \
#          {input.vcf} | \
#        bcftools annotate \
#          -a /dataref/bank/human/annotation/GRCh38/gnomad/gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz \
#          -c INFO/AF \
#          -Oz -o {output.vcf}
#
#        bcftools index {output.vcf}
#        """
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
    shell:
        r"""
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
        """

rule vaf_violin_plot_run_females:
    """
    One violin plot per run, chrX VAF, females only, blood samples only.
    Uses minor allele VAF (1-VAF when VAF>0.5) to show X-inactivation bias.
    """
    input:
        vaf_tables = lambda wc: expand(
            FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/{sample}/vaf_table.tsv",
            sample=SAMPLES
        ),
        sex = FASTQ_DIR + "/../pipeline_v0/qc/chrY_mean_expression.tsv"
    output:
        plot = FASTQ_DIR + "/../pipeline_v0/deepvariant_WES/vaf_violin_chrX_females.png"
    conda:
        PIPELINE_DIR + "/envs/python_plot_env.yml"
    shell:
        """
        python {PIPELINE_DIR}/scripts/plot_vaf_violin_by_sample.py \
            {output.plot} \
            {input.sex} \
            {input.vaf_tables}
        """

