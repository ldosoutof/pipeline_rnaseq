TOOL_VERSIONS["featurecounts"] = get_version_from_env(
    PIPELINE_DIR + "/envs/featurecounts_env.yml",
    "featureCounts -v 2>&1 | head -1"
)

rule featurecounts_gene:
    """
    Count genes using featureCounts (absolute counts, paired-end, stranded reverse)
    Uses CDS features from RefSeq GTF, paired-end, reverse-stranded.
    """
    version: TOOL_VERSIONS["featurecounts"]
    input:
        bam = lambda wc: next(
            b for r, s, b in RUN_SAMPLE
            if r == wc.run and s == wc.sample
        )
    output:
        gene = f"{PROD_ROOT}/{{run}}/pipeline_v0/featureCounts_gencode/{{sample}}/{{sample}}_gene_cds_counts.txt"
    conda:
        PIPELINE_DIR + "/envs/featurecounts_env.yml"
    params:
        gtf = GTF_REFSEQ
    benchmark:
        f"{PROD_ROOT}/{{run}}/benchmarks/featureCounts/{{sample}}.tsv"
    log:
        run_info = f"{PROD_ROOT}/{{run}}/pipeline_v0/featureCounts_gencode/{{sample}}/{{sample}}_featurecounts.log",
        time     = f"{PROD_ROOT}/{{run}}/log/featureCounts/{{sample}}/time.txt"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.gene})
        echo "start : $(date +"%d-%m-%y   %T")" > {log.time}
        featureCounts \
            -p --countReadPairs -B -s 2 -t CDS -g gene_id \
            -T {threads} \
            -a {params.gtf} \
            -o {output.gene}.tmp \
            {input.bam} \
            &> {log.run_info}
        grep -E  '^#'          {output.gene}.tmp >  {output.gene}
        grep -E  '^Geneid'     {output.gene}.tmp >> {output.gene}
        grep -Ev '^(#|Geneid)' {output.gene}.tmp >> {output.gene}
        rm {output.gene}.tmp
        echo "end : $(date +"%d-%m-%y   %T")" >> {log.time}
        """


# ── Mapping GeneID → ENSG (once, from Ensembl GTF) ───────────────────────────
# The Ensembl GTF contains: gene_id "ENSG..." + db_xref "GeneID:XXXXXX"
# This file is reusable across all future runs.
rule build_geneid_map:
    """
    Build a GeneID -> ENSG mapping table from the Ensembl GTF.
    Runs once; output is reused for all samples/runs.
    """
    input:
        gtf = GTF_ENSEMBL
    output:
        GENEID_MAP          # e.g. "geneid_to_ensg.tsv"
    conda:
        "outrider"
    log:
        "logs/build_geneid_map.log"
    shell:
        """
        mkdir -p logs
        python build_geneid_ensg_map.py \
            --gtf {input.gtf} \
            --out {output} \
            &> {log}
        """


rule map_refseq_to_ensembl:
    """
    Replace RefSeq GeneIDs with Ensembl gene IDs (ENSG) using the prebuilt map.
    """
    input:
        counts = rules.featurecounts_gene.output.gene,
    output:
        ensembl_counts = f"{PROD_ROOT}/{{run}}/pipeline_v0/featureCounts_gencode/{{sample}}/{{sample}}_gene_cds_counts_ensembl.txt"
    params:
        geneid_map = GENEID_MAP
    log:
        f"{PROD_ROOT}/{{run}}/log/featureCounts/{{sample}}/map_refseq_to_ensembl.log"
    conda:
        "outrider"
    shell:
        """
        GM="{params.geneid_map}"
        if [ ! -f "$GM" ]; then
            echo "[WARN] geneid_map not found at $GM — copying counts as-is without ENSG remapping" >> {log}
            cp {input.counts} {output.ensembl_counts}
            exit 0
        fi
        mkdir -p $(dirname {output.ensembl_counts})

        # Extract comment + header lines
        grep -E '^(#|Geneid)' {input.counts} > {output.ensembl_counts}.header

        # Extract data rows (skip comment + header)
        grep -Ev '^(#|Geneid)' {input.counts} > {output.ensembl_counts}.data

        # Map GeneID (col1) -> ENSG using the prebuilt TSV (GeneID<TAB>ENSG)
        join -1 1 -2 1 -t $'\t' \
            <(cut -f1 {output.ensembl_counts}.data | sort) \
            <(sort "$GM") \
        | awk -F'\t' '{{print $2}}' > {output.ensembl_counts}.ensg

        # Rebuild counts: replace col1 with ENSG
        paste {output.ensembl_counts}.ensg \
              <(cut -f2- {output.ensembl_counts}.data) \
        > {output.ensembl_counts}.body

        # Reassemble: header then data
        cat {output.ensembl_counts}.header \
            {output.ensembl_counts}.body \
        > {output.ensembl_counts}

        rm -f {output.ensembl_counts}.header \
              {output.ensembl_counts}.data \
              {output.ensembl_counts}.ensg \
              {output.ensembl_counts}.body

        echo "Done mapping GeneID -> ENSG for {wildcards.sample}" &> {log}
        """

# =============================================================================
# Assemble featureCounts ENSG-remapped files into a single matrix for OUTRIDER
# =============================================================================
rule matrix_featurecounts:
    """
    Aggregate all per-sample featureCounts ENSG-remapped files into one
    gene × sample count matrix — the hyper OUTRIDER input.
    Mirrors rule matrix (HTSeq) but reads featureCounts output instead.
    """
    input:
        counts = lambda wc: [
            f"{PROD_ROOT}/{r}/pipeline_v0/featureCounts_gencode/{s}/{s}_gene_cds_counts_ensembl.txt"
            for r, s, _ in RUN_SAMPLE
        ]
    output:
        matrix = FASTQ_DIR + "/../pipeline_v0/featureCounts_gencode/matrice_fc.txt"
    params:
        runs_dir     = lambda wc: PROD_ROOT,
        blacklist    = outrider_blacklist,
        run_filter   = lambda wc: RUN_NAME,
        pipeline_dir = lambda wc: config.get("pipeline_output_dir", "pipeline_v0"),
        keywords     = lambda wc: config.get("fraser_keywords", "MOINS,PUROMINS").replace(" ", ","),
        scripts      = PIPELINE_DIR,
        log_start    = lambda wc, input, threads: log_start("matrix_featurecounts", wc, threads),
        log_end      = LOG_END
    conda:
        PIPELINE_DIR + "/envs/htseq_env.yml"
    log:
        run_info = "log/matrix_featurecounts/log.txt",
        time     = "log/matrix_featurecounts/time.txt"
    benchmark:
        "benchmarks/matrix_featurecounts/benchmark.tsv"
    threads: 4
    shell:
        """
        set -euo pipefail
        {params.log_start}
        BL="{params.blacklist}"
        [ ! -f "$BL" ] && BL=""
        python {params.scripts}/scripts/create_matrice_featurecounts.py \
            --runs_dir    {params.runs_dir} \
            --output      {output.matrix} \
            --run_filter  {params.run_filter} \
            --pipeline_dir {params.pipeline_dir} \
            --blacklist   "$BL" \
            --keywords    "{params.keywords}" >> {log.run_info} 2>&1
        {params.log_end}
        """
