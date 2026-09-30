TOOL_VERSIONS["featurecounts"] = get_version_from_env("featureCounts -v 2>&1 | grep -m1 -i 'featurecounts v'")

rule featurecounts_gene:
    """
    Count genes using featureCounts (absolute counts, paired-end, stranded reverse)
    Uses CDS features from RefSeq GTF, paired-end, reverse-stranded.
    """
    version: TOOL_VERSIONS["featurecounts"]
    input:
        # Pour le run courant : prendre directement l'output de alignment_star
        # (chaîne le DAG → STAR se déclenche automatiquement, même run frais).
        # Pour les runs historiques : utiliser le BAM préexistant scanné dans RUN_SAMPLE.
        bam = lambda wc: (
            os.path.abspath(FASTQ_DIR + f"/../pipeline_v0/star/{wc.sample}_Aligned.sortedByCoord.out.bam")
            if wc.run == _CURRENT_RUN_TAG
            else next(
                b for r, s, b in RUN_SAMPLE
                if r == wc.run and s == wc.sample
            )
        )
    output:
        gene = f"{PROD_ROOT}/{{run}}/pipeline_v0/featureCounts_gencode/{{sample}}/{{sample}}_gene_cds_counts.txt"
    conda:
        PIPELINE_DIR + "/envs/featurecounts_env.yml"
    params:
        gtf       = GTF_REFSEQ,
        log_start = lambda wc, input, threads: log_start("featurecounts_gene", wc, threads),
        log_end   = LOG_END
    benchmark:
        f"{PROD_ROOT}/{{run}}/benchmarks/featureCounts/{{sample}}.tsv"
    log:
        run_info = f"{PROD_ROOT}/{{run}}/pipeline_v0/featureCounts_gencode/{{sample}}/{{sample}}_featurecounts.log",
        time     = f"{PROD_ROOT}/{{run}}/log/featureCounts/{{sample}}/time.txt"
    threads: 8
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p $(dirname {output.gene})
        featureCounts \
            -p --countReadPairs -B -s 2 -t CDS -g gene_id \
            -T {threads} \
            -a {params.gtf} \
            -o {output.gene}.tmp \
            {input.bam} \
            >> {log.run_info} 2>&1
        grep -E  '^#'          {output.gene}.tmp >  {output.gene}
        grep -E  '^Geneid'     {output.gene}.tmp >> {output.gene}
        grep -Ev '^(#|Geneid)' {output.gene}.tmp >> {output.gene}
        rm {output.gene}.tmp
        {params.log_end}
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
        gtf = GTF
    output:
        GENEID_MAP          # e.g. "geneid_to_ensg.tsv"
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    params:
        scripts   = PIPELINE_DIR,
        log_start = lambda wc, input, threads: log_start("build_geneid_map", wc, threads),
        log_end   = LOG_END
    log:
        run_info = "log/build_geneid_map/log.txt",
        time     = "log/build_geneid_map/time.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p $(dirname {output})
        python {params.scripts}/scripts/build_geneid_ensg_map.py \
            --gtf {input.gtf} \
            --out {output} >> {log.run_info} 2>&1
        {params.log_end}
        """


rule map_refseq_to_ensembl:
    """
    Replace RefSeq GeneIDs with Ensembl gene IDs (ENSG) using the prebuilt map.

    The remapping is done in Python to preserve the original featureCounts row
    order. The join/paste shell approach previously used sorted col1 before
    joining, then pasted against the unsorted data — causing systematic
    GeneID↔count misalignment. Unmapped GeneIDs are written as NA_<GeneID>
    and logged so the caller can detect coverage gaps.
    """
    input:
        counts = rules.featurecounts_gene.output.gene,
    output:
        ensembl_counts = f"{PROD_ROOT}/{{run}}/pipeline_v0/featureCounts_gencode/{{sample}}/{{sample}}_gene_cds_counts_ensembl.txt"
    params:
        geneid_map = GENEID_MAP,
        scripts    = PIPELINE_DIR,
        log_start  = lambda wc, input, threads: log_start("map_refseq_to_ensembl", wc, threads),
        log_end    = LOG_END
    conda:
        PIPELINE_DIR + "/envs/outrider_env.yml"
    log:
        run_info = f"{PROD_ROOT}/{{run}}/log/featureCounts/{{sample}}/map_refseq_to_ensembl.log",
        time     = f"{PROD_ROOT}/{{run}}/log/featureCounts/{{sample}}/map_refseq_to_ensembl_time.txt"
    shell:
        """
        set -euo pipefail
        {params.log_start}
        mkdir -p $(dirname {output.ensembl_counts})

        python - <<'PYEOF' >> {log.run_info} 2>&1
import sys

counts_path  = "{input.counts}"
map_path     = "{params.geneid_map}"
out_path     = "{output.ensembl_counts}"

# ── Load the GeneID → ENSG mapping ──────────────────────────────────────────
gid_map = {{}}
if map_path and __import__('os').path.isfile(map_path):
    with open(map_path) as fh:
        for line in fh:
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) >= 2:
                gid_map[parts[0]] = parts[1]
    print(f"[map_refseq_to_ensembl] Loaded {{len(gid_map):,}} GeneID→ENSG pairs",
          file=sys.stderr)
else:
    print(f"[WARN] geneid_map not found at {{map_path}} — copying counts as-is",
          file=sys.stderr)
    __import__('shutil').copy(counts_path, out_path)
    sys.exit(0)

# ── Remap, preserving original row order ────────────────────────────────────
# Auto-detect col1 format from first data row:
#   • Numeric GeneID  (e.g. "100287102")  → look up in gid_map (RefSeq GTF path)
#   • ENSG.version    (e.g. "ENSG00000223972.11") → strip version suffix (GENCODE GTF path)
#   • Plain ENSG      (e.g. "ENSG00000223972")     → pass through unchanged
# This handles both supported GTF types without requiring the caller to configure anything.
n_mapped   = 0
n_unmapped = 0
mode       = None   # "numeric", "ensg_versioned", "ensg_plain" — set on first data row

with open(counts_path) as fin, open(out_path, "w") as fout:
    for line in fin:
        if line.startswith("#") or line.startswith("Geneid"):
            fout.write(line)
            continue
        cols = line.rstrip("\\n").split("\\t")
        gene_id = cols[0]

        # Detect format on first data row
        if mode is None:
            if gene_id.isdigit():
                mode = "numeric"
                print(f"[map_refseq_to_ensembl] Detected numeric GeneID format "
                      f"(RefSeq GTF) — using gid_map lookup", file=sys.stderr)
            elif gene_id.startswith("ENSG") and "." in gene_id:
                mode = "ensg_versioned"
                print(f"[map_refseq_to_ensembl] Detected ENSG.version format "
                      f"(GENCODE GTF) — stripping version suffix, skipping dict lookup. "
                      f"[WARN] gid_map built from Ensembl GTF is unused in this mode. "
                      f"Verify that gtf_refseq points to the intended GTF type.",
                      file=sys.stderr)
            elif gene_id.startswith("ENSG"):
                mode = "ensg_plain"
                print(f"[map_refseq_to_ensembl] Detected plain ENSG format — "
                      f"passing through unchanged", file=sys.stderr)
            else:
                mode = "numeric"
                print(f"[WARN] Unrecognised gene_id format: '{{gene_id}}' — "
                      f"falling back to numeric GeneID lookup", file=sys.stderr)

        if mode == "numeric":
            if gene_id in gid_map:
                cols[0] = gid_map[gene_id]
                n_mapped += 1
            else:
                cols[0] = f"NA_{{gene_id}}"
                n_unmapped += 1
        elif mode == "ensg_versioned":
            cols[0] = gene_id.split(".")[0]   # strip version suffix
            n_mapped += 1
        else:   # ensg_plain — pass through
            n_mapped += 1

        fout.write("\\t".join(cols) + "\\n")

print(f"[map_refseq_to_ensembl] {{n_mapped:,}} mapped, {{n_unmapped:,}} unmapped"
      f" (written as NA_<GeneID>)", file=sys.stderr)
if n_unmapped > 0:
    print(f"[WARN] {{n_unmapped}} GeneIDs had no ENSG match — check GTF version "
          f"compatibility between featureCounts and build_geneid_ensg_map.py",
          file=sys.stderr)
PYEOF
        {params.log_end}
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
        # Attendre les fichiers ENSG-remappés du run courant (noms complets via SAMPLES).
        # Les runs historiques sont scannés directement par le script sur le disque.
        # On utilise SAMPLES (noms complets) car les chemins featureCounts_gencode/
        # sont indexés par le nom d'échantillon complet, pas l'ID court.
        counts = lambda wc: [
            f"{PROD_ROOT}/{_CURRENT_RUN_TAG}/pipeline_v0/featureCounts_gencode/{s}/{s}_gene_cds_counts_ensembl.txt"
            for s in SAMPLES
        ]
    output:
        matrix = FASTQ_DIR + "/../pipeline_v0/featureCounts_gencode/matrice_fc.txt"
    params:
        runs_dir     = lambda wc: PROD_ROOT,
        blacklist    = outrider_blacklist,
        # run_filter = run courant, pour détecter les featureCounts manquants
        run_filter   = _CURRENT_RUN_TAG,
        pipeline_dir = lambda wc: config.get("pipeline_output_dir", "pipeline_v0"),
        keywords     = lambda wc: config.get("fraser_keywords", "MOINS,PUROMOINS").replace(" ", ","),
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
            --pipeline_dir {params.pipeline_dir} \
            --blacklist   "$BL" \
            --keywords    "{params.keywords}" \
            --run_filter  "{params.run_filter}" >> {log.run_info} 2>&1
        {params.log_end}
        """
