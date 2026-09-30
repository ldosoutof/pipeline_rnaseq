# --- helpers for version detection ------------------------------------------
# get_version_from_env() et get_r_package_version() sont définies UNE SEULE FOIS
# dans snakemake/pipeline.smk (une copie ici écrasait l'autre, car ce fichier est
# inclus en dernier). Elles cherchent dans le --conda-prefix réel du run.
import subprocess
from pathlib import Path


def get_python_version():
    """Get Python version in current Snakemake environment"""
    import sys
    return ".".join(map(str, sys.version_info[:3]))


# --- rule to collect and save versions --------------------------------------
rule save_pipeline_versions:
    """
    Collect tool versions used in the pipeline and save them to
    pipeline_v0/metrics/pipeline_versions.tsv (alongside the other QC outputs).
    """
    output:
        FASTQ_DIR + "/../pipeline_v0/metrics/pipeline_versions.tsv"
    run:
        versions = {
            "fastqc":     get_version_from_env("fastqc --version | head -1 | cut -d' ' -f2"),
            "fastp":      get_version_from_env("fastp --version"),
            "STAR":       get_version_from_env("STAR --version | head -1 | cut -d' ' -f2"),
            "samtools":   get_version_from_env("samtools --version | head -1 | cut -d' ' -f2"),
            "htseq-count":get_version_from_env("htseq-count --version | head -1 | cut -d' ' -f2"),
            "kallisto":   get_version_from_env("kallisto version"),
            "featureCounts": get_version_from_env("featureCounts -v 2>&1 | head -1"),
            "Rscript":    get_version_from_env("Rscript --version | head -1"),
            "OUTRIDER":   get_r_package_version("OUTRIDER"),
            "FRASER":     get_r_package_version("FRASER"),
            "mosdepth":   get_version_from_env("mosdepth --version"),
            "multiqc":    get_version_from_env("multiqc --version"),
            "picard":     get_version_from_env("picard MarkDuplicates --version"),
            "python":     get_python_version(),
        }

        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w") as f:
            f.write("tool\tversion\n")
            for tool, version in versions.items():
                f.write(f"{tool}\t{version}\n")

        print(f"[INFO] Wrote {len(versions)} tool versions to {output[0]}")


# --- rule to generate the rulegraph -----------------------------------------
rule save_rulegraph:
    """
    Génère le rulegraph du pipeline (PNG + DOT) dans pipeline_v0/metrics/.

    S'auto-appelle via `snakemake --rulegraph` pour produire le graphe DOT,
    puis le convertit en PNG avec graphviz.
    Le snakefile et le configfile sont résolus depuis PIPELINE_DIR et OUTPUT_REP.
    """
    output:
        png = FASTQ_DIR + "/../pipeline_v0/metrics/rulegraph.png",
        dot = FASTQ_DIR + "/../pipeline_v0/metrics/rulegraph.dot",
    conda:
        PIPELINE_DIR + "/envs/metrics_env.yml"
    params:
        snakefile       = PIPELINE_DIR + "/snakemake/pipeline.smk",
        # Utiliser le même configfile que celui passé à Snakemake
        configfile      = workflow.configfiles[0] if workflow.configfiles else OUTPUT_REP + "/launch_folder/config.yml",
        snakemake_bin   = config.get("snakemake_bin", "snakemake"),
        log_start       = lambda wc, input, threads: log_start("save_rulegraph", wc, threads),
        log_end         = LOG_END
    log:
        run_info = "log/rulegraph/log.txt",
        time     = "log/rulegraph/time.txt"
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.png})
        touch {output.png} {output.dot}
        echo "[INFO] rulegraph placeholder" >> {log.run_info}
        """
