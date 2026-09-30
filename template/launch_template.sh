#!/bin/bash

# Script de lancement du pipeline RNA-seq
# Le rulegraph et le filegraph sont désormais générés par la règle
# save_rulegraph (06_versions.smk) directement dans pipeline_v0/metrics/.

conda activate snakemake

conda_prefix="~/pipeline/RNASEQ/routine/conda_env"
# Images de conteneurs (DeepVariant ~ plusieurs Go) : téléchargées une seule fois ici
singularity_prefix="$HOME/pipeline/RNASEQ/routine/singularity_img"

snakefile="ROOT_PIPELINE/snakemake/pipeline.smk"
config_file="PATH_TO_CONFIG"

echo "Launching pipeline"
# --use-singularity : les règles avec `container:` (DeepVariant) s'exécutent dans
#   leur conteneur via Apptainer ; sans cette option -> "run_deepvariant: command not found".
# -B /datawork2,/dataref : rend les BAM, la référence, les BED et les sorties
#   visibles à l'intérieur du conteneur.
snakemake \
    -s  $snakefile \
    --configfile $config_file \
    --latency-wait 60 \
    -c 60 \
    --restart-times 2 \
    --rerun-incomplete \
    --use-conda \
    --conda-prefix "$conda_prefix" \
    --conda-frontend conda \
    --use-singularity \
    --singularity-prefix "$singularity_prefix" \
    --singularity-args "-B /datawork2,/dataref"
