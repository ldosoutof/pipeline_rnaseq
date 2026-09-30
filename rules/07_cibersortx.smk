# =============================================================================
# 07_cibersortx.smk — Déconvolution cellulaire (CIBERSORTx Fractions, LM22)
# =============================================================================
# Estime les proportions de types cellulaires (22 cellules immunitaires, LM22)
# à partir de la matrice TPM gène-level du run, via le conteneur CIBERSORTx
# exécuté avec Apptainer (rootless), puis produit un graphique en barres
# empilées (façon rendu web CIBERSORTx).
#
# Chaîne : matrice_gene_tpm_gene.tsv (SYMBOLES HGNC déjà présents)
#          -> prépare la mixture (ajout en-tête GeneSymbol, agrégation doublons)
#          -> apptainer exec fractions.sif CIBERSORTxFractions (LM22)
#          -> CIBERSORTx_Results.txt (fractions par échantillon)
#          -> plot_cibersortx.py -> cibersortx_barplot.png
#
# Piloté par la section `cibersortx` du config.yml. Si non activé, la règle
# n'est pas ajoutée aux cibles finales (voir rule all).
#
# PRÉREQUIS (hors pipeline, à installer une fois) :
#   - Apptainer disponible (which apptainer -> /usr/bin/apptainer, setuid)
#   - Image CIBERSORTx fractions convertie en .SIF :
#       apptainer build cibersortx_fractions.sif docker://cibersortx/fractions
#   - Fichier de signature LM22.txt (depuis cibersortx.stanford.edu)
#   - Token + email CIBERSORTx (fichier token chmod 600, hors dépôt)
#
# NOTE IMPORTANTE (vérifié en pratique) :
#   * On part de matrice_gene_tpm_gene.tsv qui contient DÉJÀ des symboles HGNC
#     (pas des ENSG) : aucune conversion ENSG->symbole n'est nécessaire ; il
#     suffit d'ajouter l'en-tête "GeneSymbol" en tête de la 1re colonne et
#     d'agréger d'éventuels symboles en double (somme des TPM).
#   * Le token CIBERSORTx est lié à l'IP publique de sortie du serveur : il
#     doit être généré pour cette IP (sinon "Token/username/ip invalid").
#   * Le mode --absolute TRUE n'est PAS supporté dans la version conteneur
#     (renvoie des zéros). On reste en mode relatif.
# =============================================================================

rule cibersortx_prepare_mixture:
    """
    Prépare la matrice de mixture CIBERSORTx à partir de la matrice TPM
    déjà exprimée en symboles HGNC (matrice_gene_tpm_gene.tsv).
    - ajoute l'en-tête "GeneSymbol" en tête de la première colonne
    - agrège (somme) les éventuels symboles en double
    """
    input:
        matrix = rules.matrix_tpm.output.gene_gene,
    output:
        mixture = FASTQ_DIR + "/../pipeline_v0/cibersortx/mixture_cibersortx.tsv",
    conda:
        PIPELINE_DIR + "/envs/htseq_env.yml"
    params:
        out_dir   = FASTQ_DIR + "/../pipeline_v0/cibersortx",
        log_start = lambda wc, input, threads: log_start("cibersortx_prepare_mixture", wc, threads),
        log_end   = LOG_END,
    log:
        run_info = "log/cibersortx/prepare.txt",
    threads: 1
    shell:
        r"""
        set -euo pipefail
        {params.log_start}
        mkdir -p {params.out_dir}
        python - "{input.matrix}" "{output.mixture}" >> {log.run_info} 2>&1 <<'PYEOF'
import sys, pandas as pd
src, out = sys.argv[1], sys.argv[2]
df = pd.read_csv(src, sep="\t", index_col=0)
df.index.name = "GeneSymbol"
# agréger d'éventuels symboles en double (somme des TPM) — CIBERSORTx exige des symboles uniques
if df.index.duplicated().any():
    df = df.groupby(level=0).sum()
df.to_csv(out, sep="\t")
print(f"[prepare_mixture] {df.shape[0]} symboles x {df.shape[1]} echantillons")
PYEOF
        {params.log_end}
        """


rule cibersortx_fractions:
    """
    Lance CIBERSORTx Fractions (LM22) via Apptainer sur la matrice de symboles.
    Produit CIBERSORTx_Results.txt : fractions cellulaires par échantillon.
    """
    input:
        mixture = rules.cibersortx_prepare_mixture.output.mixture,
    output:
        results = FASTQ_DIR + "/../pipeline_v0/cibersortx/CIBERSORTx_Results.txt",
    params:
        out_dir     = FASTQ_DIR + "/../pipeline_v0/cibersortx",
        sif         = config.get("cibersortx", {}).get("sif", ""),
        lm22        = config.get("cibersortx", {}).get("lm22", ""),
        email       = config.get("cibersortx", {}).get("email", ""),
        token_file  = config.get("cibersortx", {}).get("token_file", ""),
        perm        = config.get("cibersortx", {}).get("perm", 100),
        qn          = config.get("cibersortx", {}).get("qn", "FALSE"),
        log_start   = lambda wc, input, threads: log_start("cibersortx_fractions", wc, threads),
        log_end     = LOG_END,
    log:
        run_info = "log/cibersortx/fractions.txt",
    threads: 4
    shell:
        """
        set -euo pipefail
        {params.log_start}

        command -v apptainer >/dev/null 2>&1 || {{ echo "[ERREUR] apptainer introuvable"; exit 1; }}
        [ -f "{params.sif}" ]        || {{ echo "[ERREUR] SIF introuvable : {params.sif}"; exit 1; }}
        [ -f "{params.lm22}" ]       || {{ echo "[ERREUR] LM22 introuvable : {params.lm22}"; exit 1; }}
        [ -f "{params.token_file}" ] || {{ echo "[ERREUR] token introuvable : {params.token_file}"; exit 1; }}

        TOKEN="$(cat {params.token_file})"
        cp {params.lm22} {params.out_dir}/LM22.txt

        apptainer exec \
            -B {params.out_dir}:/src/data \
            -B {params.out_dir}:/src/outdir \
            {params.sif} /src/CIBERSORTxFractions \
                --username {params.email} \
                --token "$TOKEN" \
                --mixture /src/data/$(basename {input.mixture}) \
                --sigmatrix /src/data/LM22.txt \
                --perm {params.perm} \
                --QN {params.qn} \
                --outdir /src/outdir \
                >> {log.run_info} 2>&1

        # ne jamais laisser le token en clair dans le log
        sed -i '/token:/d' {log.run_info} 2>/dev/null || true

        [ -f {output.results} ] || {{ echo "[ERREUR] résultats CIBERSORTx absents"; exit 1; }}
        {params.log_end}
        """


rule cibersortx_plot:
    """
    Graphique en barres empilées des fractions cellulaires (façon rendu web
    CIBERSORTx) à partir de CIBERSORTx_Results.txt.
    """
    input:
        results = rules.cibersortx_fractions.output.results,
    output:
        png = FASTQ_DIR + "/../pipeline_v0/cibersortx/cibersortx_barplot.png",
    conda:
        PIPELINE_DIR + "/envs/python_plot_env.yml"
    params:
        scripts   = PIPELINE_DIR,
        log_start = lambda wc, input, threads: log_start("cibersortx_plot", wc, threads),
        log_end   = LOG_END,
    log:
        run_info = "log/cibersortx/plot.txt",
    threads: 1
    shell:
        """
        set -euo pipefail
        {params.log_start}
        python {params.scripts}/scripts/plot_cibersortx.py \
            --results {input.results} \
            --out     {output.png} >> {log.run_info} 2>&1
        {params.log_end}
        """
