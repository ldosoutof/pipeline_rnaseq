#!/usr/bin/env bash
# =============================================================================
# run_cibersortx_all.sh — Lance CIBERSORTx Fractions (LM22) sur TOUS les runs,
# en dehors de Snakemake. Reproduit fidèlement la règle 07_cibersortx.smk.
#
# Pour chaque run ayant une matrice TPM en symboles (matrice_gene_tpm_gene.tsv) :
#   1) prépare la mixture (en-tête GeneSymbol + agrégation des doublons)
#   2) lance CIBERSORTxFractions via Apptainer (LM22, --perm 100 --QN FALSE)
#   3) écrit les résultats dans <run>/pipeline_v0/cibersortx/ (ÉCRASE si présent)
#
# Un run qui échoue n'arrête pas les autres. Log par run.
# =============================================================================
set -uo pipefail

# ---------------------------------------------------------------------------
# PARAMÈTRES — chemins confirmés (à ajuster si besoin)
# ---------------------------------------------------------------------------
RUNS_DIR="${RUNS_DIR:-/datawork2/genetique/RNASeq/diag/prod}"
SIF="${SIF:-/home/ldosoutoferreira/cibersortx_fractions.sif}"
LM22="${LM22:-/dataref/bank/human/annotation/cibersortx/LM22.txt}"
TOKEN_FILE="${TOKEN_FILE:-$HOME/secrets/cibersortx_token}"
EMAIL="${EMAIL:-laura.dosoutoferreira@chu-nantes.fr}"   # compte CIBERSORTx
PERM="${PERM:-100}"
QN="${QN:-FALSE}"
# nom de la matrice en symboles (entrée) — à ajuster si différent
MATRIX_NAME="${MATRIX_NAME:-matrice_gene_tpm_gene.tsv}"

# ---------------------------------------------------------------------------
# Vérifications préalables
# ---------------------------------------------------------------------------
command -v apptainer >/dev/null 2>&1 || { echo "[ERREUR] apptainer introuvable"; exit 1; }
[ -f "$SIF" ]        || { echo "[ERREUR] SIF introuvable : $SIF"; exit 1; }
[ -f "$LM22" ]       || { echo "[ERREUR] LM22 introuvable : $LM22"; exit 1; }
[ -f "$TOKEN_FILE" ] || { echo "[ERREUR] token introuvable : $TOKEN_FILE"; exit 1; }
[ -n "$EMAIL" ]      || { echo "[ERREUR] EMAIL non renseigné (export EMAIL=... ou éditer le script)"; exit 1; }
TOKEN="$(cat "$TOKEN_FILE")"

echo "=== CIBERSORTx sur tous les runs ==="
echo "Runs      : $RUNS_DIR"
echo "SIF       : $SIF"
echo "LM22      : $LM22"
echo "Matrice   : $MATRIX_NAME"
echo "perm=$PERM  QN=$QN"
echo ""

ok=0; fail=0; skip=0

# ---------------------------------------------------------------------------
# Boucle sur tous les runs
# ---------------------------------------------------------------------------
for run_dir in "$RUNS_DIR"/*/; do
    run="$(basename "$run_dir")"
    matrix="$run_dir/pipeline_v0/kallisto_bed/$MATRIX_NAME"
    # certains pipelines rangent la matrice ailleurs : on tente aussi pipeline_v0/
    [ -f "$matrix" ] || matrix="$run_dir/pipeline_v0/$MATRIX_NAME"

    if [ ! -f "$matrix" ]; then
        echo "[SKIP] $run : pas de $MATRIX_NAME"
        skip=$((skip+1))
        continue
    fi

    out_dir="$run_dir/pipeline_v0/cibersortx"
    mkdir -p "$out_dir"
    log="$out_dir/run_cibersortx.log"
    echo "[RUN ] $run"

    # --- 1) préparer la mixture (identique à la règle Snakemake) -------------
    mixture="$out_dir/mixture_cibersortx.tsv"
    python3 - "$matrix" "$mixture" > "$log" 2>&1 <<'PYEOF'
import sys, pandas as pd
src, out = sys.argv[1], sys.argv[2]
df = pd.read_csv(src, sep="\t", index_col=0)
df.index.name = "GeneSymbol"
# CIBERSORTx exige des symboles uniques : somme des TPM des doublons
if df.index.duplicated().any():
    df = df.groupby(level=0).sum()
df.to_csv(out, sep="\t")
print(f"[prepare_mixture] {df.shape[0]} symboles x {df.shape[1]} echantillons")
PYEOF
    if [ $? -ne 0 ]; then
        echo "[FAIL] $run : préparation mixture (voir $log)"
        fail=$((fail+1))
        continue
    fi

    # --- 2) lancer CIBERSORTx Fractions via Apptainer ------------------------
    cp "$LM22" "$out_dir/LM22.txt"
    apptainer exec \
        -B "$out_dir":/src/data \
        -B "$out_dir":/src/outdir \
        "$SIF" /src/CIBERSORTxFractions \
            --username "$EMAIL" \
            --token "$TOKEN" \
            --mixture "/src/data/$(basename "$mixture")" \
            --sigmatrix /src/data/LM22.txt \
            --perm "$PERM" \
            --QN "$QN" \
            --outdir /src/outdir \
            >> "$log" 2>&1
    rc=$?
    # ne jamais laisser le token en clair dans le log
    sed -i '/token:/d' "$log" 2>/dev/null || true

    if [ $rc -eq 0 ] && [ -f "$out_dir/CIBERSORTx_Results.txt" ]; then
        echo "[ OK ] $run -> $out_dir/CIBERSORTx_Results.txt"
        ok=$((ok+1))
    else
        echo "[FAIL] $run : CIBERSORTx (rc=$rc, voir $log)"
        fail=$((fail+1))
    fi
done

echo ""
echo "=== Terminé : $ok OK, $fail échecs, $skip sans matrice ==="
