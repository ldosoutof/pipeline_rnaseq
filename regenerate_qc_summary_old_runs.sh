#!/usr/bin/env bash
# =============================================================================
# regenerate_qc_summary_old_runs.sh
#
# Régénère metrics/qc_summary.tsv (+ qc_warnings.tsv) pour les ANCIENS runs,
# en réutilisant le script du pipeline recup_metrics.py.
#
# recup_metrics.py gère NATIVEMENT les fichiers manquants : toute métrique
# absente (OUTRIDER, FRASER, coverage, gène TPM…) est laissée VIDE (None),
# jamais inventée. Aucune donnée d'un autre run n'est empruntée.
#
# Pour chaque run :
#   1) génère matrice_gene_tpm_gene.tsv (symboles) si absente, via conversion
#      ENSG->symboles (recup_metrics a besoin d'une matrice TPM en symboles
#      valide pour lire HBA1/HBA2/HBB et les gènes DI_green).
#   2) lance recup_metrics.py -> qc_summary.tsv + qc_warnings.tsv
#
# ATTENTION : la boucle interne de recup_metrics part des fichiers
# metrics/dup/*/*_dup.txt. Un run SANS ces fichiers produira un qc_summary
# VIDE (aucun échantillon détecté) — ce n'est pas une erreur, c'est le reflet
# des données réellement disponibles.
# =============================================================================
set -uo pipefail

PROD="${PROD:-/datawork2/genetique/RNASeq/diag/prod}"
PIPE="${PIPE:-/home/ldosoutoferreira/pipeline/RNASEQ/routine/test_pipeline_v2/pipeline_corrige/pipeline_fix/pipeline}"
CONVERT_PY="${CONVERT_PY:-$PIPE/scripts/ensg_to_symbol_cibersortx.py}"
RECUP_PY="${RECUP_PY:-$PIPE/scripts/recup_metrics.py}"
GTF="${GTF:-/dataref/bank/human/annotation/GRCh38/ensembl/106/Homo_sapiens.GRCh38.106.gtf}"
DI_FILE="${DI_FILE:-/dataref/bank/human/annotation/GRCh38/current/DI_green.tsv}"

# seuils d'alerte (mêmes valeurs par défaut que la règle Snakemake)
MAX_DUP="${MAX_DUP:-40.0}"
MAX_OUT_HYPER="${MAX_OUT_HYPER:-20}"
MAX_FRASER_HYPER="${MAX_FRASER_HYPER:-20}"
MAX_HBA_TOTAL="${MAX_HBA_TOTAL:-5000.0}"

ENSG_MATRIX="${ENSG_MATRIX:-marice_gene_tpm.tsv}"
SYMBOL_MATRIX="${SYMBOL_MATRIX:-matrice_gene_tpm_gene.tsv}"
FORCE="${FORCE:-0}"   # 1 = régénérer même si qc_summary.tsv existe déjà

# ---------------------------------------------------------------------------
[ -f "$CONVERT_PY" ] || { echo "[ERREUR] conversion introuvable : $CONVERT_PY"; exit 1; }
[ -f "$RECUP_PY" ]   || { echo "[ERREUR] recup_metrics introuvable : $RECUP_PY"; exit 1; }
[ -f "$GTF" ]        || { echo "[ERREUR] GTF introuvable : $GTF"; exit 1; }
[ -f "$DI_FILE" ]    || { echo "[ATTENTION] DI_green introuvable : $DI_FILE (colonnes DI_green seront vides)"; }

echo "=== Régénération qc_summary (anciens runs) ==="
echo "Prod=$PROD"
echo "GTF=$GTF"
echo "DI_file=$DI_FILE"
echo "FORCE=$FORCE"
echo ""

ok=0; fail=0; skip=0; empty=0

for run_dir in "$PROD"/*/; do
    run="$(basename "$run_dir")"
    pv0="$run_dir/pipeline_v0"
    kb="$pv0/kallisto_bed"
    ensg="$kb/$ENSG_MATRIX"
    symbol="$kb/$SYMBOL_MATRIX"
    metrics_dir="$pv0/metrics"
    qc="$metrics_dir/qc_summary.tsv"

    # pas de pipeline_v0 -> on saute
    [ -d "$pv0" ] || { echo "[SKIP] $run : pas de pipeline_v0"; skip=$((skip+1)); continue; }
    # qc_summary déjà présent et pas de FORCE -> on saute
    if [ -f "$qc" ] && [ "$FORCE" != "1" ]; then
        echo "[SKIP] $run : qc_summary.tsv déjà présent (FORCE=1 pour refaire)"
        skip=$((skip+1)); continue
    fi

    echo "[RUN ] $run"
    mkdir -p "$metrics_dir" 2>/dev/null || { echo "[FAIL] $run : mkdir $metrics_dir (permissions ?)"; fail=$((fail+1)); continue; }
    log="$metrics_dir/regen_qc_summary.log"

    # --- 0) compat ancienne structure : recup_metrics attend
    #        metrics/dup et metrics/coverage ; les anciens runs ont dup/ et
    #        coverage/ directement sous pipeline_v0. On crée des liens
    #        symboliques (non destructif) si besoin.
    for sub in dup coverage; do
        old="$pv0/$sub"
        new="$metrics_dir/$sub"
        if [ ! -e "$new" ] && [ -d "$old" ]; then
            ln -s "../$sub" "$new" 2>/dev/null \
                && echo "        lien metrics/$sub -> ../$sub" \
                || echo "        [WARN] lien metrics/$sub non créé (permissions ?)"
        fi
    done

    # --- 1) matrice en symboles (si absente) ---------------------------------
    if [ ! -f "$symbol" ]; then
        if [ -f "$ensg" ]; then
            echo "        conversion ENSG->symboles..."
            python3 "$CONVERT_PY" --matrix "$ensg" --gtf "$GTF" --out "$symbol" > "$log" 2>&1
            if [ $? -ne 0 ] || [ ! -f "$symbol" ]; then
                echo "[FAIL] $run : conversion matrice (voir $log)"; fail=$((fail+1)); continue
            fi
        else
            echo "[FAIL] $run : ni $SYMBOL_MATRIX ni $ENSG_MATRIX -> pas de matrice TPM"
            fail=$((fail+1)); continue
        fi
    fi

    # --- 2) recup_metrics.py (gère nativement les champs manquants) ----------
    python3 "$RECUP_PY" \
        --run_path        "$run_dir" \
        --tpm_file        "$symbol" \
        --di_file         "$DI_FILE" \
        --output          "$metrics_dir/qc_summary.tsv" \
        --warnings_output "$metrics_dir/qc_warnings.tsv" \
        --max_dup         "$MAX_DUP" \
        --max_out_hyper   "$MAX_OUT_HYPER" \
        --max_fraser_hyper "$MAX_FRASER_HYPER" \
        --max_hba_total   "$MAX_HBA_TOTAL" >> "$log" 2>&1
    rc=$?

    if [ $rc -ne 0 ]; then
        echo "[FAIL] $run : recup_metrics (rc=$rc, voir $log)"; fail=$((fail+1)); continue
    fi
    if [ ! -s "$metrics_dir/qc_summary.tsv" ]; then
        echo "[VIDE] $run : qc_summary généré mais vide (pas de fichiers dup ?)"
        empty=$((empty+1)); continue
    fi
    # nb de lignes de données (hors en-tête)
    n=$(( $(wc -l < "$metrics_dir/qc_summary.tsv") - 1 ))
    echo "[ OK ] $run : qc_summary.tsv ($n échantillons)"
    ok=$((ok+1))
done

echo ""
echo "=== Terminé : $ok OK, $empty vides, $fail échecs, $skip ignorés ==="
