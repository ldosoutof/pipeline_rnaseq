#!/usr/bin/env bash
# =============================================================================
# make_pilot_env.sh — prépare un run pilote isolé de la production
#
# Rejoue un run déjà traité en production avec le code de la branche dev, pour la
# porte de validation (procédure de gestion des versions, § 5). Le run pilote :
#   - voit la même cohorte de référence que la production (liens symboliques
#     vers les runs de prod/, en lecture seule) ;
#   - recalcule tout à partir des FASTQ du run (liens vers ceux de prod/) ;
#   - écrit UNIQUEMENT dans DEV_ROOT : sorties du run, caches FRASER, table des
#     identifiants de gènes. Rien n'est jamais écrit dans prod/ ;
#   - n'envoie rien vers sitatst (sync_sitatst désactivé).
#
# Usage :
#   bash tools/make_pilot_env.sh <run_tag> [--copy-fraser-cache] [--all-runs]
#   ex. bash tools/make_pilot_env.sh 20260819_RUN53_NextSeq_High_16RNASEQ
#
# Variables (valeurs par défaut entre crochets) :
#   SITE_PATHS  site_paths.yml de production [<dépôt>/template/site_paths.yml]
#   DEV_ROOT    racine de développement      [/datawork2/genetique/RNASeq/diag/dev]
#
# --copy-fraser-cache : copie les caches de comptage FRASER de production dans
#   DEV_ROOT (plus rapide). Sans cette option, FRASER recompte toute la cohorte
#   à partir des BAM : plus long, mais totalement indépendant de la production.
# --all-runs : lie aussi les runs POSTÉRIEURS au run rejoué. Par défaut, seuls les
#   runs antérieurs (date du nom de dossier) sont liés, pour retrouver la cohorte
#   avec laquelle la production a calculé ce run : OUTRIDER et FRASER dépendent de
#   la cohorte, et des runs plus récents créeraient des écarts étrangers au code.
# =============================================================================
set -euo pipefail

RUN="${1:-}"
COPY_CACHE=0; ALL_RUNS=0
for opt in "${@:2}"; do
    case "$opt" in
        --copy-fraser-cache) COPY_CACHE=1 ;;
        --all-runs)          ALL_RUNS=1 ;;
        *) echo "[ERROR] option inconnue : $opt"; exit 1 ;;
    esac
done
if [ -z "$RUN" ]; then
    sed -n '2,26p' "$0"; exit 1
fi

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SITE_PATHS="${SITE_PATHS:-$REPO/template/site_paths.yml}"
DEV_ROOT="${DEV_ROOT:-/datawork2/genetique/RNASeq/diag/dev}"
python3 -c "import yaml" 2>/dev/null || { echo "[ERROR] python3 sans PyYAML : activer l'env conda snakemake."; exit 1; }

get() { python3 -c "import yaml,sys; v=(yaml.safe_load(open('$SITE_PATHS')) or {}).get('$1',''); print(v or '')"; }
PROD="$(get runs_dir)";            PROD="${PROD%/}"
FC_PROD="$(get fraser_count)";     FC_PROD="${FC_PROD%/}"
FCH_PROD="$(get fraser_count_hyper)"; FCH_PROD="${FCH_PROD%/}"
GMAP_PROD="$(get geneid_map)"

[ -n "$PROD" ] && [ -d "$PROD" ] || { echo "[ERROR] runs_dir introuvable dans $SITE_PATHS"; exit 1; }
[ -d "$PROD/$RUN/fastq" ]       || { echo "[ERROR] $PROD/$RUN/fastq introuvable"; exit 1; }
case "$DEV_ROOT" in "$PROD"|"$PROD"/*) echo "[ERROR] DEV_ROOT ne doit pas être dans $PROD"; exit 1;; esac

COHORT="$DEV_ROOT/cohort"
PILOT="$COHORT/$RUN"
[ -e "$PILOT" ] && { echo "[ERROR] $PILOT existe déjà : le supprimer pour repartir de zéro"; exit 1; }
mkdir -p "$COHORT"

echo "[1/5] Cohorte : liens vers les runs de production (sauf $RUN)"
n=0; skipped=0
for d in "$PROD"/20*_RUN*/; do
    name="$(basename "$d")"
    [ "$name" = "$RUN" ] && continue
    if [ "$ALL_RUNS" = 0 ] && [[ "$name" > "$RUN" ]]; then skipped=$((skipped+1)); continue; fi
    ln -sfn "${d%/}" "$COHORT/$name"; n=$((n+1))
done
echo "      $n runs liés dans $COHORT ; $skipped runs postérieurs exclus (--all-runs pour les inclure)"

echo "[2/5] Run pilote : dossier réel, FASTQ liés depuis la production"
mkdir -p "$PILOT/fastq"
for f in "$PROD/$RUN"/fastq/*.fastq.gz; do ln -s "$f" "$PILOT/fastq/"; done
echo "      $(ls "$PILOT/fastq" | wc -l) FASTQ liés"

echo "[3/5] Caches FRASER de développement"
mkdir -p "$DEV_ROOT/Fraser2_2" "$DEV_ROOT/Fraser2_2_hyper"
if [ "$COPY_CACHE" = 1 ]; then
    for pair in "$FC_PROD:$DEV_ROOT/Fraser2_2" "$FCH_PROD:$DEV_ROOT/Fraser2_2_hyper"; do
        src="${pair%%:*}"; dst="${pair##*:}"
        [ -d "$src" ] && { echo "      copie de $src ($(du -sh "$src" | cut -f1))"; cp -a "$src"/. "$dst"/; }
    done
else
    echo "      caches vides : FRASER recomptera la cohorte (plus long)"
fi

echo "[4/5] Table geneid -> ENSG"
if [ -n "$GMAP_PROD" ] && [ -f "$GMAP_PROD" ]; then
    cp "$GMAP_PROD" "$DEV_ROOT/geneid_to_ensg.tsv"
else
    echo "      absente en production : elle sera construite dans DEV_ROOT"
fi

echo "[5/5] site_paths de développement"
python3 - "$SITE_PATHS" "$DEV_ROOT" <<'PYEOF'
import sys, yaml
src, dev = sys.argv[1], sys.argv[2].rstrip("/")
d = yaml.safe_load(open(src)) or {}
d["runs_dir"]            = f"{dev}/cohort/"
d["fraser_count"]        = f"{dev}/Fraser2_2/"
d["fraser_count_hyper"]  = f"{dev}/Fraser2_2_hyper/"
d["geneid_map"]          = f"{dev}/geneid_to_ensg.tsv"
d["cohort_follow_links"] = True
d.setdefault("sync_sitatst", {})
d["sync_sitatst"]["enabled"] = False
out = f"{dev}/site_paths_dev.yml"
with open(out, "w") as f:
    f.write("# Généré par tools/make_pilot_env.sh — run pilote, NE PAS utiliser en production\n")
    yaml.safe_dump(d, f, allow_unicode=True, sort_keys=False)
print(f"      {out}")
PYEOF

cat <<EOF

Environnement prêt. Depuis le clone de développement (branche dev) :

  cd $REPO
  python3 scripts/preprocess.py --path $REPO/ --workDir $REPO/$RUN/ \\
      --dataDir $PILOT/ --site-paths $DEV_ROOT/site_paths_dev.yml
  cd $REPO/$RUN && bash launch_folder/launch.sh

Sorties du run pilote : $PILOT/pipeline_v0/
À comparer avec       : $PROD/$RUN/pipeline_v0/  (non-régression)
EOF
