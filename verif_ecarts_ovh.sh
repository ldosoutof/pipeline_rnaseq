#!/usr/bin/env bash
# =============================================================================
# verif_ecarts_ovh.sh — Vérification des écarts de la qualification initiale
# (QUAL-INIT-RNASEQ-001, § 7) côté OVH. LECTURE SEULE : ne modifie rien.
# Usage : bash verif_ecarts_ovh.sh | tee verif_ecarts_ovh_$(date +%F).txt
# =============================================================================
set -u
PIPE="${PIPE:-/home/ldosoutoferreira/pipeline/RNASEQ/routine/test_pipeline_v2/pipeline_corrige/pipeline_fix/pipeline}"
PROD="${PROD:-/datawork2/genetique/RNASeq/diag/prod}"
RUN_DV="${RUN_DV:-RUN54}"        # run relancé avec DeepVariant (NC-02)
RUN_ANNOT="${RUN_ANNOT:-RUN53}"  # run réannoté (NC-04, NC-05)

ok()   { echo "  [OK]      $*"; }
ko()   { echo "  [KO]      $*"; }
voir() { echo "  [À VOIR]  $*"; }
echo "Vérification des écarts — $(hostname) — $(date '+%F %T')"
echo "PIPE=$PIPE"; echo "PROD=$PROD"; echo

# ---------------------------------------------------------------- NC-01 / QI-03
echo "NC-01 — pipeline_versions.tsv complet (aucun 'not found')"
f=$(find "$PROD" -maxdepth 5 -name pipeline_versions.tsv -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1 | cut -d' ' -f2-)
if [ -z "$f" ]; then voir "aucun pipeline_versions.tsv trouvé sous $PROD"
else
  echo "  fichier le plus récent : $f ($(date -r "$f" '+%F %T'))"
  n=$(grep -ci "not found" "$f")
  if [ "$n" -eq 0 ]; then ok "aucune version manquante ($(($(wc -l < "$f") - 1)) outils)"
  else ko "$n outil(s) 'not found' :"; grep -i "not found" "$f" | sed 's/^/            /'; fi
  echo "  -> vérifier que ce fichier a été produit APRÈS le déploiement du correctif."
fi
echo

# ---------------------------------------------------------------- NC-02 / QI-04
echo "NC-02 — DeepVariant exécuté en conteneur ($RUN_DV)"
if grep -q -- "--use-singularity\|--sdm .*apptainer" "$PIPE/template/launch_template.sh" 2>/dev/null
then ok "launch_template.sh contient l'exécution en conteneur"
else ko "launch_template.sh sans --use-singularity"; fi
R=$(ls -d "$PROD"/*"$RUN_DV"_* 2>/dev/null | head -1)
if [ -z "$R" ]; then voir "dossier $RUN_DV introuvable sous $PROD"
else
  nb=$(ls "$R"/pipeline_v0/star/*_Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l)
  nv=$(find "$R/pipeline_v0/deepvariant_WES" -name "*_chrX.vcf.gz" 2>/dev/null | wc -l)
  echo "  BAM : $nb   VCF chrX : $nv"
  if [ "$nv" -gt 0 ] && [ "$nv" -eq "$nb" ]; then ok "un VCF chrX par échantillon"
  elif [ "$nv" -gt 0 ]; then voir "VCF présents mais nombre différent des BAM (à expliquer)"
  else ko "aucun VCF chrX : DeepVariant n'a pas abouti"; fi
  g=$(grep -rl "command not found" "$PIPE"/*"$RUN_DV"_*/log/deepvariant_WES/ 2>/dev/null | wc -l)
  [ "$g" -gt 0 ] && voir "$g log(s) DeepVariant contiennent encore 'command not found' (anciens essais ?)"
fi
echo

# ---------------------------------------------------------------- NC-03 / QO-01
echo "NC-03 — FASTQ déposés dans prod/<run>/fastq (et plus dans prod/fastq)"
if [ -d "$PROD/fastq" ]; then
  nf=$(find "$PROD/fastq" -type f 2>/dev/null | wc -l)
  if [ "$nf" -eq 0 ]; then ok "prod/fastq existe mais est vide"
  else voir "prod/fastq contient encore $nf fichier(s) à déplacer vers leur run"
       ls -lt "$PROD/fastq" | head -5 | sed 's/^/            /'; fi
else ok "prod/fastq n'existe plus"; fi
last=$(ls -dt "$PROD"/20*_RUN*/ 2>/dev/null | head -1)
if [ -n "$last" ]; then
  nq=$(ls "$last"fastq/*.fastq.gz 2>/dev/null | wc -l)
  echo "  run le plus récent : $(basename "$last") — $nq FASTQ dans son dossier fastq/"
  [ "$nq" -gt 0 ] && ok "le dernier run a ses FASTQ au bon endroit" || voir "aucun FASTQ dans $(basename "$last")/fastq"
  echo "  -> à confirmer sur un run transféré APRÈS la correction du watcher sitatst."
fi
echo

# ---------------------------------------------------------------- NC-04 / QO-03
echo "NC-04 — fichiers annotés sans gène non testé (p-value vide)"
touched=0
for a in "$PROD"/*/pipeline_v0/outrider*/outrider_htseq_annot.tsv; do
  [ -f "$a" ] || continue
  n=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="pValue") p=i; next} p && ($p=="" || $p=="NA"){c++} END{print c+0}' "$a")
  if [ "$n" -gt 0 ]; then
    touched=$((touched+1))
    echo "  [KO]      $n lignes non testées : ${a#$PROD/}"
  fi
done
[ "$touched" -eq 0 ] && ok "aucun fichier annoté touché" || echo "  -> $touched fichier(s) à réannoter (forcerun de la règle d'annotation du run)"
echo

# ---------------------------------------------------------------- NC-05 / QO-06
echo "NC-05 — durée de l'annotation optimisée ($RUN_ANNOT)"
L=$(ls -d "$PIPE"/*"$RUN_ANNOT"_* 2>/dev/null | head -1)
d=$(grep -rh "\[DONE\].*significant events" "$L"/log 2>/dev/null | tail -2)
if [ -n "$d" ]; then echo "$d" | sed 's/^/  /'; ok "durée journalisée (à comparer aux 21 h CPU initiales)"
else voir "ligne [DONE] introuvable dans $L/log (annotation pas encore rejouée ?)"; fi
if grep -q "sets_by_key" "$PIPE/scripts/annotation_outrider_hits_bis2.py" 2>/dev/null
then ok "version optimisée déployée"; else ko "ancienne version du script d'annotation"; fi
echo

# ---------------------------------------------------------------- NC-06 / QI-02
echo "NC-06 — chemins de référence chrX renseignés dans le template"
T="$PIPE/template/config_template.yml"
grep -nE "^(chrx_bed|g1000_chrx|gnomad_chrx):" "$T" | sed 's/^/  template: /'
if grep -E "^(chrx_bed|g1000_chrx|gnomad_chrx):" "$T" | grep -q "/path/to"
then ko "chemins fictifs encore présents"
     echo "  valeurs trouvées dans les config.yml de runs (candidates) :"
     grep -hE "^(chrx_bed|g1000_chrx|gnomad_chrx):" "$PIPE"/*/launch_folder/config.yml 2>/dev/null \
       | grep -v "/path/to" | sort | uniq -c | sed 's/^/            /'
else ok "aucun chemin fictif"; fi
echo

# ---------------------------------------------------------------- NC-08 / QI-07
echo "NC-08 — pipeline de production sous git, sur un tag"
if git -C "$PIPE" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  v=$(git -C "$PIPE" describe --tags --dirty 2>/dev/null || echo "aucun tag")
  echo "  version : $v"
  case "$v" in *dirty*) ko "modifications locales non commitées" ;; "aucun tag") voir "dépôt sans tag" ;; *) ok "positionné sur $v" ;; esac
else ko "$PIPE n'est pas un dépôt git"; fi
echo
echo "Fin. Coller cette sortie pour mise à jour du dossier de qualification (v1.1)."
