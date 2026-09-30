# Journal des modifications – pipeline RNA-seq diagnostique

Une section par version publiée (tag git). `[Non publié]` regroupe les modifications présentes sur `dev`,
en attente de la porte de validation (procédure de gestion des versions, § 5). Chaque entrée indique son
impact sur les résultats.

## [Non publié]

## [v2.0.0] – 2026-09-30

État de la production au 30/09/2026, versionné a posteriori. Entre juillet et septembre 2026, le code de
production a évolué hors de git : cette version l'enregistre tel qu'il a été exécuté. Elle sert de version
de référence à la qualification initiale (QUAL-INIT-RNASEQ-001).

**Impact sur les résultats : oui** (scripts OUTRIDER, FRASER, matrices et métriques modifiés) → version majeure.

### Ajouté
- Déconvolution cellulaire CIBERSORTx (signature LM22) : `rules/07_cibersortx.smk`, `scripts/ensg_to_symbol_cibersortx.py`, `scripts/plot_cibersortx.py`, `scripts/run_cibersortx_all.sh`, `docs/CIBERSORTX_INTEGRATION.md`.
- Journalisation structurée des règles : `rules/logging.smk`.
- Configuration propre au site séparée du template : `template/site_paths.yml` (fusionné par `preprocess.py` au-dessus de `config_template.yml`), `site_paths_template.yml`, `blacklist_template.tsv`.
- Blacklists FRASER et OUTRIDER versionnées (identifiants courts d'échantillons, sans nom) et outils associés : `convert_blacklist.py`, `merge_blacklists.py`.
- Watcher OVH (`watcher_ovh.py`) et son script de déploiement (`deploy_watcher_ovh.sh`).
- Scripts : `annotation_coverage.py`, `build_geneid_ensg_map.py`, `update_db.py`, `backfill_db.py`, `backfill_qc_summary.py`, `regenerate_qc_summary_old_runs.sh`, `verif_ecarts_ovh.sh`, `audit_dead_code.py`, `make_release.sh`.

### Modifié
- Règles 01 à 06, `snakemake/pipeline.smk`, templates de configuration et de lancement.
- Scripts OUTRIDER et FRASER (`outrider_newVersion.R`, `fraser.R`, `fraser_newVer.R`, `outrider_1file.py`, `fraser_1file.py`, `fraser_config_create.py`), construction des matrices, métriques (`recup_metrics.py`), analyse par échantillon.
- Module biais d'inactivation X : règle déplacée dans `rules/experimental/07_variantCalling_bis.smk`, activée par défaut.
- Correctifs de septembre 2026, déjà déployés en production :
  - Annotation OUTRIDER : le filtre p < 0,05 exclut désormais les gènes non testés (p-value vide) ; réécriture vectorisée pour les performances. Écarts NC-04 et NC-05. Livrables par échantillon non affectés.
  - Détection des versions d'outils dans le `--conda-prefix` réel du run (écart NC-01).
  - DeepVariant exécuté en conteneur : `--use-singularity` et montages dans `launch_template.sh` (écart NC-02).
  - `qc_db` vidé dans `site_paths.yml` : fin des écritures dans une base QC parallèle sur OVH.

### Retiré
- Scripts et règles obsolètes absents de la production : `rules/05_test.smk`, `rules/07_variantCalling.smk`, anciennes variantes d'annotation, de boîtes à moustaches, de matrices et de requêtes (`annotation_outrider_v2_rare.py`, `annotation_test_fraser.py`, `bam_config.py`, `boxplots3.R`, `boxplots4.R`, `create_matrice.R`, `create_matrice_tpm_gene.R`, `fraser_config_wblacklist.py`, `matrice_of_target.py`, `mean_kallisto.py`, `outrider_new.R`, `request_datalake2_listGenes.py`, `request_fraser_datalake2_listGenes_v4_multithreads.py`, `script_recup_metrics.py`, `taux_dup.py`, `tx2gene_abund.R`, `volcano2.R`).
- `scripts/concatenation.sh` : utilisé par le watcher de sitatst, qui en possède sa propre copie (à versionner avec l'outillage de sitatst).

### Limites connues
- `pipeline_versions.tsv` : versions de featureCounts et de MultiQC non détectées (voir `[Non publié]`).
- `monitoring/` est conservé dans cette version (appelé par `generate_metrics`) ; son retrait est prévu sur `dev`.

## Historique antérieur (reconstitué depuis git)
- **Mars à juillet 2026 (non étiqueté)** : règle d'appel de variants chrX et graphique de biais X, mises à jour des règles, scripts, templates et environnements, tests unitaires, téléchargement des fichiers PanelApp, HPO et pLI.
- **[v1.0.1] – 2025-11-13** : correctif d'affichage (seuil 10 TPM des points du run).
- **[v1.0] – 2025-11-10** : versions des outils collectées et fusionnées, notification par e-mail, nettoyage des scripts.
- **Septembre à octobre 2025** : mise en place initiale du dépôt, corrections OUTRIDER et FRASER (blacklist), nouvelles métriques, README.
