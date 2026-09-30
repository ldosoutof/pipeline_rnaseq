# Journal des modifications – pipeline RNA-seq diagnostique

Une section par version publiée (tag git). `[Non publié]` regroupe les modifications présentes sur `dev`,
en attente de la porte de validation (procédure de gestion des versions, § 5). Chaque entrée indique son
impact sur les résultats.

## [Non publié]

### Ajouté
- `pipeline_versions.tsv` : ligne `pipeline` avec la version git du pipeline exécuté (`git describe --tags --dirty`). Chaque run est rattaché à une version identifiée ; le suffixe `-dirty` signale une modification hors git (procédure de gestion des versions, § 3 et § 9). Impact sur les résultats : aucun.

### Modifié
- **Un seul chemin d'ingestion des métriques** : le watcher de sitatst, au retour des résultats. Le pipeline ne pousse plus vers le Pushgateway (`push_metrics.py`) et n'écrit plus dans aucune base QC : ni localement (`qc_db`, règle `generate_metrics`), ni à distance (`update_db.py` lancé après la synchronisation). Le transfert des résultats vers sitatst et le contrôle d'intégrité sont conservés. Impact sur les résultats : aucun (fichiers produits identiques).
- Configuration : clés `qc_db`, `prometheus_*`, `sync_sitatst.update_db`, `remote_db` et `remote_update_db_script` retirées de `config_template.yml`, `site_paths.yml` et `site_paths_template.yml`.

### Corrigé
- `pipeline_versions.tsv` : version de featureCounts détectée (`featureCounts -v` affiche une ligne vide avant la version) ; si un outil ne renvoie pas sa version à l'exécution (cas de MultiQC), elle est lue dans les métadonnées du paquet conda et marquée « via conda-meta ». Clôt le reliquat de l'écart NC-01. Impact sur les résultats : aucun.
- Journal d'événements des règles (`rules/logging.smk`) : les événements de début et de fin étaient écrits dans un fichier littéralement nommé `{log.run_info}` (un paramètre Snakemake n'est pas formaté une seconde fois), et le code de sortie enregistré valait toujours 0. Ils sont désormais écrits dans `log/events/<règle>.jsonl` ; l'événement de fin est posé par un `trap EXIT`, avec le vrai code de sortie, en cas de succès comme d'échec. `collect_failed_logs` lit ce nouvel emplacement, signale aussi les tâches interrompues et joint la fin du `log.txt` de la règle. Tests mis à jour. Impact sur les résultats : aucun.

### Sécurité
- Synchronisation vers sitatst : `sshpass -e` (mot de passe lu dans la variable `SSHPASS`) au lieu de `sshpass -p <mot de passe>`, qui l'exposait à tous les utilisateurs du serveur dans la liste des processus. Utilisé seulement si `pass_file` est renseigné ; la clé SSH dédiée reste la méthode recommandée.

### Documentation
- `config_template.yml` : commentaire expliquant que les chemins chrX réels viennent de `site_paths.yml` (écart NC-06 levé : aucun run n'utilise les chemins fictifs).

### Retiré
- `monitoring/` (Pushgateway, Prometheus, Grafana), `envs/dashboard_env.yml`, `scripts/update_db.py`, `scripts/backfill_db.py` : ils reçoivent ou affichent les métriques et rejoignent le dépôt du dashboard.

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
- Tests unitaires : 80 réussis, 3 en échec dans `tests/test_fraser_config_create.py`. Les tests attendent l'ancien format à 5 colonnes (avec `gene`) ; `fraser_config_create.py` écrit désormais 4 colonnes (`sampleID`, `bamFile`, `group`, `pairedEnd`). Format à confirmer, puis tests à mettre à jour.

## Historique antérieur (reconstitué depuis git)
- **Mars à juillet 2026 (non étiqueté)** : règle d'appel de variants chrX et graphique de biais X, mises à jour des règles, scripts, templates et environnements, tests unitaires, téléchargement des fichiers PanelApp, HPO et pLI.
- **[v1.0.1] – 2025-11-13** : correctif d'affichage (seuil 10 TPM des points du run).
- **[v1.0] – 2025-11-10** : versions des outils collectées et fusionnées, notification par e-mail, nettoyage des scripts.
- **Septembre à octobre 2025** : mise en place initiale du dépôt, corrections OUTRIDER et FRASER (blacklist), nouvelles métriques, README.
