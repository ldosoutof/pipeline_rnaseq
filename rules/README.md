# rules/ — Règles Snakemake actives et archivées

## Règles actives (incluses dans `snakemake/pipeline.smk`)

| Fichier | Étapes couvertes |
|---|---|
| `logging.smk` | Helpers JSON structurés (`log_start`, `LOG_END`, `collect_failed_logs`) |
| `01_trim_fastqc.smk` | FastQC brut · fastp trimming · FastQC post-trim |
| `02_alignement.smk` | Alignement STAR 2-pass · indexation BAM samtools |
| `03_comptage.smk` | HTSeq · Kallisto · tx2gene (TPM) · matrice multi-run |
| `03bis_featureCount.smk` | featureCounts CDS · GeneID→ENSG remap · matrice hyper |
| `04_outrider_fraser.smk` | OUTRIDER · FRASER (pipeline normal HTSeq) · annotations |
| `04_outrider_fraser_hyper.smk` | OUTRIDER hyper · FRASER hyper (pipeline featureCounts CDS) |
| `04_make_zip.smk` | Bundling ZIP par run · analyse par échantillon |
| `05_metrics.smk` | Picard MarkDup · bam_stats · MultiQC · volcano · boxplots · notebook PCA |
| `06_versions.smk` | Capture des versions d'outils → `pipeline_versions.tsv` |

## Module optionnel

| Fichier | Statut |
|---|---|
| `experimental/07_variantCalling_bis.smk` | Activable via `run_variant_calling: true` dans le config |

## Scripts CLI / ops (appelés manuellement, pas depuis les règles)

| Script | Usage |
|---|---|
| `scripts/annotation_coverage.py` | QC couverture annotation panel DI vs GTF (voir F24) |
| `scripts/backfill_db.py` | Backfill historique base de données monitoring |
| `scripts/backfill_qc_summary.py` | Backfill résumé QC historique |
| `scripts/preprocess.py` | Prétraitement ad hoc (usage manuel) |

## Fichiers archivés (`rules/archive/`)

Ces fichiers sont des variantes ou versions antérieures conservées à titre d'archive.
Ils ne sont pas inclus dans `pipeline.smk` et ne sont pas exécutés.

Pour les déplacer depuis la racine `rules/` vers `archive/`, utiliser `git mv` afin
de conserver l'historique git. Le script `audit_dead_code.py` peut le faire avec `--move`.

> Ne pas éditer ces fichiers. Toute évolution doit se faire dans les règles actives.

## Fichiers à supprimer

- `rules/.zip` — archive ad hoc non documentée, ne pas versionner
  (`git rm rules/.zip` puis `git commit`)

## Utilitaire d'hygiène

```bash
# Lister les .smk non atteignables depuis pipeline.smk
python audit_dead_code.py --pipeline snakemake/pipeline.smk

# Déplacer automatiquement les archives
python audit_dead_code.py --pipeline snakemake/pipeline.smk --move
```
