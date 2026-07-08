# Pipeline RNA-seq diagnostique — CHU Nantes

Pipeline Snakemake de détection d'anomalies d'expression (**OUTRIDER**) et
d'épissage (**FRASER**) à partir de RNA-seq, pour le diagnostic en génétique
clinique. Ré-implémentation interne inspirée de DROP.

Le pipeline traite un run de séquençage en modélisant ses échantillons contre
une cohorte de référence (run courant + historique), produit les listes
d'événements aberrants annotées par échantillon, ainsi qu'un ensemble de
métriques qualité.

---

## Vue d'ensemble

Le pipeline enchaîne les étapes suivantes (détail par fichier de règles dans
[`rules/README.md`](rules/README.md)) :

| Étape | Règles | Sortie principale |
|---|---|---|
| Contrôle qualité lecture + trimming | `01_trim_fastqc.smk` | FastQC, fastp |
| Alignement | `02_alignement.smk` | BAM STAR (2-pass) indexés |
| Comptage | `03_comptage.smk` | comptages HTSeq, TPM Kallisto, matrices multi-run |
| Comptage CDS (branche « hyper ») | `03bis_featureCount.smk` | matrice featureCounts cohorte |
| Expression / épissage aberrants | `04_outrider_fraser.smk` | OUTRIDER + FRASER (branche HTSeq) |
| Idem branche hyper | `04_outrider_fraser_hyper.smk` | OUTRIDER + FRASER (branche featureCounts CDS) |
| Bundling livrables | `04_make_zip.smk` | ZIP par run, analyse par échantillon |
| Métriques | `05_metrics.smk` | MarkDup, bam_stats, MultiQC, volcano, boxplots, notebook PCA |
| Versions d'outils | `06_versions.smk` | `pipeline_versions.tsv` |

Deux branches d'analyse coexistent :
- **branche normale (HTSeq)** : comptages gène-level HTSeq → OUTRIDER/FRASER ;
- **branche « hyper » (featureCounts CDS)** : comptages sur CDS → OUTRIDER/FRASER
  hyper. Les fichiers par échantillon sont produits dans les répertoires
  `outrider(_hyper)/filesbysample/` et `fraser(_hyper)/filesbysample/`.

Un module de variant calling est disponible en option
(`experimental/07_variantCalling_bis.smk`), activable via
`run_variant_calling: true` dans le fichier de configuration.

---

## Prérequis

- **Snakemake** et **conda/mamba** (les environnements par règle sont résolus
  via `--use-conda`).
- Les environnements conda sont décrits dans `envs/` (dont `outrider_env.yml`
  et `fraser_env.yml`, tous deux sur `r-base=4.4` / Bioconductor 3.20).
- Accès aux fichiers de référence (génome, index STAR, index Kallisto, GTF,
  fichiers d'annotation) listés dans le fichier de configuration (voir plus bas).

> **Note environnements.** OUTRIDER 1.26.3 et FRASER 2.2.0 utilisent tous deux
> Bioconductor 3.20 (R 4.4). Pour les paquets Bioconductor, épingler uniquement
> `r-base` + le paquet et laisser le solveur dériver les dépendances
> (DESeq2, SummarizedExperiment, …).

---

## Ressources d'annotation (pLI, PanelApp, HPO)

Trois ressources d'annotation à récupérer/régénérer et à référencer dans le
`config.yml` (`pli`, `panelapp`, `hpo`, `pheno`). Toutes en GRCh38. Pour chaque
ressource : consigner **version + date de tirage + md5** (traçabilité ISO 15189),
et rejouer la vérification après tout changement de version.

### pLI / LOEUF (gnomAD v4.1.1, GRCh38)

Fichier de contrainte v4.1.1 déjà coordonné GRCh38, distribué par Ensembl.

```bash
wget https://ftp.ensembl.org/pub/current/variation/LOEUF/loeuf_dataset_grch38.tsv.gz
wget https://ftp.ensembl.org/pub/current/variation/LOEUF/loeuf_dataset_grch38.tsv.gz.tbi

# Extraction pLI/LOEUF, une ligne par gène (MANE Select Ensembl)
{
  printf 'gene\tgene_id\ttranscript\tpLI\tLOEUF\tLOEUF_decile\n'
  zcat loeuf_dataset_grch38.tsv.gz \
    | awk -F'\t' '$5=="true" && $3 ~ /^ENST/ {print $1"\t"$2"\t"$3"\t"$111"\t"$94"\t"$97}'
} > pli_by_gene.grch38.v4.1.1.tsv

md5sum pli_by_gene.grch38.v4.1.1.tsv | tee pli_by_gene.grch38.v4.1.1.md5
```

Notes : pLI = colonne **111** (`lof.pLI`), LOEUF = **94**, flag MANE `true` en
minuscules, filtrer `^ENST` (MANE en double RefSeq/Ensembl). **Ne jamais convertir
les `NA` en 0** (contrainte non calculable ≠ pLI bas). Seuils v4 : pLI ≥ 0,9,
LOEUF < 0,6. Vérifier le build : `MED13` doit ressortir à ~61,9 Mb (GRCh38).

### Panel PanelApp — Intellectual disability (Genomics England, panel 285)

Liste verte (niveau diagnostique, `confidence_level=3`) via l'API REST, paginée.
Nécessite `jq`.

```bash
# version courante
curl -s "https://panelapp.genomicsengland.co.uk/api/v1/panels/285/" | grep -o '"version":"[^"]*"' | head -1

# extraction de la liste verte (coords sous ensembl_genes.GRch38.<version>.location)
BASE="https://panelapp.genomicsengland.co.uk/api/v1/panels/285/genes/?confidence_level=3"
url="$BASE"; : > green_grch38.raw.tsv
while [ -n "$url" ] && [ "$url" != "null" ]; do
  page="$(curl -s "$url")"
  echo "$page" | jq -r '
    .results[] | .gene_data as $g
    | ($g.ensembl_genes.GRch38 // {}) as $b | ($b | keys_unsorted[-1]) as $ver
    | select($ver != null)
    | [ $g.gene_symbol, $g.hgnc_id, $g.ensembl_genes.GRch38[$ver].ensembl_id,
        $b[$ver].location, .mode_of_inheritance ] | @tsv' >> green_grch38.raw.tsv
  url="$(echo "$page" | jq -r '.next')"
done

# contrôle de complétude : doit correspondre au count de l'API (~1513 en v10.47)
wc -l green_grch38.raw.tsv ; curl -s "$BASE" | jq -r '.count'

# conversion en BED trié GRCh38
awk -F'\t' '{split($4,a,"[:-]"); if(a[1]!="") print a[1]"\t"a[2]"\t"a[3]"\t"$1"\t"$2"\t"$3"\t"$5}' \
  green_grch38.raw.tsv | sort -k1,1 -k2,2n > DI_green_$(date +%d%m%y)_grch38.bed

md5sum green_grch38.raw.tsv DI_green_*_grch38.bed | tee panelapp_285.md5
```

### Annotations phénotypiques HPO

Fichiers agnostiques au génome (aucune coordonnée). Utiliser les PURL persistantes ;
télécharger les deux le **même jour** (la version de `genes_to_phenotype.txt` est
celle du `.hpoa`).

```bash
wget http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa
wget http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt
# ontologie, si le pipeline interprète les termes (libellés, hiérarchie) :
wget http://purl.obolibrary.org/obo/hp.obo

grep -E '^#(version|hpo-version|date)' phenotype.hpoa      # relever la version
md5sum phenotype.hpoa genes_to_phenotype.txt | tee hpo_annotations.md5
```

> Procédures détaillées (contrôles de build, pièges de colonnes, obsolescence de
> termes HP) : voir `README_annotation_resources.md`.

---

## Installation

```bash
# 1. cloner le dépôt
git clone <URL_DU_DEPOT> pipeline
cd pipeline

# 2. créer les environnements conda (résolus automatiquement au 1er run via --use-conda)
#    ou les pré-construire dans un préfixe partagé :
#    (adapter le chemin conda_env à votre installation)
```

> Le chemin exact du dépôt distant et la procédure de pré-construction des
> environnements dépendent de votre installation locale — à compléter selon
> votre infrastructure.

---

## Prétraitement (`preprocess.py`) — à lancer avant le pipeline

`scripts/preprocess.py` prépare, pour un run donné, les deux fichiers dont
Snakemake a besoin : `config.yml` et `launch.sh`. Il n'est appelé par aucune
règle (usage manuel) et constitue la **première étape opérationnelle** du flux.

### Ce qu'il fait

1. **Détecte automatiquement les échantillons** en scannant le sous-dossier
   `fastq/` du run : tout fichier `*_R1.fastq.gz` ne commençant pas par
   `Undetermined` est retenu ; le nom d'échantillon est le nom de fichier privé
   de son suffixe `_R1` (p. ex. `26D0643-…-PUROMOINS_R1.fastq.gz`
   → `26D0643-…-PUROMOINS`).
2. **Génère `config.yml`** en fusionnant, par priorité croissante :
   `template/config_template.yml` (valeurs par défaut documentées) →
   `template/site_paths.yml` (chemins réels de l'installation) → les valeurs du
   run auto-détectées (`samples`, `outputDir`, `fastq_dir`, `pipeline_dir`).
   Les valeurs vides d'un niveau supérieur n'écrasent pas celles du template.
3. **Génère `launch.sh`** depuis `template/launch_template.sh` (substitution des
   jetons `PATH_TO_CONFIG` et `ROOT_PIPELINE`), rendu exécutable.

Les deux fichiers sont écrits dans `<workDir>/launch_folder/`.

### Usage

```bash
python scripts/preprocess.py \
    --path    /chemin/vers/pipeline/ \
    --workDir /chemin/vers/répertoire de travail du run/      \
    --dataDir /chemin/vers/run/         # doit contenir un sous-dossier fastq/

# site_paths hors de template/ :
python scripts/preprocess.py -p … -w … -d … \
    --site-paths /chemin/vers/site_paths.yml
```

| Argument | Rôle |
|---|---|
| `-p / --path` | racine du pipeline |
| `-w / --workDir` | répertoire de travail du run (résultats) |
| `-d / --dataDir` | répertoire du run (contient `fastq/`) |
| `--site-paths` | chemins de l'installation (défaut : `template/site_paths.yml`) |

### Points de vigilance

- **Convention FASTQ en dur** : la détection ne reconnaît que le motif
  `*_R1.fastq.gz`. Un autre schéma (`*_R1_001.fastq.gz`, `.fq.gz`, …) ne sera pas
  détecté correctement.
- Si `site_paths.yml` est absent, le script **n'échoue pas** mais avertit et
  génère un config avec les chemins génériques `/path/to/…` du template ; il
  liste alors les clés restées génériques à compléter avant lancement.

### Enchaînement

Une fois `config.yml` et `launch.sh` générés, le run se lance simplement avec
`bash <workDir>/launch_folder/launch.sh` — voir
[Lancement › Lancement simple](#lancement-simple-via-launchsh).

---

## Lancement

### Lancement simple (via `launch.sh`)

Voie normale : après le prétraitement, `preprocess.py` a généré un `launch.sh`
prêt à l'emploi (qui appelle Snakemake avec le `config.yml` du run). Il suffit
alors de :

```bash
# 1. préparer config.yml + launch.sh (voir section Prétraitement)
python scripts/preprocess.py --path … --workDir … --dataDir …

# 2. lancer
bash <workDir>/launch_folder/launch.sh
```

Les sous-sections suivantes décrivent l'appel **manuel** de Snakemake, utile pour
un run partiel, un dry-run, ou pour ajuster les options (`--rerun-triggers`, etc.).

### Run complet (appel Snakemake manuel)

```bash
snakemake -s snakemake/pipeline.smk --configfile config.yml \
    --use-conda --conda-prefix "$CONDA_PREFIX_DIR" \
    -c 60 \
    --rerun-triggers mtime
```

- `-c 60` : nombre de cœurs (adapter à la machine).
- `--rerun-triggers mtime` : ne redéclenche les règles que sur la date des
  fichiers. **Recommandé** : sans ce drapeau, un changement de code/params peut
  déclencher une reconstruction en cascade (jusqu'à STAR).
- `--rerun-incomplete` : à ajouter pour reprendre après une interruption
  (fichiers laissés incomplets par un job arrêté).

> Un seul run à la fois par machine : deux exécutions concurrentes corrompent
> les fichiers intermédiaires partagés.

### Cible partielle : OUTRIDER + branche hyper uniquement

Pré-requis : les matrices HTSeq (`matrice.txt`) et featureCounts
(`matrice_fc.txt`) doivent déjà exister (produites par un run complet
précédent).

```bash
snakemake -s snakemake/pipeline.smk --configfile config.yml \
    --use-conda -c 60 --rerun-triggers mtime \
    outrider_and_hyper
```

### Dry-run (fortement recommandé avant tout run)

```bash
snakemake -s snakemake/pipeline.smk --configfile config.yml \
    -n -r --rerun-triggers mtime --rerun-incomplete
```

Le dry-run liste les jobs qui seront exécutés sans rien lancer : vérifier que
seules les étapes attendues s'y trouvent (et pas, par exemple, un
ré-alignement STAR non désiré).

---

## Sorties

Sous `<outputDir>/pipeline_v0/` :

| Répertoire | Contenu |
|---|---|
| `htseq/`, `featureCounts_gencode/` | comptages et matrices |
| `kallisto_bed/` | TPM gène-level et `matrice_gene_tpm.tsv` |
| `outrider/`, `outrider_hyper/` | résultats OUTRIDER + `filesbysample/` par échantillon |
| `fraser/`, `fraser_hyper/` | résultats FRASER + `filesbysample/` par échantillon |
| `analysis_output*/` | ZIP de livrables par run |
| `metrics/` | `qc_summary.tsv`, MultiQC, volcano, boxplots, notebook PCA |

---

## Dépannage rapide

- **Un correctif ne semble pas pris en compte** : vérifier que le fichier a bien
  été déployé à l'emplacement d'exécution (`pipeline_dir`), pas seulement dans
  une copie. Un `grep` du changement sur le fichier réellement exécuté avant de
  relancer évite un run inutile.
- **Reconstruction en cascade inattendue (STAR, fastp…)** : ajouter
  `--rerun-triggers mtime`.
- **`IncompleteFilesException`** : relancer avec `--rerun-incomplete`.
- **Régénérer une étape précise** : supprimer son fichier de sortie (Snakemake
  reconstruit la règle et son aval) ou utiliser `-R <règle>` (toujours avec
  `--rerun-triggers mtime` pour éviter la cascade).

Pour l'historique détaillé des corrections apportées au pipeline, voir
`CHANGELOG_CORRECTIONS.md` (à la racine du dépôt corrigé).

---

## Structure du dépôt

```
pipeline/
├── snakemake/pipeline.smk      # Snakefile principal (include des règles, rule all)
├── rules/                      # règles Snakemake (voir rules/README.md)
│   ├── 01_trim_fastqc.smk … 06_versions.smk
│   └── experimental/           # module variant calling optionnel
├── scripts/                    # scripts Python / R appelés par les règles
├── envs/                       # environnements conda par outil
├── template/                   # gabarits (rapports, notebooks)
├── monitoring/                 # stack de monitoring (voir monitoring/README.md)
└── tests/                      # tests
```

---

## Configuration

Le pipeline lit un fichier `config.yml` (passé via `--configfile`). Clés
attendues (extrait — voir le fichier de configuration d'exemple du dépôt pour la
liste complète) :

```yaml
runs_dir: "/chemin/vers/les/runs/"
configuration:
  fastq_dir:    "/chemin/vers/run/fastq/"
  outputDir:    "/chemin/vers/run/"
  pipeline_dir: "/chemin/vers/pipeline/"     # dossier CONTENANT envs/ scripts/ rules/
  samples:
    - "SAMPLE_ID_1"
    - "SAMPLE_ID_2"

# références
genome:       "/ref/GRCh38/genome.fa"
star_genome:  "/ref/GRCh38/STAR/"
kallisto_idx: "/ref/GRCh38/kallisto/transcripts.idx"
gtfFile:      "/ref/GRCh38/annotation/Homo_sapiens.GRCh38.106.gtf"
gtf_refseq:   "/ref/GRCh38/annotation/gencode.v49.basic.annotation.nochr.gtf"

# annotation / panels
panelapp:  "/ref/annotation/panelapp_DI_green.tsv"
hpo:       "/ref/annotation/phenotype.hpoa"
pheno:     "/ref/annotation/genes_to_phenotype.txt"
pli:       "/ref/annotation/PLI.vcf"
gnomad:    "/ref/annotation/gnomad_v4_by_gene.tsv"

# paramètres
fraser_keywords: "MOINS PUROMOINS"   # mots-clés d'inclusion des échantillons sang/lympho
run_variant_calling: false
```

> **Important.** `pipeline_dir` doit pointer vers le dossier qui **contient**
> `envs/`, `scripts/` et `rules/` — pas vers un sous-dossier. Les scripts sont
> appelés via `{params.scripts}/scripts/...`.
