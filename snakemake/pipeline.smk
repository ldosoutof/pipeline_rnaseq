shell.executable('bash')
import os
import re
import sys
import csv
import glob as _glob
import smtplib
import subprocess
from os.path import join
from time import strftime, localtime
from datetime import datetime
from email.mime.text import MIMEText
from pathlib import Path

##─────────────────────────────────────────────────────────────────────────────
## Pipeline RNA-seq — CHU Nantes
## Auteur  : Laura DO SOUTO FERREIRA (dosoutoferreira.laura@chu-nantes.fr)
## Mise à jour : 2026
##─────────────────────────────────────────────────────────────────────────────


##─────────────────────────────────────────────────────────────────────────────
## Constantes
##─────────────────────────────────────────────────────────────────────────────

RUNS_DIR     = config['runs_dir'].rstrip("/")
PROD_ROOT    = RUNS_DIR                                # alias explicite
FASTQ_DIR    = config['configuration']['fastq_dir']
OUTPUT_REP   = config['configuration']['outputDir']
PIPELINE_DIR = config['configuration']['pipeline_dir']
SAMPLES      = config['configuration']['samples']
GENOME       = config['genome']
DI_BED       = config['di_bed']
PANELAPP     = config['panelapp']
HPO          = config['hpo']
PHENO        = config['pheno']
PLI          = config['pli']
BED          = config['bed']
RSEQ_BED     = config['rseq_bed']
STAR_GENOME  = config['star_genome']
RSEM_GENOME  = config['rsem_genome']
KALLISTO_IDX = config['kallisto_idx']
PADDED       = config['padded']
GTF          = config['gtfFile']
MATRICES     = config['matrices']
TPM          = config['matrice_tpm']

GTF_REFSEQ   = config.get('gtf_refseq', '')
GENEID_MAP   = config.get('geneid_map', PROD_ROOT + '/geneid_to_ensg.tsv')
GNOMAD       = config.get('gnomad', '')
MENDELIOME   = config.get('mendeliome', '')

fraser_count       = config['fraser_count']
fraser_count_hyper = config.get('fraser_count_hyper',
                                config['fraser_count'].rstrip('/') + '_hyper/')

# Blacklist unifiée (TSV : sample_id <TAB> tool <TAB> reason)
# tool : outrider | fraser | pca | all
BLACKLIST_FILE = config.get('blacklist', '')

# Tag du run courant — défini ici pour être disponible dans tous les includes
_CURRENT_RUN_TAG = os.path.basename(os.path.dirname(FASTQ_DIR))


##─────────────────────────────────────────────────────────────────────────────
## Wildcard constraints globaux
## Empêchent Snakemake de résoudre un wildcard avec une valeur non conforme
## et lèvent une AmbiguousRuleException lisible à la place d'un match silencieux.
##─────────────────────────────────────────────────────────────────────────────

wildcard_constraints:
    # Identifiant complet d'échantillon : ex. 25D2693-STEMC-PUROMOINS-AVITI
    # Autorise alphanum, tiret et underscore ; interdit le séparateur de chemin /
    sample     = r"[A-Za-z0-9][A-Za-z0-9_-]+",
    # Identifiant court (7 premiers caractères) : ex. 25D2693
    samples_id = r"[A-Za-z0-9]{5,10}",
    # Tag de run : commence par une année (20xx) suivie de caractères quelconques
    # ex. 20251015_RUN41_NS2000_xxx
    run        = r"20[0-9]{6}[A-Za-z0-9_-]*"

##─────────────────────────────────────────────────────────────────────────────
## Détection des versions d'outils (traçabilité ISO 15189)
##─────────────────────────────────────────────────────────────────────────────
## Les environnements conda sont créés par Snakemake dans le dossier passé via
## --conda-prefix (launch_template.sh : ~/pipeline/RNASEQ/routine/conda_env),
## et NON dans PIPELINE_DIR/conda_env (ancienne hypothèse -> "not found").
## Ordre de recherche des préfixes :
##   1. config['configuration']['conda_prefix'] (optionnel, explicite)
##   2. préfixe conda de la session Snakemake en cours (--conda-prefix)
##   3. PIPELINE_DIR/conda_env (ancien emplacement)
##   4. ~/pipeline/RNASEQ/routine/conda_env (défaut de launch_template.sh)

def _conda_prefixes():
    cands = []
    cfg = config.get('configuration', {}).get('conda_prefix')
    if cfg:
        cands.append(cfg)
    try:                                    # Snakemake >= 8
        cands.append(workflow.deployment_settings.conda_prefix)
    except Exception:
        pass
    try:                                    # Snakemake 7
        cands.append(workflow.conda_prefix)
    except Exception:
        pass
    cands.append(os.path.join(PIPELINE_DIR, 'conda_env'))
    cands.append('~/pipeline/RNASEQ/routine/conda_env')
    seen, out = set(), []
    for c in cands:
        if not c:
            continue
        p = Path(os.path.expanduser(str(c))).resolve()
        if p.is_dir() and p not in seen:
            seen.add(p)
            out.append(p)
    return out


def _envs_with(rel_path):
    """Dossiers d'environnement contenant rel_path (ex. 'bin/STAR').
    Recherche directe <prefix>/<env>/<rel_path> (pas de glob récursif : rapide)."""
    hits = []
    for prefix in _conda_prefixes():
        for env_dir in sorted(prefix.iterdir()):
            if env_dir.is_dir() and (env_dir / rel_path).exists():
                hits.append(env_dir)
    return hits


_VERSION_CACHE = {}

# Paquet conda qui fournit chaque binaire (repli sur conda-meta quand l'outil ne
# renvoie pas sa version à l'exécution : lenteur au démarrage, locale, format).
_CONDA_PKG = {
    'featureCounts': 'subread', 'htseq-count': 'htseq', 'STAR': 'star',
    'Rscript': 'r-base', 'fastqc': 'fastqc', 'fastp': 'fastp',
    'samtools': 'samtools', 'kallisto': 'kallisto', 'mosdepth': 'mosdepth',
    'multiqc': 'multiqc', 'picard': 'picard',
}


def _conda_meta_version(env_dir, pkg):
    """Version d'un paquet lue dans <env>/conda-meta/<paquet>-<version>-<build>.json."""
    import json
    for meta in sorted((env_dir / 'conda-meta').glob(f'{pkg}-*.json')):
        try:
            d = json.loads(meta.read_text())
        except Exception:
            continue
        if d.get('name') == pkg and d.get('version'):
            return d['version']
    return ''

def get_version_from_env(cmd, env_yml=None):
    """
    Version d'un outil : exécute `cmd` dans l'environnement conda qui contient
    son binaire (bin/ de l'env ajouté en tête du PATH, pour que les wrappers
    java/perl/python et les pipes fonctionnent).
    Si plusieurs environnements contiennent l'outil avec des versions
    DIFFÉRENTES, toutes sont rapportées (ambiguïté explicite plutôt qu'un choix
    arbitraire). env_yml : conservé pour compatibilité, non utilisé.
    """
    if cmd in _VERSION_CACHE:
        return _VERSION_CACHE[cmd]
    bin_name = cmd.split()[0]
    envs = _envs_with(f'bin/{bin_name}')
    if not envs:
        prefixes = ', '.join(str(p) for p in _conda_prefixes()) or 'aucun préfixe existant'
        print(f'[WARN] {bin_name} introuvable dans les environnements conda ({prefixes})')
        res = f'{bin_name}: not found'
        _VERSION_CACHE[cmd] = res
        return res
    versions = []
    for env_dir in envs:
        try:
            r = subprocess.run(
                f'export PATH="{env_dir}/bin:$PATH"; {cmd}',
                shell=True, capture_output=True, text=True,
                executable='/bin/bash', timeout=120,
            )
            v = (r.stdout + r.stderr).strip().split('\n')[0].strip()
        except Exception as e:
            v = ''
            print(f'[WARN] échec de `{cmd}` dans {env_dir.name} : {e}')
        if v and v not in versions:
            versions.append(v)
    if not versions:
        # repli : métadonnées du paquet conda (source indiquée pour la traçabilité)
        pkg = _CONDA_PKG.get(bin_name, bin_name.lower())
        for env_dir in envs:
            v = _conda_meta_version(env_dir, pkg)
            if v and f'{v} (via conda-meta)' not in versions:
                versions.append(f'{v} (via conda-meta)')
    if not versions:
        res = f'{bin_name}: version not found'
    elif len(versions) == 1:
        res = versions[0]
    else:
        print(f'[WARN] {bin_name} : versions différentes selon les environnements : {versions}')
        res = ' | '.join(versions)
    _VERSION_CACHE[cmd] = res
    return res


def get_r_package_version(pkg):
    """Version d'un package R lue dans DESCRIPTION de la bibliothèque R des
    environnements conda (pas d'exécution de R nécessaire)."""
    key = f'Rpkg:{pkg}'
    if key in _VERSION_CACHE:
        return _VERSION_CACHE[key]
    versions = []
    for env_dir in _envs_with(f'lib/R/library/{pkg}/DESCRIPTION'):
        with open(env_dir / f'lib/R/library/{pkg}/DESCRIPTION') as fh:
            for line in fh:
                if line.startswith('Version:'):
                    v = line.split(':', 1)[1].strip()
                    if v not in versions:
                        versions.append(v)
                    break
    if not versions:
        print(f'[WARN] package R {pkg} introuvable dans les environnements conda')
        res = f'{pkg}: not found'
    elif len(versions) == 1:
        res = versions[0]
    else:
        print(f'[WARN] {pkg} : versions différentes selon les environnements : {versions}')
        res = ' | '.join(versions)
    _VERSION_CACHE[key] = res
    return res


##─────────────────────────────────────────────────────────────────────────────
## Helpers — sample ID extraction
##─────────────────────────────────────────────────────────────────────────────

def _extract_short_id(sample_name):
    """
    Extract the core sample ID from any naming convention:
      25D2693-STEMC-PUROMOINS-AVITI  →  25D2693
      26D0198-MOINS                   →  26D0198
      26D0198_MOINS                   →  26D0198
    """
    return re.split(r'[-_]', sample_name)[0]


SAMPLES_ID = [_extract_short_id(s) for s in SAMPLES]

# FRASER / OUTRIDER : uniquement les échantillons sang
FRASER_KEYWORDS = tuple(
    kw.upper() for kw in
    config.get('fraser_keywords', 'MOINS PUROMOINS').split()
)

def _get_full_sample_names(samples, fastq_dir):
    """
    Si le nom de l'échantillon contient déjà un keyword sang, on l'utilise
    tel quel. Sinon on scanne fastq_dir/*_R1.fastq.gz pour retrouver le nom
    complet (supporte NextSeq 26D0198-MOINS et AVITI 25D2693-STEMC-PUROMOINS-AVITI).
    """
    full = []
    for s in samples:
        if any(kw in s.upper() for kw in FRASER_KEYWORDS):
            full.append(s)
        else:
            hits = sorted(_glob.glob(os.path.join(fastq_dir, f'{s}*_R1.fastq.gz')))
            if hits:
                full.append(os.path.basename(hits[0]).replace('_R1.fastq.gz', ''))
            else:
                full.append(s)
    return full

SAMPLES_FULL     = _get_full_sample_names(SAMPLES, FASTQ_DIR)
SAMPLES_BLOOD    = [s for s in SAMPLES_FULL if any(kw in s.upper() for kw in FRASER_KEYWORDS)]
SAMPLES_ID_BLOOD = [_extract_short_id(s) for s in SAMPLES_BLOOD]

print(f'[INFO] FRASER_KEYWORDS    : {FRASER_KEYWORDS}',                       file=sys.stderr)
print(f'[INFO] SAMPLES       ({len(SAMPLES):>3}) : {SAMPLES[:3]}...',         file=sys.stderr)
print(f'[INFO] SAMPLES_BLOOD ({len(SAMPLES_BLOOD):>3}) : {SAMPLES_BLOOD[:3]}...', file=sys.stderr)
print(f'[INFO] SAMPLES_ID_BLOOD ({len(SAMPLES_ID_BLOOD):>3}) : {SAMPLES_ID_BLOOD[:3]}...', file=sys.stderr)


##─────────────────────────────────────────────────────────────────────────────
## Blacklist unifiée
##─────────────────────────────────────────────────────────────────────────────

def _load_blacklist(tool):
    """
    Charge la blacklist unifiée (TSV : sample_id, tool, reason).
    tool : 'outrider' | 'fraser' | 'pca'
    Les lignes avec tool='all' s'appliquent à tous les outils.
    """
    excluded = set()
    if not BLACKLIST_FILE or not os.path.exists(BLACKLIST_FILE):
        return excluded
    with open(BLACKLIST_FILE) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            t   = row.get('tool', '').strip().lower()
            sid = row.get('sample_id', '').strip()
            if sid and not sid.startswith('#') and (t == tool or t == 'all'):
                excluded.add(sid)
    return excluded

def active_samples(tool, samples):
    """Retourne les samples non blacklistés pour un outil donné."""
    excluded = _load_blacklist(tool)
    if excluded:
        print(f'[blacklist] {tool} : {len(excluded)} exclus : {sorted(excluded)}')
    return [s for s in samples if s not in excluded]

# Conservés pour compatibilité avec les règles existantes
outrider_blacklist = BLACKLIST_FILE
fraser_blacklist   = BLACKLIST_FILE
PCA_blacklist      = BLACKLIST_FILE

ACTIVE_FRASER   = active_samples('fraser',   samples=SAMPLES_ID_BLOOD)
ACTIVE_OUTRIDER = active_samples('outrider', samples=SAMPLES_ID_BLOOD)
ACTIVE_PCA      = active_samples('pca',      samples=SAMPLES_ID)


##─────────────────────────────────────────────────────────────────────────────
## RUN_SAMPLE — scan des BAMs historiques
##─────────────────────────────────────────────────────────────────────────────

def _build_run_sample(prod_root, samples_id):
    """
    Scanne prod_root pour les BAMs alignés et retourne
    [(run_tag, sample_id, bam_path), ...] pour chaque échantillon actif.
    """
    triples = []
    pattern = os.path.join(prod_root, '20*', 'pipeline_v0', 'star',
                           '*_Aligned.sortedByCoord.out.bam')
    for bam in sorted(_glob.glob(pattern)):
        parts   = bam.split(os.sep)
        run_tag = next((p for p in parts if p.startswith('20')), None)
        if run_tag is None:
            continue
        basename  = os.path.basename(bam)
        sample_id = basename.split('-')[0] if '-' in basename else basename.split('_')[0]
        if sample_id in samples_id:
            triples.append((run_tag, sample_id, bam))
    return triples


##─────────────────────────────────────────────────────────────────────────────
## Includes
##─────────────────────────────────────────────────────────────────────────────

include: '../rules/logging.smk'
include: '../rules/01_trim_fastqc.smk'
include: '../rules/02_alignement.smk'
include: '../rules/03_comptage.smk'
include: '../rules/03bis_featureCount.smk'
include: '../rules/04_outrider_fraser.smk'
include: '../rules/04_outrider_fraser_hyper.smk'
include: '../rules/04_make_zip.smk'
include: '../rules/05_metrics.smk'
include: '../rules/06_versions.smk'

if config.get('run_variant_calling', True):
    include: '../rules/experimental/07_variantCalling_bis.smk'

if config.get('cibersortx', {}).get('enabled', False):
    include: '../rules/07_cibersortx.smk'

workdir: OUTPUT_REP


##─────────────────────────────────────────────────────────────────────────────
## RUN_SAMPLE (construit après les includes pour que les wildcards soient connus)
##─────────────────────────────────────────────────────────────────────────────

RUN_SAMPLE         = _build_run_sample(PROD_ROOT, set(SAMPLES_ID))
RUN_SAMPLE_CURRENT = [(r, s, b) for r, s, b in RUN_SAMPLE if r == _CURRENT_RUN_TAG]

# Dict pré-calculé : run_tag → liste des sample_ids (utilisé dans 04_make_zip.smk)
# Défini ici (après RUN_SAMPLE) pour être disponible dans tous les includes
_RUN_TO_SAMPLES = {}
for _r, _s, _ in RUN_SAMPLE:
    _RUN_TO_SAMPLES.setdefault(_r, []).append(_s)

# Validation : tous les échantillons du run courant devraient avoir un BAM.
# _build_run_sample scanne les BAMs existants. Au PREMIER lancement d'un run,
# STAR n'a pas encore tourné, donc aucun BAM n'existe — c'est normal.
# On AVERTIT (non bloquant) : Snakemake construit le DAG et déclenche STAR,
# puis featurecounts_gene une fois les BAMs produits.
_current_samples_found   = {s for r, s, b in RUN_SAMPLE_CURRENT}
_current_samples_missing = [s for s in SAMPLES_ID if s not in _current_samples_found]
if _current_samples_missing:
    if len(_current_samples_found) == 0:
        # Run frais : aucun BAM encore — STAR sera déclenché par le DAG. Normal.
        print(
            f'[INFO] Run courant ({_CURRENT_RUN_TAG}) : aucun BAM encore présent '
            f'({len(_current_samples_missing)} échantillon(s)). '
            f'STAR sera exécuté par le pipeline ; featureCounts suivra.',
            file=sys.stderr
        )
    else:
        # Run partiellement aligné : certains BAMs manquent. À signaler.
        print(
            f'[WARN] Run courant ({_CURRENT_RUN_TAG}) : '
            f'{len(_current_samples_missing)} échantillon(s) sans BAM, '
            f'{len(_current_samples_found)} présent(s).\n'
            f'  Sans BAM : {_current_samples_missing}\n'
            f'  → STAR sera (re)déclenché pour les manquants si les FASTQ existent.',
            file=sys.stderr
        )

# Log non bloquant : résumé par run et signalement des runs historiques
# dont certains échantillons attendus sont absents.
# Le run courant est déjà couvert par la WorkflowError ci-dessus.
for _run in sorted(_RUN_TO_SAMPLES):
    _found    = _RUN_TO_SAMPLES[_run]
    _expected = SAMPLES_ID if _run == _CURRENT_RUN_TAG else [
        s for s in SAMPLES_ID if s in {s2 for _, s2, _ in RUN_SAMPLE}
    ]
    print(
        f'[RUN_SAMPLE] {_run}: {len(_found)} BAM(s) trouvé(s)',
        file=sys.stderr
    )

# Signaler les échantillons du config absents de tout run historique
_all_found = {s for _, s, _ in RUN_SAMPLE}
_never_found = sorted(s for s in SAMPLES_ID if s not in _all_found)
if _never_found:
    print(
        f'[WARN] {len(_never_found)} échantillon(s) du config sans BAM dans aucun run '
        f'historique : {_never_found}\n'
        f'  → featurecounts_gene et les règles dépendantes ne seront pas déclenchés '
        f'pour ces échantillons.',
        file=sys.stderr
    )


##─────────────────────────────────────────────────────────────────────────────
## Rule all
##─────────────────────────────────────────────────────────────────────────────

rule all:
    input:
        # ── FASTQ ──────────────────────────────────────────────────────────
        expand(FASTQ_DIR + '/{sample}_R1.fastq.gz', sample=SAMPLES),
        expand(FASTQ_DIR + '/{sample}_R2.fastq.gz', sample=SAMPLES),
        # ── QC & trimming ──────────────────────────────────────────────────
        expand(rules.fastqc_report.output,      sample=SAMPLES),
        expand(rules.fastp.output.html,         sample=SAMPLES),
        expand(rules.fastqc_trim_report.output, sample=SAMPLES),
        # ── Alignement ─────────────────────────────────────────────────────
        expand(rules.alignment_star.output.bam, sample=SAMPLES),
        expand(rules.index_bam.output.bai,      sample=SAMPLES),
        # ── Comptage ───────────────────────────────────────────────────────
        expand(rules.htseq_gene.output.gene,          sample=SAMPLES),
        expand(rules.matrix.output,                   sample=SAMPLES),
        expand(rules.matrix_tpm.output.gene,          sample=SAMPLES),
        expand(rules.kallistoBed.output.h5,           sample=SAMPLES),
        expand(rules.kallisto2gene.output,            sample=SAMPLES),
        # ── featureCounts (run courant uniquement) ──────────────────────────
        [rules.featurecounts_gene.output.gene.format(run=_CURRENT_RUN_TAG, sample=s)
         for s in SAMPLES],
        [rules.map_refseq_to_ensembl.output.ensembl_counts.format(run=_CURRENT_RUN_TAG, sample=s)
         for s in SAMPLES],
        rules.matrix_featurecounts.output.matrix,
        # ── OUTRIDER / FRASER ───────────────────────────────────────────────
        expand(rules.outrider.output.out_file,    sample=SAMPLES),
        rules.annotation_outrider.output.annot,
        expand(rules.fraser_config.output,        sample=SAMPLES),
        rules.fraser.output.fraser,
        expand(rules.volcano.output,              samples_id=ACTIVE_OUTRIDER),
        expand(rules.boxplot.output.filt,         samples_id=ACTIVE_OUTRIDER),
        expand(rules.fraser_boxplot.output.filt,  samples_id=ACTIVE_FRASER),
        # ── Pipeline hyper ──────────────────────────────────────────────────
        rules.annotation_outrider_hyper.output.annot,
        rules.annotation_fraser_hyper.output.annot_fraser,
        rules.rnaseq_per_sample_hyper.output.zip_file,
        # ── Métriques & QC ─────────────────────────────────────────────────
        expand(rules.bam_stats.output.on_target, sample=SAMPLES),
        expand(rules.rseqc.output,               sample=SAMPLES),
        rules.multiqc.output,
        rules.generate_metrics.output.metrics,
        rules.generate_metrics.output.warnings,
        # ── Bundling & rapports ─────────────────────────────────────────────
        # Pour le run COURANT, cleanup_final supprime analysis_input/,
        # analysis_output/ et analysis_output_current/ en fin de run : les exiger
        # ici faisait échouer « all » APRÈS le nettoyage (fichiers absents), donc
        # tout run complet finissait en erreur et onsuccess (courriel de succès,
        # synchronisation vers sitatst) ne s'exécutait jamais. Ces zips
        # intermédiaires ne sont donc exigés que pour les AUTRES runs de
        # RUN_SAMPLE, dont les dossiers ne sont pas nettoyés. Pour le run courant,
        # le zip hyper reste produit (entrée de run_rnaseq_analysis_hyper).
        [rules.make_analysis_zip.output.zip.format(run=r)
         for r in sorted({r for r, s, _ in RUN_SAMPLE} - {_CURRENT_RUN_TAG})],
        [rules.run_rnaseq_analysis.output.result_zip.format(run=r)
         for r in sorted({r for r, s, _ in RUN_SAMPLE} - {_CURRENT_RUN_TAG})],
        [rules.run_rnaseq_analysis_current_run.output.result_zip.format(run=r)
         for r in sorted({r for r, s, _ in RUN_SAMPLE} - {_CURRENT_RUN_TAG})],
        # ── Bundling & rapports — pipeline hyper ───────────────────────────
        [rules.make_analysis_zip_hyper.output.zip.format(run=r)
         for r in sorted({r for r, s, _ in RUN_SAMPLE} - {_CURRENT_RUN_TAG})],
        [rules.run_rnaseq_analysis_hyper.output.result_zip.format(run=r)
         for r in sorted({_CURRENT_RUN_TAG} | {r for r, s, _ in RUN_SAMPLE})],
        rules.generate_and_run_param_notebook.output.executed_nb,
        # ── Versions & rulegraph ────────────────────────────────────────────
        rules.save_pipeline_versions.output,
        rules.save_rulegraph.output.png,
        # ── Biais d'inactivation X (activé par défaut) ─────────────────────
        *(  [rules.mean_chrY_expression.output.tsv,
             rules.vaf_violin_plot_run_females.output.plot]
            if config.get('run_variant_calling', True) else []  ),
        # ── Déconvolution cellulaire CIBERSORTx (optionnel) ────────────────
        *(  [rules.cibersortx_plot.output.png]
            if config.get('cibersortx', {}).get('enabled', False) else []  ),
        # ── Nettoyage final (dernier) : retire analysis_input/output/current ─
        # Dépend des livrables hyper -> s'exécute en tout dernier.
        rules.cleanup_final.output.marker,


##─────────────────────────────────────────────────────────────────────────────
## Cible partielle : relancer uniquement OUTRIDER + pipeline hyper
##
## Usage :
##   snakemake -s pipeline.smk --configfile config.yml \
##       --use-conda -c 60 --target-rules outrider_and_hyper
##
## Pré-requis : matrice HTSeq (matrice.txt) et matrice featureCounts
## doivent déjà exister (produites par un run complet précédent).
##─────────────────────────────────────────────────────────────────────────────

rule outrider_and_hyper:
    input:
        expand(rules.outrider.output.out_file, sample=SAMPLES),
        rules.annotation_outrider.output.annot,
        rules.annotation_outrider.output.files,
        rules.annotation_outrider_hyper.output.annot,
        rules.annotation_fraser_hyper.output.annot_fraser,
        rules.rnaseq_per_sample_hyper.output.zip_file,


##─────────────────────────────────────────────────────────────────────────────
## Notifications mail
##─────────────────────────────────────────────────────────────────────────────

_SMTP_HOST = config.get('smtp_host', 'localhost')
_SMTP_PORT = int(config.get('smtp_port', 25))
_MAIL_FROM = config.get('mail_from', '')
_MAIL_TO   = config.get('mail_to',   [])

def send_email(subject, body):
    try:
        msg            = MIMEText(body)
        msg['Subject'] = subject
        msg['From']    = _MAIL_FROM
        msg['To']      = ', '.join(_MAIL_TO)
        with smtplib.SMTP(_SMTP_HOST, _SMTP_PORT, timeout=10) as server:
            server.send_message(msg)
        print(f'[INFO] Email envoyé à {_MAIL_TO}')
    except Exception as e:
        print(f'[WARN] Envoi email impossible : {e}')

# Le retour des résultats vers sitatst est fait par le watcher OVH
# (scripts/send_results_to_sitatst.py), pas par le pipeline.


onsuccess:
    send_email(
        subject=f'✅ [{_CURRENT_RUN_TAG}] Pipeline RNA-seq terminé avec succès',
        body=(
            f'Le pipeline a terminé sans erreur le '
            f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}.\n'
            f'Run : {_CURRENT_RUN_TAG}'
        )
    )

onerror:
    now      = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    failures = collect_failed_logs(log_dir='log')
    if failures:
        parts = []
        for f in failures:
            parts.append(
                f"Règle      : {f['rule']}\n"
                f"Exit code  : {f['exit_code']}\n"
                f"Timestamp  : {f['timestamp']}\n"
                f"Log file   : {f['log_path']}\n"
                f"--- 30 dernières lignes ---\n{f['tail']}"
            )
        body = (
            f'Erreur dans le pipeline le {now}.\nRun : {_CURRENT_RUN_TAG}\n\n'
            + ('\n' + '='*60 + '\n').join(parts)
        )
    else:
        body = (
            f'Erreur dans le pipeline le {now}.\nRun : {_CURRENT_RUN_TAG}\n'
            'Aucun log structuré trouvé — consultez la sortie Snakemake.\n'
            'Répertoire de logs : log/'
        )
    send_email(subject=f'❌ [{_CURRENT_RUN_TAG}] Pipeline RNA-seq — ERREUR', body=body)
