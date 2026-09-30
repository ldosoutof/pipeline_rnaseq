#!/usr/bin/env python3
"""
RNA-Seq Analysis Pipeline - Per Sample Output
Processes FRASER2 and OUTRIDER outputs and generates one file per sample per tool

Stratégie de parallélisation :
  Étape 1 — Chargement des fichiers  : ThreadPoolExecutor  (I/O-bound, libère le GIL)
  Étape 2 — Annotation + écriture    : ProcessPoolExecutor (CPU-bound par sample,
                                        contourne le GIL pour pandas)
  FRASER et OUTRIDER sont traités dans des pools distincts mais peuvent tourner
  simultanément si les ressources le permettent.
"""

import argparse
import json
import os
import zipfile
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
import logging
import sys

import re
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

DEFAULT_WORKERS = min(os.cpu_count() or 4, 16)


# =============================================================================
# Fonctions module-level (picklables -> ProcessPoolExecutor)
# =============================================================================

def _process_and_save_sample(args):
    """
    Traite et sauvegarde les donnees d'UN sample.
    Recoit les donnees du sample deja filtrees (en dict pour la picklabilite),
    applique les annotations et ecrit le fichier TSV.
    Execute dans un process worker -> vraie parallelisation CPU.
    """
    (
        sample_full,
        sample_dict,
        tool_name,
        files_dir,
        gtf_dict,
        gnomad_dict,
        mendeliome_dict,
        gene_col,
    ) = args

    # Extraire l'ID : premier segment avant '-' ou '.' (ex: 25D2507-HOL-Hay → 25D2507)
    sample_short = re.split(r'[-_.]', sample_full)[0]
    filepath = Path(files_dir) / f"{sample_short}.{tool_name}.tab"

    df = pd.DataFrame(sample_dict)

    # -- Annotation GTF -------------------------------------------------------
    if tool_name == 'outrider' and gtf_dict:
        by_gene = gtf_dict.get('by_gene', gtf_dict)  # compat ancien format

        def _gtf_field(gid, field):
            if pd.isna(gid):
                return None
            return by_gene.get(str(gid).split('.')[0], {}).get(field)

        for field in ('gene_name', 'chrom', 'start', 'end', 'strand'):
            df[field] = df['geneID'].apply(lambda x, f=field: _gtf_field(x, f))

        # Récupérer le gene_id complet depuis le GTF
        df['gene_id'] = df['geneID'].apply(
            lambda gid: _gtf_field(gid, 'gene_id') or gid
        )

        if 'gene_name' in df.columns and df['gene_name'].notna().any():
            gene_col = 'gene_name'

    elif tool_name == 'fraser':
        # chrom = seqnames (copie directe)
        df['chrom'] = df['seqnames'].astype(str) if 'seqnames' in df.columns else None

        # hgncSymbol → gene_name
        if 'hgncSymbol' in df.columns and df['hgncSymbol'].notna().any():
            df['gene_name'] = df['hgncSymbol']
            gene_col = 'gene_name'
        elif 'gene_name' in df.columns and df['gene_name'].notna().any():
            gene_col = 'gene_name'
        else:
            gene_col = None

        # Compléter gene_id à partir de gene_name via le GTF (by_gene_name).
        # Le bloc hgncSymbol ci-dessus remplit gene_name sans gene_id : on récupère
        # l'ENSG correspondant au symbole HGNC pour que les deux soient renseignés.
        if gene_col == 'gene_name' and gtf_dict and 'by_gene_name' in gtf_dict:
            by_name = gtf_dict['by_gene_name']
            if 'gene_id' not in df.columns:
                df['gene_id'] = None
            missing_id = df['gene_id'].isna() | df['gene_id'].astype(str).isin(('', 'nan', 'None'))
            df.loc[missing_id, 'gene_id'] = df.loc[missing_id, 'gene_name'].map(
                lambda g: by_name.get(str(g), {}).get('gene_id')
            )

        # Coordinate-based GTF lookup for rows missing gene_name
        if gtf_dict and 'by_interval' in gtf_dict and gtf_dict['by_interval']:
            interval_index = gtf_dict['by_interval']
            gnames = []
            gids   = []
            for _, row in df.iterrows():
                chrom = str(row.get('seqnames', ''))
                try:
                    s = int(row.get('start', 0))
                    e = int(row.get('end',   0))
                except (ValueError, TypeError):
                    gnames.append(None); gids.append(None); continue

                # Check if gene_name already populated
                existing = str(row.get('gene_name', '')) if gene_col else ''
                if existing and existing not in ('', 'nan', 'None'):
                    gnames.append(existing)
                    gids.append(str(row.get('gene_id', '')) or None)
                    continue

                if chrom not in interval_index:
                    gnames.append(None); gids.append(None); continue

                mode, tree, entries = interval_index[chrom]
                found_name = found_id = None
                if mode == 'ncls':
                    hits = list(tree.find_overlap(s, e))
                    if hits:
                        found_name = entries[hits[0][2]].get('gene_name')
                        found_id   = entries[hits[0][2]].get('gene_id')
                else:
                    for entry in entries:
                        if entry['start'] <= e and entry['end'] >= s:
                            found_name = entry.get('gene_name')
                            found_id   = entry.get('gene_id')
                            break
                gnames.append(found_name)
                gids.append(found_id)

            df['gene_name'] = gnames
            df['gene_id']   = gids
            gene_col = 'gene_name'

        elif gene_col and gtf_dict and 'by_gene_name' in gtf_dict:
            by_name = gtf_dict['by_gene_name']
            df['gene_id'] = df[gene_col].map(
                lambda g: by_name.get(str(g), {}).get('gene_id')
            )

    # -- Annotation gnomAD ----------------------------------------------------
    if gnomad_dict and gene_col in df.columns:
        for metric in ('pLI', 'oe_lof', 'lof_z', 'mis_z', 'syn_z',
                       'constraint_flag', 'oe_mis', 'oe_syn'):
            df[metric] = df[gene_col].map(
                lambda g, m=metric: gnomad_dict.get(str(g), {}).get(m)
            )

    # -- Annotation Mendeliome ------------------------------------------------
    if mendeliome_dict and gene_col in df.columns:
        for col in ('confidence_level', 'Mode_Of_Inheritance', 'Phenotypes'):
            df[col] = df[gene_col].map(
                lambda g, c=col: mendeliome_dict.get(str(g), {}).get(c)
            )

    # -- Selection et ordre des colonnes de sortie ----------------------------
    # Colonnes supprimées directement à l'écriture
    DROP_FRASER   = {'hgncSymbol'}
    DROP_OUTRIDER = {'gene_id'}

    if tool_name == 'fraser':
        ordered = [
            'seqnames', 'start', 'end', 'width', 'strand', 'sampleID',
            'type', 'pValue', 'padjust', 'psiValue', 'deltaPsi', 'counts',
            'totalCounts', 'meanCounts', 'meanTotalCounts', 'nonsplitCounts',
            'nonsplitProportion', 'nonsplitProportion_99quantile',
            'gene_name', 'gene_id', 'chrom',
            'pLI', 'oe_lof', 'lof_z', 'mis_z',
            'confidence_level', 'Mode_Of_Inheritance', 'Phenotypes',
        ]
        drop = DROP_FRASER
    else:
        ordered = [
            'geneID', 'sampleID', 'pValue', 'padjust', 'zScore', 'l2fc',
            'rawcounts', 'rawCounts', 'meanRawcounts', 'normcounts', 'meanCorrected',
            'theta', 'aberrant', 'AberrantBySample', 'AberrantByGene',
            'padj_rank', 'hgnc_symbol', 'gene_name', 'chrom', 'start', 'end', 'strand',
            'pLI', 'oe_lof', 'lof_z', 'mis_z',
            'confidence_level', 'Mode_Of_Inheritance', 'Phenotypes',
        ]
        drop = DROP_OUTRIDER

    out_cols = [c for c in ordered if c in df.columns and c not in drop]
    df[out_cols].to_csv(filepath, sep='\t', index=False)
    return filepath, len(df), sample_short


# =============================================================================
# Helpers de conversion dict (picklables)
# =============================================================================

def _gtf_to_dict(gtf_df):
    if gtf_df is None:
        return {}
    gtf_df = gtf_df.copy()
    gtf_df['gene_id_clean'] = gtf_df['gene_id'].str.split('.').str[0]

    # Index par gene_id_clean (pour OUTRIDER)
    by_gene = (
        gtf_df.set_index('gene_id_clean')
        [['gene_id', 'gene_name', 'chrom', 'start', 'end', 'strand']]
        .to_dict('index')
    )

    # Index par gene_name/symbole HGNC (pour FRASER : hgncSymbol → ENSG)
    gtf_nodup = gtf_df.drop_duplicates(subset='gene_name', keep='first')
    by_gene_name = (
        gtf_nodup.set_index('gene_name')
        [['gene_id', 'gene_id_clean', 'chrom', 'start', 'end', 'strand']]
        .to_dict('index')
    )
    # Interval index for coordinate-based FRASER annotation
    # Uses ncls if available, falls back to sorted list for linear scan
    by_interval = {}
    try:
        import numpy as np
        from ncls import NCLS
        for chrom, grp in gtf_df.groupby('chrom'):
            chrom = str(chrom)
            entries = grp[['start','end','gene_id','gene_name']].to_dict('records')
            starts = np.array([e['start'] for e in entries], dtype=np.int64)
            ends   = np.array([e['end']   for e in entries], dtype=np.int64)
            ids    = np.arange(len(entries), dtype=np.int64)
            by_interval[chrom] = ('list', None, entries)  # Force list mode for picklability
    except ImportError:
        # Linear scan fallback — group entries by chromosome
        for chrom, grp in gtf_df.groupby('chrom'):
            chrom = str(chrom)
            entries = grp[['start','end','gene_id','gene_name']].to_dict('records')
            by_interval[chrom] = ('list', None, entries)

    return {'by_gene': by_gene, 'by_gene_name': by_gene_name, 'by_interval': by_interval}


def _gnomad_to_dict(gnomad_df):
    """
    Convertit gnomAD en dict picklable : gene -> {pLI, oe_lof, lof_z, ...}
    Le renommage v4 (syn.z_score -> syn_z) est deja fait dans load_gnomad.
    drop_duplicates en filet de securite si canonical non filtre.
    """
    if gnomad_df is None:
        return {}
    wanted = ["pLI", "oe_lof", "lof_z", "mis_z", "syn_z",
              "constraint_flag", "oe_mis", "oe_syn"]
    cols = [c for c in wanted if c in gnomad_df.columns]
    df = gnomad_df[["gene"] + cols].copy()
    if "pLI" in df.columns:
        df = df.sort_values("pLI", ascending=False)
    df = df.drop_duplicates(subset="gene", keep="first")
    return df.set_index("gene")[cols].to_dict("index")


def _mendeliome_to_dict(mendel_df):
    if mendel_df is None:
        return {}
    cols = [c for c in ('confidence_level', 'Mode_Of_Inheritance', 'Phenotypes')
            if c in mendel_df.columns]
    return mendel_df.set_index('gene_symbol')[cols].to_dict('index')


# =============================================================================
# Processeur principal
# =============================================================================

class RNASeqProcessorPerSample:
    """Traite et annote les resultats FRASER2/OUTRIDER avec une sortie par sample."""

    def __init__(self, fraser_file, outrider_file, samples_file, gtf_file,
                 output_dir, gnomad_file=None, mendeliome_file=None,
                 mode='samples', pvalue_filter=None, create_zip=True,
                 partial_match=False, workers=DEFAULT_WORKERS):

        self.fraser_file     = Path(fraser_file)     if fraser_file     else None
        self.outrider_file   = Path(outrider_file)   if outrider_file   else None
        self.samples_file    = Path(samples_file)    if samples_file    else None
        self.gtf_file        = Path(gtf_file)
        self.gnomad_file     = Path(gnomad_file)     if gnomad_file     else None
        self.mendeliome_file = Path(mendeliome_file) if mendeliome_file else None

        self.output_dir    = Path(output_dir)
        self.mode          = mode
        self.pvalue_filter = pvalue_filter
        self.create_zip    = create_zip
        self.partial_match = partial_match
        self.workers       = workers

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.files_dir = self.output_dir / "per_sample_files"
        self.files_dir.mkdir(parents=True, exist_ok=True)

        self.samples         = None
        self.fraser_data     = None
        self.outrider_data   = None
        self.gtf_data        = None
        self.gnomad_data     = None
        self.mendeliome_data = None

        # Dicts picklables (calcules une seule fois avant le pool de process)
        self._gtf_dict        = None
        self._gnomad_dict     = None
        self._mendeliome_dict = None

    # -------------------------------------------------------------------------
    # Chargement
    # -------------------------------------------------------------------------

    def load_samples(self):
        if self.mode == 'all':
            logger.info("Mode : TOUS les samples")
            self.samples = None
            return None
        logger.info(f"Chargement des samples depuis {self.samples_file}")
        with open(self.samples_file) as f:
            self.samples = [l.strip() for l in f if l.strip()]
        logger.info(f"{len(self.samples)} samples charges")
        return self.samples

    def load_fraser(self):
        """
        Charge FRASER en normalisant les noms de colonnes selon la version du pipeline R.
        Aliases connus :
          padjValue  -> padjust
          hgncSymbol -> hgncSymbol  (déjà correct, utilisé tel quel)
        """
        if self.fraser_file is None:
            return None
        logger.info(f"Chargement FRASER : {self.fraser_file}")
        self.fraser_data = pd.read_csv(self.fraser_file, sep='\t', low_memory=False)

        col_aliases = {
            "padjValue": "padjust",
            "padj":      "padjust",
        }
        rename_map = {c: col_aliases[c] for c in self.fraser_data.columns
                      if c in col_aliases}
        if rename_map:
            self.fraser_data = self.fraser_data.rename(columns=rename_map)
            logger.info(f"  colonnes renommées : {rename_map}")

        logger.info(f"  -> {len(self.fraser_data):,} enregistrements")
        return self.fraser_data

    def load_outrider(self):
        """
        Charge OUTRIDER en détectant automatiquement la colonne sampleID.
        Certains fichiers R exportent une colonne d'index sans nom en position 0
        qui peut masquer ou déplacer sampleID — on la supprime si nécessaire.
        """
        if self.outrider_file is None:
            return None
        logger.info(f"Chargement OUTRIDER : {self.outrider_file}")

        # Lecture sans index_col pour ne rien masquer
        self.outrider_data = pd.read_csv(self.outrider_file, sep='\t', low_memory=False)

        # Supprimer la colonne d'index R (Unnamed: 0) si présente
        unnamed = [c for c in self.outrider_data.columns if c.startswith("Unnamed")]
        if unnamed:
            self.outrider_data = self.outrider_data.drop(columns=unnamed)
            logger.info(f"  colonne(s) d'index R supprimée(s) : {unnamed}")

        # Normaliser la casse de sampleID
        col_map = {c: "sampleID" for c in self.outrider_data.columns
                   if c.lower() == "sampleid" and c != "sampleID"}
        if col_map:
            self.outrider_data = self.outrider_data.rename(columns=col_map)
            logger.info(f"  colonne renommée : {col_map}")

        if "sampleID" not in self.outrider_data.columns:
            raise ValueError(
                f"Colonne sampleID introuvable. Colonnes disponibles : "
                f"{list(self.outrider_data.columns[:10])}"
            )

        # Normaliser les noms de colonnes OUTRIDER (variants selon version R/pipeline)
        col_aliases = {
            "padjValue":   "padjust",   # alias fréquent de padjust
            "padj":        "padjust",
            "rawCounts":   "rawcounts",
            "hgnc_symbol": "gene_name", # symbole HGNC → gene_name utilisé pour annotations
        }
        rename_map = {c: col_aliases[c] for c in self.outrider_data.columns
                      if c in col_aliases}
        if rename_map:
            self.outrider_data = self.outrider_data.rename(columns=rename_map)
            logger.info(f"  colonnes renommées : {rename_map}")

        logger.info(f"  -> {len(self.outrider_data):,} enregistrements")
        return self.outrider_data

    def load_gtf(self):
        logger.info(f"Chargement GTF : {self.gtf_file}")
        records = []
        with open(self.gtf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene':
                    continue
                attrs = {}
                for attr in fields[8].split(';'):
                    attr = attr.strip()
                    if not attr:
                        continue
                    kv = attr.split(' ', 1)
                    if len(kv) == 2:
                        attrs[kv[0]] = kv[1].strip('"')
                records.append({
                    'chrom':     fields[0],
                    'start':     int(fields[3]),
                    'end':       int(fields[4]),
                    'strand':    fields[6],
                    'gene_id':   attrs.get('gene_id'),
                    'gene_name': attrs.get('gene_name'),
                })
        self.gtf_data = pd.DataFrame(records)
        logger.info(f"  -> {len(self.gtf_data):,} genes")
        return self.gtf_data

    def load_gnomad(self):
        """
        Charge gnomAD v2 ou v4.
        - v4 : filtre sur canonical=True (une ligne par transcrit → une par gene)
               et renomme syn.z_score → syn_z pour uniformiser avec v2.
        - v2 : deja une ligne par gene, aucun traitement special.
        """
        if self.gnomad_file is None:
            return None
        logger.info(f"Chargement gnomAD : {self.gnomad_file}")
        self.gnomad_data = pd.read_csv(self.gnomad_file, sep='\t', low_memory=False)

        # gnomAD v4 : filtre canonical
        if "canonical" in self.gnomad_data.columns:
            n_before = len(self.gnomad_data)
            self.gnomad_data = self.gnomad_data[
                self.gnomad_data["canonical"] == True
            ].copy()
            logger.info(f"  gnomAD v4 detecte — filtre canonical : {n_before:,} -> {len(self.gnomad_data):,}")
            # Renommer syn.z_score -> syn_z pour uniformiser
            self.gnomad_data = self.gnomad_data.rename(columns={"syn.z_score": "syn_z"})
        else:
            logger.info(f"  gnomAD v2 detecte")

        logger.info(f"  -> {len(self.gnomad_data):,} genes")
        return self.gnomad_data

    def load_mendeliome(self):
        if self.mendeliome_file is None:
            return None
        if not self.mendeliome_file.exists():
            logger.warning(f"Mendeliome introuvable : {self.mendeliome_file}")
            return None
        logger.info(f"Chargement Mendeliome JSON : {self.mendeliome_file}")
        with open(self.mendeliome_file) as f:
            payload = json.load(f)
        version = payload.get("version", "?")
        records = []
        for entry in payload.get("genes", []):
            gd     = entry.get("gene_data", {})
            symbol = gd.get("gene_symbol") or gd.get("hgnc_symbol")
            if not symbol:
                continue
            phenotypes = entry.get("phenotypes", [])
            records.append({
                "gene_symbol":         symbol,
                "confidence_level":    str(entry.get("confidence_level", "")),
                "Mode_Of_Inheritance": entry.get("mode_of_inheritance", ""),
                "Phenotypes":          " | ".join(p for p in phenotypes if p),
            })
        self.mendeliome_data = (
            pd.DataFrame(records)
            .drop_duplicates(subset="gene_symbol")
        )
        n_green = (self.mendeliome_data["confidence_level"] == "3").sum()
        logger.info(
            f"  -> {len(self.mendeliome_data):,} genes uniques "
            f"(v{version}, {n_green} verts)"
        )
        return self.mendeliome_data

    def load_all_data(self):
        """Chargement parallele via ThreadPoolExecutor (I/O-bound)."""
        logger.info("Chargement parallele des fichiers...")
        tasks = {
            "GTF":        self.load_gtf,
            "gnomAD":     self.load_gnomad,
            "Mendeliome": self.load_mendeliome,
            "FRASER":     self.load_fraser,
            "OUTRIDER":   self.load_outrider,
        }
        errors = {}
        with ThreadPoolExecutor(max_workers=len(tasks)) as ex:
            futures = {ex.submit(fn): name for name, fn in tasks.items()}
            for future in as_completed(futures):
                name = futures[future]
                try:
                    future.result()
                    logger.info(f"  OK {name}")
                except Exception as e:
                    errors[name] = e
                    logger.error(f"  ERREUR {name} : {e}")

        if "GTF" in errors:
            raise RuntimeError(f"GTF indisponible : {errors['GTF']}")
        if errors:
            logger.warning(f"Donnees non chargees : {list(errors.keys())}")

        # Conversion en dicts picklables une seule fois pour tous les workers
        logger.info("Construction des dicts picklables pour ProcessPool...")
        self._gtf_dict        = _gtf_to_dict(self.gtf_data)
        self._gnomad_dict     = _gnomad_to_dict(self.gnomad_data)
        self._mendeliome_dict = _mendeliome_to_dict(self.mendeliome_data)
        logger.info(
            f"  GTF : {len(self._gtf_dict):,} | "
            f"gnomAD : {len(self._gnomad_dict):,} | "
            f"Mendeliome : {len(self._mendeliome_dict):,}"
        )

    # -------------------------------------------------------------------------
    # Filtrage samples
    # -------------------------------------------------------------------------

    @staticmethod
    def _short_id(s):
        """
        Extrait l'ID court (ex: 26D0643) depuis n'importe quel format :
          26D0643-MOINS                  -> 26D0643
          25D2693-STEMC-PUROMOINS-AVITI  -> 25D2693
          26D0643.MOINS                  -> 26D0643
          26D0643                        -> 26D0643
        """
        import re as _re
        return _re.split(r'[-_.]', str(s))[0]

    def _get_matched_samples(self, data_samples):
        """
        Résolution en deux passes :
          1. Match exact (noms identiques)
          2. Match normalisé : compare le short ID des deux côtés
             pour gérer le cas où samples.txt contient "26D0643-MOINS"
             mais sampleID dans le fichier annoté est "26D0643"
             (normalisé par annotation_outrider_hits_bis2.py).
        """
        if self.mode == 'all':
            return list(data_samples)

        # Dictionnaire short_id -> [full_name_in_data]
        data_by_short = {}
        for ds in data_samples:
            data_by_short.setdefault(self._short_id(ds), []).append(ds)

        # Short IDs des samples demandés
        requested_shorts = {self._short_id(s): s for s in self.samples}

        matched = []
        found_requested = set()

        for short, data_full_list in data_by_short.items():
            if short in requested_shorts:
                matched.extend(data_full_list)
                found_requested.add(requested_shorts[short])

        # Conserver aussi un match exact pour les cas où partial_match est activé
        if self.partial_match:
            for ds in data_samples:
                if ds not in matched:
                    if any(ls in ds for ls in self.samples):
                        matched.append(ds)

        logger.info(
            f"{len(matched)} samples correspondants trouvés "
            f"({len(found_requested)}/{len(self.samples or [])} demandés)"
        )

        not_found = set(self.samples or []) - found_requested
        if not_found:
            logger.warning(
                f"Samples non trouvés dans les données "
                f"(vérifier nom/normalisation) : {sorted(not_found)}"
            )
        return matched

    def _filter_data(self, data, label):
        if self.mode != 'all':
            matched = self._get_matched_samples(data['sampleID'].unique())
            data = data[data['sampleID'].isin(matched)].copy()
        if self.pvalue_filter is not None:
            # Chercher la colonne padjust (ou alias)
            padj_col = next(
                (c for c in ('padjust', 'padjValue', 'padj') if c in data.columns),
                None
            )
            if padj_col:
                n_before = len(data)
                data = data[data[padj_col] < self.pvalue_filter].copy()
                logger.info(f"Filtre p-value ({label}, col={padj_col}) : {n_before:,} -> {len(data):,}")
            else:
                logger.warning(f"Colonne padjust introuvable pour {label}, filtre p-value ignoré")
        logger.info(f"{label} filtre : {len(data):,} enregistrements")
        return data

    # -------------------------------------------------------------------------
    # Traitement parallele par sample (coeur du pipeline)
    # -------------------------------------------------------------------------

    def _run_tool_parallel(self, data, tool_name, gene_col):
        """
        Repartit les samples sur self.workers process.
        Chaque process recoit les donnees d'UN sample + les dicts de reference
        et effectue annotation + ecriture de facon totalement independante.
        -> Vraie parallelisation CPU : pas de GIL, pas de partage memoire.
        """
        unique_samples = data['sampleID'].unique()
        logger.info(
            f"Traitement {tool_name} : {len(unique_samples)} samples "
            f"sur {self.workers} workers"
        )

        tasks = [
            (
                sample,
                data[data['sampleID'] == sample].to_dict('list'),
                tool_name,
                str(self.files_dir),
                self._gtf_dict,
                self._gnomad_dict,
                self._mendeliome_dict,
                gene_col,
            )
            for sample in unique_samples
        ]

        saved_files = []
        with ProcessPoolExecutor(max_workers=self.workers) as executor:
            futures = {executor.submit(_process_and_save_sample, t): t[0]
                       for t in tasks}
            for future in as_completed(futures):
                sample = futures[future]
                try:
                    filepath, n_records, short = future.result()
                    saved_files.append(filepath)
                    logger.info(f"  OK {filepath.name} ({n_records} lignes)")
                except Exception as exc:
                    logger.error(f"  ERREUR {sample} : {exc}")

        logger.info(f"{len(saved_files)}/{len(unique_samples)} fichiers {tool_name} crees")
        return saved_files

    def process_fraser(self):
        if self.fraser_data is None:
            return []
        logger.info("--- FRASER ---")
        data = self._filter_data(self.fraser_data, 'FRASER')
        # Priorité : hgncSymbol (FRASER natif) > gene_name (après résolution GTF)
        gene_col = next(
            (c for c in ('hgncSymbol', 'gene_name') if c in data.columns),
            None
        )
        if gene_col is None:
            logger.warning("Aucune colonne gene trouvée pour FRASER, annotations désactivées")
            gene_col = 'hgncSymbol'  # fallback, sera absent → pas d'annotation
        logger.info(f"  gene_col FRASER : {gene_col}")
        return self._run_tool_parallel(data, 'fraser', gene_col)

    def process_outrider(self):
        if self.outrider_data is None:
            return []
        logger.info("--- OUTRIDER ---")
        data = self._filter_data(self.outrider_data, 'OUTRIDER')
        # Priorité : gene_name (ex hgnc_symbol renommé) > geneID
        gene_col = 'gene_name' if 'gene_name' in data.columns else 'geneID'
        logger.info(f"  gene_col OUTRIDER : {gene_col}")
        return self._run_tool_parallel(data, 'outrider', gene_col)

    # -------------------------------------------------------------------------
    # ZIP
    # -------------------------------------------------------------------------

    def create_zip_archive(self, files_list):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        zip_path  = self.output_dir / f"run_{timestamp}.zip"
        logger.info(f"Creation ZIP : {zip_path.name} ({len(files_list)} fichiers)")
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
            for fp in files_list:
                zf.write(fp, Path(fp).name)
        logger.info(f"ZIP cree : {zip_path}")
        return zip_path

    # -------------------------------------------------------------------------
    # Point d'entree
    # -------------------------------------------------------------------------

    def run(self):
        logger.info("=" * 60)
        logger.info("Pipeline RNA-Seq — sortie par sample")
        logger.info(f"Workers : {self.workers} | Mode : {self.mode.upper()}")
        logger.info("=" * 60)

        # Etape 1 : chargement parallele (threads, I/O-bound)
        self.load_samples()
        self.load_all_data()

        # Etape 2 : annotation + ecriture par sample (process, CPU-bound)
        fraser_files   = self.process_fraser()
        outrider_files = self.process_outrider()

        all_files = fraser_files + outrider_files

        if self.create_zip and all_files:
            zip_path = self.create_zip_archive(all_files)
            logger.info("=" * 60)
            logger.info("Pipeline termine !")
            logger.info(f"ZIP : {zip_path}")
            logger.info("=" * 60)
            return zip_path

        logger.info("=" * 60)
        logger.info("Pipeline termine !")
        logger.info(f"Fichiers : {self.files_dir}")
        logger.info("=" * 60)
        return all_files


# =============================================================================
# CLI autonome
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Pipeline RNA-Seq — sortie par sample (FRASER2 + OUTRIDER)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Exemples :

  # Samples specifiques + ZIP
  %(prog)s --fraser fraser.tab --outrider outrider.tab --samples samples.txt \\
           --gtf genes.gtf --output results/

  # Tous les samples, filtre p-value, 8 workers
  %(prog)s --fraser fraser.tab --outrider outrider.tab --gtf genes.gtf \\
           --output results/ --mode all --pvalue 0.05 --workers 8

  # Avec gnomAD et Mendeliome Australia (JSON)
  %(prog)s --fraser fraser.tab --outrider outrider.tab --samples samples.txt \\
           --gtf genes.gtf --gnomad gnomad.txt \\
           --mendeliome references/mendeliome_australia.json \\
           --partial-match --output results/

Workers par defaut : {DEFAULT_WORKERS} (base sur les CPU disponibles)
        """
    )

    parser.add_argument('--fraser',   help='Fichier FRASER2 (TSV)')
    parser.add_argument('--outrider', help='Fichier OUTRIDER (TSV)')
    parser.add_argument('--gtf',    required=True, help='Fichier GTF GENCODE')
    parser.add_argument('--output', required=True, help='Dossier de sortie')

    parser.add_argument('--mode', choices=['samples', 'all'], default='samples')
    parser.add_argument('--samples', help='Fichier liste de samples')
    parser.add_argument('--partial-match', action='store_true')
    parser.add_argument('--pvalue', type=float)
    parser.add_argument('--gnomad')
    parser.add_argument('--mendeliome')
    parser.add_argument('--no-zip', action='store_true')
    parser.add_argument('--workers', type=int, default=DEFAULT_WORKERS)
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()

    if args.mode == 'samples' and not args.samples:
        parser.error("--samples est requis avec --mode samples")
    if not args.fraser and not args.outrider:
        parser.error("Au moins --fraser ou --outrider est requis")

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    processor = RNASeqProcessorPerSample(
        fraser_file=args.fraser,
        outrider_file=args.outrider,
        samples_file=args.samples,
        gtf_file=args.gtf,
        output_dir=args.output,
        gnomad_file=args.gnomad,
        mendeliome_file=args.mendeliome,
        mode=args.mode,
        pvalue_filter=args.pvalue,
        create_zip=not args.no_zip,
        partial_match=args.partial_match,
        workers=args.workers,
    )

    try:
        processor.run()
        sys.exit(0)
    except Exception as e:
        logger.error(f"Echec : {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
