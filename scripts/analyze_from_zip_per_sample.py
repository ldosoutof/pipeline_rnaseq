#!/usr/bin/env python3
"""
Auto-detect and process RNA-Seq files from a ZIP archive
Generates one output file per sample per tool
Includes automatic download/update of reference files (GENCODE GTF, gnomAD, Mendeliome Australia)

Parallel strategy:
  - Reference downloads   → ThreadPoolExecutor  (I/O réseau, pas de GIL)
  - Mendeliome pagination → ThreadPoolExecutor  (I/O réseau)
  - Sample processing     → ProcessPoolExecutor (CPU-bound, contourne le GIL)
"""

import argparse
import zipfile
import tempfile
import shutil
from pathlib import Path
import logging
import sys
import re
import json
import os
import urllib.request
import urllib.error
import subprocess
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Nombre de workers par défaut = nombre de CPU logiques
DEFAULT_WORKERS = min(os.cpu_count() or 4, 8)

# ─────────────────────────────────────────────────────────────
# REFERENCE DOWNLOAD
# ─────────────────────────────────────────────────────────────

REFERENCES_DIR = Path("references")

GENCODE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
    "/release_44/gencode.v44.annotation.gtf.gz"
)
# ← UPDATED: gnomAD v4.1 (plain TSV, no decompression needed)
GNOMAD_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad"
    "/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
)
PANELAPP_AUS_PANEL_URL = "https://panelapp-aus.org/api/v1/panels/137/"
PANELAPP_AUS_GENES_URL = "https://panelapp-aus.org/api/v1/panels/137/genes/"

GTF_FILENAME    = "gencode.v44.annotation.gtf"
# ← UPDATED: gnomAD v4.1 filename
GNOMAD_FILENAME = "gnomad.v4.1.constraint_metrics.tsv"
MENDEL_FILENAME = "mendeliome_australia.json"


def _progress_hook(block_num, block_size, total_size):
    downloaded = block_num * block_size
    if total_size > 0:
        pct = min(downloaded / total_size * 100, 100)
        mb = downloaded / 1_048_576
        print(f"\r   {pct:5.1f}%  ({mb:.1f} MB)", end="", flush=True)


def _wget_download(url: str, dest: Path):
    """Download with wget if available, else urllib fallback."""
    try:
        subprocess.run(
            ["wget", "-q", "--show-progress", "-O", str(dest), url],
            check=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.info("   wget non disponible, utilisation de urllib...")
        urllib.request.urlretrieve(url, dest, reporthook=_progress_hook)
        print()


def _gunzip(src: Path, dst: Path):
    """Decompress a .gz file."""
    import gzip
    logger.info(f"   Décompression de {src.name}...")
    with gzip.open(src, 'rb') as f_in, open(dst, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    src.unlink()


def download_gencode(ref_dir: Path) -> Path:
    gtf = ref_dir / GTF_FILENAME
    if gtf.exists():
        logger.info("✓ GTF GENCODE déjà présent")
        return gtf

    gz = ref_dir / (GTF_FILENAME + ".gz")
    logger.info("   Téléchargement GTF GENCODE v44 (~50 MB compressé)...")
    _wget_download(GENCODE_URL, gz)
    _gunzip(gz, gtf)
    logger.info(f"✓ GTF téléchargé ({gtf})")
    return gtf


# ← UPDATED: gnomAD v4.1 is a plain TSV — no gunzip step needed
def download_gnomad(ref_dir: Path) -> Path:
    tsv = ref_dir / GNOMAD_FILENAME
    if tsv.exists():
        logger.info("✓ gnomAD déjà présent")
        return tsv

    logger.info("   Téléchargement gnomAD v4.1 (~10 MB)...")
    _wget_download(GNOMAD_URL, tsv)
    logger.info(f"✓ gnomAD téléchargé ({tsv})")
    return tsv


def _api_get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode())


def download_mendeliome(ref_dir: Path) -> Path:
    local = ref_dir / MENDEL_FILENAME

    try:
        panel_info = _api_get(PANELAPP_AUS_PANEL_URL)
        remote_version = panel_info.get("version", "")
    except Exception as e:
        logger.warning(f"   Impossible de contacter PanelApp Australia: {e}")
        if local.exists():
            logger.warning("   Utilisation du fichier local existant.")
            return local
        else:
            raise RuntimeError("Pas de fichier Mendeliome local et API inaccessible.")

    logger.info(f"   Version distante Mendeliome : {remote_version}")

    local_version = None
    if local.exists():
        try:
            with open(local) as f:
                local_version = json.load(f).get("version")
        except Exception:
            pass
    logger.info(f"   Version locale  Mendeliome : {local_version or 'aucune'}")

    if remote_version and remote_version == local_version:
        logger.info("✓ Mendeliome déjà à jour")
        return local

    logger.info(f"   Nouvelle version ({local_version or '—'} → {remote_version}). Téléchargement...")

    logger.info("   Collecte des URLs de pages...")
    page_urls = []
    url = PANELAPP_AUS_GENES_URL
    while url:
        page_urls.append(url)
        data = _api_get(url)
        url = data.get("next")

    logger.info(f"   {len(page_urls)} pages à télécharger en parallèle...")

    pages_results = {}
    with ThreadPoolExecutor(max_workers=min(len(page_urls), 8)) as executor:
        future_to_idx = {
            executor.submit(_api_get, page_url): i
            for i, page_url in enumerate(page_urls)
        }
        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            pages_results[idx] = future.result()["results"]
            print(f"\r   {len(pages_results)}/{len(page_urls)} pages récupérées...", end="", flush=True)
    print()

    genes = []
    for i in sorted(pages_results):
        genes.extend(pages_results[i])

    payload = {
        "panel_id": panel_info["id"],
        "name":     panel_info["name"],
        "version":  remote_version,
        "stats":    panel_info.get("stats", {}),
        "genes":    genes,
    }
    with open(local, "w") as f:
        json.dump(payload, f, indent=2)

    logger.info(f"✓ Mendeliome téléchargé — {len(genes)} gènes (v{remote_version})")
    return local


def setup_references(ref_dir: Path = REFERENCES_DIR) -> dict:
    """
    Télécharge (si nécessaire) GTF GENCODE, gnomAD et le Mendeliome Australia.
    GTF et gnomAD sont téléchargés en parallèle (threads I/O réseau).
    Retourne un dict avec les chemins vers chaque fichier.
    """
    ref_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Vérification / téléchargement des fichiers de référence")
    logger.info("=" * 60)

    results = {}
    errors = {}

    download_tasks = {
        "gtf":    (download_gencode, ref_dir),
        "gnomad": (download_gnomad,  ref_dir),
    }

    logger.info("\n[1-2/3] GTF GENCODE + gnomAD (téléchargement parallèle)...")
    with ThreadPoolExecutor(max_workers=2) as executor:
        future_to_key = {
            executor.submit(fn, ref_dir): key
            for key, (fn, ref_dir) in download_tasks.items()
        }
        for future in as_completed(future_to_key):
            key = future_to_key[future]
            try:
                results[key] = future.result()
            except Exception as e:
                errors[key] = e
                logger.error(f"Erreur téléchargement {key}: {e}")

    if errors:
        raise RuntimeError(f"Échec téléchargement : {errors}")

    logger.info("\n[3/3] Mendeliome Australia (PanelApp)")
    results["mendeliome"] = download_mendeliome(ref_dir)

    logger.info("\n" + "=" * 60)
    logger.info("Références prêtes :")
    for k, v in results.items():
        logger.info(f"  {k:12s} → {v}")
    logger.info("=" * 60)

    return results


# ─────────────────────────────────────────────────────────────
# ZIP ANALYZER
# ─────────────────────────────────────────────────────────────

class ZipAnalyzer:
    """Analyze and extract RNA-Seq files from ZIP archives"""

    def __init__(self, zip_path):
        self.zip_path = Path(zip_path)
        self.temp_dir = None
        self.fraser_file = None
        self.outrider_file = None

    def extract_zip(self):
        logger.info(f"Extraction de l'archive ZIP : {self.zip_path}")
        self.temp_dir = Path(tempfile.mkdtemp(prefix="rnaseq_zip_"))
        with zipfile.ZipFile(self.zip_path, 'r') as zip_ref:
            zip_ref.extractall(self.temp_dir)
        logger.info(f"Fichiers extraits dans : {self.temp_dir}")
        return self.temp_dir

    def detect_fraser_file(self):
        logger.info("Détection du fichier FRASER...")
        patterns = [
            r'fraser.*\.tab$', r'fraser.*\.tsv$',
            r'.*fraser.*\.txt$', r'fraser[_-]results.*',
        ]
        for root, dirs, files in self.temp_dir.walk():
            for file in files:
                for pattern in patterns:
                    if re.search(pattern, file.lower(), re.IGNORECASE):
                        candidate = root / file
                        if self._verify_fraser_format(candidate):
                            logger.info(f"✓ FRASER : {candidate.relative_to(self.temp_dir)}")
                            self.fraser_file = candidate
                            return candidate
        logger.warning("Fichier FRASER non détecté automatiquement")
        return None

    def detect_outrider_file(self):
        logger.info("Détection du fichier OUTRIDER...")
        patterns = [
            r'outrider.*\.tab$', r'outrider.*\.tsv$',
            r'.*outrider.*\.txt$', r'outrider[_-]results.*', r'.*htseq.*',
        ]
        for root, dirs, files in self.temp_dir.walk():
            for file in files:
                for pattern in patterns:
                    if re.search(pattern, file.lower(), re.IGNORECASE):
                        candidate = root / file
                        if self._verify_outrider_format(candidate):
                            logger.info(f"✓ OUTRIDER : {candidate.relative_to(self.temp_dir)}")
                            self.outrider_file = candidate
                            return candidate
        logger.warning("Fichier OUTRIDER non détecté automatiquement")
        return None

    def _verify_fraser_format(self, file_path):
        try:
            with open(file_path) as f:
                header = f.readline().strip().lower()
                return all(c in header for c in ['sampleid', 'hgncsymbol', 'pvalue', 'deltapsi'])
        except Exception as e:
            logger.debug(f"Erreur vérification FRASER {file_path}: {e}")
            return False

    def _verify_outrider_format(self, file_path):
        try:
            with open(file_path) as f:
                header = f.readline().strip().lower()
                return all(c in header for c in ['geneid', 'sampleid', 'zscore', 'pvalue'])
        except Exception as e:
            logger.debug(f"Erreur vérification OUTRIDER {file_path}: {e}")
            return False

    def list_all_files(self):
        logger.info("Fichiers présents dans le ZIP :")
        files = []
        for root, dirs, filenames in self.temp_dir.walk():
            for filename in filenames:
                file_path = root / filename
                files.append((file_path.relative_to(self.temp_dir), file_path.stat().st_size))
        files.sort(key=lambda x: x[1], reverse=True)
        for rel_path, size in files[:20]:
            logger.info(f"  {rel_path} ({size / 1_048_576:.2f} MB)")
        if len(files) > 20:
            logger.info(f"  ... et {len(files) - 20} autres fichiers")
        return files

    def manual_file_selection(self):
        data_files = [f for f in self.temp_dir.rglob('*.*')
                      if f.suffix.lower() in ['.tsv', '.tab', '.txt']]
        if not data_files:
            logger.error("Aucun fichier de données trouvé dans le ZIP")
            return False

        print("\n" + "=" * 60)
        print("Sélection manuelle des fichiers")
        print("=" * 60)
        for i, file in enumerate(data_files, 1):
            size = file.stat().st_size / 1_048_576
            print(f"{i:2d}. {file.relative_to(self.temp_dir)} ({size:.2f} MB)")

        for attr, label in [('fraser_file', 'FRASER'), ('outrider_file', 'OUTRIDER')]:
            if not getattr(self, attr):
                print(f"\nQuel fichier est la sortie {label} ?")
                try:
                    choice = int(input("Numéro (0 pour ignorer) : "))
                    if 1 <= choice <= len(data_files):
                        setattr(self, attr, data_files[choice - 1])
                        logger.info(f"Fichier {label} : {getattr(self, attr)}")
                except (ValueError, KeyboardInterrupt):
                    logger.warning(f"Sélection {label} ignorée")

        return bool(self.fraser_file or self.outrider_file)

    def cleanup(self):
        if self.temp_dir and self.temp_dir.exists():
            logger.info(f"Nettoyage du répertoire temporaire : {self.temp_dir}")
            shutil.rmtree(self.temp_dir)


# ─────────────────────────────────────────────────────────────
# PIPELINE
# ─────────────────────────────────────────────────────────────

def run_pipeline(fraser_file, outrider_file, samples_file, gtf_file, output_dir,
                 gnomad_file=None, mendeliome_file=None,
                 mode='samples', pvalue_filter=None,
                 verbose=False, create_zip=True, partial_match=False,
                 workers=1):
    """Lance le pipeline d'analyse RNA-Seq avec sortie par sample."""
    logger.info("Lancement du pipeline RNA-Seq (sortie par sample)...")

    import rnaseq_analysis_per_sample

    processor = rnaseq_analysis_per_sample.RNASeqProcessorPerSample(
        fraser_file=fraser_file,
        outrider_file=outrider_file,
        samples_file=samples_file,
        gtf_file=gtf_file,
        output_dir=output_dir,
        gnomad_file=gnomad_file,
        mendeliome_file=mendeliome_file,
        mode=mode,
        pvalue_filter=pvalue_filter,
        create_zip=create_zip,
        partial_match=partial_match,
        workers=workers,
    )

    return processor.run()


# ─────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Analyse RNA-Seq depuis une archive ZIP (sortie par sample).\n"
            "Télécharge automatiquement les références si nécessaire."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemples :

  # Télécharger / mettre à jour les références uniquement
  %(prog)s --download-refs

  # Tous les samples, références auto-téléchargées
  %(prog)s --zip results.zip --mode all --output results/

  # Samples spécifiques depuis un fichier
  %(prog)s --zip results.zip --samples samples.txt --output results/

  # Références déjà présentes dans un dossier custom
  %(prog)s --zip results.zip --mode all --output results/ --refs-dir /data/references

  # Avec seuil de p-value et correspondance partielle des noms
  %(prog)s --zip results.zip --samples samples.txt --output results/ \\
           --pvalue 0.05 --partial-match
        """
    )

    io_group = parser.add_argument_group("Entrée / Sortie")
    io_group.add_argument('--zip',
                          help='Archive ZIP contenant les fichiers FRASER et/ou OUTRIDER')
    io_group.add_argument('--output',
                          help='Dossier de sortie')

    ref_group = parser.add_argument_group("Références")
    ref_group.add_argument('--download-refs', action='store_true',
                           help='Télécharger / mettre à jour les références puis quitter')
    ref_group.add_argument('--refs-dir', default='references',
                           help='Dossier des fichiers de référence (défaut : ./references)')
    ref_group.add_argument('--gtf',
                           help='Chemin GTF custom (désactive le téléchargement auto du GTF)')
    ref_group.add_argument('--gnomad',
                           help='Chemin gnomAD custom (désactive le téléchargement auto)')
    ref_group.add_argument('--mendeliome',
                           help='Chemin Mendeliome JSON custom (désactive le téléchargement auto)')

    proc_group = parser.add_argument_group("Mode de traitement")
    proc_group.add_argument('--mode', choices=['samples', 'all'], default='samples',
                            help=(
                                '"samples" : traiter uniquement les samples du fichier --samples  '
                                '(défaut) | "all" : traiter tous les samples du ZIP'
                            ))
    proc_group.add_argument('--samples',
                            help='Fichier listant les samples à traiter (requis si --mode samples)')
    proc_group.add_argument('--partial-match', action='store_true',
                            help='Correspondance partielle des noms de samples '
                                 '(ex. "23D1192" → "23D1192.HOL.Hay")')

    filter_group = parser.add_argument_group("Filtres")
    filter_group.add_argument('--pvalue', type=float,
                              help='Seuil de p-value ajustée (ex. 0.05)')

    adv_group = parser.add_argument_group("Options avancées")
    adv_group.add_argument('--no-zip', action='store_true',
                           help='Ne pas créer de ZIP de sortie (conserver les fichiers individuels)')
    adv_group.add_argument('--interactive', action='store_true',
                           help='Mode interactif si la détection automatique échoue')
    adv_group.add_argument('--verbose', action='store_true',
                           help='Logging détaillé')
    adv_group.add_argument('--workers', type=int, default=1,
                           metavar='N',
                           help=(
                               f'Nombre de workers parallèles pour le traitement des samples '
                               f'(défaut : 1 = séquentiel). '
                               f'Utilise des process réels (ProcessPoolExecutor) pour contourner '
                               f'le GIL Python. Recommandé : nombre de CPU - 1 '
                               f'(détecté : {DEFAULT_WORKERS})'
                           ))

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    refs_dir = Path(args.refs_dir)

    if args.download_refs:
        setup_references(refs_dir)
        logger.info("Références à jour. Sortie.")
        sys.exit(0)

    if not args.zip:
        parser.error("--zip est requis (sauf avec --download-refs)")
    if not args.output:
        parser.error("--output est requis")
    if args.mode == 'samples' and not args.samples:
        parser.error("--samples est requis avec --mode samples")

    need_auto = not (args.gtf and args.gnomad and args.mendeliome)

    if need_auto:
        refs = setup_references(refs_dir)
    else:
        refs = {}

    gtf_file        = Path(args.gtf)        if args.gtf        else refs["gtf"]
    gnomad_file     = Path(args.gnomad)     if args.gnomad     else refs["gnomad"]
    mendeliome_file = Path(args.mendeliome) if args.mendeliome else refs["mendeliome"]

    logger.info("=" * 60)
    logger.info("Pipeline RNA-Seq — Auto-détection depuis ZIP")
    logger.info(f"Mode          : {args.mode.upper()}")
    if args.mode == 'samples':
        logger.info(f"Fichier samples: {args.samples}")
    logger.info(f"GTF           : {gtf_file}")
    logger.info(f"gnomAD        : {gnomad_file}")
    logger.info(f"Mendeliome    : {mendeliome_file}")
    if args.pvalue:
        logger.info(f"Filtre p-value: < {args.pvalue}")
    logger.info(f"Workers       : {args.workers} {'(séquentiel)' if args.workers == 1 else '(parallèle — ProcessPoolExecutor)'}")
    logger.info("=" * 60)

    analyzer = ZipAnalyzer(args.zip)

    try:
        analyzer.extract_zip()
        analyzer.list_all_files()
        analyzer.detect_fraser_file()
        analyzer.detect_outrider_file()

        if args.interactive and not (analyzer.fraser_file and analyzer.outrider_file):
            logger.info("Détection incomplète — passage en mode interactif...")
            if not analyzer.manual_file_selection():
                logger.error("Aucun fichier sélectionné. Abandon.")
                sys.exit(1)

        if not analyzer.fraser_file and not analyzer.outrider_file:
            logger.error("Impossible de détecter les fichiers FRASER ou OUTRIDER.")
            logger.error("Utilisez --interactive pour la sélection manuelle.")
            sys.exit(1)

        result = run_pipeline(
            fraser_file=analyzer.fraser_file,
            outrider_file=analyzer.outrider_file,
            samples_file=args.samples,
            gtf_file=gtf_file,
            output_dir=args.output,
            gnomad_file=gnomad_file,
            mendeliome_file=mendeliome_file,
            mode=args.mode,
            pvalue_filter=args.pvalue,
            verbose=args.verbose,
            create_zip=not args.no_zip,
            partial_match=args.partial_match,
            workers=args.workers,
        )

        logger.info("=" * 60)
        logger.info("✓ Pipeline terminé avec succès !")
        if isinstance(result, Path):
            logger.info(f"✓ ZIP de sortie : {result}")
        else:
            logger.info(f"✓ {len(result)} fichiers individuels créés")
        logger.info("=" * 60)
        sys.exit(0)

    except Exception as e:
        logger.error(f"Échec du pipeline : {e}", exc_info=True)
        sys.exit(1)

    finally:
        analyzer.cleanup()


if __name__ == '__main__':
    main()
