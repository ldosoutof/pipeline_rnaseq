#!/usr/bin/env python3
"""
Watcher OVH — détecte l'arrivée complète des FASTQ (fichier sentinelle déposé
par le watcher sitatst) et lance le pipeline RNA-seq automatiquement.

Flux global :
  [sitatst] watcher_nextseq.py : concatène + rsync FASTQ -> OVH, puis dépose
            un fichier sentinelle TRANSFER_COMPLETE dans le dossier du run.
  [OVH]     CE watcher : détecte TRANSFER_COMPLETE -> preprocess.py -> launch.sh,
            puis, si le pipeline a réussi, renvoie pipeline_v0/ vers sitatst
            (scripts/send_results_to_sitatst.py) et y dépose RESULTS_COMPLETE.

Points clés :
  - Détection par SENTINELLE (et non présence de fichiers) : garantit que le
    transfert est terminé (le sentinelle est déposé en dernier, après rsync +
    contrôle d'intégrité côté sitatst).
  - SÉRIALISATION : un seul run à la fois (les intermédiaires du pipeline sont
    partagés ; deux runs concurrents les corrompent). Un verrou global empêche
    tout lancement parallèle.
  - Idempotence : un run déjà traité (marqueur .launched) n'est pas relancé.
  - Retour des résultats : marqueurs .pipeline_ok / .pipeline_failed (issue du
    pipeline) et .results_synced (transfert vérifié). Un run réussi mais non
    transféré (sitatst injoignable…) est retenté au plus toutes les
    RESULTS_RETRY_DELAY secondes.
"""

import os
import sys
import time
import subprocess
import logging
from pathlib import Path
from datetime import datetime

# ──────────────────────────────────────────────
# CONFIGURATION (à adapter au serveur OVH)
# ──────────────────────────────────────────────
# Répertoire où le watcher sitatst dépose les runs (FASTQ + sentinelle).
# NB : vérifier datawork vs datawork2 selon l'installation réelle.
WATCH_DIR       = "/datawork2/genetique/RNASeq/diag/prod"

# Racine du pipeline (dossier contenant scripts/, rules/, snakemake/…)
PIPELINE_DIR    = "/home/ldosoutoferreira/pipeline/RNASEQ/routine/test_pipeline_v2/pipeline_corrige/pipeline_fix/pipeline"

# Scripts du pipeline
PREPROCESS      = f"{PIPELINE_DIR}/scripts/preprocess.py"
# site_paths.yml rempli (chemins réels de l'installation OVH)
SITE_PATHS      = f"{PIPELINE_DIR}/template/site_paths.yml"

# Nom du fichier sentinelle (doit correspondre à SENTINEL_NAME du watcher sitatst)
SENTINEL_NAME   = "TRANSFER_COMPLETE"
# Marqueur déposé par CE watcher une fois un run lancé (évite un double lancement)
LAUNCHED_MARKER = ".launched"

POLL_INTERVAL   = 300           # secondes entre deux scans (5 min)

# Retour des résultats vers sitatst (paramètres : section sync_sitatst de SITE_PATHS)
SEND_RESULTS        = f"{PIPELINE_DIR}/scripts/send_results_to_sitatst.py"
PIPELINE_OK         = ".pipeline_ok"       # pipeline terminé avec le code 0
PIPELINE_FAILED     = ".pipeline_failed"   # pipeline en échec (code dans le fichier)
RESULTS_MARKER      = ".results_synced"    # écrit par send_results_to_sitatst.py
RESULTS_RETRY_DELAY = 1800                 # s entre deux tentatives de retour (30 min)

LOG_FILE        = "/home/ldosoutoferreira/logs/watcher_ovh.log"
# ──────────────────────────────────────────────

os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(LOG_FILE)],
)
log = logging.getLogger(__name__)


def pipeline_is_running() -> bool:
    """
    Vrai si un pipeline Snakemake tourne déjà (sérialisation : un seul run à la
    fois). On cherche un processus snakemake actif. Robuste et sans dépendance :
    lit /proc plutôt que d'appeler pgrep (évite un faux positif sur soi-même).
    """
    me = os.getpid()
    for pid_dir in Path("/proc").iterdir():
        if not pid_dir.name.isdigit():
            continue
        pid = int(pid_dir.name)
        if pid == me:
            continue
        try:
            cmdline = (pid_dir / "cmdline").read_bytes().replace(b"\0", b" ").decode(errors="replace")
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
        # snakemake lancé via launch.sh -> la ligne de commande contient "snakemake"
        if "snakemake" in cmdline and "pipeline.smk" in cmdline:
            return True
    return False


def find_ready_runs(watch_path: Path) -> list:
    """
    Renvoie la liste des runs PRÊTS à être lancés :
      - contient le sentinelle SENTINEL_NAME (transfert complet)
      - ne contient PAS le marqueur LAUNCHED_MARKER (pas déjà lancé)
    Triés par ordre de nom (chronologique via l'horodatage du nom de run).
    """
    ready = []
    for run_dir in sorted(watch_path.iterdir()):
        if not run_dir.is_dir():
            continue
        sentinel = run_dir / SENTINEL_NAME
        launched = run_dir / LAUNCHED_MARKER
        if sentinel.exists() and not launched.exists():
            ready.append(run_dir)
    return ready


def run_preprocess(run_dir: Path) -> bool:
    """
    Lance preprocess.py pour générer config.yml + launch.sh.
      --path    = racine pipeline
      --workDir = dossier du run (résultats y seront écrits)
      --dataDir = dossier du run (contient fastq/)
    """
    cmd = [
        "python", PREPROCESS,
        "--path",    PIPELINE_DIR,
        "--workDir", str(run_dir),
        "--dataDir", str(run_dir),
    ]
    if Path(SITE_PATHS).is_file():
        cmd += ["--site-paths", SITE_PATHS]

    log.info(f"preprocess.py : {' '.join(cmd)}")
    log_out = run_dir / "preprocess.log"
    try:
        with open(log_out, "w") as fout:
            ret = subprocess.run(cmd, stdout=fout, stderr=subprocess.STDOUT, timeout=1800)
        if ret.returncode != 0:
            log.error(f"preprocess.py a échoué (code {ret.returncode}). Voir {log_out}")
            return False
        # Vérifier que le config généré ne contient plus de placeholders /path/to/
        config = run_dir / "launch_folder" / "config.yml"
        if config.is_file():
            txt = config.read_text()
            if "/path/to/" in txt:
                log.error(f"config.yml contient encore des /path/to/ — "
                          f"site_paths.yml incomplet. Lancement annulé.")
                return False
        else:
            log.error(f"config.yml introuvable après preprocess ({config}).")
            return False
        log.info("preprocess.py OK (config.yml + launch.sh générés).")
        return True
    except Exception as e:
        log.error(f"Erreur preprocess.py : {e}")
        return False


def run_pipeline(run_dir: Path) -> None:
    """
    Lance launch.sh (le pipeline). BLOQUANT : cette fonction attend la fin du
    pipeline (elle est appelée dans un contexte déjà sérialisé). Le run est
    marqué .launched AVANT le lancement pour éviter tout double départ même si
    le pipeline est long.
    """
    launch = run_dir / "launch_folder" / "launch.sh"
    if not launch.is_file():
        log.error(f"launch.sh introuvable ({launch}) — run ignoré.")
        return

    # Marquer AVANT de lancer (idempotence : pas de second départ)
    marker = run_dir / LAUNCHED_MARKER
    marker.write_text(datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")

    log_out = run_dir / "pipeline_run.log"
    log.info(f"Lancement pipeline : bash {launch}  (log → {log_out})")
    try:
        with open(log_out, "w") as fout:
            ret = subprocess.run(
                ["bash", str(launch)],
                stdout=fout, stderr=subprocess.STDOUT,
                cwd=PIPELINE_DIR,
            )
        stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if ret.returncode == 0:
            log.info(f"Pipeline terminé avec succès pour {run_dir.name}.")
            (run_dir / PIPELINE_OK).write_text(stamp + "\n")
            send_results(run_dir)
        else:
            log.error(f"Pipeline a échoué (code {ret.returncode}) pour "
                      f"{run_dir.name}. Voir {log_out} — résultats NON renvoyés vers sitatst.")
            (run_dir / PIPELINE_FAILED).write_text(f"{stamp} code={ret.returncode}\n")
    except Exception as e:
        log.error(f"Erreur lancement pipeline : {e}")


def send_results(run_dir: Path) -> bool:
    """
    Renvoie pipeline_v0/ vers sitatst via scripts/send_results_to_sitatst.py
    (transfert, vérification, sentinelle RESULTS_COMPLETE, marqueur local).
    Codes du script : 0 = fait (ou déjà fait), 2 = désactivé, autre = échec.
    """
    log_out = run_dir / "send_results.log"
    cmd = [sys.executable, SEND_RESULTS, str(run_dir), "--site-paths", SITE_PATHS]
    log.info(f"Retour des résultats vers sitatst : {run_dir.name} (log → {log_out})")
    try:
        with open(log_out, "a") as fout:
            fout.write(f"===== {datetime.now():%Y-%m-%d %H:%M:%S}\n"); fout.flush()
            ret = subprocess.run(cmd, stdout=fout, stderr=subprocess.STDOUT, timeout=86400)
    except Exception as e:
        log.error(f"Retour des résultats impossible pour {run_dir.name} : {e}")
        return False
    if ret.returncode == 0:
        log.info(f"Résultats de {run_dir.name} transférés et vérifiés (RESULTS_COMPLETE déposé).")
        return True
    if ret.returncode == 2:
        log.info("Retour des résultats désactivé (sync_sitatst.enabled) — rien envoyé.")
        return False
    log.error(f"Retour des résultats en échec pour {run_dir.name} (code {ret.returncode}). "
              f"Nouvelle tentative dans {RESULTS_RETRY_DELAY // 60} min. Voir {log_out}")
    return False


def find_unsent_runs(watch_path: Path) -> list:
    """Runs réussis (.pipeline_ok) pas encore transférés (.results_synced),
    dont la dernière tentative date de plus de RESULTS_RETRY_DELAY."""
    now = time.time()
    pending = []
    for run_dir in sorted(watch_path.iterdir()):
        if not run_dir.is_dir():
            continue
        if not (run_dir / PIPELINE_OK).exists() or (run_dir / RESULTS_MARKER).exists():
            continue
        last = run_dir / "send_results.log"
        if last.exists() and now - last.stat().st_mtime < RESULTS_RETRY_DELAY:
            continue
        pending.append(run_dir)
    return pending


def process_ready_run(run_dir: Path) -> None:
    """preprocess -> pipeline pour un run prêt."""
    name = run_dir.name
    log.info("=" * 60)
    log.info(f"Run prêt détecté (sentinelle présente) : {name}")

    if not run_preprocess(run_dir):
        log.error(f"Run {name} : preprocess échoué — pipeline non lancé. "
                  f"(le run sera re-tenté au prochain scan si non marqué .launched)")
        return

    run_pipeline(run_dir)


def main():
    watch_path = Path(WATCH_DIR)
    if not watch_path.is_dir():
        log.critical(f"Répertoire surveillé introuvable : {WATCH_DIR}")
        sys.exit(1)

    log.info(f"Watcher OVH démarré — surveillance de {WATCH_DIR} "
             f"(intervalle {POLL_INTERVAL}s, sentinelle '{SENTINEL_NAME}')")

    while True:
        try:
            # Sérialisation : ne rien lancer si un pipeline tourne déjà.
            if pipeline_is_running():
                log.info("Un pipeline tourne déjà — attente (pas de lancement).")
            else:
                ready = find_ready_runs(watch_path)
                if ready:
                    # On traite UN SEUL run par cycle (le plus ancien), puis on
                    # laisse la boucle re-vérifier — garantit la sérialisation.
                    run_dir = ready[0]
                    process_ready_run(run_dir)
                else:
                    # Rien à lancer : retenter les retours de résultats en attente
                    unsent = find_unsent_runs(watch_path)
                    if unsent:
                        send_results(unsent[0])
        except Exception as e:
            log.error(f"Erreur dans la boucle principale : {e}", exc_info=True)

        time.sleep(POLL_INTERVAL)


if __name__ == "__main__":
    main()

