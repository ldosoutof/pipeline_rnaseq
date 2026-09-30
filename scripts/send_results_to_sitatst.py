#!/usr/bin/env python3
"""
send_results_to_sitatst.py — renvoie les résultats d'un run (pipeline_v0/) d'OVH
vers sitatst (baie Isilon S:), puis dépose une sentinelle RESULTS_COMPLETE.

Symétrique du transfert aller (watcher_nextseq.py : FASTQ -> OVH + sentinelle
TRANSFER_COMPLETE). Appelé par watcher_ovh.py après un pipeline réussi ;
peut aussi être lancé à la main pour un run donné.

    python3 scripts/send_results_to_sitatst.py /datawork2/.../prod/<run>
    python3 scripts/send_results_to_sitatst.py <run_dir> --dry-run   # affiche seulement
    python3 scripts/send_results_to_sitatst.py <run_dir> --force     # renvoie même si déjà fait

Paramètres : section `sync_sitatst` de site_paths.yml
    enabled, host, user, dest_path, ssh_key (recommandé), pass_file (option).
Destination : <dest_path>/<run>/pipeline_v0/ ; sentinelle : <dest_path>/<run>/RESULTS_COMPLETE

Étapes (chacune bloquante) :
  1. rsync de pipeline_v0/ avec des options compatibles CIFS (la destination sur
     sitatst est un partage CIFS : ni permissions, ni propriétaire, ni liens
     symboliques, ni dates de dossiers) ;
  2. vérification : un second rsync à blanc (--dry-run --itemize-changes) ne
     doit plus trouver aucun fichier à transférer ;
  3. dépôt de la sentinelle distante RESULTS_COMPLETE (en dernier : sa présence
     garantit que le transfert est complet et vérifié) ;
  4. marqueur local <run>/.results_synced (idempotence).
Code de sortie : 0 = transféré et vérifié (ou déjà fait), 1 = échec, 2 = désactivé.
Mode local de test : host vide -> copie locale vers dest_path.
"""
import argparse
import os
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import yaml

SENTINEL = "RESULTS_COMPLETE"
MARKER = ".results_synced"
RSYNC_OPTS = [
    "-rtz", "--checksum", "--partial",
    "--no-perms", "--no-owner", "--no-group", "--omit-dir-times",   # destination CIFS
    "--skip-compress=gz/bam/bai/zip/png/pdf/xz/bgz",
]


def load_sync(site_paths: Path) -> dict:
    d = yaml.safe_load(site_paths.read_text()) or {}
    return d.get("sync_sitatst", {}) or {}


def build(sync: dict):
    """Renvoie (préfixe distant, commande ssh de base, env) selon la config."""
    host, user = sync.get("host", ""), sync.get("user", "")
    ssh_key, pass_file = sync.get("ssh_key", ""), sync.get("pass_file", "")
    env = dict(os.environ)
    ssh = ["ssh", "-o", "StrictHostKeyChecking=no"]
    if ssh_key:
        ssh += ["-i", os.path.expanduser(ssh_key)]
    wrap = []
    if pass_file and os.path.isfile(os.path.expanduser(pass_file)):
        env["SSHPASS"] = Path(os.path.expanduser(pass_file)).read_text().strip()
        wrap = ["sshpass", "-e"]            # mot de passe via SSHPASS, jamais visible dans ps
    else:
        ssh += ["-o", "BatchMode=yes"]      # clé seule : échoue au lieu d'attendre un mot de passe
    remote = f"{user}@{host}:" if host else ""
    return remote, ssh, wrap, env, bool(host)


def _exec(cmd, env, dry, capture=False):
    print("  $ " + " ".join(shlex.quote(c) for c in cmd if not c.startswith("SSHPASS")))
    if dry:
        return subprocess.CompletedProcess(cmd, 0, "", "")
    return subprocess.run(cmd, env=env, capture_output=capture, text=True, timeout=86400)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="dossier du run sur OVH (contient pipeline_v0/)")
    ap.add_argument("--site-paths", default=str(Path(__file__).resolve().parent.parent / "template" / "site_paths.yml"))
    ap.add_argument("--dry-run", action="store_true", help="affiche les commandes sans rien transférer")
    ap.add_argument("--force", action="store_true", help="renvoie même si le run est marqué .results_synced")
    a = ap.parse_args()

    run_dir = Path(a.run_dir).resolve()
    run = run_dir.name
    src = run_dir / "pipeline_v0"
    marker = run_dir / MARKER
    if not src.is_dir():
        print(f"[ERROR] {src} introuvable"); return 1
    if marker.exists() and not a.force:
        print(f"[INFO] {run} déjà transféré ({marker.read_text().strip()}) — rien à faire (--force pour renvoyer)")
        return 0

    sync = load_sync(Path(a.site_paths))
    if not sync.get("enabled", False):
        print(f"[INFO] sync_sitatst.enabled est faux dans {a.site_paths} — pas de transfert"); return 2
    dest_root = str(sync.get("dest_path", "")).rstrip("/")
    if not dest_root:
        print("[ERROR] sync_sitatst.dest_path manquant"); return 1
    remote, ssh, wrap, env, is_remote = build(sync)
    dest_run = f"{dest_root}/{run}"
    e_opt = ["-e", " ".join(shlex.quote(c) for c in ssh)] if is_remote else []
    print(f"[INFO] {run} : {src}/ -> {remote}{dest_run}/pipeline_v0/  ({datetime.now():%Y-%m-%d %H:%M:%S})")

    # 0. dossier de destination (rsync ne crée que le dernier niveau)
    mk = (wrap + ssh + [remote.rstrip(":"), f"mkdir -p {shlex.quote(dest_run + '/pipeline_v0')}"]) if is_remote \
        else ["mkdir", "-p", f"{dest_run}/pipeline_v0"]
    if _exec(mk, env, a.dry_run).returncode != 0:
        print("[ERROR] création du dossier distant impossible"); return 1

    # 1. transfert
    r = _exec(wrap + ["rsync"] + RSYNC_OPTS + e_opt + [f"{src}/", f"{remote}{dest_run}/pipeline_v0/"], env, a.dry_run)
    if r.returncode != 0:
        print(f"[ERROR] rsync a échoué (code {r.returncode}) — sentinelle NON déposée"); return 1

    # 2. vérification : plus aucun fichier à transférer (taille + date)
    v = _exec(wrap + ["rsync", "-rtn", "--itemize-changes", "--no-perms", "--no-owner", "--no-group",
                    "--omit-dir-times"] + e_opt + [f"{src}/", f"{remote}{dest_run}/pipeline_v0/"],
            env, a.dry_run, capture=True)
    pending = [l for l in (v.stdout or "").splitlines() if l[:2] in (">f", "<f", "cf")]
    if v.returncode != 0 or pending:
        print(f"[ERROR] vérification KO (code {v.returncode}, {len(pending)} fichier(s) différent(s)) — sentinelle NON déposée")
        for l in pending[:10]:
            print("    " + l)
        return 1
    print("[INFO] vérification OK : destination identique à la source")

    # 3. sentinelle distante, en dernier
    stamp = f"{datetime.now():%Y-%m-%d %H:%M:%S} {os.uname().nodename}"
    sent_cmd = f"echo {shlex.quote(stamp)} > {shlex.quote(dest_run + '/' + SENTINEL)}"
    s = _exec((wrap + ssh + [remote.rstrip(":"), sent_cmd]) if is_remote else ["bash", "-c", sent_cmd], env, a.dry_run)
    if s.returncode != 0:
        print("[ERROR] dépôt de la sentinelle impossible"); return 1

    # 4. marqueur local
    if not a.dry_run:
        marker.write_text(f"{stamp} -> {remote}{dest_run}\n")
    print(f"[INFO] {run} : transfert terminé, vérifié, sentinelle {SENTINEL} déposée")
    return 0


if __name__ == "__main__":
    sys.exit(main())
