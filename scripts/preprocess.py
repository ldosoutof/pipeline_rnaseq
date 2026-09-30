#!/usr/bin/env python3
"""
preprocess.py
=============
Script utilitaire à lancer AVANT le pipeline Snakemake.

Pour chaque run, il :
  1. Scanne fastq/ pour détecter automatiquement les échantillons
  2. Génère un config.yml prêt à l'emploi en fusionnant :
       - config_template.yml  (modèle documenté, ne pas modifier)
       - site_paths.yml       (chemins réels de votre installation)
       - les valeurs du run   (fastq_dir, outputDir, samples — auto-détectés)
  3. Génère un launch.sh depuis launch_template.sh

Usage :
  python preprocess.py \\
      --path    /path/to/pipeline/ \\
      --workDir /path/to/run/      \\
      --dataDir /path/to/run/

  # Avec un fichier site_paths ailleurs que dans template/
  python preprocess.py \\
      --path      /path/to/pipeline/ \\
      --workDir   /path/to/run/      \\
      --dataDir   /path/to/run/      \\
      --site-paths /path/to/my_site_paths.yml

Arguments :
  -p / --path        Répertoire racine du pipeline
  -w / --workDir     Répertoire de travail du run (où iront les résultats)
  -d / --dataDir     Répertoire du run (doit contenir un sous-dossier fastq/)
  --site-paths       Fichier de chemins de l'installation (défaut : template/site_paths.yml)
"""

import argparse
import codecs
import os
import sys
import yaml
from pathlib import Path


# ── Helpers ───────────────────────────────────────────────────────────────────

def deep_merge(base, override):
    """
    Fusionne deux dicts YAML de manière récursive.
    Les valeurs non vides de override écrasent celles de base.
    Les valeurs vides ("", [], None) dans override sont ignorées
    pour ne pas écraser les valeurs par défaut du template.
    """
    result = base.copy()
    for key, val in override.items():
        if val in ("", [], None):
            continue                         # clé vide → garder la valeur template
        if key in result and isinstance(result[key], dict) and isinstance(val, dict):
            result[key] = deep_merge(result[key], val)
        else:
            result[key] = val
    return result


# ── Arguments ─────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser(
    prog="preprocess.py",
    description="Génère le config.yml et launch.sh pour un run RNA-seq",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
parser.add_argument("-p",  "--path",       required=True, help="Répertoire racine du pipeline")
parser.add_argument("-w",  "--workDir",    required=True, help="Répertoire de travail du run")
parser.add_argument("-d",  "--dataDir",    required=True, help="Répertoire du run (contient fastq/)")
parser.add_argument("--site-paths", default="",
                    help="Fichier de chemins de l'installation (défaut : template/site_paths.yml)")
args = parser.parse_args()

path_pipeline = Path(args.path).resolve()
path_data     = Path(args.dataDir).resolve()
path_workdir  = Path(args.workDir).resolve()
workdir       = path_workdir.name


# ── Validation ────────────────────────────────────────────────────────────────

fastq_dir = path_data / "fastq"
if not fastq_dir.is_dir():
    print(f"[ERROR] Dossier fastq/ introuvable : {fastq_dir}", file=sys.stderr)
    sys.exit(1)

template_config = path_pipeline / "template" / "config_template.yml"
if not template_config.exists():
    print(f"[ERROR] config_template.yml introuvable : {template_config}", file=sys.stderr)
    sys.exit(1)

template_launch = path_pipeline / "template" / "launch_template.sh"
if not template_launch.exists():
    print(f"[ERROR] launch_template.sh introuvable : {template_launch}", file=sys.stderr)
    sys.exit(1)

# Localiser site_paths.yml
site_paths_path = Path(args.site_paths).resolve() if args.site_paths else \
                  path_pipeline / "template" / "site_paths.yml"

if not site_paths_path.exists():
    print(f"[WARN] site_paths.yml introuvable : {site_paths_path}", file=sys.stderr)
    print(f"       Le config généré contiendra les chemins génériques du template.", file=sys.stderr)
    print(f"       → Copiez template/site_paths.yml, remplissez vos chemins réels,",  file=sys.stderr)
    print(f"         puis relancez avec --site-paths /chemin/vers/site_paths.yml",     file=sys.stderr)
    site_paths = {}
else:
    site_paths = yaml.safe_load(site_paths_path.read_text()) or {}
    print(f"[INFO] Chemins de l'installation chargés depuis : {site_paths_path}")


# ── Dossier de sortie ─────────────────────────────────────────────────────────

launch_folder = path_pipeline / workdir / "launch_folder"
launch_folder.mkdir(parents=True, exist_ok=True)
print(f"[INFO] Dossier de lancement : {launch_folder}")


# ── Détection des échantillons ────────────────────────────────────────────────

sample_ids = []
for f in sorted(os.listdir(fastq_dir)):
    if f.endswith("_R1.fastq.gz") and not f.startswith("Undetermined"):
        sample_ids.append(f.rsplit("_", 1)[0])

if not sample_ids:
    print(f"[ERROR] Aucun fichier *_R1.fastq.gz trouvé dans {fastq_dir}", file=sys.stderr)
    sys.exit(1)

print(f"[INFO] {len(sample_ids)} échantillon(s) détecté(s) :")
for s in sample_ids:
    print(f"       {s}")


# ── Construction du config ────────────────────────────────────────────────────
# Ordre de priorité (du moins au plus prioritaire) :
#   1. config_template.yml  — valeurs par défaut documentées
#   2. site_paths.yml       — chemins réels de l'installation
#   3. valeurs du run       — auto-détectées (écrasent tout)

base_config = yaml.safe_load(template_config.read_text()) or {}
config      = deep_merge(base_config, site_paths)

# Valeurs spécifiques au run — toujours écrasées en dernier
config.setdefault("configuration", {})
config["configuration"]["samples"]      = sample_ids
config["configuration"]["outputDir"]    = str(path_workdir)
config["configuration"]["fastq_dir"]    = str(fastq_dir)
config["configuration"]["pipeline_dir"] = str(path_pipeline)


# ── Écriture du config.yml ────────────────────────────────────────────────────

new_config = launch_folder / "config.yml"
with open(new_config, "w") as fp:
    yaml.dump(config, fp, default_flow_style=False, allow_unicode=True, sort_keys=False)

# Vérifier si des chemins génériques /path/to restent dans le config
remaining_placeholders = [
    k for k, v in config.items()
    if isinstance(v, str) and "/path/to/" in v
]
if remaining_placeholders:
    print(f"\n[WARN] Ces clés contiennent encore des chemins génériques :")
    for k in remaining_placeholders:
        print(f"       {k}: {config[k]}")
    print(f"       → Complétez {new_config} ou remplissez site_paths.yml")
else:
    print(f"\n[INFO] Tous les chemins sont renseignés ✅")

print(f"[INFO] Config généré : {new_config}")


# ── Génération du launch.sh ───────────────────────────────────────────────────

with codecs.open(template_launch, encoding="utf-8", errors="ignore") as fp:
    launch_content = fp.read()

launch_content = launch_content.replace("PATH_TO_CONFIG", str(new_config))
launch_content = launch_content.replace("ROOT_PIPELINE",  str(path_pipeline))

launch_sh = launch_folder / "launch.sh"
with open(launch_sh, "w") as fp:
    fp.write(launch_content)
launch_sh.chmod(0o755)

print(f"[INFO] Script de lancement généré : {launch_sh}")
print()
print("Pour lancer le pipeline :")
print(f"  bash {launch_sh}")
