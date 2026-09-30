#!/usr/bin/env bash
# =============================================================
# deploy_watcher_ovh.sh
# Déploiement du watcher OVH (lance le pipeline à l'arrivée des FASTQ)
# en service systemd, sous l'utilisateur ldosoutoferreira.
# Usage : bash deploy_watcher_ovh.sh
# =============================================================
set -uo pipefail

# ──────────────────────────────────────────────
# CONFIGURATION — À AJUSTER selon l'installation OVH
# ──────────────────────────────────────────────
# Chemin de l'env conda 'snakemake' sur OVH.
# Base conda OVH : /appli/miniconda3 (confirmé via whereis conda).
# Confirmer le chemin EXACT de l'env avec :  conda env list | grep -i snakemake
CONDA_ENV_PATH="/appli/miniconda3/envs/snakemake"   # <-- vérifier avec conda env list

# Répertoire des scripts / logs sur OVH
SCRIPTS_DIR="/home/ldosoutoferreira/scripts"
LOGS_DIR="/home/ldosoutoferreira/logs"

WATCHER_SRC_NAME="watcher_ovh.py"              # le script (à côté de ce deploy)
WATCHER_SCRIPT="${SCRIPTS_DIR}/watcher_ovh.py" # destination
SERVICE_NAME="watcher_ovh"
SERVICE_FILE="/etc/systemd/system/${SERVICE_NAME}.service"
RUN_USER="ldosoutoferreira"
# ──────────────────────────────────────────────

echo "======================================================"
echo " Déploiement watcher OVH (${SERVICE_NAME})"
echo "======================================================"

# ── 1. Vérifier le python de l'env conda ────────────────────
PYTHON="${CONDA_ENV_PATH}/bin/python"
echo "[1/5] Vérification de l'env conda snakemake..."
if [ ! -x "${PYTHON}" ]; then
    echo "  ERREUR : python introuvable dans l'env conda : ${PYTHON}"
    echo "  Vérifiez CONDA_ENV_PATH (conda env list | grep -i snakemake)."
    exit 1
fi
# Vérifier que snakemake est bien présent dans cet env
if [ ! -x "${CONDA_ENV_PATH}/bin/snakemake" ]; then
    echo "  ATTENTION : snakemake introuvable dans ${CONDA_ENV_PATH}/bin/"
    echo "  Le watcher lance launch.sh qui appelle snakemake — vérifiez l'env."
fi
echo "  Python : ${PYTHON}"

# ── 2. Dossiers ─────────────────────────────────────────────
echo "[2/5] Création des répertoires..."
mkdir -p "${SCRIPTS_DIR}" "${LOGS_DIR}"

# ── 3. Copie du script ──────────────────────────────────────
echo "[3/5] Installation de ${WATCHER_SRC_NAME}..."
SCRIPT_SRC="$(cd "$(dirname "$0")" && pwd)/${WATCHER_SRC_NAME}"
if [ ! -f "${SCRIPT_SRC}" ]; then
    echo "  ERREUR : ${WATCHER_SRC_NAME} introuvable à côté de ce script."
    exit 1
fi
if [ "${SCRIPT_SRC}" != "${WATCHER_SCRIPT}" ]; then
    cp "${SCRIPT_SRC}" "${WATCHER_SCRIPT}"
    echo "  Copié : ${SCRIPT_SRC} -> ${WATCHER_SCRIPT}"
else
    echo "  Déjà en place."
fi
chmod 750 "${WATCHER_SCRIPT}"

# ── 4. Service systemd ──────────────────────────────────────
echo "[4/5] Installation du service systemd ${SERVICE_NAME}..."
# On appelle directement le python de l'env conda (pas de 'conda activate' :
# non disponible dans un shell systemd non interactif). Le PATH inclut le bin
# de l'env pour que launch.sh -> snakemake trouve les bons exécutables.
sudo tee "${SERVICE_FILE}" > /dev/null << UNIT
[Unit]
Description=Watcher OVH — lance le pipeline RNA-seq à l'arrivée des FASTQ
After=network.target

[Service]
Type=simple
User=${RUN_USER}
ExecStart=${PYTHON} ${WATCHER_SCRIPT}
Restart=on-failure
RestartSec=30
StandardOutput=append:${LOGS_DIR}/${SERVICE_NAME}.log
StandardError=append:${LOGS_DIR}/${SERVICE_NAME}.log
Environment=PATH=${CONDA_ENV_PATH}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
Environment=CONDA_PREFIX=${CONDA_ENV_PATH}

[Install]
WantedBy=multi-user.target
UNIT

sudo systemctl daemon-reload
sudo systemctl enable "${SERVICE_NAME}"

# ── 5. Démarrage ────────────────────────────────────────────
echo "[5/5] Démarrage du service..."
sudo systemctl restart "${SERVICE_NAME}"
sleep 2
if sudo systemctl is-active --quiet "${SERVICE_NAME}"; then
    echo ""
    echo "======================================================"
    echo " Service ${SERVICE_NAME} démarré avec succès."
    echo " Logs   : tail -f ${LOGS_DIR}/${SERVICE_NAME}.log"
    echo " Statut : sudo systemctl status ${SERVICE_NAME}"
    echo "======================================================"
else
    echo ""
    echo "======================================================"
    echo " ATTENTION : le service n'est pas actif. Diagnostic :"
    echo "   sudo systemctl status ${SERVICE_NAME}"
    echo "   tail -n 50 ${LOGS_DIR}/${SERVICE_NAME}.log"
    echo "======================================================"
fi

