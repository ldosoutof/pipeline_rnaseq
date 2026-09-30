#!/usr/bin/env bash
# make_release.sh — empaquetage sécurisé du pipeline RNA-seq
#
# Remplace l'usage de `zip -r` qui ignore le .gitignore.
# Usage :
#   bash make_release.sh [version]
#   bash make_release.sh 2.1.0
#
# Prérequis : git (utilisé pour lire .gitignore via git archive),
#             pre_archive_check.sh dans le même répertoire ou dans PATH.
#
# Sortie : pipeline-<version>.zip dans le répertoire courant.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

VERSION="${1:-$(date +%Y%m%d)}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="pipeline-${VERSION}.zip"
GUARD="${REPO_ROOT}/pre_archive_check.sh"

# ── 1. Refuser si le dépôt n'est pas propre ──────────────────────────────────
# git archive ne prend que les fichiers commités — les modifications non
# commitées ne sont pas incluses dans l'archive. Si un fichier sensible a été
# commité par erreur, il faut le retirer de l'historique (git filter-repo)
# avant de lancer ce script.
if ! git -C "$REPO_ROOT" rev-parse --git-dir > /dev/null 2>&1; then
    echo "[ERROR] Ce répertoire n'est pas un dépôt git."
    echo "        git archive est requis pour que .gitignore soit respecté."
    exit 1
fi

# ── 2. Lancer le garde de sécurité sur l'arbre de travail ───────────────────
# On contrôle l'arbre de travail (pas l'archive) pour attraper les fichiers
# non encore commités qui pourraient être ajoutés manuellement à l'archive.
if [ -x "$GUARD" ]; then
    echo "[INFO] Contrôle de sécurité pré-empaquetage..."
    if ! bash "$GUARD" "$REPO_ROOT"; then
        echo "[ERROR] Contenu sensible détecté — empaquetage annulé."
        echo "        Corrigez les problèmes signalés, puis relancez."
        exit 1
    fi
else
    echo "[WARN] pre_archive_check.sh introuvable à $GUARD — contrôle ignoré."
fi

# ── 3. Créer l'archive via git archive ──────────────────────────────────────
# git archive respecte .gitignore ET .gitattributes (export-ignore).
# Tous les fichiers listés dans .gitignore sont automatiquement exclus.
echo "[INFO] Création de $OUT depuis HEAD..."

git -C "$REPO_ROOT" archive \
    --format=zip \
    --prefix="pipeline-${VERSION}/" \
    -o "$REPO_ROOT/$OUT" \
    HEAD

# ── 4. Vérification post-archive ────────────────────────────────────────────
# Extraire dans un tmpdir et re-passer le garde pour confirmer l'absence de
# fichiers sensibles dans l'archive finale.
TMPDIR_CHECK="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_CHECK"' EXIT

unzip -q "$REPO_ROOT/$OUT" -d "$TMPDIR_CHECK"

if [ -x "$GUARD" ]; then
    echo "[INFO] Contrôle de sécurité post-archive..."
    if ! bash "$GUARD" "$TMPDIR_CHECK"; then
        rm -f "$REPO_ROOT/$OUT"
        echo "[ERROR] L'archive contient du contenu sensible — fichier supprimé."
        echo "        Vérifiez que .gitignore est à jour et que les fichiers"
        echo "        sensibles n'ont jamais été commités."
        exit 1
    fi
fi

# ── 5. Résumé ────────────────────────────────────────────────────────────────
FILE_COUNT="$(unzip -l "$REPO_ROOT/$OUT" | tail -1 | awk '{print $2}')"
echo ""
echo "✅  Archive créée : $OUT"
echo "    Version : $VERSION"
echo "    Fichiers : $FILE_COUNT"
echo "    SHA256  : $(sha256sum "$REPO_ROOT/$OUT" | cut -d' ' -f1)"
