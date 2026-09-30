#!/usr/bin/env python3
"""
plot_cibersortx.py

Génère un graphique en barres empilées des fractions cellulaires estimées par
CIBERSORTx, à la manière du rendu HTML du site en ligne :
un échantillon par barre, chaque segment = un type cellulaire LM22.

Usage :
    python plot_cibersortx.py \
        --results CIBERSORTx_Results.txt \
        --out cibersortx_barplot.png

Options :
    --top N        n'afficher que les N types cellulaires les plus abondants
                   (les autres regroupés en "Other") — défaut : tous (22)
    --html FICHIER produire aussi une version HTML interactive (plotly)
"""
import argparse
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Les 3 dernières colonnes du fichier CIBERSORTx ne sont pas des types cellulaires
META_COLS = ["P-value", "Correlation", "RMSE"]


def load_results(path):
    df = pd.read_csv(path, sep="\t")
    # première colonne = identifiant échantillon (souvent "Mixture")
    sample_col = df.columns[0]
    df = df.set_index(sample_col)
    # séparer types cellulaires (tout sauf les métadonnées finales)
    cell_cols = [c for c in df.columns if c not in META_COLS]
    fractions = df[cell_cols].astype(float)
    meta = df[[c for c in META_COLS if c in df.columns]]
    return fractions, meta


def plot_stacked(fractions, out_png, top=None):
    frac = fractions.copy()

    # optionnel : ne garder que les top N types, regrouper le reste en "Other"
    if top is not None and top < frac.shape[1]:
        means = frac.mean(axis=0).sort_values(ascending=False)
        keep = means.head(top).index.tolist()
        other = frac.drop(columns=keep).sum(axis=1)
        frac = frac[keep]
        frac["Other"] = other

    # couleurs distinctes (palette tab20, adaptée à ~22 catégories)
    n = frac.shape[1]
    cmap = plt.colormaps["tab20"].resampled(max(n, 3))
    colors = [cmap(i) for i in range(n)]

    fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(frac) + 4), 7))

    bottom = np.zeros(len(frac))
    x = np.arange(len(frac))
    for i, col in enumerate(frac.columns):
        ax.bar(x, frac[col].values, bottom=bottom, color=colors[i],
               width=0.8, label=col, edgecolor="white", linewidth=0.3)
        bottom += frac[col].values

    ax.set_xticks(x)
    ax.set_xticklabels(frac.index, rotation=90, fontsize=8)
    ax.set_ylabel("Fraction cellulaire estimée", fontsize=11)
    ax.set_ylim(0, 1)
    ax.set_title("CIBERSORTx — composition cellulaire relative (LM22)", fontsize=12)
    # légende à droite, hors du plot
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=7, ncol=1, frameon=False)
    plt.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    sys.stderr.write(f"[plot_cibersortx] barplot écrit : {out_png}\n")


def plot_html(fractions, out_html):
    """Version interactive optionnelle (plotly)."""
    try:
        import plotly.graph_objects as go
    except ImportError:
        sys.stderr.write("[plot_cibersortx] plotly non installé, HTML ignoré.\n")
        return
    frac = fractions
    fig = go.Figure()
    for col in frac.columns:
        fig.add_bar(x=list(frac.index), y=frac[col].values, name=col)
    fig.update_layout(
        barmode="stack",
        title="CIBERSORTx — composition cellulaire relative (LM22)",
        yaxis_title="Fraction cellulaire estimée",
        xaxis_tickangle=-90,
        legend_title="Type cellulaire",
        height=650,
    )
    fig.write_html(out_html)
    sys.stderr.write(f"[plot_cibersortx] HTML interactif écrit : {out_html}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="CIBERSORTx_Results.txt")
    ap.add_argument("--out", required=True, help="image PNG de sortie")
    ap.add_argument("--top", type=int, default=None,
                    help="n'afficher que les N types les plus abondants")
    ap.add_argument("--html", default=None, help="produire aussi un HTML interactif")
    args = ap.parse_args()

    fractions, meta = load_results(args.results)
    sys.stderr.write(
        f"[plot_cibersortx] {fractions.shape[0]} échantillons, "
        f"{fractions.shape[1]} types cellulaires.\n"
    )
    plot_stacked(fractions, args.out, top=args.top)
    if args.html:
        plot_html(fractions, args.html)


if __name__ == "__main__":
    main()
