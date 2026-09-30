#!/usr/bin/env python3
"""
plot_vaf_violin_by_sample.py

Violin plot du biais d'inactivation du chromosome X (VAF repliée) pour les
échantillons féminins d'un run. Une colonne par échantillon.

Contrôle positif (optionnel) : un échantillon de référence au biais connu
(p. ex. biais quasi-complet) peut être superposé pour servir de repère visuel,
via --control <vaf_table.tsv> <label>. Il est tracé dans une couleur distincte
et placé en dernière position.

Usage :
  plot_vaf_violin_by_sample.py <output.png> <chrY_mean_expression.tsv> \
      [--control <vaf_table.tsv> <label>] \
      <vaf_table1.tsv> [vaf_table2.tsv ...]
"""

import argparse
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

VAF_COLS = ["CHROM", "POS", "REF", "ALT", "DP", "VAF", "folded"]
CONTROL_COLOR = "#8B0000"   # dark red — distinct des violons standards


def load_vaf_table(path, sample_label):
    """Charge une vaf_table.tsv, filtre chrX, etiquette avec sample_label."""
    df = pd.read_csv(path, sep="\t", header=None, names=VAF_COLS)
    required = {"CHROM", "POS", "VAF", "folded"}
    if not required.issubset(df.columns):
        sys.exit(f"{path} : colonnes manquantes {required}, "
                 f"trouve {list(df.columns)}")
    df = df[df["CHROM"].isin(["X", "chrX"])].copy()
    df["sample"] = sample_label
    return df


def main():
    ap = argparse.ArgumentParser(description="Violin plot biais X (VAF repliee)")
    ap.add_argument("output_png", help="chemin de l'image de sortie")
    ap.add_argument("sex_table", help="chrY_mean_expression.tsv (inference sexe)")
    ap.add_argument("vaf_tables", nargs="+", help="vaf_table.tsv des echantillons")
    ap.add_argument("--control", nargs=2, metavar=("VAF_TABLE", "LABEL"),
                    default=None,
                    help="echantillon controle (fichier VAF + label affiche)")
    args = ap.parse_args()

    # -- Inference du sexe : ne garder que les femmes (chrY faible) --
    sex_df = pd.read_csv(args.sex_table, sep="\t")
    if not {"sample", "mean_chrY_TPM"}.issubset(sex_df.columns):
        sys.exit(f"Table sexe : colonnes 'sample','mean_chrY_TPM' requises, "
                 f"trouve {list(sex_df.columns)}")
    female_samples = sex_df.loc[sex_df["mean_chrY_TPM"] < 1, "sample"].tolist()
    if not female_samples:
        sys.exit("Aucun echantillon feminin detecte (mean_chrY_TPM < 1)")

    # -- Charger les VAF des echantillons feminins --
    dfs = []
    for vaf_file in args.vaf_tables:
        sample = Path(vaf_file).parent.name.replace("-", ".")
        if sample not in female_samples:
            continue
        dfs.append(load_vaf_table(vaf_file, sample))
    if not dfs:
        sys.exit("Aucune donnee VAF pour les echantillons feminins")
    vaf_df = pd.concat(dfs, ignore_index=True)

    sample_order = sorted(vaf_df["sample"].unique())

    # -- Charger le controle (optionnel) --
    control_df = None
    control_label = None
    if args.control:
        ctrl_path, control_label = args.control
        if not Path(ctrl_path).is_file():
            sys.exit(f"Controle introuvable : {ctrl_path}")
        control_df = load_vaf_table(ctrl_path, control_label)
        if control_df.empty:
            sys.exit(f"Controle {ctrl_path} : aucune donnee chrX")
        sample_order = sample_order + [control_label]
        plot_df = pd.concat([vaf_df, control_df], ignore_index=True)
    else:
        plot_df = vaf_df

    # -- Plot --
    plt.figure(figsize=(max(10, len(sample_order) * 0.6), 6))

    palette = {s: "#4C72B0" for s in sample_order}
    if control_label is not None:
        palette[control_label] = CONTROL_COLOR

    sns.violinplot(data=plot_df, x="sample", y="folded", order=sample_order,
                   hue="sample", palette=palette, legend=False,
                   inner=None, cut=0, linewidth=1)
    sns.stripplot(data=plot_df, x="sample", y="folded", order=sample_order,
                  color="black", size=2, alpha=0.4)

    medians = plot_df.groupby("sample")["folded"].median()
    for i, sample in enumerate(sample_order):
        if sample not in medians.index:
            continue
        median = medians[sample]
        mcol = CONTROL_COLOR if sample == control_label else "red"
        plt.plot(i, median, "_", color=mcol, markersize=12, markeredgewidth=2)
        plt.text(i, median + 0.01, f"{median:.2f}",
                 ha="center", va="bottom", fontsize=8, color=mcol)

    plt.axhline(0.5, color="green",  linestyle="--", linewidth=0.8,
                label="Balanced (median=0.5)")
    plt.axhline(0.3, color="orange", linestyle="--", linewidth=0.8,
                label="Moderate skew (0.3)")
    plt.axhline(0.1, color="red",    linestyle="--", linewidth=0.8,
                label="Extreme skew (0.1)")

    if control_label is not None:
        ctrl_idx = sample_order.index(control_label)
        plt.axvspan(ctrl_idx - 0.5, ctrl_idx + 0.5, color=CONTROL_COLOR, alpha=0.06)

    plt.ylim(0.0, 0.5)
    plt.ylabel("Folded VAF  [ min(VAF, 1-VAF) ]")
    plt.xlabel("Sample (females only)"
               + ("  |  dernier = controle" if control_label else ""))
    plt.title("chrX X-inactivation bias\n"
              "Median 0.5 = balanced  |  Median -> 0.0 = complete skew")
    plt.legend(fontsize=8)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(args.output_png, dpi=300)
    plt.close()

    msg = f"[OK] Violin plot ecrit : {args.output_png}"
    if control_label:
        msg += f" (avec controle '{control_label}')"
    print(msg)


if __name__ == "__main__":
    main()
