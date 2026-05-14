#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from pathlib import Path
import plotly.express as px

def save_interactive_html(df, x, y, color, html_path, xlab, ylab, title):
    """
    Save an interactive volcano plot as HTML.
    Hover shows gene_name, pValue, zScore, l2fc, Model_Of_Inheritance (if present).
    """
    hover_cols = ["gene_name", "pValue", "zScore", "l2fc"]
    if "Model_Of_Inheritance" in df.columns:
        hover_cols.append("Model_Of_Inheritance")

    fig = px.scatter(
        df,
        x=x,
        y=y,
        color=color,
        color_discrete_map={
            "Up-regulated": "red",
            "Down-regulated": "blue",
            "Unchanged": "grey",
        },
        hover_data=hover_cols,
        title=title,
    )
    fig.update_layout(
        xaxis_title=xlab,
        yaxis_title=ylab,
        template="simple_white",
        width=900,
        height=650,
    )
    fig.write_html(html_path, include_plotlyjs="cdn")
    print(f"🌐 interactive plot saved to {html_path}")


def main(outrider_file: str, output_png: str):
    out1 = Path(output_png)
    out1.parent.mkdir(parents=True, exist_ok=True)

    # ---------- Read & clean ----------
    df = pd.read_csv(outrider_file, sep="\t", dtype=str)

    # Early exit for samples with no aberrant events (empty or header-only file)
    if df.empty or len(df.columns) < 3:
        print(f"[INFO] No aberrant events for {Path(outrider_file).stem} — skipping volcano plot")
        # Create empty placeholder files so Snakemake output check passes
        for path in [out1,
                     out1.with_suffix(".html"),
                     out1.with_name(out1.stem + "_fc.png"),
                     out1.with_name(out1.stem + "_fc.html")]:
            path.touch()
        return

    # Numeric conversions (replace commas with dots first)
    for col in ["pValue", "zScore", "l2fc"]:
        df[col] = (
            df[col].str.replace(",", ".", regex=False).astype(float)
            if col in df.columns else np.nan
        )
    df["fc"] = 2 ** df["l2fc"]

    # Ensure gene_name column exists
    if "gene_name" not in df.columns:
        raise ValueError("Input must contain a 'gene_name' column")

    # Keep rows with a gene name
    top = df[df["gene_name"].notna() & (df["gene_name"] != "")].copy()

    # -------------------------------------------------------------------------
    # 1) Volcano Plot by zScore
    # -------------------------------------------------------------------------
    top["Expression"] = np.select(
        [
            (top["zScore"] > 3) & (top["pValue"] <= 0.01),
            (top["zScore"] < -3) & (top["pValue"] <= 0.01),
        ],
        ["Up-regulated", "Down-regulated"],
        default="Unchanged",
    )

    plt.figure(figsize=(10, 7))
    sns.scatterplot(
        data=top,
        x="zScore",
        y=-np.log10(top["pValue"]),
        hue="Expression",
        palette={"Up-regulated": "red", "Down-regulated": "blue", "Unchanged": "grey"},
        edgecolor=None,
        s=25,
    )

    plt.axhline(y=2, linestyle="dotted", color="black")
    plt.axvline(x=3, linestyle="dashed", color="black")
    plt.axvline(x=-3, linestyle="dashed", color="black")
    plt.xlim(-10, 10)
    plt.xlabel("zScore")
    plt.ylabel("-log10(pValue)")
    plt.title("Volcano Plot (zScore)")

    # Label significant genes with Model_Of_Inheritance
    texts = []
    if "Model_Of_Inheritance" in top.columns:
        label_rows = top[
            top["Expression"].isin(["Up-regulated", "Down-regulated"])
            & top["Model_Of_Inheritance"].notna()
            & (top["Model_Of_Inheritance"] != "")
        ]
        for _, r in label_rows.iterrows():
            texts.append(
                plt.text(
                    r["zScore"],
                    -np.log10(r["pValue"]),
                    r["gene_name"],
                    fontsize=8,
                    color="black",
                )
            )
        adjust_text(texts, arrowprops=dict(arrowstyle="->", color="gray", lw=0.5),expand_points=(1.2, 1.2),force_text=0.5)

    out1 = Path(output_png)
    plt.tight_layout()
    plt.savefig(out1, dpi=300)
    plt.close()
    print(f"✅ zScore volcano saved to {out1}")

    out1_html = out1.with_suffix(".html")
    save_interactive_html(
        top,
        x="zScore",
        y=-np.log10(top["pValue"]),
        color="Expression",
        html_path=out1_html,
        xlab="zScore",
        ylab="-log10(pValue)",
        title="Volcano Plot (zScore)",
    )



    # -------------------------------------------------------------------------
    # 2) Volcano Plot by log2 Fold Change (l2fc)
    # -------------------------------------------------------------------------
    top["Expression_fc"] = np.select(
        [
            (top["l2fc"] > 1) & (top["pValue"] <= 0.01),
            (top["l2fc"] < -1) & (top["pValue"] <= 0.01),
        ],
        ["Up-regulated", "Down-regulated"],
        default="Unchanged",
    )

    plt.figure(figsize=(10, 7))
    sns.scatterplot(
        data=top,
        x="l2fc",
        y=-np.log10(top["pValue"]),
        hue="Expression_fc",
        palette={"Up-regulated": "red", "Down-regulated": "blue", "Unchanged": "grey"},
        edgecolor=None,
        s=25,
    )

    plt.axhline(y=2, linestyle="dotted", color="black")
    plt.axvline(x=1, linestyle="dashed", color="black")
    plt.axvline(x=-1, linestyle="dashed", color="black")
    plt.xlim(-5, 5)
    plt.xlabel("log2 Fold Change")
    plt.ylabel("-log10(pValue)")
    plt.title("Volcano Plot (log2 Fold Change)")

    texts = []
    if "Model_Of_Inheritance" in top.columns:
        label_rows = top[
            top["Expression_fc"].isin(["Up-regulated", "Down-regulated"])
            & top["Model_Of_Inheritance"].notna()
            & (top["Model_Of_Inheritance"] != "")
        ]
        for _, r in label_rows.iterrows():
            texts.append(
                plt.text(
                    r["l2fc"],
                    -np.log10(r["pValue"]),
                    r["gene_name"],
                    fontsize=8,
                    color="black",
                )
            )
        adjust_text(texts, arrowprops=dict(arrowstyle="->", color="gray", lw=0.5),expand_points=(1.2, 1.2),force_text=0.5)

    out2 = out1.with_name(out1.stem + "_fc.png")
    plt.tight_layout()
    plt.savefig(out2, dpi=300)
    plt.close()
    print(f"✅ log2FC volcano saved to {out2}")

    out2_html = out2.with_suffix(".html")
    save_interactive_html(
        top,
        x="l2fc",
        y=-np.log10(top["pValue"]),
        color="Expression_fc",
        html_path=out2_html,
        xlab="log2 Fold Change",
        ylab="-log10(pValue)",
        title="Volcano Plot (log2 Fold Change)",
    )

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python volcano2.py <outrider_file.tsv> <output_plot.png>")
    main(sys.argv[1], sys.argv[2])

