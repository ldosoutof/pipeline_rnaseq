#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import re


# ----------------- Helper Functions -----------------
def read_fraser_files(directory: Path) -> pd.DataFrame:
    files = list(directory.glob("*.fraser.tab"))
    print("Detected files:", [str(f) for f in files])
    if not files:
        raise FileNotFoundError("No .fraser.tab files found")

    dfs = []
    for f in files:
        try:
            dfs.append(pd.read_csv(f, sep="\t", dtype=str))
        except Exception as e:
            print(f"Error reading {f}: {e}")
    if not dfs:
        raise FileNotFoundError("No valid .fraser.tab data read")
    return pd.concat(dfs, ignore_index=True)


def convert_columns(df: pd.DataFrame) -> pd.DataFrame:
    for col in ["padjust", "l2fc", "theta"]:
        if col in df.columns:
            df[col] = (
                df[col].astype(str)
                .str.replace(",", ".", regex=False)
                .apply(lambda x: pd.to_numeric(x, errors="coerce"))
            )
    return df


def shorten_sample_id(s: str) -> str:
    s = str(s)
    if s.startswith("X"):
        s = s[1:]
    return re.split(r"[._-]", s)[0]


def plot_png(count_df, title, outfile):
    """Static PNG boxplot with outliers labeled."""
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(9, 6))
    order = ["Older Runs", "Actual Run"]

    ax = sns.boxplot(
        x="RunCategory",
        y="count",
        data=count_df,
        order=order,
        hue="RunCategory",
        palette={"Actual Run": "lightblue", "Older Runs": "lightgreen"},
        legend=False,
    )

    # Outliers
    for category, group in count_df.groupby("RunCategory"):
        q1 = group["count"].quantile(0.25)
        q3 = group["count"].quantile(0.75)
        iqr = q3 - q1
        upper = q3 + 1.5 * iqr
        outliers = group[group["count"] > upper]

        for _, row in outliers.iterrows():
            ax.text(
                x={"Older Runs": 0, "Actual Run": 1}[category],
                y=row["count"] + 0.05 * row["count"],
                s=row["sampleID_short"],
                color="red",
                ha="center",
                fontsize=8,
                rotation=45,
            )

    ax.set_title(title)
    ax.set_xlabel("Sample Group")
    ax.set_ylabel("Event Count")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"PNG saved to {outfile}")


def plot_html(count_df, title, outfile_html):
    """Interactive HTML boxplot fully aligned with Seaborn (with hover stats)."""
    fig = go.Figure()

    categories = ["Older Runs", "Actual Run"]
    colors = {"Older Runs": "lightgreen", "Actual Run": "lightblue"}

    for i, cat in enumerate(categories):
        group = count_df[count_df["RunCategory"] == cat]
        if group.empty:
            continue

        y = group["count"].values
        q1 = group["count"].quantile(0.25)
        q3 = group["count"].quantile(0.75)
        median = group["count"].median()
        iqr = q3 - q1

        # Whiskers (like seaborn: closest points inside fences)
        lower_whisker = group["count"][group["count"] >= (q1 - 1.5 * iqr)].min()
        upper_whisker = group["count"][group["count"] <= (q3 + 1.5 * iqr)].max()

        # --- Draw box ---
        fig.add_shape(type="rect",
            x0=i-0.25, x1=i+0.25, y0=q1, y1=q3,
            line=dict(color="black"),
            fillcolor=colors[cat],
            layer="below"
        )

        # --- Median line ---
        fig.add_shape(type="line",
            x0=i-0.25, x1=i+0.25, y0=median, y1=median,
            line=dict(color="black", width=2)
        )

        # --- Whiskers ---
        fig.add_shape(type="line", x0=i, x1=i, y0=lower_whisker, y1=q1, line=dict(color="black"))
        fig.add_shape(type="line", x0=i, x1=i, y0=q3, y1=upper_whisker, line=dict(color="black"))
        fig.add_shape(type="line", x0=i-0.15, x1=i+0.15, y0=lower_whisker, y1=lower_whisker, line=dict(color="black"))
        fig.add_shape(type="line", x0=i-0.15, x1=i+0.15, y0=upper_whisker, y1=upper_whisker, line=dict(color="black"))

        # --- Overlay sample points ---
        for _, row in group.iterrows():
            color = "red" if row["count"] < lower_whisker or row["count"] > upper_whisker else colors[cat]
            fig.add_trace(go.Scatter(
                x=[i],
                y=[row["count"]],
                mode="markers",
                marker=dict(color=color, size=6, line=dict(width=0.5, color="black")),
                hovertext=row["sampleID_short"],
                hoverinfo="y+text",
                showlegend=False
            ))

        # --- Add invisible hover markers for summary stats ---
        stats_hover = [
            f"{cat} Q1: {q1:.2f}",
            f"{cat} Median: {median:.2f}",
            f"{cat} Q3: {q3:.2f}",
            f"{cat} Lower whisker: {lower_whisker:.2f}",
            f"{cat} Upper whisker: {upper_whisker:.2f}"
        ]
        stats_y = [q1, median, q3, lower_whisker, upper_whisker]

        fig.add_trace(go.Scatter(
            x=[i]*len(stats_y),
            y=stats_y,
            mode="markers",
            marker=dict(color="rgba(0,0,0,0)", size=10, symbol="circle-open"),
            hovertext=stats_hover,
            hoverinfo="text",
            showlegend=False
        ))

    fig.update_xaxes(
        tickmode="array",
        tickvals=list(range(len(categories))),
        ticktext=categories,
        title="Sample Group"
    )
    fig.update_yaxes(title="Event Count")

    fig.update_layout(
        title=title,
        template="simple_white",
        width=900,
        height=600
    )

    fig.write_html(outfile_html, include_plotlyjs="cdn")
    print(f"HTML saved to {outfile_html}")


# ----------------- Main Script -----------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate FRASER event count boxplots (PNG + HTML)."
    )
    parser.add_argument("directory", help="Directory with .fraser.tab files")
    parser.add_argument("output_png", help="Output PNG for overall boxplot")
    parser.add_argument("output_png_filtered", help="Output PNG for filtered data boxplot")
    parser.add_argument("actual_samples", help="Comma-separated actual sample IDs")
    args = parser.parse_args()

    print("Raw actual_samples:", args.actual_samples)

    # Read and preprocess data
    df = convert_columns(read_fraser_files(Path(args.directory)))
    if "sampleID" not in df.columns:
        raise ValueError("Missing 'sampleID' column in input data.")

    # Normalize sample IDs
    df["sampleID_short"] = df["sampleID"].apply(shorten_sample_id)
    actual_samples_short = set(shorten_sample_id(s) for s in args.actual_samples.split(","))
    print("Normalized actual_samples_short:", actual_samples_short)

    # Assign RunCategory
    df["RunCategory"] = df["sampleID_short"].apply(
        lambda x: "Actual Run" if x in actual_samples_short else "Older Runs"
    )

    # Aggregate counts per sample
    count_df = df.groupby(["sampleID_short", "RunCategory"]).size().reset_index(name="count")
    print("Sample counts by RunCategory:\n", count_df["RunCategory"].value_counts())

    # ---- Overall plot ----
    out_png = Path(args.output_png)
    plot_png(count_df, "FRASER Events per Sample (All Data)", out_png)
    plot_html(count_df, "FRASER Events per Sample (All Data)", out_png.with_suffix(".html"))

    # ---- Filtered plot (padjust < 0.05) ----
    if "padjust" in df.columns:
        filt = df[df["padjust"].notna() & (df["padjust"] < 0.05)]
        count_filt = filt.groupby(["sampleID_short", "RunCategory"]).size().reset_index(name="count")
        out_png_f = Path(args.output_png_filtered)
        plot_png(count_filt, "FRASER Events per Sample (padj < 0.05)", out_png_f)
        plot_html(count_filt, "FRASER Events per Sample (padj < 0.05)", out_png_f.with_suffix(".html"))
    else:
        print("Column 'padjust' missing—skipping filtered plot.")


if __name__ == "__main__":
    main()

