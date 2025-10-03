#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import re

# ----------------- Helper Functions -----------------
def read_outrider_files(directory: Path) -> pd.DataFrame:
    all_files = list(directory.glob("*.outrider.tab"))
    if not all_files:
        raise FileNotFoundError(f"No .outrider.tab files found in {directory}")
    dfs = []
    for f in all_files:
        try:
            dfs.append(pd.read_csv(f, sep="\t", dtype=str))
        except Exception as e:
            print(f"Error reading {f}: {e}")
    if not dfs:
        raise FileNotFoundError("No valid .outrider.tab files could be read")
    return pd.concat(dfs, ignore_index=True)


def convert_columns(df: pd.DataFrame) -> pd.DataFrame:
    for col in ["padjust", "l2fc", "theta"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.replace(",", ".", regex=False).apply(
                lambda x: pd.to_numeric(x, errors="coerce")
            )
    return df


def shorten_sample_id(s: str) -> str:
    s = str(s)
    if s.startswith("X"):
        s = s[1:]  # remove only first leading X
    return re.split(r"[._-]", s)[0]


def compute_box_stats(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["is_outlier"] = False
    stats = []
    for cat, group in df.groupby("RunCategory"):
        q1 = group["count"].quantile(0.25)
        q3 = group["count"].quantile(0.75)
        median = group["count"].median()
        iqr = q3 - q1
        lower_whisker = group[group["count"] >= (q1 - 1.5 * iqr)]["count"].min()
        upper_whisker = group[group["count"] <= (q3 + 1.5 * iqr)]["count"].max()

        group["is_outlier"] = (group["count"] < lower_whisker) | (group["count"] > upper_whisker)
        df.loc[group.index, "is_outlier"] = group["is_outlier"]

        stats.append({
            "RunCategory": cat,
            "q1": q1,
            "q3": q3,
            "median": median,
            "lower_whisker": lower_whisker,
            "upper_whisker": upper_whisker,
        })
    return pd.DataFrame(stats), df


def plot_png(count_df, stats_df, title, outfile):
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(8, 6))
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

    outliers = count_df[count_df["is_outlier"]]
    x_map = {cat: i for i, cat in enumerate(order)}
    for _, row in outliers.iterrows():
        ax.text(
            x_map[row["RunCategory"]],
            row["count"] + 0.05 * row["count"],
            row["sampleID_short"],
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


def plot_html(count_df, stats_df, title, outfile_html):
    """
    Interactive HTML boxplot aligned with PNG using precomputed stats.
    """
    fig = go.Figure()
    categories = ["Older Runs", "Actual Run"]
    colors = {"Older Runs": "lightgreen", "Actual Run": "lightblue"}

    # Map RunCategory to stats
    stats_map = {row["RunCategory"]: row for _, row in stats_df.iterrows()}

    for i, cat in enumerate(categories):
        group = count_df[count_df["RunCategory"] == cat]
        if group.empty or cat not in stats_map:
            continue

        stats = stats_map[cat]
        q1 = stats["q1"]
        q3 = stats["q3"]
        median = stats["median"]
        lower_whisker = stats["lower_whisker"]
        upper_whisker = stats["upper_whisker"]

        # Draw box
        fig.add_shape(type="rect", x0=i-0.25, x1=i+0.25, y0=q1, y1=q3,
                      line=dict(color="black"), fillcolor=colors[cat], layer="below")
        # Median
        fig.add_shape(type="line", x0=i-0.25, x1=i+0.25, y0=median, y1=median,
                      line=dict(color="black", width=2))
        # Whiskers
        fig.add_shape(type="line", x0=i, x1=i, y0=lower_whisker, y1=q1, line=dict(color="black"))
        fig.add_shape(type="line", x0=i, x1=i, y0=q3, y1=upper_whisker, line=dict(color="black"))
        fig.add_shape(type="line", x0=i-0.15, x1=i+0.15, y0=lower_whisker, y1=lower_whisker, line=dict(color="black"))
        fig.add_shape(type="line", x0=i-0.15, x1=i+0.15, y0=upper_whisker, y1=upper_whisker, line=dict(color="black"))

        # Sample points
        for _, row in group.iterrows():
            color = "red" if row["is_outlier"] else colors[cat]
            fig.add_trace(go.Scatter(
                x=[i], y=[row["count"]],
                mode="markers",
                marker=dict(color=color, size=6, line=dict(width=0.5, color="black")),
                hovertext=row["sampleID_short"],
                hoverinfo="y+text",
                showlegend=False
            ))

    fig.update_xaxes(tickmode="array", tickvals=list(range(len(categories))), ticktext=categories, title="Sample Group")
    fig.update_yaxes(title="Event Count")
    fig.update_layout(title=title, template="simple_white", width=900, height=600)
    fig.write_html(outfile_html, include_plotlyjs="cdn")
    print(f"HTML saved to {outfile_html}")


# ----------------- Main Script -----------------
def main():
    parser = argparse.ArgumentParser(description="Generate Outrider boxplots (PNG + HTML).")
    parser.add_argument("directory", help="Directory containing .outrider.tab files")
    parser.add_argument("output_file", help="Output PNG for overall boxplot")
    parser.add_argument("output_file_filter", help="Output PNG for filtered data boxplot")
    parser.add_argument("actual_samples", help="Comma-separated actual sample IDs")
    args = parser.parse_args()

    df = convert_columns(read_outrider_files(Path(args.directory)))
    if "sampleID" not in df.columns:
        raise ValueError("'sampleID' column missing in input data")

    df["sampleID_short"] = df["sampleID"].apply(shorten_sample_id)
    actual_samples_short = set(shorten_sample_id(s) for s in args.actual_samples.split(","))
    df["RunCategory"] = df["sampleID_short"].apply(
        lambda x: "Actual Run" if x in actual_samples_short else "Older Runs"
    )

    # ---- Overall ----
    count_df = df.groupby(["sampleID", "RunCategory"]).size().reset_index(name="count")
    count_df = count_df.merge(df[["sampleID", "sampleID_short"]].drop_duplicates(),
                              on="sampleID", how="left")
    print("Sample counts by RunCategory:\n", count_df["RunCategory"].value_counts())
    stats_df, count_df = compute_box_stats(count_df)
    plot_png(count_df, stats_df, "Outrider Events per Sample (All Data)", args.output_file)
    plot_html(count_df, stats_df, "Outrider Events per Sample (All Data)", Path(args.output_file).with_suffix(".html"))

    # ---- Filtered ----
    if "padjust" in df.columns:
        filt = df[df["padjust"].notna() & (df["padjust"] < 0.05)]
        count_filt = filt.groupby(["sampleID", "RunCategory"]).size().reset_index(name="count")
        count_filt = count_filt.merge(df[["sampleID", "sampleID_short"]].drop_duplicates(),
                                      on="sampleID", how="left")
        stats_df_f, count_filt = compute_box_stats(count_filt)
        plot_png(count_filt, stats_df_f, "Outrider Events per Sample (padj < 0.05)", args.output_file_filter)
        plot_html(count_filt, stats_df_f, "Outrider Events per Sample (padj < 0.05)",
                  Path(args.output_file_filter).with_suffix(".html"))
    else:
        print("Column 'padjust' missing—skipping filtered plot.")


if __name__ == "__main__":
    main()

