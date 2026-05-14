#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np
import re


# ----------------- Helper Functions -----------------
def read_outrider_files(directory: Path) -> pd.DataFrame:
    all_files = list(directory.glob("*.outrider.tab"))
    if not all_files:
        return pd.DataFrame()
    dfs = []
    for f in all_files:
        try:
            df = pd.read_csv(f, sep="\t", dtype=str)
            if not df.empty:
                dfs.append(df)
        except Exception as e:
            print(f"Error reading {f}: {e}")
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def convert_columns(df: pd.DataFrame) -> pd.DataFrame:
    # Normalise column names — outrider_newVersion.R writes padjValue/rawCounts
    if 'padjValue' in df.columns and 'padjust' not in df.columns:
        df = df.rename(columns={'padjValue': 'padjust'})
    if 'rawCounts' in df.columns and 'rawcounts' not in df.columns:
        df = df.rename(columns={'rawCounts': 'rawcounts'})
    for col in ["padjust", "l2fc", "theta"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.replace(",", ".", regex=False).apply(
                lambda x: pd.to_numeric(x, errors="coerce")
            )
    return df


def shorten_sample_id(s: str) -> str:
    s = str(s)
    if s.startswith("X"):
        s = s[1:]
    return re.split(r"[._-]", s)[0]


def compute_box_stats(df: pd.DataFrame):
    df = df.copy()
    df["is_zero"] = df["count"] == 0
    df["is_outlier"] = False

    stats = []
    for cat, group in df.groupby("RunCategory"):
        q1 = group["count"].quantile(0.25)
        q3 = group["count"].quantile(0.75)
        median = group["count"].median()
        iqr = q3 - q1
        lower = group[group["count"] >= (q1 - 1.5 * iqr)]["count"].min()
        upper = group[group["count"] <= (q3 + 1.5 * iqr)]["count"].max()

        group["is_outlier"] = (group["count"] < lower) | (group["count"] > upper)
        df.loc[group.index, "is_outlier"] = group["is_outlier"]

        stats.append({
            "RunCategory": cat,
            "q1": q1,
            "q3": q3,
            "median": median,
            "lower_whisker": lower,
            "upper_whisker": upper,
        })

    return pd.DataFrame(stats), df


# ----------------- PNG Plot (boxplot + stripplot) -----------------
def plot_png(count_df, stats_df, title, outfile):
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(8, 6))
    order = ["Older Runs", "Actual Run"]

    # Base boxplot
    sns.boxplot(
        x="RunCategory",
        y="count",
        data=count_df,
        order=order,
        palette={"Actual Run": "lightblue", "Older Runs": "lightgreen"},
        showfliers=False
    )

    # Color column for stripplot
    count_df["point_color"] = count_df.apply(
        lambda r: "red"
        if (r["is_outlier"] or r["is_zero"])
        else ("lightblue" if r["RunCategory"] == "Actual Run" else "lightgreen"),
        axis=1
    )

    # Stripplot on top
    sns.stripplot(
        x="RunCategory",
        y="count",
        data=count_df,
        order=order,
        jitter=True,
        size=6,
        alpha=0.9,
        palette={"Actual Run": "lightblue", "Older Runs": "lightgreen"},
        dodge=False
    )

    # Overlay red points
    red_points = count_df[count_df["point_color"] == "red"]
    for _, row in red_points.iterrows():
        plt.text(
            order.index(row["RunCategory"]) - 0.05,
            row["count"] + 0.02 * (1 + abs(row["count"])),
            row["sampleID_short"],
            color="red",
            fontsize=7,
            rotation=45,
            ha="center"
        )

    plt.title(title)
    plt.xlabel("Sample Group")
    plt.ylabel("Event Count")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"PNG saved to {outfile}")


# ----------------- HTML Plot (boxplot + jittered scatter) -----------------
def plot_html(count_df, stats_df, title, outfile_html):
    fig = go.Figure()
    categories = ["Older Runs", "Actual Run"]
    colors = {"Older Runs": "lightgreen", "Actual Run": "lightblue"}

    stats_map = {row["RunCategory"]: row for _, row in stats_df.iterrows()}

    for i, cat in enumerate(categories):
        group = count_df[count_df["RunCategory"] == cat]

        if cat not in stats_map:
            continue

        s = stats_map[cat]

        # Box
        fig.add_trace(go.Box(
            y=group["count"],
            x=[cat] * len(group),
            name=cat,
            boxpoints=False,
            marker_color=colors[cat]
        ))

        # Scatter with jitter
        jitter_width = 0.18
        x_positions = [i + np.random.uniform(-jitter_width, jitter_width) for _ in range(len(group))]

        for (idx, row), xj in zip(group.iterrows(), x_positions):
            color = "red" if (row["is_outlier"] or row["is_zero"]) else colors[cat]

            fig.add_trace(go.Scatter(
                x=[xj],
                y=[row["count"]],
                mode="markers",
                marker=dict(size=7, color=color, line=dict(width=0.5, color="black")),
                hovertext=row["sampleID_short"],
                hoverinfo="text+y",
                showlegend=False
            ))

    fig.update_xaxes(
        tickmode="array",
        tickvals=list(range(len(categories))),
        ticktext=categories,
        title="Sample Group"
    )
    fig.update_yaxes(title="Event Count")
    fig.update_layout(title=title, width=900, height=600, template="simple_white")

    fig.write_html(outfile_html, include_plotlyjs="cdn")
    print(f"HTML saved to {outfile_html}")


# ----------------- Main Script -----------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory")
    parser.add_argument("output_file")
    parser.add_argument("output_file_filter")
    parser.add_argument("actual_samples")
    args = parser.parse_args()

    df = convert_columns(read_outrider_files(Path(args.directory)))

    if df.empty:
        print("[INFO] No OUTRIDER events found — creating empty placeholder plots")
        Path(args.output_box).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output_box).touch()
        Path(args.output_filt).touch()
        return

    df["sampleID_short"] = df["sampleID"].apply(shorten_sample_id)

    actual_samples_short = set(shorten_sample_id(s) for s in args.actual_samples.split(","))

    df["RunCategory"] = df["sampleID_short"].apply(
        lambda x: "Actual Run" if x in actual_samples_short else "Older Runs"
    )

    # ---- Overall ----
    count_df = df.groupby(["sampleID", "RunCategory"]).size().reset_index(name="count")
    count_df = count_df.merge(
        df[["sampleID", "sampleID_short"]].drop_duplicates(),
        on="sampleID"
    )

    # Add missing samples with count=0
    all_short = set(count_df["sampleID_short"]) | actual_samples_short
    id_map = df[["sampleID_short", "sampleID"]].drop_duplicates()
    full = pd.DataFrame({"sampleID_short": list(all_short)})
    full = full.merge(id_map, on="sampleID_short", how="left")
    full["sampleID"] = full["sampleID"].fillna(full["sampleID_short"])

    full["RunCategory"] = full["sampleID_short"].apply(
        lambda x: "Actual Run" if x in actual_samples_short else "Older Runs"
    )

    count_df = full.merge(count_df, on=["sampleID", "sampleID_short", "RunCategory"], how="left")
    count_df["count"] = count_df["count"].fillna(0).astype(int)

    stats_df, count_df = compute_box_stats(count_df)

    plot_png(count_df, stats_df, "Outrider Events per Sample (All Data)", args.output_file)
    plot_html(count_df, stats_df, "Outrider Events per Sample (All Data)",
              Path(args.output_file).with_suffix(".html"))

    # ---- Filtered ----
    if "padjust" in df.columns:
        fdf = df[df["padjust"] < 0.05]
        count_f = fdf.groupby(["sampleID", "RunCategory"]).size().reset_index(name="count")
        count_f = count_f.merge(
            df[["sampleID", "sampleID_short"]].drop_duplicates(),
            on="sampleID"
        )

        count_f = full.merge(count_f, on=["sampleID", "sampleID_short", "RunCategory"], how="left")
        count_f["count"] = count_f["count"].fillna(0).astype(int)

        stats_df_f, count_f = compute_box_stats(count_f)

        plot_png(count_f, stats_df_f,
                 "Outrider Events per Sample (padj < 0.05)", args.output_file_filter)
        plot_html(count_f, stats_df_f,
                  "Outrider Events per Sample (padj < 0.05)",
                  Path(args.output_file_filter).with_suffix(".html"))

    else:
        print("No padjust column — skipping filtered plot.")


if __name__ == "__main__":
    main()

