#!/usr/bin/env python3
"""Generate barcode- and gene-level QC summaries from an all.poolcount table."""

from __future__ import annotations
import argparse
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    """Define and parse command-line options for the QC workflow."""
    parser = argparse.ArgumentParser(
        description="Create CSV and PNG QC outputs for a tab-delimited barcode count table."
    )
    parser.add_argument("-i", "--input", required=True, type=Path, help="all.poolcount TSV")
    parser.add_argument("-o", "--output", required=True, type=Path, help="output directory")
    parser.add_argument(
        "-m",
        "--name-map",
        dest="name_map",
        type=Path,
        help="optional CSV with old_colnames and new_colnames columns",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="write CSV reports only (matplotlib is not required)",
    )
    parser.add_argument(
        "--annotation-columns",
        type=int,
        default=7,
        help="number of leading non-sample columns (default: 7)",
    )
    parser.add_argument("--t0-prefix", default="T0", help="T0 sample prefix (default: T0)")
    parser.add_argument("--top-n", type=int, default=10, help="top features per sample (default: 10)")
    parser.add_argument(
        "--low-count-threshold",
        type=int,
        default=3,
        help="count values in 1..threshold-1 are low (default: 3)",
    )
    return parser.parse_args()


def require_file(path: Path, label: str) -> None:
    """Fail early with a clear message when a required file is missing."""
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")


def load_data(args: argparse.Namespace) -> tuple[pd.DataFrame, list[str]]:
    """Load the count table, optionally rename columns, and validate counts."""
    require_file(args.input, "Input file")
    allpool = pd.read_csv(args.input, sep="\t", low_memory=False)

    # The name-map file translates original experiment names into shorter,
    # human-readable sample names used throughout the reports and plots.
    if args.name_map:
        require_file(args.name_map, "Name-map file")
        name_map = pd.read_csv(args.name_map, sep=",")
        needed = {"old_colnames", "new_colnames"}
        if not needed.issubset(name_map.columns):
            raise ValueError(f"Name-map must contain columns: {sorted(needed)}")
        mapping = dict(zip(name_map["old_colnames"], name_map["new_colnames"]))
        allpool = allpool.rename(columns=mapping)

    if args.annotation_columns < 1 or args.annotation_columns >= allpool.shape[1]:
        raise ValueError("--annotation-columns must leave at least one sample column")
    if "barcode" not in allpool.columns:
        raise ValueError("Input must contain a 'barcode' column")

    # Everything after the leading annotation fields is treated as a sample.
    sample_cols = allpool.columns[args.annotation_columns:].tolist()

    # Coerce sample columns to numbers and reject text or negative counts.
    # Missing count values are interpreted as zero reads.
    converted = allpool[sample_cols].apply(pd.to_numeric, errors="coerce")
    invalid = converted.isna() & allpool[sample_cols].notna()
    if invalid.any().any():
        bad_cols = invalid.any()[invalid.any()].index.tolist()
        raise ValueError(f"Non-numeric values found in sample columns: {bad_cols}")
    if (converted.fillna(0) < 0).any().any():
        raise ValueError("Sample counts must be non-negative")
    allpool.loc[:, sample_cols] = converted.fillna(0)
    return allpool, sample_cols


def plot_top_features(
    data: pd.DataFrame,
    samples: list[str],
    label_col: str,
    output_file: Path,
    title: str,
) -> None:
    """Plot the most abundant barcodes or genes separately for each sample."""
    # Arrange panels for a fixed 16:9 canvas so the complete figure fills a
    # standard widescreen presentation slide without becoming tall or square.
    ncols = max(1, math.ceil(math.sqrt(len(samples) * 16 / 9)))
    nrows = math.ceil(len(samples) / ncols)
    dense_grid = ncols >= 6
    if len(samples) >= 24:
        y_label_size = 2.5
    elif len(samples) >= 16:
        y_label_size = 3.5
    else:
        y_label_size = 5
    panel_title_size = 8 if dense_grid else 10
    axis_label_size = 6 if dense_grid else 8
    x_tick_size = 6 if dense_grid else 7
    value_label_size = 5 if dense_grid else 6
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(16, 9),
        squeeze=False,
    )
    flat_axes = axes.ravel()
    for ax, sample in zip(flat_axes, samples):
        subset = data.loc[data["sample"] == sample].sort_values("population_in_pool_pct")
        bars = ax.barh(subset[label_col].astype(str), subset["population_in_pool_pct"])
        for bar in bars:
            value = bar.get_width()
            ax.text(
                value,
                bar.get_y() + bar.get_height() / 2,
                f" {value:.2f}%",
                va="center",
                fontsize=value_label_size,
            )
        ax.set_title(str(sample), fontsize=panel_title_size)
        ax.set_xlabel("Percent of reads in sample", fontsize=axis_label_size)
        ax.tick_params(axis="x", labelsize=x_tick_size)
        ax.tick_params(axis="y", labelsize=y_label_size)
        ax.margins(x=0.2)
    for ax in flat_axes[len(samples):]:
        ax.axis("off")
    fig.suptitle(title, fontsize=14, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(output_file, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_sample_qc(
    qc: pd.DataFrame,
    output_file: Path,
    low_count_threshold: int,
) -> None:
    """Plot sample QC, combining related count and percentage measurements."""
    samples = qc.index.astype(str)
    # A fixed 16:9 canvas fits a standard widescreen presentation slide. Reduce
    # tick-label size when an experiment contains many samples instead of making
    # the figure wider than the slide.
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), squeeze=False)
    tick_font_size = max(5, min(7, 10 - len(samples) // 5))
    percent_label_size = max(4, min(6, 8 - len(samples) // 6))
    flat_axes = axes.ravel()

    def add_bar_labels(ax, bars, decimals: int = 0) -> None:
        labels = [
            f"{bar.get_height():,.{decimals}f}"
            if math.isfinite(bar.get_height()) else ""
            for bar in bars
        ]
        ax.bar_label(bars, labels=labels, padding=3, fontsize=7, rotation=90)
        ax.margins(y=0.18)

    def plot_single_metric(
        ax, metric: str, title: str, ylabel: str, decimals: int = 0
    ) -> None:
        values = pd.to_numeric(qc[metric], errors="coerce")
        bars = ax.bar(samples, values, color="#4C78A8")
        add_bar_labels(ax, bars, decimals)
        ax.set_title(title)
        ax.set_ylabel(ylabel)

    def plot_count_and_percent(
        ax, count_metric: str, percent_metric: str, title: str
    ) -> None:
        count_values = pd.to_numeric(qc[count_metric], errors="coerce")
        percent_values = pd.to_numeric(qc[percent_metric], errors="coerce")
        bars = ax.bar(samples, count_values, color="#4C78A8", label="Count")
        count_labels = [
            f"{bar.get_height():,.0f}"
            if math.isfinite(bar.get_height()) else ""
            for bar in bars
        ]
        # Keep count labels above the bars, matching the other QC panels. Extra
        # headroom reserves a clear strip for the labels and combined legend.
        ax.bar_label(
            bars,
            labels=count_labels,
            label_type="edge",
            padding=5,
            fontsize=7,
            rotation=90,
            color="black",
        )
        finite_counts = count_values.dropna()
        count_upper = finite_counts.max() if not finite_counts.empty else 1.0
        ax.set_ylim(0, max(1.0, count_upper * 1.3))
        ax.set_title(title)
        ax.set_ylabel("Count")

        percent_ax = ax.twinx()
        percent_ax.plot(
            samples, percent_values, color="#E45756", marker="o",
            linewidth=1.8, label="%",
        )
        for sample, value in zip(samples, percent_values):
            if pd.notna(value):
                percent_ax.annotate(
                    f"{value:.1f}%", (sample, value), xytext=(0, -10),
                    textcoords="offset points", ha="center", va="top",
                    fontsize=percent_label_size,
                    color="#B23A39",
                    bbox={
                        "boxstyle": "round,pad=0.08",
                        "facecolor": "white",
                        "edgecolor": "none",
                        "alpha": 0.85,
                    },
                )
        percent_ax.set_ylabel("%")
        percent_ax.set_ylim(0, max(140, percent_values.max(skipna=True) * 1.45))
        count_handles, count_legend_labels = ax.get_legend_handles_labels()
        percent_handles, percent_legend_labels = percent_ax.get_legend_handles_labels()
        ax.legend(
            count_handles + percent_handles,
            count_legend_labels + percent_legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.99),
            ncol=2,
            fontsize=8,
        )

    plot_single_metric(flat_axes[0], "total_reads", "Total Reads", "Reads")
    plot_single_metric(
        flat_axes[1], "mean_reads_per_detected_barcode",
        "Mean Reads per Detected Barcode", "Mean reads", decimals=1,
    )
    plot_count_and_percent(
        flat_axes[2], "detected_barcodes", "detected_barcodes_pct",
        "Detected Barcodes: Count and %",
    )
    plot_count_and_percent(
        flat_axes[3],
        f"barcodes_with_1_to_{low_count_threshold - 1}_reads",
        f"barcodes_with_1_to_{low_count_threshold - 1}_reads_pct",
        f"Barcodes with 1–{low_count_threshold - 1} Reads: Count and %",
    )
    plot_single_metric(
        flat_axes[4], "t0_barcode_retention_pct", "T0 Barcode Retention",
        "%", decimals=1,
    )
    flat_axes[5].axis("off")

    for ax in flat_axes[:5]:
        ax.tick_params(axis="x", labelrotation=60, labelsize=tick_font_size)
        ax.grid(axis="y", linestyle=":", alpha=0.4)

    fig.suptitle("Sample QC Summary", fontsize=16, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(output_file, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_barcode_count_distribution(
    counts: pd.DataFrame,
    output_file: Path,
) -> None:
    """Plot clear percentile ranges for nonzero barcode counts."""
    # For each sample, the dot shows the median, the thick line shows the middle
    # 50%, and the thin line shows the 5th-95th percentile range. Zero counts are
    # reported in sample_qc.csv and omitted because a log axis cannot display 0.
    positive_counts = counts.where(counts > 0)
    percentiles = positive_counts.quantile([0.05, 0.25, 0.50, 0.75, 0.95])
    positions = list(range(1, len(counts.columns) + 1))
    figure_width = max(10, len(counts.columns) * 0.55)
    fig, ax = plt.subplots(figsize=(figure_width, 6))

    ax.vlines(
        positions, percentiles.loc[0.05], percentiles.loc[0.95],
        color="#9ECAE1", linewidth=2, label="5th-95th percentile",
    )
    ax.vlines(
        positions, percentiles.loc[0.25], percentiles.loc[0.75],
        color="#3182BD", linewidth=8, label="25th-75th percentile",
    )
    ax.scatter(
        positions, percentiles.loc[0.50], color="#B22222", s=28,
        zorder=3, label="Median",
    )
    for position, median in zip(positions, percentiles.loc[0.50]):
        if pd.notna(median):
            ax.annotate(
                f"{median:,.0f}", (position, median), xytext=(5, 0),
                textcoords="offset points", fontsize=7, va="center",
            )

    ax.set_yscale("log")
    ax.set_xticks(list(positions), counts.columns.astype(str), rotation=60, ha="right")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Read count per barcode (detected barcodes only; log scale)")
    ax.set_title("Barcode Count Distribution by Percentile")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(axis="y", which="both", linestyle=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(output_file, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_barcode_count_violin(
    counts: pd.DataFrame,
    output_file: Path,
) -> None:
    """Plot log-transformed count distributions for detected barcodes."""
    sample_names = []
    log_counts = []

    for sample in counts.columns:
        detected_counts = counts.loc[counts[sample] > 0, sample].to_numpy()
        if detected_counts.size == 0:
            continue
        sample_names.append(str(sample))
        log_counts.append(np.log10(detected_counts))

    if not log_counts:
        return

    positions = list(range(1, len(sample_names) + 1))
    figure_width = max(10, len(sample_names) * 0.55)
    fig, ax = plt.subplots(figsize=(figure_width, 6))
    violin_parts = ax.violinplot(
        log_counts,
        positions=positions,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        widths=0.8,
    )
    for body in violin_parts["bodies"]:
        body.set_facecolor("#4C78A8")
        body.set_edgecolor("#2F4B6C")
        body.set_alpha(0.75)

    # Overlay narrow boxplots to show the median, interquartile range, and
    # 1.5-IQR whiskers without obscuring the violin distributions.
    ax.boxplot(
        log_counts,
        positions=positions,
        widths=0.14,
        patch_artist=True,
        showfliers=False,
        boxprops={"facecolor": "white", "edgecolor": "#222222", "linewidth": 1},
        whiskerprops={"color": "#222222", "linewidth": 1},
        capprops={"color": "#222222", "linewidth": 1},
        medianprops={"color": "#B22222", "linewidth": 1.8},
    )

    ax.set_xticks(positions, sample_names, rotation=60, ha="right")
    ax.set_xlabel("Sample")
    ax.set_ylabel("log10(read count per detected barcode)")
    ax.set_title("Barcode Count Distribution — Violin Plot")
    ax.text(
        0.01,
        0.99,
        "Zeros excluded; white boxes show IQR and red lines show medians",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(output_file, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_gene_insertion_distribution(
    gene_insertions: pd.DataFrame,
    output_file: Path,
) -> None:
    """Plot the insertion-count distribution across genes for each condition."""
    conditions = gene_insertions.columns.astype(str).tolist()
    insertion_counts = [
        gene_insertions[condition].to_numpy(dtype=float)
        for condition in gene_insertions.columns
    ]

    # A violin density cannot be estimated when every gene has the same count.
    # Add a negligible spread for density estimation only; the boxplots below
    # always use the unmodified insertion counts.
    violin_values = []
    for values in insertion_counts:
        if values.size < 2 or np.allclose(values, values[0]):
            center = values[0] if values.size else 0.0
            violin_values.append(np.array([center - 0.001, center, center + 0.001]))
        else:
            violin_values.append(values)

    positions = list(range(1, len(conditions) + 1))
    figure_width = max(10, len(conditions) * 0.55)
    fig, ax = plt.subplots(figsize=(figure_width, 6))
    violin_parts = ax.violinplot(
        violin_values,
        positions=positions,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        widths=0.8,
    )
    for body in violin_parts["bodies"]:
        body.set_facecolor("#59A14F")
        body.set_edgecolor("#2F5D2A")
        body.set_alpha(0.75)

    ax.boxplot(
        insertion_counts,
        positions=positions,
        widths=0.14,
        patch_artist=True,
        showfliers=False,
        boxprops={"facecolor": "white", "edgecolor": "#222222", "linewidth": 1},
        whiskerprops={"color": "#222222", "linewidth": 1},
        capprops={"color": "#222222", "linewidth": 1},
        medianprops={"color": "#B22222", "linewidth": 1.8},
    )
    ax.set_xticks(positions, conditions, rotation=60, ha="right")
    ax.set_xlabel("Condition")
    ax.set_ylabel("Detected insertions per gene")
    ax.set_title("Insertion Number per Gene")
    ax.set_ylim(bottom=0)
    ax.text(
        0.01,
        0.99,
        "White boxes show IQR; red lines show medians; outliers hidden",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(output_file, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.low_count_threshold < 2:
        raise ValueError("--low-count-threshold must be at least 2")
    if not args.skip_plots and plt is None:
        raise SystemExit("matplotlib is required for plots; install it or use --skip-plots")
    args.output.mkdir(parents=True, exist_ok=True)
    plot_dir = args.output / "plots"
    plot_dir.mkdir(exist_ok=True)

    allpool, sample_cols = load_data(args)

    # Calculate barcode-level coverage statistics independently per sample.
    counts = allpool[sample_cols].astype(float)
    totals = counts.sum()
    detected = counts.gt(0).sum()
    low = counts.gt(0) & counts.lt(args.low_count_threshold)

    # T0 columns represent the starting pool; all remaining columns represent
    # later selections or conditions.
    t0_cols = [c for c in sample_cols if str(c).startswith(args.t0_prefix)]
    selection_cols = [c for c in sample_cols if c not in t0_cols]

    # Define the reference barcode set as barcodes detected in every T0 sample.
    # Retention is the percentage of that reference set detected after selection.
    retained = pd.Series(index=sample_cols, dtype=float)
    if t0_cols:
        present_in_all_t0 = counts[t0_cols].gt(0).all(axis=1)
        retained.loc[t0_cols] = 100.0
        if selection_cols:
            retained.loc[selection_cols] = counts.loc[present_in_all_t0, selection_cols].gt(0).mean().mul(100)

    # Include the exact low-count range in the CSV headers. With the default
    # threshold of 3, low counts are barcode counts of 1 or 2 reads.
    low_count_name = f"barcodes_with_1_to_{args.low_count_threshold - 1}_reads"
    qc = pd.DataFrame(
        {
            "total_reads": totals,
            "detected_barcodes": detected,
            "mean_reads_per_detected_barcode": totals.div(detected.replace(0, pd.NA)),
            "detected_barcodes_pct": counts.gt(0).mean().mul(100),
            low_count_name: low.sum(),
            f"{low_count_name}_pct": low.mean().mul(100),
            "t0_barcode_retention_pct": retained,
        }
    )
    qc.index.name = "sample"
    qc.round(3).to_csv(args.output / "sample_qc.csv")
    if not args.skip_plots:
        plot_sample_qc(
            qc,
            plot_dir / "sample_qc.png",
            args.low_count_threshold,
        )

    counts.describe(percentiles=[0.5, 0.75, 0.9, 0.95, 0.99]).to_csv(
        args.output / "barcode_count_distribution.csv"
    )
    if not args.skip_plots:
        plot_barcode_count_distribution(
            counts,
            plot_dir / "barcode_count_distribution.png",
        )
        plot_barcode_count_violin(
            counts,
            plot_dir / "barcode_count_distribution_violin.png",
        )

    # Identify the most abundant individual barcodes in each sample. Location
    # fields are included when available to help trace each barcode to a gene.
    top_barcodes = []
    location_cols = [c for c in ("barcode", "locusId", "scaffold", "pos") if c in allpool.columns]
    for sample in sample_cols:
        top = allpool.nlargest(args.top_n, sample)[location_cols + [sample]].copy()
        top = top.rename(columns={sample: "hits"})
        top["sample"] = sample
        top["population_in_pool_pct"] = (top["hits"] / totals[sample] * 100) if totals[sample] else 0.0
        top_barcodes.append(top)
    top_barcode_df = pd.concat(top_barcodes, ignore_index=True)

    # Give intergenic barcodes a stable plotting label based on genomic location
    # when no locusId annotation is available.
    if "locusId" in top_barcode_df:
        fallback = "intergenic"
        if {"scaffold", "pos"}.issubset(top_barcode_df.columns):
            fallback = (
                "intergenic_"
                + top_barcode_df["scaffold"].fillna("unknown").astype(str)
                + "_"
                + top_barcode_df["pos"].fillna("unknown").astype(str)
            )
        top_barcode_df["locusId"] = top_barcode_df["locusId"].fillna(fallback)
        # Keep complete barcode sequences in the CSV, but abbreviate them in
        # plots so long sequences do not squeeze the bar-chart panels.
        barcode_text = top_barcode_df["barcode"].astype(str)
        barcode_plot_text = barcode_text.where(
            barcode_text.str.len() <= 15,
            barcode_text.str[:8] + "…" + barcode_text.str[-4:],
        )
        top_barcode_df["plot_label"] = (
            top_barcode_df["locusId"].astype(str) + "\n" + barcode_plot_text
        )
    else:
        top_barcode_df["plot_label"] = top_barcode_df["barcode"].astype(str)
    top_barcode_df.to_csv(args.output / "top_barcodes.csv", index=False)
    if not args.skip_plots:
        plot_top_features(top_barcode_df, sample_cols, "plot_label", plot_dir / "top_barcodes.png", "Top barcodes")

    # Exclude barcodes that have zero reads in every sample before feature-level
    # summaries. Genic barcodes are summed by locus; intergenic barcodes remain
    # individual features because they have no locusId.
    detected_rows = counts.gt(0).any(axis=1)
    filtered = allpool.loc[detected_rows].copy()
    if "locusId" in filtered.columns:
        genic = filtered.loc[filtered["locusId"].notna()].copy()
        intergenic = filtered.loc[filtered["locusId"].isna()].copy()
        gene_hits = genic.groupby("locusId", dropna=False)[sample_cols].sum()

        # Count independently detected barcode insertions per gene and condition.
        gene_insertions = (
            genic[sample_cols]
            .gt(0)
            .groupby(genic["locusId"])
            .sum()
            .astype(int)
        )
        gene_insertions.index.name = "locusId"
        gene_insertions.to_csv(args.output / "gene_insertion_counts.csv")
        if not args.skip_plots:
            plot_gene_insertion_distribution(
                gene_insertions,
                plot_dir / "gene_insertion_counts_violin.png",
            )

        # Report the distribution of total barcode reads assigned to each gene.
        gene_hits.describe(percentiles=[0.5, 0.75, 0.9, 0.95, 0.99]).to_csv(
            args.output / "gene_count_distribution.csv"
        )
        feature_qc = pd.DataFrame(
            {
                "detected_genes": gene_hits.gt(0).sum(),
                "detected_intergenic_barcodes": intergenic[sample_cols].gt(0).sum(),
                "mean_insertions_per_gene": gene_insertions.sum().div(
                    gene_insertions.gt(0).sum().replace(0, pd.NA)
                ),
            }
        )
        feature_qc["total_detected_features"] = (
            feature_qc["detected_genes"]
            + feature_qc["detected_intergenic_barcodes"]
        )
        feature_qc.index.name = "sample"
        feature_qc.to_csv(args.output / "gene_qc.csv")
        if not args.skip_plots:
            ax = feature_qc[["detected_genes", "detected_intergenic_barcodes"]].plot(
                kind="bar", stacked=True, figsize=(max(8, len(sample_cols) * 0.45), 5)
            )
            ax.set_ylabel("Count")

            # Show the actual feature counts inside each stacked bar segment.
            for container in ax.containers:
                labels = [
                    f"{int(bar.get_height())}" if bar.get_height() > 0 else ""
                    for bar in container
                ]
                ax.bar_label(
                    container,
                    labels=labels,
                    label_type="center",
                    fontsize=7,
                    color="white",
                    fontweight="bold",
                )

            # Mean insertions per gene use a much smaller scale than feature
            # counts, so show them as a line on a separate right-hand y-axis.
            mean_ax = ax.twinx()
            mean_insertions = pd.to_numeric(
                feature_qc["mean_insertions_per_gene"], errors="coerce"
            )
            positions = np.arange(len(feature_qc))
            mean_line = mean_ax.plot(
                positions,
                mean_insertions,
                color="#E45756",
                marker="o",
                linewidth=2,
                label="Mean insertions per gene",
            )[0]
            for position, value in zip(positions, mean_insertions):
                if pd.notna(value):
                    mean_ax.annotate(
                        f"{value:.2f}",
                        (position, value),
                        xytext=(0, 4),
                        textcoords="offset points",
                        ha="center",
                        fontsize=7,
                        color="#B23A39",
                    )
            mean_ax.set_ylabel("Mean insertions per detected gene")
            finite_means = mean_insertions.dropna()
            mean_upper = finite_means.max() if not finite_means.empty else 1.0
            mean_ax.set_ylim(0, max(1.0, mean_upper * 1.2))

            bar_handles, bar_labels = ax.get_legend_handles_labels()
            ax.legend(
                bar_handles + [mean_line],
                bar_labels + [mean_line.get_label()],
                loc="lower center",
                bbox_to_anchor=(0.5, 1.04),
                ncol=3,
                fontsize=8,
                frameon=False,
            )
            mean_ax.legend().remove()

            # Reserve space above the axes so the legend never covers bars or
            # the mean-insertion line, even when values approach the axis limit.
            ax.figure.tight_layout(rect=(0, 0, 1, 0.9))
            ax.figure.savefig(
                plot_dir / "gene_intergenic_counts.png",
                dpi=200,
                bbox_inches="tight",
            )
            plt.close(ax.figure)

        # Rank genes by their summed barcode counts in each sample. The two
        # percentages distinguish abundance among genes from abundance in the
        # complete pool, which also contains intergenic barcodes.
        top_genes = []
        for sample in sample_cols:
            top = gene_hits.nlargest(args.top_n, sample)[[sample]].reset_index()
            top = top.rename(columns={sample: "hits"})
            top["sample"] = sample
            gene_total = gene_hits[sample].sum()
            top["population_in_genes_pct"] = (top["hits"] / gene_total * 100) if gene_total else 0.0
            top["population_in_pool_pct"] = (top["hits"] / totals[sample] * 100) if totals[sample] else 0.0
            top_genes.append(top)
        top_gene_df = pd.concat(top_genes, ignore_index=True)
        top_gene_df.to_csv(args.output / "top_genes.csv", index=False)
        if not args.skip_plots:
            plot_top_features(top_gene_df, sample_cols, "locusId", plot_dir / "top_genes.png", "Top genes")

    print(f"QC complete: {args.output.resolve()}")


if __name__ == "__main__":
    main()
