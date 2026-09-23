#!/usr/bin/env python3
"""Summarize and visualize FEBA gene-fitness results.

The script reads fit_logratios.tab and fit_t.tab, renames experiment columns
with a sample-name map, and retains gene-condition measurements whose absolute t-like
statistic exceeds the selected threshold.
"""

import argparse
import json
import logging
import math
import sys
import time
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde, pearsonr, spearmanr
from sklearn.decomposition import PCA


ID_COLUMNS = ["locusId", "sysName", "desc"]
SENSITIVITY_THRESHOLDS = (0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 4.0)
COMBINED_HEATMAP_THRESHOLDS = (1.5, 2.0, 2.5)
LOGGER = logging.getLogger("fit_analysis")

# Apply one typeface consistently to titles, axes, ticks, legends,
# annotations, color bars, and any math-text labels in every figure.
plt.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
})


def setup_logging(output_path):
    """Log each run to both the terminal and the output directory."""
    LOGGER.setLevel(logging.INFO)
    LOGGER.handlers.clear()
    LOGGER.propagate = False
    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    terminal_handler = logging.StreamHandler(sys.stderr)
    terminal_handler.setFormatter(formatter)
    file_handler = logging.FileHandler(
        output_path / "fit_analysis.log", mode="w", encoding="utf-8"
    )
    file_handler.setFormatter(formatter)
    LOGGER.addHandler(terminal_handler)
    LOGGER.addHandler(file_handler)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="Directory containing fit_logratios.tab and fit_t.tab")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Directory for output tables and plots")
    parser.add_argument("-m", "--name-map", dest="name_map", type=Path, required=True,
                        help="CSV containing a new_colnames column")
    parser.add_argument("-g", "--groups", type=Path, required=True,
                        help="JSON mapping conditions to replicate samples")
    parser.add_argument("-t", "--threshold", type=float, required=True,
                        help="Retain values with absolute t-like statistic above this threshold")
    return parser.parse_args()


def load_inputs(args):
    """Read inputs, validate matching tables, and apply readable sample names."""
    fit = pd.read_csv(args.input / "fit_logratios.tab", sep="\t")
    fit_t = pd.read_csv(args.input / "fit_t.tab", sep="\t")
    name_map = pd.read_csv(args.name_map)
    with args.groups.open(encoding="utf-8") as handle:
        groups = json.load(handle)

    if fit.columns.tolist() != fit_t.columns.tolist():
        raise ValueError("fit_logratios.tab and fit_t.tab must have identical columns")
    if "new_colnames" not in name_map:
        raise ValueError("Name-map must contain a 'new_colnames' column")

    old_names = fit.columns[len(ID_COLUMNS):].tolist()
    name_map = name_map.dropna(subset=["new_colnames"]).copy()
    new_names_all = name_map["new_colnames"].astype(str).tolist()
    if len(new_names_all) != len(old_names):
        raise ValueError(
            f"Column mismatch: {len(new_names_all)} new names for {len(old_names)} fitness columns"
        )

    # FEBA retains Time0 columns for quality control. They are not experimental
    # fitness measurements, so remove each baseline column and its paired
    # name-map entry before downstream filtering and comparisons.
    paired_names = list(zip(old_names, new_names_all))
    experimental_pairs = [
        (old, new) for old, new in paired_names
        if new.split("-", 1)[0] != "T0"
    ]
    experimental_old_names = [old for old, _ in experimental_pairs]
    experimental_new_names = [new for _, new in experimental_pairs]
    rename_map = dict(experimental_pairs)
    keep_columns = ID_COLUMNS + experimental_old_names
    sample_key = pd.DataFrame({
        "feba_column": old_names,
        "analysis_sample": new_names_all,
    })
    sample_key["role"] = np.where(
        sample_key["analysis_sample"].str.split("-", n=1).str[0].eq("T0"),
        "Time0 baseline", "experimental",
    )
    sample_key["included_in_fit_analysis"] = sample_key["role"].eq("experimental")

    return (
        fit.loc[:, keep_columns].rename(columns=rename_map),
        fit_t.loc[:, keep_columns].rename(columns=rename_map),
        groups,
        experimental_new_names,
        sample_key,
    )


def write_results_guide(output_path, csv_path, sample_key, groups, consensus, threshold):
    """Write a compact guide and condition-level table for reviewing results."""
    sample_key = sample_key.copy()
    sample_key["replicate_group"] = sample_key["analysis_sample"].map(
        lambda sample: ";".join(
            condition for condition, samples in groups.items() if sample in samples
        )
    )
    sample_key.to_csv(csv_path / "sample_key.csv", index=False)

    condition_rows = []
    if consensus is not None and not consensus.empty:
        for condition, data in consensus.groupby("condition", sort=True):
            condition_rows.append({
                "condition": condition,
                "replicate_samples": data["samples"].iloc[0],
                "replicates_found": int(data["replicates_found"].iloc[0]),
                "genes_testable": int(data["replicates_evaluable"].gt(0).sum()),
                "high_confidence_calls": int(data["confidence_tier"].eq("high").sum()),
                "moderate_confidence_calls": int(data["confidence_tier"].eq("moderate").sum()),
                "exploratory_calls": int(data["confidence_tier"].eq("low").sum()),
                "negative_high_confidence": int(
                    (data["confidence_tier"].eq("high") & data["mean_fitness"].lt(0)).sum()
                ),
                "positive_high_confidence": int(
                    (data["confidence_tier"].eq("high") & data["mean_fitness"].gt(0)).sum()
                ),
            })
    condition_summary = pd.DataFrame(condition_rows)
    condition_summary.to_csv(csv_path / "condition_summary.csv", index=False)

    high_calls = 0 if consensus is None else int(consensus["confidence_tier"].eq("high").sum())
    lines = [
        "# Fitness-analysis results guide", "",
        f"Primary evidence threshold: **|t-like| > {threshold:g}**.", "",
        "## Start here", "",
        "- `csv/condition_summary.csv`: number of supported gene-condition calls per condition.",
        "- `csv/high_confidence_genes.csv`: main results; use this for genes supported by at least two evaluable, direction-consistent replicates.",
        "- `csv/moderate_confidence_genes.csv`: weaker but still replicated evidence.",
        "- `csv/replicate_consensus_all.csv`: complete gene-by-condition evidence table, including genes without support.",
        "- `csv/sample_key.csv`: mapping from FEBA columns to readable sample names and replicate groups.",
        "- `csv/replicate_pair_agreement.csv`: replicate agreement statistics.",
        "", "## What is included", "",
        f"- Experimental samples analyzed: {int(sample_key['included_in_fit_analysis'].sum())}",
        f"- Time0 baseline samples excluded from downstream fitness comparisons: {int((sample_key['role'] == 'Time0 baseline').sum())}",
        f"- Configured replicate groups: {len(groups)}",
        f"- High-confidence gene-condition calls: {high_calls}",
        "", "## Confidence labels", "",
        "- **High:** at least two evaluable replicates, consistent direction, |mean fitness| > 1, and at least two replicates with |t-like| > 2.5.",
        "- **Moderate:** at least two evaluable, direction-consistent replicates, |mean fitness| > 1, and at least one replicate with |t-like| > 2.",
        "- **Exploratory:** |mean fitness| > 1 with at least one replicate having |t-like| > 1.5.",
        "", "## Plots and their underlying tables", "",
        "| Plot | Main source table | Purpose |",
        "|---|---|---|",
        "| `plots/confidence_tier_counts.png` | `csv/replicate_consensus_all.csv` | Count high, moderate, exploratory, and unsupported calls by condition. |",
        "| `plots/candidate_counts_by_threshold.png` | `csv/threshold_sensitivity_by_condition.csv` | Show how candidate counts change with the t-like cutoff. |",
        "| `plots/threshold_sensitivity.png` | `csv/threshold_sensitivity_by_sample.csv`, `csv/threshold_sensitivity_overall.csv` | Show the effect of threshold choice across samples and overall. |",
        "| `plots/heatmap_filtered.png` | `csv/fit_filtered.csv` | Display gene-sample values that pass the selected primary threshold. |",
        "| `plots/heatmap_top10.png` | `csv/top_positive.csv`, `csv/top_negative.csv` | Show the ten strongest mean positive and negative sample-level effects. |",
        "| `plots/top_bottom_heatmap_t*.png` | `csv/fitness_results_long.csv` | Display genes with fitness above 1 or below −1 at each shown t-like cutoff. |",
        "| `plots/fitness_t_heatmap_t*.png` | `csv/fitness_results_long.csv` | Display thresholded fitness values with their t-like statistics. |",
        "| `plots/rep_cor_scatter.png`, `plots/rep_cor_heatmap.png`, `plots/replicate_bland_altman.png` | `csv/replicate_pair_agreement.csv` plus per-gene values in `csv/fitness_results_long.csv` | Replicate quality control. |",
        "| `plots/fitness_hist.png`, `plots/t-like_hist.png`, `plots/fitness_vs_t_scatter.png`, `plots/simple_pca_plot.png`, `plots/fitness_heatmap.png` | `csv/fitness_results_long.csv` | Sample-level distribution and pattern quality control. |",
    ]
    (output_path / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    LOGGER.info("Wrote results guide: %s", output_path / "README.md")
    LOGGER.info("Wrote sample key: %s", csv_path / "sample_key.csv")
    LOGGER.info("Wrote condition summary: %s", csv_path / "condition_summary.csv")


def save_figure(path, use_tight_layout=True):
    if use_tight_layout:
        plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    LOGGER.info("Wrote figure: %s", path)


def make_distribution_plot(data, output_file):
    """Plot comparable histograms and KDE curves for every sample."""
    n_plots = len(data.columns)
    n_cols = math.ceil(math.sqrt(n_plots))
    n_rows = math.ceil(n_plots / n_cols)
    _, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows),
                           sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    x_values = np.linspace(data.min().min(), data.max().max(), 300)

    for axis, column in zip(axes, data.columns):
        values = data[column].dropna()
        axis.hist(values, bins=30, density=True, alpha=0.6, edgecolor="black")
        if len(values) > 1 and values.nunique() > 1:
            kde = gaussian_kde(values)
            axis.plot(x_values, kde(x_values), color="red", linewidth=2)
        axis.set_title(column)
    for axis in axes[n_plots:]:
        axis.set_visible(False)
    save_figure(output_file)


def make_fitness_vs_t_plot(numeric_fit, numeric_t, threshold, output_file):
    """Plot fitness against its t-like statistic in a panel for every sample."""
    n_plots = len(numeric_fit.columns)
    n_cols = math.ceil(math.sqrt(n_plots))
    n_rows = math.ceil(n_plots / n_cols)
    _, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows),
                           sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()

    for axis, sample in zip(axes, numeric_fit.columns):
        paired = pd.concat(
            [numeric_fit[sample].rename("fitness"), numeric_t[sample].rename("t_like")],
            axis=1,
        ).dropna()
        significant = paired["t_like"].abs().gt(threshold)
        axis.scatter(paired.loc[~significant, "fitness"],
                     paired.loc[~significant, "t_like"],
                     s=7, alpha=0.25, color="0.55", linewidths=0,
                     label="Below threshold")
        axis.scatter(paired.loc[significant, "fitness"],
                     paired.loc[significant, "t_like"],
                     s=9, alpha=0.55, color="tab:red", linewidths=0,
                     label="Above threshold")
        axis.axhline(threshold, color="black", linestyle="--", linewidth=0.8,
                     label=f"t-like threshold: ±{threshold:g}")
        axis.axhline(-threshold, color="black", linestyle="--", linewidth=0.8,
                     label="_nolegend_")
        axis.axvline(1, color="tab:blue", linestyle=":", linewidth=1.0,
                     label="Fitness reference: ±1")
        axis.axvline(-1, color="tab:blue", linestyle=":", linewidth=1.0,
                     label="_nolegend_")
        axis.axvline(0, color="black", linewidth=0.5, alpha=0.5)
        x_ticks = np.unique(np.append(axis.get_xticks(), [-1.0, 0.0, 1.0]))
        y_ticks = np.unique(np.append(axis.get_yticks(), [-threshold, 0.0, threshold]))
        axis.set_xticks(x_ticks)
        axis.set_yticks(y_ticks)
        axis.set(
            title=sample,
            xlabel="Gene fitness (reference lines at −1 and 1)",
            ylabel=f"t-like statistic (threshold lines at ±{threshold:g})",
        )
        axis.text(0.03, 0.97, f"{significant.sum():,}/{len(paired):,} pass",
                  transform=axis.transAxes, ha="left", va="top", fontsize=8)

    for axis in axes[n_plots:]:
        axis.set_visible(False)
    if n_plots:
        axes[0].legend(loc="lower right", fontsize=7, frameon=False)
    save_figure(output_file)


def make_threshold_sensitivity_outputs(numeric_fit, numeric_t, csv_path, plot_path):
    """Summarize how retained gene-fitness results vary across t-like cutoffs."""
    per_sample_rows = []
    overall_rows = []

    for threshold in SENSITIVITY_THRESHOLDS:
        mask = numeric_t.abs().gt(threshold) & numeric_fit.notna()
        retained_fitness = numeric_fit.where(mask)
        for sample in numeric_fit.columns:
            sample_mask = mask[sample]
            sample_fitness = numeric_fit.loc[sample_mask, sample]
            per_sample_rows.append({
                "threshold": threshold,
                "sample": sample,
                "values_tested": int((numeric_t[sample].notna()
                                      & numeric_fit[sample].notna()).sum()),
                "values_retained": int(sample_mask.sum()),
                "positive_fitness": int(sample_fitness.gt(0).sum()),
                "negative_fitness": int(sample_fitness.lt(0).sum()),
            })

        genes_retained = retained_fitness.notna().any(axis=1)
        overall_rows.append({
            "threshold": threshold,
            "values_tested": int((numeric_t.notna() & numeric_fit.notna()).sum().sum()),
            "values_retained": int(mask.sum().sum()),
            "genes_retained": int(genes_retained.sum()),
            "positive_fitness": int((retained_fitness > 0).sum().sum()),
            "negative_fitness": int((retained_fitness < 0).sum().sum()),
        })

    per_sample = pd.DataFrame(per_sample_rows)
    per_sample["percent_retained"] = np.where(
        per_sample["values_tested"].gt(0),
        100 * per_sample["values_retained"] / per_sample["values_tested"],
        np.nan,
    )
    overall = pd.DataFrame(overall_rows)
    overall["percent_retained"] = np.where(
        overall["values_tested"].gt(0),
        100 * overall["values_retained"] / overall["values_tested"],
        np.nan,
    )
    per_sample.to_csv(csv_path / "threshold_sensitivity_by_sample.csv", index=False)
    overall.to_csv(csv_path / "threshold_sensitivity_overall.csv", index=False)
    LOGGER.info("Threshold-sensitivity summary:")
    for row in overall.itertuples(index=False):
        LOGGER.info(
            "  |t-like| > %-3g: %s values (%.2f%%), %s genes, %s positive, %s negative",
            row.threshold, f"{row.values_retained:,}", row.percent_retained,
            f"{row.genes_retained:,}", f"{row.positive_fitness:,}",
            f"{row.negative_fitness:,}",
        )
    LOGGER.info("Wrote table: %s", csv_path / "threshold_sensitivity_by_sample.csv")
    LOGGER.info("Wrote table: %s", csv_path / "threshold_sensitivity_overall.csv")

    _, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    sns.lineplot(data=per_sample, x="threshold", y="percent_retained",
                 hue="sample", marker="o", ax=axes[0])
    axes[0].set(title="Retained values by sample", xlabel="|t-like| threshold",
                ylabel="Values retained (%)")
    axes[0].legend(title="Sample", bbox_to_anchor=(1.02, 1), loc="upper left",
                   fontsize=7, title_fontsize=8)
    sns.lineplot(data=overall, x="threshold", y="genes_retained",
                 marker="o", color="tab:blue", ax=axes[1])
    axes[1].set(title="Genes with at least one retained value",
                xlabel="|t-like| threshold", ylabel="Genes retained")
    axes[1].set_xticks(SENSITIVITY_THRESHOLDS)
    save_figure(plot_path / "threshold_sensitivity.png")


def make_fitness_t_heatmaps(fit, numeric_fit, numeric_t, output_path):
    """Make threshold-specific heatmaps with fitness colors and fitness/t labels."""
    for threshold in COMBINED_HEATMAP_THRESHOLDS:
        mask = numeric_t.abs().gt(threshold) & numeric_fit.notna()
        retained_rows = mask.any(axis=1)
        if not retained_rows.any():
            LOGGER.warning(
                "Skipping combined heatmap at |t-like| > %g: no values pass", threshold
            )
            continue

        heat = numeric_fit.loc[retained_rows].where(mask.loc[retained_rows])
        t_values = numeric_t.loc[retained_rows].where(mask.loc[retained_rows])
        row_labels = (
            fit.loc[retained_rows, "locusId"].astype(str)
            + "|"
            + fit.loc[retained_rows, "desc"].fillna("").astype(str)
        )
        heat.index = row_labels
        t_values.index = row_labels

        # Put the highest mean fitness at the top and the lowest at the bottom,
        # matching the high-to-low direction of the fitness color scale.
        row_order = heat.mean(axis=1).sort_values(ascending=False).index
        heat = heat.loc[row_order]
        t_values = t_values.loc[row_order]
        annotations = pd.DataFrame("", index=heat.index, columns=heat.columns)
        for sample in heat.columns:
            present = heat[sample].notna() & t_values[sample].notna()
            annotations.loc[present, sample] = [
                f"F {fitness_value:.2f}\nT {t_value:.2f}"
                for fitness_value, t_value in zip(
                    heat.loc[present, sample], t_values.loc[present, sample]
                )
            ]

        # Scale with the number of retained genes, with an upper bound that
        # avoids impractically large figures at permissive thresholds.
        figure_height = min(40, max(8, 0.28 * len(heat)))
        figure_width = max(10, 1.5 * len(heat.columns))
        figure, axis = plt.subplots(figsize=(figure_width, figure_height))
        figure.subplots_adjust(top=0.93, bottom=0.12, left=0.22, right=0.98)
        colorbar_axis = figure.add_axes([0.81, 0.965, 0.16, 0.008])
        sns.heatmap(
            heat, ax=axis, center=0, cmap="RdYlGn", annot=annotations, fmt="",
            linewidths=0.25, linecolor="white",
            annot_kws={"fontsize": 5},
            cbar_ax=colorbar_axis,
            cbar_kws={"orientation": "horizontal"},
        )
        colorbar_axis.set_title("Gene fitness", fontsize=7, pad=2)
        colorbar_axis.tick_params(labelsize=6, length=2)
        plt.title(f"Fitness and t-like values (|t-like| > {threshold:g})")
        plt.xlabel("Sample")
        plt.ylabel("Gene")
        plt.xticks(rotation=45, ha="right", fontsize=8)
        plt.yticks(fontsize=5)
        threshold_label = str(threshold).replace(".", "p")
        save_figure(output_path / f"fitness_t_heatmap_t{threshold_label}.png",
                    use_tight_layout=False)


def make_threshold_top_bottom_heatmaps(fit, numeric_fit, numeric_t, output_path):
    """Plot all genes with fitness > 1 or < -1 at each selected t-like cutoff."""
    for threshold in COMBINED_HEATMAP_THRESHOLDS:
        mask = numeric_t.abs().gt(threshold) & numeric_fit.notna()
        retained = numeric_fit.where(mask).dropna(how="all")
        if retained.empty:
            LOGGER.warning(
                "Skipping top/bottom heatmap at |t-like| > %g: no values pass", threshold
            )
            continue

        # Match the top_positive.csv and top_negative.csv selection rule: a
        # gene qualifies when at least one value passing this t-like cutoff has
        # a fitness effect above 1 or below -1. Genes can occur in both panels.
        top = retained.loc[retained.gt(1).any(axis=1)].copy()
        bottom = retained.loc[retained.lt(-1).any(axis=1)].copy()

        for table in (top, bottom):
            table.index = (
                fit.loc[table.index, "locusId"].astype(str)
                + "|"
                + fit.loc[table.index, "desc"].fillna("").astype(str)
            )

        # Apply the same high-to-low ordering to both panels.
        top = top.loc[top.mean(axis=1).sort_values(ascending=False).index]
        bottom = bottom.loc[bottom.mean(axis=1).sort_values(ascending=False).index]
        plotted_values = pd.concat([top.stack(), bottom.stack()])
        color_limit = plotted_values.abs().max()
        if not np.isfinite(color_limit) or color_limit == 0:
            color_limit = 1.0

        figure_height = min(40, max(7, 0.28 * max(len(top), len(bottom))))
        figure, axes = plt.subplots(
            1, 2,
            figsize=(max(14, 1.5 * len(retained.columns)), figure_height),
        )
        figure.subplots_adjust(top=0.84, bottom=0.20, left=0.15, right=0.98,
                               wspace=0.45)
        colorbar_panel = 0 if not top.empty else 1
        colorbar_axis = figure.add_axes([0.84, 0.93, 0.12, 0.012])
        for panel_index, (axis, table, title) in enumerate(zip(
            axes,
            (top, bottom),
            ("Fitness > 1 in at least one retained value",
             "Fitness < -1 in at least one retained value"),
        )):
            if table.empty:
                axis.set_axis_off()
                axis.set_title(f"{title}\nNo qualifying genes")
                continue
            sns.heatmap(
                table, ax=axis, center=0, cmap="RdYlGn", annot=True, fmt=".2f",
                vmin=-color_limit, vmax=color_limit,
                annot_kws={"fontsize": 6},
                cbar=panel_index == colorbar_panel,
                cbar_ax=colorbar_axis if panel_index == colorbar_panel else None,
                cbar_kws={"orientation": "horizontal"},
            )
            axis.set(title=title, xlabel="Sample", ylabel="Gene")
            axis.tick_params(axis="x", labelrotation=45, labelsize=7)
            axis.tick_params(axis="y", labelsize=6)
            for label in axis.get_xticklabels():
                label.set_horizontalalignment("right")

        colorbar_axis.set_title("Gene fitness", fontsize=7, pad=2)
        colorbar_axis.tick_params(labelsize=6, length=2)

        plt.suptitle(
            f"Genes with strong fitness effects (|t-like| > {threshold:g})", y=1.02
        )
        LOGGER.info(
            "Top/bottom heatmap at |t-like| > %g includes %s positive and %s "
            "negative genes",
            threshold, f"{len(top):,}", f"{len(bottom):,}",
        )
        threshold_label = str(threshold).replace(".", "p")
        save_figure(output_path / f"top_bottom_heatmap_t{threshold_label}.png",
                    use_tight_layout=False)


def threshold_column_name(prefix, threshold):
    """Return a stable column suffix such as t1_5 or t2_0."""
    return f"{prefix}_t{threshold:.1f}".replace(".", "_")


def highest_threshold_passed(abs_t_like):
    """Return the strictest configured threshold exceeded by each value."""
    highest = pd.Series(np.nan, index=abs_t_like.index, dtype=float)
    for threshold in SENSITIVITY_THRESHOLDS:
        highest = highest.mask(abs_t_like.gt(threshold), threshold)
    return highest


def classify_confidence_tiers(fitness, t_like):
    """Classify genes using the confidence rules shared by tables and plots."""
    evaluable = fitness.notna() & t_like.notna()
    n_evaluable = evaluable.sum(axis=1)
    positive = (fitness.gt(0) & evaluable).sum(axis=1)
    negative = (fitness.lt(0) & evaluable).sum(axis=1)
    direction_consistent = (
        n_evaluable.gt(0)
        & (positive.eq(n_evaluable) | negative.eq(n_evaluable))
    )
    mean_fitness = fitness.where(evaluable).mean(axis=1)
    high = (
        n_evaluable.ge(2)
        & direction_consistent
        & mean_fitness.abs().gt(1)
        & (t_like.abs().gt(2.5) & evaluable).sum(axis=1).ge(2)
    )
    moderate = (
        n_evaluable.ge(2)
        & direction_consistent
        & mean_fitness.abs().gt(1)
        & (t_like.abs().gt(2.0) & evaluable).sum(axis=1).ge(1)
    )
    exploratory = (
        mean_fitness.abs().gt(1)
        & (t_like.abs().gt(1.5) & evaluable).sum(axis=1).ge(1)
    )
    return pd.Series(
        np.select(
            [high, moderate, exploratory],
            ["high", "moderate", "low"],
            default="not_supported",
        ),
        index=fitness.index,
        name="confidence_tier",
    )


def make_sample_qc_summary(numeric_fit, numeric_t, primary_threshold, csv_path):
    """Write one row of measurement and threshold QC metrics per sample."""
    rows = []
    for sample in numeric_fit.columns:
        fitness = numeric_fit[sample]
        t_like = numeric_t[sample]
        evaluable = fitness.notna() & t_like.notna()
        row = {
            "sample": sample,
            "genes_total": len(fitness),
            "genes_with_fitness": int(fitness.notna().sum()),
            "genes_with_t_like": int(t_like.notna().sum()),
            "complete_pairs": int(evaluable.sum()),
            "missing_pairs": int((~evaluable).sum()),
            "percent_evaluable": 100 * evaluable.mean(),
            "median_abs_fitness": fitness[evaluable].abs().median(),
            "median_abs_t_like": t_like[evaluable].abs().median(),
            "max_abs_fitness": fitness[evaluable].abs().max(),
            "max_abs_t_like": t_like[evaluable].abs().max(),
        }
        for threshold in SENSITIVITY_THRESHOLDS:
            passes = evaluable & t_like.abs().gt(threshold)
            row[threshold_column_name("passing", threshold)] = int(passes.sum())
        primary_passes = evaluable & t_like.abs().gt(primary_threshold)
        row["passing_primary_threshold"] = int(primary_passes.sum())
        row["primary_positive"] = int((primary_passes & fitness.gt(0)).sum())
        row["primary_negative"] = int((primary_passes & fitness.lt(0)).sum())
        row["primary_abs_fitness_gt_1"] = int(
            (primary_passes & fitness.abs().gt(1)).sum()
        )
        rows.append(row)

    summary = pd.DataFrame(rows)
    output_file = csv_path / "sample_qc_summary.csv"
    summary.to_csv(output_file, index=False)
    LOGGER.info("Wrote sample QC for %s samples: %s", f"{len(summary):,}", output_file)


def make_long_results(fit, numeric_fit, numeric_t, groups, primary_threshold, csv_path):
    """Write all gene-sample measurements and non-destructive review flags."""
    sample_groups = {
        sample: ";".join(condition for condition, samples in groups.items() if sample in samples)
        for sample in numeric_fit.columns
    }
    frames = []
    for sample in numeric_fit.columns:
        frame = fit[ID_COLUMNS].copy()
        frame["sample"] = sample
        frame["condition"] = sample_groups[sample]
        frame["fitness"] = numeric_fit[sample]
        frame["t_like"] = numeric_t[sample]
        frame["abs_fitness"] = frame["fitness"].abs()
        frame["abs_t_like"] = frame["t_like"].abs()
        frame["evaluable"] = frame["fitness"].notna() & frame["t_like"].notna()
        frame["direction"] = np.select(
            [frame["fitness"].gt(0), frame["fitness"].lt(0)],
            ["positive", "negative"], default="zero_or_missing",
        )
        for threshold in SENSITIVITY_THRESHOLDS:
            frame[threshold_column_name("passes", threshold)] = (
                frame["evaluable"] & frame["abs_t_like"].gt(threshold)
            )
        frame["passes_primary_threshold"] = (
            frame["evaluable"] & frame["abs_t_like"].gt(primary_threshold)
        )
        frame["fitness_gt_0_5"] = frame["abs_fitness"].gt(0.5)
        frame["fitness_gt_1"] = frame["abs_fitness"].gt(1.0)
        frame["fitness_gt_2"] = frame["abs_fitness"].gt(2.0)
        frame["highest_tested_t_threshold_passed"] = highest_threshold_passed(
            frame["abs_t_like"]
        )
        frames.append(frame)

    long_results = pd.concat(frames, ignore_index=True)
    output_file = csv_path / "fitness_results_long.csv"
    long_results.to_csv(output_file, index=False)
    LOGGER.info("Wrote %s unfiltered gene-sample records: %s",
                f"{len(long_results):,}", output_file)


def make_replicate_consensus(fit, numeric_fit, numeric_t, groups,
                             primary_threshold, csv_path):
    """Summarize unfiltered measurements across each configured replicate group."""
    frames = []
    for condition, configured_samples in groups.items():
        samples = [sample for sample in configured_samples if sample in numeric_fit]
        if not samples:
            LOGGER.warning("No recognized samples for replicate group: %s", condition)
            continue
        fitness = numeric_fit[samples]
        t_like = numeric_t[samples]
        evaluable = fitness.notna() & t_like.notna()
        n_evaluable = evaluable.sum(axis=1)
        positive = (fitness.gt(0) & evaluable).sum(axis=1)
        negative = (fitness.lt(0) & evaluable).sum(axis=1)
        consensus = fit[ID_COLUMNS].copy()
        consensus["condition"] = condition
        consensus["samples"] = ";".join(samples)
        consensus["replicates_configured"] = len(configured_samples)
        consensus["replicates_found"] = len(samples)
        consensus["replicates_evaluable"] = n_evaluable
        consensus["mean_fitness"] = fitness.where(evaluable).mean(axis=1)
        consensus["fitness_sd"] = fitness.where(evaluable).std(axis=1, ddof=1)
        consensus["min_fitness"] = fitness.where(evaluable).min(axis=1)
        consensus["max_fitness"] = fitness.where(evaluable).max(axis=1)
        consensus["mean_abs_t_like"] = t_like.where(evaluable).abs().mean(axis=1)
        consensus["max_abs_t_like"] = t_like.where(evaluable).abs().max(axis=1)
        consensus["positive_replicates"] = positive
        consensus["negative_replicates"] = negative
        consensus["direction_consistent"] = (
            n_evaluable.gt(0) & (positive.eq(n_evaluable) | negative.eq(n_evaluable))
        )
        for threshold in SENSITIVITY_THRESHOLDS:
            pass_count = (t_like.abs().gt(threshold) & evaluable).sum(axis=1)
            consensus[threshold_column_name("replicates_passing", threshold)] = pass_count
            consensus[threshold_column_name("any_passes", threshold)] = pass_count.gt(0)
            consensus[threshold_column_name("all_evaluable_pass", threshold)] = (
                n_evaluable.gt(0) & pass_count.eq(n_evaluable)
            )
        primary_count = (t_like.abs().gt(primary_threshold) & evaluable).sum(axis=1)
        consensus["replicates_passing_primary"] = primary_count
        consensus["any_replicate_passes_primary"] = primary_count.gt(0)
        consensus["all_evaluable_replicates_pass_primary"] = (
            n_evaluable.gt(0) & primary_count.eq(n_evaluable)
        )
        consensus["abs_mean_fitness_gt_1"] = consensus["mean_fitness"].abs().gt(1)
        consensus["coverage_status"] = np.select(
            [n_evaluable.eq(0), n_evaluable.eq(1), n_evaluable.ge(2)],
            ["not_testable", "limited", "adequate"],
            default="not_testable",
        )
        consensus["confidence_tier"] = classify_confidence_tiers(fitness, t_like)
        frames.append(consensus)

    if not frames:
        LOGGER.warning("Skipping replicate consensus: no group samples were recognized")
        return None
    consensus = pd.concat(frames, ignore_index=True)
    confidence_order = {"high": 1, "moderate": 2, "low": 3, "not_supported": 4}
    consensus["confidence_rank"] = consensus["confidence_tier"].map(confidence_order)
    consensus["evidence_summary"] = consensus.apply(
        lambda row: (
            f"{int(row['replicates_evaluable'])}/{int(row['replicates_found'])} replicates "
            f"evaluable; direction_consistent={row['direction_consistent']}; "
            f"|mean fitness|={abs(row['mean_fitness']):.3g}; "
            f"max |t|={row['max_abs_t_like']:.3g}; "
            f"primary passes={int(row['replicates_passing_primary'])}"
        ) if row["replicates_evaluable"] else "No evaluable fitness/t-like pair",
        axis=1,
    )
    consensus = consensus.sort_values(
        ["condition", "confidence_rank", "abs_mean_fitness_gt_1",
         "mean_abs_t_like", "mean_fitness"],
        ascending=[True, True, False, False, True],
    )
    consensus["candidate_rank_within_condition"] = (
        consensus.groupby("condition").cumcount() + 1
    )
    all_file = csv_path / "replicate_consensus_all.csv"
    high_file = csv_path / "high_confidence_genes.csv"
    moderate_file = csv_path / "moderate_confidence_genes.csv"
    exploratory_file = csv_path / "exploratory_candidates.csv"
    unresolved_file = csv_path / "unresolved_genes.csv"
    consensus.to_csv(all_file, index=False)
    consensus[consensus["confidence_tier"].eq("high")].to_csv(high_file, index=False)
    consensus[consensus["confidence_tier"].eq("moderate")].to_csv(
        moderate_file, index=False
    )
    consensus[consensus["confidence_tier"].eq("low")].to_csv(
        exploratory_file, index=False
    )
    unresolved = consensus[consensus["coverage_status"].ne("adequate")].copy()
    unresolved["unresolved_reason"] = np.where(
        unresolved["replicates_evaluable"].eq(0),
        "no evaluable fitness/t-like pairs",
        "only one evaluable replicate",
    )
    unresolved.to_csv(unresolved_file, index=False)
    LOGGER.info("Wrote replicate consensus: %s", all_file)
    LOGGER.info("Wrote %s high-confidence condition-gene rows: %s",
                f"{consensus['confidence_tier'].eq('high').sum():,}", high_file)
    LOGGER.info("Wrote %s moderate-confidence condition-gene rows: %s",
                f"{consensus['confidence_tier'].eq('moderate').sum():,}", moderate_file)
    LOGGER.info("Wrote %s exploratory condition-gene rows: %s",
                f"{consensus['confidence_tier'].eq('low').sum():,}", exploratory_file)
    LOGGER.info("Wrote %s unresolved condition-gene rows: %s",
                f"{len(unresolved):,}", unresolved_file)
    return consensus


def make_condition_candidate_sensitivity(consensus, csv_path, plot_path):
    """Summarize threshold effects after accounting for replicate support."""
    if consensus is None or consensus.empty:
        LOGGER.warning("Skipping condition candidate sensitivity: no consensus data")
        return

    rows = []
    for condition, condition_data in consensus.groupby("condition", sort=False):
        for threshold in SENSITIVITY_THRESHOLDS:
            pass_column = threshold_column_name("replicates_passing", threshold)
            any_pass = condition_data[pass_column].gt(0)
            rows.append({
                "condition": condition,
                "threshold": threshold,
                "genes_total": len(condition_data),
                "genes_testable": int(condition_data["replicates_evaluable"].gt(0).sum()),
                "genes_with_any_replicate_passing": int(any_pass.sum()),
                "genes_passing_and_abs_mean_fitness_gt_1": int(
                    (any_pass & condition_data["abs_mean_fitness_gt_1"]).sum()
                ),
                "genes_passing_and_direction_consistent": int(
                    (any_pass & condition_data["direction_consistent"]).sum()
                ),
                "genes_passing_with_2plus_replicates": int(
                    (condition_data[pass_column].ge(2)
                     & condition_data["direction_consistent"]).sum()
                ),
            })

    summary = pd.DataFrame(rows)
    output_file = csv_path / "threshold_sensitivity_by_condition.csv"
    summary.to_csv(output_file, index=False)
    LOGGER.info("Wrote condition-level threshold sensitivity: %s", output_file)

    plot_data = summary.melt(
        id_vars=["condition", "threshold"],
        value_vars=[
            "genes_with_any_replicate_passing",
            "genes_passing_and_abs_mean_fitness_gt_1",
            "genes_passing_with_2plus_replicates",
        ],
        var_name="candidate_rule", value_name="genes",
    )
    labels = {
        "genes_with_any_replicate_passing": "Any replicate passes",
        "genes_passing_and_abs_mean_fitness_gt_1": "Also |mean fitness| > 1",
        "genes_passing_with_2plus_replicates": "2+ passing, consistent direction",
    }
    plot_data["candidate_rule"] = plot_data["candidate_rule"].map(labels)
    grid = sns.relplot(
        data=plot_data, x="threshold", y="genes", hue="candidate_rule",
        col="condition", col_wrap=3, kind="line", marker="o",
        facet_kws={"sharey": False}, height=3.2, aspect=1.15,
    )
    grid.set_axis_labels("|t-like| threshold", "Genes retained")
    grid.set_titles("{col_name}")
    grid.fig.subplots_adjust(top=0.90)
    grid.fig.suptitle("Candidate sensitivity after replicate and effect-size checks")
    grid.savefig(plot_path / "candidate_counts_by_threshold.png",
                 dpi=300, bbox_inches="tight")
    plt.close(grid.fig)
    LOGGER.info("Wrote figure: %s", plot_path / "candidate_counts_by_threshold.png")


def make_confidence_tier_plot(consensus, output_file):
    """Plot confidence-tier counts for every replicate group."""
    if consensus is None or consensus.empty:
        return
    counts = (
        consensus.groupby(["condition", "confidence_tier"]).size()
        .rename("genes").reset_index()
    )
    if counts.empty:
        LOGGER.warning("Skipping confidence-tier plot: no supported candidates")
        return
    plt.figure(figsize=(max(8, 1.1 * counts["condition"].nunique()), 5))
    ax = sns.barplot(
        data=counts, x="condition", y="genes", hue="confidence_tier",
        hue_order=["high", "moderate", "low", "not_supported"],
        palette={
            "high": "#2166ac",
            "moderate": "#fdae61",
            "low": "#7b4ab5",
            "not_supported": "#bdbdbd",
        },
    )
    for container in ax.containers:
        ax.bar_label(
            container,
            labels=[f"{int(bar.get_height())}" if bar.get_height() > 0 else ""
                    for bar in container],
            padding=2,
            fontsize=8,
        )
    ax.margins(y=0.12)
    plt.title("Gene-condition results by confidence tier")
    plt.xlabel("Condition")
    plt.ylabel("Gene-condition results")
    plt.xticks(rotation=45, ha="right")
    plt.legend(
        title="Confidence tier",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=True,
    )
    plt.subplots_adjust(right=0.80)
    save_figure(output_file)


def make_replicate_outputs(numeric_fit, numeric_t, groups, threshold,
                           csv_path, plot_path):
    """Create correlation and Bland-Altman outputs from unfiltered replicates."""
    valid_groups = {
        condition: [sample for sample in samples if sample in numeric_fit]
        for condition, samples in groups.items()
    }
    valid_groups = {key: value for key, value in valid_groups.items() if len(value) >= 2}
    replicate_pairs = [
        (condition, sample_1, sample_2)
        for condition, samples in valid_groups.items()
        for sample_1, sample_2 in combinations(samples, 2)
    ]
    results = []
    valid_pairs = []
    for condition, sample_1, sample_2 in replicate_pairs:
        paired = numeric_fit[[sample_1, sample_2]].dropna()
        if len(paired) < 2:
            LOGGER.warning("Skipping replicate pair %s/%s: fewer than two genes",
                           sample_1, sample_2)
            continue
        difference = paired[sample_1] - paired[sample_2]
        difference_sd = difference.std(ddof=1)
        results.append({
            "condition": condition, "replicate_1": sample_1, "replicate_2": sample_2,
            "genes": len(paired),
            "pearson": pearsonr(paired[sample_1], paired[sample_2]).statistic,
            "spearman": spearmanr(paired[sample_1], paired[sample_2]).statistic,
            "mean_difference": difference.mean(),
            "difference_sd": difference_sd,
            "lower_95_agreement": difference.mean() - 1.96 * difference_sd,
            "upper_95_agreement": difference.mean() + 1.96 * difference_sd,
        })
        valid_pairs.append((condition, sample_1, sample_2))
    agreement_file = csv_path / "replicate_pair_agreement.csv"
    pd.DataFrame(results).to_csv(agreement_file, index=False)
    LOGGER.info("Wrote replicate pair agreement: %s", agreement_file)

    if valid_pairs:
        n_plots = len(valid_pairs)
        n_cols = math.ceil(math.sqrt(n_plots))
        n_rows = math.ceil(n_plots / n_cols)
        _, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
        axes = np.atleast_1d(axes).ravel()
        for axis, (condition, sample_1, sample_2) in zip(axes, valid_pairs):
            x, y = numeric_fit[sample_1], numeric_fit[sample_2]
            axis.scatter(x, y, s=5, alpha=0.3)
            axis.axline((0, 0), slope=1, linewidth=1)
            axis.text(0.05, 0.95, f"Spearman r = {x.corr(y, method='spearman'):.2f}",
                      transform=axis.transAxes, ha="left", va="top")
            axis.set(xlabel=sample_1, ylabel=sample_2, title=condition)
        for axis in axes[n_plots:]:
            axis.set_visible(False)
        save_figure(plot_path / "rep_cor_scatter.png")

        _, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
        axes = np.atleast_1d(axes).ravel()
        for axis, (condition, sample_1, sample_2) in zip(axes, valid_pairs):
            group_samples = valid_groups[condition]
            confidence_tier = classify_confidence_tiers(
                numeric_fit[group_samples], numeric_t[group_samples]
            )
            paired = pd.concat(
                [numeric_fit[sample_1].rename("fit_1"),
                 numeric_fit[sample_2].rename("fit_2"),
                 numeric_t[sample_1].rename("t_1"),
                 numeric_t[sample_2].rename("t_2")], axis=1,
            ).dropna(subset=["fit_1", "fit_2"])
            paired["mean"] = paired[["fit_1", "fit_2"]].mean(axis=1)
            paired["difference"] = paired["fit_1"] - paired["fit_2"]
            paired["confidence_tier"] = confidence_tier.reindex(paired.index).fillna(
                "not_supported"
            )
            tier_styles = {
                "not_supported": ("0.72", 7, 0.22),
                "low": ("tab:purple", 9, 0.42),
                "moderate": ("tab:orange", 11, 0.55),
                "high": ("tab:blue", 13, 0.65),
            }
            for tier in ("not_supported", "low", "moderate", "high"):
                color, size, alpha = tier_styles[tier]
                selected = paired["confidence_tier"].eq(tier)
                axis.scatter(
                    paired.loc[selected, "mean"], paired.loc[selected, "difference"],
                    s=size, alpha=alpha, color=color, linewidths=0,
                )
            bias = paired["difference"].mean()
            difference_sd = paired["difference"].std(ddof=1)
            axis.axhline(bias, color="black", linewidth=1)
            axis.axhline(bias + 1.96 * difference_sd, color="black", linestyle="--",
                         linewidth=0.8)
            axis.axhline(bias - 1.96 * difference_sd, color="black", linestyle="--",
                         linewidth=0.8)
            axis.set(title=condition, xlabel="Mean replicate fitness",
                     ylabel=f"Fitness difference\n{sample_1} - {sample_2}")
        legend_handles = [
            Line2D([0], [0], marker="o", linestyle="none", color="tab:blue",
                   markersize=4, label="High confidence"),
            Line2D([0], [0], marker="o", linestyle="none", color="tab:orange",
                   markersize=4, label="Moderate confidence"),
            Line2D([0], [0], marker="o", linestyle="none", color="tab:purple",
                   markersize=4, label="Exploratory (low confidence)"),
            Line2D([0], [0], marker="o", linestyle="none", color="0.72",
                   markersize=4, label="Not supported"),
            Line2D([0], [0], color="black", linewidth=1,
                   label="Solid line: mean replicate difference (bias)"),
            Line2D([0], [0], color="black", linestyle="--", linewidth=0.8,
                   label="Dashed lines: bias ± 1.96 SD (95% agreement limits)"),
        ]
        axes[0].legend(handles=legend_handles, fontsize=6, frameon=True, loc="best")
        plt.suptitle(
            "Bland–Altman replicate agreement; colors use the confidence tiers "
            "reported in the CSV outputs"
        )
        for axis in axes[n_plots:]:
            axis.set_visible(False)
        save_figure(plot_path / "replicate_bland_altman.png")
    else:
        LOGGER.warning("No replicate pairs were available for replicate plots")

    plt.figure(figsize=(8, 6))
    sns.heatmap(numeric_fit.corr(method="spearman"), annot=True, fmt=".2f",
                annot_kws={"fontsize": 6}, vmin=-1, vmax=1, center=0, cmap="vlag")
    plt.title("Spearman correlation of fitness values")
    save_figure(plot_path / "rep_cor_heatmap.png")


def make_pca_plot(numeric_fit, output_file):
    """Plot samples in the first two principal-component dimensions."""
    matrix = numeric_fit.T.dropna(axis=1)
    if min(matrix.shape) < 2:
        LOGGER.warning("Skipping PCA: at least two samples and two complete genes are required")
        return
    pca = PCA(n_components=2)
    components = pca.fit_transform(matrix)
    plt.figure(figsize=(8, 6))
    plt.scatter(components[:, 0], components[:, 1])
    for index, sample in enumerate(matrix.index):
        plt.annotate(sample, components[index, :2])
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    plt.title("PCA of gene fitness")
    save_figure(output_file)


def make_heatmaps(fit, numeric_fit, filtered, top_negative, top_positive, output_path):
    """Create complete, thresholded, and extreme-fitness heatmaps."""
    heat = numeric_fit.copy()
    heat.index = fit["locusId"]
    heat = heat.loc[heat.mean(axis=1).sort_values(ascending=False).index]
    # Hierarchical clustering cannot calculate distances from missing values.
    # Keep the missing cells visible, but disable column clustering rather than
    # aborting the complete analysis for a partially covered library.
    cluster_columns = not heat.isna().any().any() and heat.shape[1] > 1
    if not cluster_columns and heat.isna().any().any():
        LOGGER.warning(
            "fitness_heatmap.png contains missing values; displaying columns in input "
            "order because hierarchical column clustering requires complete data"
        )
    grid = sns.clustermap(heat, cmap="RdYlGn", center=0, figsize=(8, 10),
                          row_cluster=False, col_cluster=cluster_columns, xticklabels=True,
                          yticklabels=False, linewidths=0,
                          dendrogram_ratio=(0.15, 0.15),
                          cbar_kws={"label": "Gene fitness"})
    plt.setp(grid.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=10)
    grid.savefig(output_path / "fitness_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(grid.fig)

    filtered_heat = filtered.copy()
    filtered_heat.index = fit.loc[filtered.index, "locusId"]
    filtered_heat = filtered_heat.loc[filtered_heat.mean(axis=1).sort_values(ascending=False).index]
    if not filtered_heat.empty:
        plt.figure(figsize=(12, 10))
        sns.heatmap(filtered_heat, center=0, cmap="RdYlGn", annot=True,
                    annot_kws={"fontsize": 5})
        plt.xlabel("Condition")
        plt.ylabel("Gene")
        plt.yticks(fontsize=5)
        save_figure(output_path / "heatmap_filtered.png")

    _, axes = plt.subplots(1, 2, figsize=(14, 5))
    for axis, table, title in zip(
        axes, [top_negative, top_positive], ["Bottom 10 Mean Fitness", "Top 10 Mean Fitness"]
    ):
        values = table[numeric_fit.columns]
        values.index = table["locusId"].astype(str) + "|" + table["desc"].fillna("")
        sns.heatmap(values, ax=axis, center=0, cmap="RdYlGn", annot=True,
                    fmt=".2f", annot_kws={"fontsize": 4})
        axis.set(title=title, xlabel="Condition", ylabel="Gene")
        axis.tick_params(axis="both", labelsize=5)
    save_figure(output_path / "heatmap_top10.png")


def save_barcode_outputs(input_path, fit, filtered_full, output_path):
    """Extract barcode counts for genes retained by the t-like filter."""
    all_pool_file = input_path / "all.poolcount"
    if not all_pool_file.exists():
        # Compatibility for the previous layout where fit tables were stored
        # one directory beneath all.poolcount.
        all_pool_file = input_path.parent / "all.poolcount"
    if not all_pool_file.exists():
        LOGGER.warning("Skipping barcode outputs: %s was not found", all_pool_file)
        return
    all_pool = pd.read_csv(all_pool_file, sep="\t")
    selected = all_pool[all_pool["locusId"].isin(filtered_full["locusId"])].copy()
    count_columns = selected.columns[7:]
    selected = selected[selected[count_columns].gt(0).any(axis=1)].copy()
    selected.to_csv(output_path / "selected_barcodes.csv", index=False)
    selected["nTot"] = selected[count_columns].sum(axis=1)
    distribution = selected.groupby("locusId")["nTot"].describe()
    distribution["desc"] = fit.groupby("locusId")["desc"].first()
    distribution.sort_values("max", ascending=False).to_csv(
        output_path / "barcode_distribution.csv"
    )
    LOGGER.info("Selected %s barcodes across %s genes", f"{len(selected):,}",
                f"{selected['locusId'].nunique():,}")
    LOGGER.info("Wrote table: %s", output_path / "selected_barcodes.csv")
    LOGGER.info("Wrote table: %s", output_path / "barcode_distribution.csv")


def main():
    start_time = time.monotonic()
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    csv_path = args.output / "csv"
    plot_path = args.output / "plots"
    csv_path.mkdir(parents=True, exist_ok=True)
    plot_path.mkdir(parents=True, exist_ok=True)
    setup_logging(args.output)
    LOGGER.info("Starting FEBA fitness analysis")
    LOGGER.info("Input directory: %s", args.input.resolve())
    LOGGER.info("Output directory: %s", args.output.resolve())
    LOGGER.info("CSV directory: %s", csv_path.resolve())
    LOGGER.info("Plot directory: %s", plot_path.resolve())
    LOGGER.info("Name-map file: %s", args.name_map.resolve())
    LOGGER.info("Groups file: %s", args.groups.resolve())
    LOGGER.info("Primary threshold: |t-like| > %g", args.threshold)
    fit, fit_t, groups, condition_columns, sample_key = load_inputs(args)
    numeric_fit = fit[condition_columns].apply(pd.to_numeric, errors="coerce")
    numeric_t = fit_t[condition_columns].apply(pd.to_numeric, errors="coerce")
    paired_values = numeric_fit.notna() & numeric_t.notna()
    LOGGER.info("Loaded %s genes, %s samples, and %s replicate groups",
                f"{len(fit):,}", f"{len(condition_columns):,}", f"{len(groups):,}")
    LOGGER.info("Samples: %s", ", ".join(condition_columns))
    LOGGER.info("Complete fitness/t-like pairs: %s of %s (%.2f%%)",
                f"{paired_values.sum().sum():,}", f"{paired_values.size:,}",
                100 * paired_values.sum().sum() / paired_values.size)
    LOGGER.info("Fitness range: %.4g to %.4g", numeric_fit.min().min(),
                numeric_fit.max().max())
    LOGGER.info("T-like range: %.4g to %.4g", numeric_t.min().min(),
                numeric_t.max().max())

    # Filtering is performed per gene-condition value. A gene remains in the
    # table if at least one condition passes the absolute t-like threshold.
    significant_mask = numeric_t.abs().gt(args.threshold) & numeric_fit.notna()
    filtered = numeric_fit.where(significant_mask).dropna(how="all")
    filtered_full = pd.concat(
        [fit.loc[filtered.index, ["locusId", "desc"]], filtered], axis=1
    )
    filtered_full["mean_fit"] = filtered.mean(axis=1)
    filtered_full = filtered_full.sort_values("mean_fit")
    filtered_full.to_csv(csv_path / "fit_filtered.csv", index=False)

    # Export every gene with a strong fitness effect in at least one value that
    # also passes the run's primary t-like filter. A gene can appear in both
    # files if it has fitness > 1 in one condition and fitness < -1 in another.
    positive_candidates = filtered_full.loc[
        filtered_full[condition_columns].gt(1).any(axis=1)
    ].copy()
    positive_candidates = positive_candidates.loc[
        positive_candidates[condition_columns].max(axis=1)
        .sort_values(ascending=False).index
    ]
    negative_candidates = filtered_full.loc[
        filtered_full[condition_columns].lt(-1).any(axis=1)
    ].copy()
    negative_candidates = negative_candidates.loc[
        negative_candidates[condition_columns].min(axis=1)
        .sort_values(ascending=True).index
    ]
    positive_candidates.to_csv(csv_path / "top_positive.csv", index=False)
    negative_candidates.to_csv(csv_path / "top_negative.csv", index=False)

    # Keep the existing ten-gene summary heatmap compact. The CSV exports above
    # contain the complete positive and negative candidate sets.
    top_positive = filtered_full[filtered_full["mean_fit"] > 0].nlargest(10, "mean_fit")
    top_negative = filtered_full[filtered_full["mean_fit"] < 0].nsmallest(10, "mean_fit")

    LOGGER.info("Primary filter retained %s gene-sample values in %s genes",
                f"{significant_mask.sum().sum():,}", f"{len(filtered):,}")
    LOGGER.info("Primary-filter direction: %s positive and %s negative values",
                f"{numeric_fit.where(significant_mask).gt(0).sum().sum():,}",
                f"{numeric_fit.where(significant_mask).lt(0).sum().sum():,}")
    LOGGER.info("Wrote table: %s", csv_path / "fit_filtered.csv")
    LOGGER.info(
        "Wrote %s genes with at least one primary-supported fitness value > 1: %s",
        f"{len(positive_candidates):,}", csv_path / "top_positive.csv"
    )
    LOGGER.info(
        "Wrote %s genes with at least one primary-supported fitness value < -1: %s",
        f"{len(negative_candidates):,}", csv_path / "top_negative.csv"
    )

    make_sample_qc_summary(numeric_fit, numeric_t, args.threshold, csv_path)
    make_long_results(fit, numeric_fit, numeric_t, groups, args.threshold, csv_path)
    consensus = make_replicate_consensus(
        fit, numeric_fit, numeric_t, groups, args.threshold, csv_path
    )
    write_results_guide(
        args.output, csv_path, sample_key, groups, consensus, args.threshold
    )
    make_condition_candidate_sensitivity(consensus, csv_path, plot_path)
    make_confidence_tier_plot(consensus, plot_path / "confidence_tier_counts.png")
    save_barcode_outputs(args.input, fit, filtered_full, csv_path)
    make_distribution_plot(numeric_fit, plot_path / "fitness_hist.png")
    make_distribution_plot(numeric_t, plot_path / "t-like_hist.png")
    make_fitness_vs_t_plot(
        numeric_fit, numeric_t, args.threshold, plot_path / "fitness_vs_t_scatter.png"
    )
    make_threshold_sensitivity_outputs(numeric_fit, numeric_t, csv_path, plot_path)
    make_fitness_t_heatmaps(fit, numeric_fit, numeric_t, plot_path)
    make_threshold_top_bottom_heatmaps(fit, numeric_fit, numeric_t, plot_path)
    make_replicate_outputs(
        numeric_fit, numeric_t, groups, args.threshold, csv_path, plot_path
    )
    make_pca_plot(numeric_fit, plot_path / "simple_pca_plot.png")
    make_heatmaps(fit, numeric_fit, filtered, top_negative, top_positive, plot_path)
    LOGGER.info("Analysis completed successfully in %.1f seconds", time.monotonic() - start_time)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        LOGGER.exception("Analysis failed")
        raise
