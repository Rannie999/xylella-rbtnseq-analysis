#!/usr/bin/env python3
"""Run pairwise Welch tests across all valid FEBA gene-fitness measurements.

FDR correction is performed separately within each condition comparison over
every valid gene-level test. FEBA t-like values are joined afterward as
within-condition evidence; they never prefilter the hypotheses used to
calculate p- or q-values.
"""

import argparse
import json
import logging
import math
import re
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


ID_COLUMNS = ["locusId", "sysName", "desc"]
RESULT_COLUMNS = [
    "locusId", "sysName", "desc", "condition_1", "condition_2",
    "replicates_1", "replicates_2", "mean_1", "mean_2",
    "fitness_difference", "t_statistic", "p_value", "q_value", "fdr_test_count",
]
LOGGER = logging.getLogger("pairwise_fitness")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="FEBA output directory containing fit_logratios.tab and fit_t.tab")
    parser.add_argument("-f", "--filtered", type=Path,
                        help="Optional fit_filtered.csv for a clearly labelled exploratory analysis")
    parser.add_argument("-m", "--name-map", dest="name_map", type=Path, required=True,
                        help="CSV containing the new_colnames column")
    parser.add_argument("-g", "--groups", type=Path, required=True,
                        help="JSON mapping conditions to replicate samples")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output directory")
    return parser.parse_args()


def setup_logging(output_path):
    LOGGER.setLevel(logging.INFO)
    LOGGER.handlers.clear()
    LOGGER.propagate = False
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s",
                                  datefmt="%Y-%m-%d %H:%M:%S")
    handlers = [
        logging.StreamHandler(sys.stderr),
        logging.FileHandler(output_path / "pairwise_ttest.log", mode="w", encoding="utf-8"),
    ]
    for handler in handlers:
        handler.setFormatter(formatter)
        LOGGER.addHandler(handler)


def load_inputs(args):
    raw_path = args.input / "fit_logratios.tab"
    t_like_path = args.input / "fit_t.tab"
    if not raw_path.is_file() or not t_like_path.is_file():
        raise FileNotFoundError(
            "--input must be a FEBA output directory containing "
            "fit_logratios.tab and fit_t.tab"
        )
    raw = pd.read_csv(raw_path, sep="\t")
    t_like = pd.read_csv(t_like_path, sep="\t")
    name_map = pd.read_csv(args.name_map)
    with args.groups.open(encoding="utf-8") as handle:
        groups = json.load(handle)

    if not set(ID_COLUMNS).issubset(raw.columns):
        raise ValueError(f"Unfiltered table must contain {ID_COLUMNS}")
    if raw.columns.tolist() != t_like.columns.tolist():
        raise ValueError("fit_logratios.tab and fit_t.tab must have identical columns")
    if "new_colnames" not in name_map:
        raise ValueError("Name-map must contain a 'new_colnames' column")
    old_names = raw.columns[len(ID_COLUMNS):].tolist()
    name_map = name_map.dropna(subset=["new_colnames"]).copy()
    new_names_all = name_map["new_colnames"].astype(str).tolist()
    if len(new_names_all) != len(old_names):
        raise ValueError(
            f"Column mismatch: {len(new_names_all)} new names for {len(old_names)} fitness columns"
        )
    pairs = list(zip(old_names, new_names_all))
    experimental_pairs = [
        (old, new) for old, new in pairs if new.split("-", 1)[0] != "T0"
    ]
    selected_old = [old for old, _ in experimental_pairs]
    rename_map = dict(experimental_pairs)
    raw = raw.loc[:, ID_COLUMNS + selected_old].rename(columns=rename_map)
    t_like = t_like.loc[:, ID_COLUMNS + selected_old].rename(columns=rename_map)
    filtered = None
    if args.filtered:
        filtered = pd.read_csv(args.filtered)
        if not {"locusId", "desc"}.issubset(filtered.columns):
            raise ValueError("Filtered table must contain locusId and desc columns")
        if "sysName" not in filtered:
            filtered = filtered.merge(raw[["locusId", "sysName"]], on="locusId", how="left")
    return raw, t_like, filtered, groups, t_like_path


def classify_condition_confidence(fitness, t_like):
    """Match the confidence rules used by fit_analysis.py."""
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
    confidence = pd.Series(
        np.select(
            [high, moderate, exploratory],
            ["high", "moderate", "low"],
            default="not_supported",
        ),
        index=fitness.index,
    )
    return confidence, evaluable, direction_consistent


def summarize_condition_evidence(raw, t_like, groups):
    """Summarize original fitness/t-like evidence for every gene-condition."""
    frames = []
    for condition, configured_samples in groups.items():
        samples = [sample for sample in configured_samples
                   if sample in raw.columns and sample in t_like.columns]
        if not samples:
            continue
        fitness = raw[samples].apply(pd.to_numeric, errors="coerce")
        t_values = t_like[samples].apply(pd.to_numeric, errors="coerce")
        confidence, evaluable, direction_consistent = classify_condition_confidence(
            fitness, t_values
        )
        frame = raw[ID_COLUMNS].copy()
        frame["condition"] = condition
        frame["samples"] = ";".join(samples)
        frame["replicates_configured"] = len(configured_samples)
        frame["replicates_evaluable"] = evaluable.sum(axis=1)
        frame["mean_fitness"] = fitness.where(evaluable).mean(axis=1)
        frame["fitness_sd"] = fitness.where(evaluable).std(axis=1, ddof=1)
        frame["mean_t_like"] = t_values.where(evaluable).mean(axis=1)
        frame["max_abs_t_like"] = t_values.where(evaluable).abs().max(axis=1)
        frame["direction_consistent"] = direction_consistent
        frame["confidence_tier"] = confidence
        frames.append(frame)
    if not frames:
        raise ValueError("No group samples were found in both original FEBA tables")
    return pd.concat(frames, ignore_index=True)


def add_original_evidence(results, evidence):
    """Attach within-condition FEBA evidence to each pairwise comparison."""
    annotated = results.copy()
    evidence_columns = [
        column for column in evidence.columns
        if column not in {"sysName", "desc", "condition"}
    ]
    for side in (1, 2):
        side_evidence = evidence[["condition", *evidence_columns]].copy()
        rename = {
            column: f"condition_{side}_{column}"
            for column in evidence_columns if column != "locusId"
        }
        side_evidence = side_evidence.rename(
            columns={"condition": f"condition_{side}", **rename}
        )
        annotated = annotated.merge(
            side_evidence,
            on=["locusId", f"condition_{side}"],
            how="left",
            validate="many_to_one",
        )
    annotated["phenotype_support_summary"] = (
        annotated["condition_1_confidence_tier"].fillna("not_available")
        + " vs "
        + annotated["condition_2_confidence_tier"].fillna("not_available")
    )
    condition_1_lower = annotated["fitness_difference"].lt(0)
    condition_2_lower = annotated["fitness_difference"].gt(0)
    annotated["lower_fitness_condition"] = np.select(
        [condition_1_lower, condition_2_lower],
        [annotated["condition_1"], annotated["condition_2"]],
        default="equal_or_undefined",
    )
    annotated["lower_fitness_confidence_tier"] = np.select(
        [condition_1_lower, condition_2_lower],
        [annotated["condition_1_confidence_tier"],
         annotated["condition_2_confidence_tier"]],
        default="not_available",
    )
    return annotated


def run_pairwise_tests(data, groups):
    """Return all genes per valid comparison, retaining untestable rows."""
    valid_groups = {
        condition: [sample for sample in samples if sample in data.columns]
        for condition, samples in groups.items()
    }
    for condition, samples in groups.items():
        missing = [sample for sample in samples if sample not in data.columns]
        if missing:
            LOGGER.warning("%s: samples absent from table: %s", condition, ", ".join(missing))
    valid_groups = {name: samples for name, samples in valid_groups.items()
                    if len(samples) >= 2}
    if len(valid_groups) < 2:
        raise ValueError("At least two groups with two recognized replicates are required")

    rows = []
    for condition_1, condition_2 in combinations(valid_groups, 2):
        samples_1, samples_2 = valid_groups[condition_1], valid_groups[condition_2]
        for _, gene in data.iterrows():
            values_1 = pd.to_numeric(gene[samples_1], errors="coerce").dropna()
            values_2 = pd.to_numeric(gene[samples_2], errors="coerce").dropna()
            t_statistic = p_value = np.nan
            if len(values_1) >= 2 and len(values_2) >= 2:
                test = ttest_ind(values_1, values_2, equal_var=False, nan_policy="omit")
                t_statistic, p_value = test.statistic, test.pvalue
            rows.append({
                "locusId": gene["locusId"], "sysName": gene.get("sysName", np.nan),
                "desc": gene.get("desc", np.nan), "condition_1": condition_1,
                "condition_2": condition_2, "replicates_1": len(values_1),
                "replicates_2": len(values_2),
                "mean_1": values_1.mean() if len(values_1) else np.nan,
                "mean_2": values_2.mean() if len(values_2) else np.nan,
                "fitness_difference": (values_1.mean() - values_2.mean()
                                       if len(values_1) and len(values_2) else np.nan),
                "t_statistic": t_statistic, "p_value": p_value, "q_value": np.nan,
                "fdr_test_count": np.nan,
            })

    results = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    for indices in results.groupby(["condition_1", "condition_2"]).groups.values():
        valid_indices = results.loc[list(indices), "p_value"].dropna().index
        if len(valid_indices):
            results.loc[valid_indices, "q_value"] = multipletests(
                results.loc[valid_indices, "p_value"], method="fdr_bh"
            )[1]
            results.loc[valid_indices, "fdr_test_count"] = len(valid_indices)
    return results


def save_results(results, prefix, output_path):
    all_file = output_path / f"{prefix}_ttest_results.csv"
    significant_file = output_path / f"{prefix}_ttest_q_le_0p10.csv"
    results.to_csv(all_file, index=False)
    significant = results[results["q_value"].le(0.10)].copy()
    significant["abs_fitness_difference"] = significant["fitness_difference"].abs()
    significant = significant.sort_values(
        ["q_value", "abs_fitness_difference"],
        ascending=[True, False],
        kind="stable",
    ).drop(columns="abs_fitness_difference")
    significant.to_csv(significant_file, index=False)
    LOGGER.info("%s: %s rows tested, %s with q <= 0.10", prefix,
                f"{results['q_value'].notna().sum():,}",
                f"{results['q_value'].le(0.10).sum():,}")
    LOGGER.info("Wrote: %s", all_file)
    LOGGER.info("Wrote: %s", significant_file)
    LOGGER.info("%s lower-fitness confidence tiers: %s", prefix,
                results["lower_fitness_confidence_tier"].value_counts().to_dict())
    return results


def safe_filename(value):
    """Convert a comparison label into a portable filename component."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_")


def save_figure(path):
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    LOGGER.info("Wrote figure: %s", path)


def make_comparison_plots(results, prefix, plot_path):
    """Create consolidated volcano, p-value, and mean-difference figures."""
    exploratory = " — candidate-focused exploratory analysis" if prefix == "filtered" else ""
    panels = []
    for (condition_1, condition_2), group in results.groupby(
        ["condition_1", "condition_2"], sort=False
    ):
        tested = group.dropna(subset=["p_value", "q_value", "fitness_difference"]).copy()
        if tested.empty:
            LOGGER.warning("Skipping plots for %s vs %s: no valid tests",
                           condition_1, condition_2)
            continue
        comparison = f"{condition_1}_vs_{condition_2}"
        filename = safe_filename(comparison) + ".png"

        positive_q = tested.loc[tested["q_value"].gt(0), "q_value"]
        floor = positive_q.min() / 10 if len(positive_q) else 1e-300
        tested["minus_log10_q"] = -np.log10(tested["q_value"].clip(lower=floor))
        tested["significance"] = np.select(
            [tested["q_value"].lt(0.05), tested["q_value"].lt(0.10)],
            ["q < 0.05", "0.05 ≤ q < 0.10"], default="q ≥ 0.10",
        )
        tested["mean_fitness"] = tested[["mean_1", "mean_2"]].mean(axis=1)
        panels.append((condition_1, condition_2, tested))

    if panels:
        palette = {"q ≥ 0.10": "0.65", "0.05 ≤ q < 0.10": "tab:orange",
                   "q < 0.05": "tab:red"}
        hue_order = list(palette)
        n_panels = len(panels)
        n_cols = min(3, n_panels)
        n_rows = math.ceil(n_panels / n_cols)

        # All volcano comparisons in one figure.
        figure, axes = plt.subplots(
            n_rows, n_cols, figsize=(5.2 * n_cols + 1.5, 4.2 * n_rows),
            sharex=False, sharey=False,
        )
        axes = np.atleast_1d(axes).ravel()
        legend_handles = legend_labels = None
        for axis, (condition_1, condition_2, tested) in zip(axes, panels):
            sns.scatterplot(
                data=tested, x="fitness_difference", y="minus_log10_q",
                hue="significance", hue_order=hue_order, palette=palette,
                s=14, alpha=0.65, linewidth=0, ax=axis,
            )
            axis.axhline(-np.log10(0.05), color="black", linestyle="--",
                         linewidth=0.8)
            axis.axvline(0, color="black", linewidth=0.6)
            # Limit labels to the strongest q-values and effect-size extremes.
            label_indices = set(tested.nsmallest(3, "q_value").index)
            label_indices.update(tested.nlargest(2, "fitness_difference").index)
            label_indices.update(tested.nsmallest(2, "fitness_difference").index)
            labeled = tested.loc[list(label_indices)].sort_values(
                "minus_log10_q", ascending=False
            )
            side_number = {"left": 0, "right": 0}
            for _, row in labeled.iterrows():
                side = "right" if row["fitness_difference"] >= 0 else "left"
                number = side_number[side]
                side_number[side] += 1
                vertical_offset = (6 + 5 * (number // 2)) * (1 if number % 2 == 0 else -1)
                horizontal_offset = 5 if side == "right" else -5
                axis.annotate(
                    str(row["locusId"]),
                    (row["fitness_difference"], row["minus_log10_q"]),
                    xytext=(horizontal_offset, vertical_offset),
                    textcoords="offset points", fontsize=4.5,
                    ha="left" if side == "right" else "right", va="center",
                    arrowprops={"arrowstyle": "-", "color": "0.45", "lw": 0.35},
                    bbox={"boxstyle": "round,pad=0.08", "fc": "white",
                          "ec": "none", "alpha": 0.72},
                )
            axis.set_title(f"{condition_1} vs {condition_2}", fontsize=9)
            axis.set_xlabel(f"Fitness difference\n({condition_1} − {condition_2})")
            axis.set_ylabel("−log10(q-value)")
            if legend_handles is None:
                legend_handles, legend_labels = axis.get_legend_handles_labels()
            if axis.get_legend() is not None:
                axis.get_legend().remove()
        for axis in axes[n_panels:]:
            axis.set_visible(False)
        figure.suptitle(f"Volcano plots for all {prefix} comparisons{exploratory}")
        figure.legend(legend_handles, legend_labels, title="Significance",
                      bbox_to_anchor=(0.995, 0.98), loc="upper right",
                      fontsize=7, title_fontsize=8, frameon=False)
        figure.tight_layout(rect=(0, 0, 0.88, 0.95))
        save_figure(plot_path / f"{prefix}_volcano_plots.png")

        # All mean-fitness versus difference comparisons in one figure.
        figure, axes = plt.subplots(
            n_rows, n_cols, figsize=(5.2 * n_cols + 1.5, 4.2 * n_rows),
            sharex=False, sharey=False,
        )
        axes = np.atleast_1d(axes).ravel()
        legend_handles = legend_labels = None
        for axis, (condition_1, condition_2, tested) in zip(axes, panels):
            sns.scatterplot(
                data=tested, x="mean_fitness", y="fitness_difference",
                hue="significance", hue_order=hue_order, palette=palette,
                s=14, alpha=0.65, linewidth=0, ax=axis,
            )
            axis.axhline(0, color="black", linewidth=0.6)
            axis.set_title(f"{condition_1} vs {condition_2}", fontsize=9)
            axis.set_xlabel("Mean fitness across both conditions")
            axis.set_ylabel(f"Fitness difference\n({condition_1} − {condition_2})")
            if legend_handles is None:
                legend_handles, legend_labels = axis.get_legend_handles_labels()
            if axis.get_legend() is not None:
                axis.get_legend().remove()
        for axis in axes[n_panels:]:
            axis.set_visible(False)
        figure.suptitle(
            f"Mean fitness versus difference for all {prefix} comparisons{exploratory}"
        )
        figure.legend(legend_handles, legend_labels, title="Significance",
                      bbox_to_anchor=(0.995, 0.98), loc="upper right",
                      fontsize=7, title_fontsize=8, frameon=False)
        figure.tight_layout(rect=(0, 0, 0.88, 0.95))
        save_figure(plot_path / f"{prefix}_mean_difference_plots.png")

        # All raw p-value distributions in one figure.
        figure, axes = plt.subplots(
            n_rows, n_cols, figsize=(5 * n_cols, 3.8 * n_rows),
            sharex=True, sharey=False,
        )
        axes = np.atleast_1d(axes).ravel()
        bins = np.linspace(0, 1, 21)
        for axis, (condition_1, condition_2, tested) in zip(axes, panels):
            axis.hist(tested["p_value"], bins=bins, color="steelblue",
                      edgecolor="white")
            axis.axvline(0.05, color="tab:red", linestyle="--", linewidth=0.8)
            axis.set_title(f"{condition_1} vs {condition_2}", fontsize=9)
            axis.set_xlabel("Raw p-value")
            axis.set_ylabel("Genes")
        for axis in axes[n_panels:]:
            axis.set_visible(False)
        figure.suptitle(f"P-value distributions for all {prefix} comparisons")
        figure.tight_layout(rect=(0, 0, 1, 0.96))
        save_figure(plot_path / f"{prefix}_pvalue_histograms.png")


def make_difference_heatmap(results, prefix, plot_path):
    """Heatmap genes with q < 0.10 in any comparison, or the best 30 if none."""
    tested = results.dropna(subset=["q_value", "fitness_difference"]).copy()
    if tested.empty:
        LOGGER.warning("Skipping %s difference heatmap: no valid tests", prefix)
        return
    selected_genes = tested.loc[tested["q_value"].lt(0.10), "locusId"].unique()
    selection_note = "q < 0.10 in at least one comparison"
    if not len(selected_genes):
        selected_genes = (
            tested.groupby("locusId")["q_value"].min().nsmallest(30).index.to_numpy()
        )
        selection_note = "30 genes with smallest minimum q-values"
    selected = tested[tested["locusId"].isin(selected_genes)].copy()
    selected["comparison"] = selected["condition_1"] + " vs " + selected["condition_2"]
    descriptions = selected.groupby("locusId")["desc"].first().fillna("")
    matrix = selected.pivot_table(index="locusId", columns="comparison",
                                  values="fitness_difference", aggfunc="first")
    mean_1 = selected.pivot_table(index="locusId", columns="comparison",
                                  values="mean_1", aggfunc="first")
    mean_2 = selected.pivot_table(index="locusId", columns="comparison",
                                  values="mean_2", aggfunc="first")
    row_order = matrix.abs().max(axis=1).sort_values(ascending=False).index
    matrix = matrix.loc[row_order]
    mean_1 = mean_1.reindex(index=matrix.index, columns=matrix.columns)
    mean_2 = mean_2.reindex(index=matrix.index, columns=matrix.columns)
    annotations = pd.DataFrame("", index=matrix.index, columns=matrix.columns)
    for row in matrix.index:
        for column in matrix.columns:
            if pd.notna(mean_1.loc[row, column]) and pd.notna(mean_2.loc[row, column]):
                annotations.loc[row, column] = (
                    f"C1 {mean_1.loc[row, column]:.2f}\n"
                    f"C2 {mean_2.loc[row, column]:.2f}"
                )
    display_labels = [f"{gene}|{descriptions.get(gene, '')}" for gene in matrix.index]
    matrix.index = display_labels
    annotations.index = display_labels
    color_limit = matrix.abs().max().max()
    if not np.isfinite(color_limit) or color_limit == 0:
        color_limit = 1.0
    figure_width = max(12, 0.32 * len(matrix.columns))
    figure_height = max(7, min(30, 0.30 * len(matrix)))
    figure, axis = plt.subplots(figsize=(figure_width, figure_height))
    figure.subplots_adjust(top=0.88, bottom=0.25, left=0.25, right=0.80)
    colorbar_axis = figure.add_axes([0.86, 0.28, 0.016, 0.46])
    sns.heatmap(matrix, ax=axis, cmap="RdBu_r", center=0,
                vmin=-color_limit, vmax=color_limit, mask=matrix.isna(),
                annot=annotations, fmt="", annot_kws={"fontsize": 5},
                cbar_ax=colorbar_axis, cbar_kws={"orientation": "vertical"})
    colorbar_axis.set_ylabel("Fitness difference", fontsize=7, labelpad=5)
    colorbar_axis.tick_params(labelsize=6, length=2)
    exploratory = " — exploratory" if prefix == "filtered" else ""
    axis.set_title(
        f"Pairwise fitness differences ({selection_note}){exploratory}\n"
        "Cell text: mean fitness in condition 1 (C1) and condition 2 (C2)"
    )
    axis.set_xlabel("Condition comparison")
    axis.set_ylabel("Gene")
    axis.tick_params(axis="x", rotation=45, labelsize=7)
    axis.tick_params(axis="y", labelsize=6)
    for label in axis.get_xticklabels():
        label.set_horizontalalignment("right")
    save_figure(plot_path / f"{prefix}_difference_heatmap.png")


def make_comparison_summary(results, csv_path, plot_path):
    """Summarize all-gene test and discovery counts for every comparison."""
    frames = []
    for (condition_1, condition_2), group in results.groupby(
        ["condition_1", "condition_2"], sort=False
    ):
        q_significant = group["q_value"].lt(0.05)
        frames.append({
            "condition_1": condition_1, "condition_2": condition_2,
            "genes_in_input": len(group),
            "fdr_test_count": int(group["q_value"].notna().sum()),
            "raw_p_lt_0p05": int(group["p_value"].lt(0.05).sum()),
            "q_lt_0p10": int(group["q_value"].lt(0.10).sum()),
            "q_lt_0p05": int(q_significant.sum()),
            "q_lt_0p05_positive_difference": int(
                (q_significant & group["fitness_difference"].gt(0)).sum()
            ),
            "q_lt_0p05_negative_difference": int(
                (q_significant & group["fitness_difference"].lt(0)).sum()
            ),
        })
    summary = pd.DataFrame(frames)
    output_file = csv_path / "comparison_summary.csv"
    summary.to_csv(output_file, index=False)
    LOGGER.info("Wrote comparison summary: %s", output_file)
    summary["comparison"] = summary["condition_1"] + " vs " + summary["condition_2"]
    figure_height = max(7, 0.28 * summary["comparison"].nunique())
    plt.figure(figsize=(9, figure_height))
    sns.barplot(data=summary, y="comparison", x="q_lt_0p05", color="tab:blue")
    plt.title("Pairwise discoveries after all-gene FDR correction")
    plt.xlabel("Genes with q < 0.05")
    plt.ylabel("Condition comparison")
    plt.tight_layout()
    save_figure(plot_path / "comparison_summary.png")


def write_results_guide(output_path, includes_filtered):
    """Document the all-gene FDR scope and the primary output files."""
    lines = [
        "# Pairwise fitness-test results", "",
        "Each condition pair is tested gene by gene with Welch's t-test. Benjamini–Hochberg q-values are calculated separately for each condition pair across **all valid gene tests** in that pair. Neither the fit-analysis t-like filter nor an effect-size filter is applied before q-value calculation.",
        "", "## Start here", "",
        "- `csv/all_genes_ttest_q_le_0p10.csv`: genes with q ≤ 0.10, ordered by q-value then absolute fitness difference.",
        "- `csv/comparison_summary.csv`: the number of FDR tests and discoveries for every condition pair.",
        "- `csv/all_genes_ttest_results.csv`: all tests, including non-significant and untestable genes. `fdr_test_count` shows the correction scope for that comparison.",
        "- `plots/all_genes_volcano_plots.png`: effect size against q-value for each comparison.",
        "- `plots/all_genes_difference_heatmap.png`: fitness differences for genes with q < 0.10 in at least one comparison.",
        "", "## Interpretation", "",
        "A q-value answers whether fitness differs between the two conditions after accounting for all genes tested in that comparison. The `condition_1_*` and `condition_2_*` columns add FEBA within-condition evidence and confidence tiers; they do not affect the pairwise q-value.",
        "", "## Plots and their underlying tables", "",
        "| Plot | Main source table | Purpose |",
        "|---|---|---|",
        "| `plots/all_genes_volcano_plots.png` | `csv/all_genes_ttest_results.csv` | Fitness difference versus q-value for each condition comparison. |",
        "| `plots/all_genes_mean_difference_plots.png` | `csv/all_genes_ttest_results.csv` | Mean fitness versus the difference between conditions. |",
        "| `plots/all_genes_pvalue_histograms.png` | `csv/all_genes_ttest_results.csv` | Raw p-value distributions; a test-quality diagnostic. |",
        "| `plots/all_genes_difference_heatmap.png` | `csv/all_genes_ttest_results.csv` | Fitness differences for genes with q < 0.10 in at least one comparison. |",
        "| `plots/comparison_summary.png` | `csv/comparison_summary.csv` | Number of q < 0.05 discoveries per comparison after all-gene FDR correction. |",
    ]
    if includes_filtered:
        lines.extend([
            "", "## Exploratory filtered results", "",
            "`filtered_exploratory_*` files repeat the comparisons only for genes retained by the fit-analysis t-like filter. Their q-values are conditional on that preselection, so use them for curiosity or candidate review only; do not treat them as the formal FDR-controlled result.",
            "", "| Exploratory plot | Main source table | Purpose |",
            "|---|---|---|",
            "| `plots/filtered_exploratory_volcano_plots.png` | `csv/filtered_exploratory_ttest_results.csv` | Candidate-only fitness difference versus conditional q-value. |",
            "| `plots/filtered_exploratory_mean_difference_plots.png` | `csv/filtered_exploratory_ttest_results.csv` | Candidate-only mean fitness versus difference. |",
            "| `plots/filtered_exploratory_pvalue_histograms.png` | `csv/filtered_exploratory_ttest_results.csv` | Candidate-only raw p-value diagnostic. |",
            "| `plots/filtered_exploratory_difference_heatmap.png` | `csv/filtered_exploratory_ttest_results.csv` | Candidate-only fitness differences for q < 0.10 in at least one comparison. |",
        ])
    (output_path / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    LOGGER.info("Wrote results guide: %s", output_path / "README.md")


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    csv_path = args.output / "csv"
    plot_path = args.output / "plots"
    csv_path.mkdir(parents=True, exist_ok=True)
    plot_path.mkdir(parents=True, exist_ok=True)
    setup_logging(args.output)
    LOGGER.info("Starting all-gene pairwise Welch tests")
    raw, t_like, filtered, groups, t_like_path = load_inputs(args)
    LOGGER.info("Loaded %s genes across %s groups; Time0 columns excluded",
                f"{len(raw):,}", f"{len(groups):,}")
    LOGGER.info("Original t-like table: %s", t_like_path)
    evidence = summarize_condition_evidence(raw, t_like, groups)
    evidence_file = csv_path / "original_gene_condition_confidence.csv"
    evidence.to_csv(evidence_file, index=False)
    LOGGER.info("Wrote original gene-condition confidence evidence: %s", evidence_file)
    results = save_results(
        add_original_evidence(run_pairwise_tests(raw, groups), evidence),
        "all_genes", csv_path
    )
    make_comparison_plots(results, "all_genes", plot_path)
    make_difference_heatmap(results, "all_genes", plot_path)
    make_comparison_summary(results, csv_path, plot_path)
    if filtered is not None:
        exploratory = save_results(
            add_original_evidence(run_pairwise_tests(filtered, groups), evidence),
            "filtered_exploratory", csv_path,
        )
        make_comparison_plots(exploratory, "filtered_exploratory", plot_path)
        make_difference_heatmap(exploratory, "filtered_exploratory", plot_path)
    write_results_guide(args.output, filtered is not None)
    LOGGER.info("Pairwise analysis completed successfully")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        LOGGER.exception("Pairwise analysis failed")
        raise
