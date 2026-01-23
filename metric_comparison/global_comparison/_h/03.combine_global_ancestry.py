#!/usr/bin/env python3
"""
Combine per-chromosome phased global ancestry results and generate summary plots.

This script collects per-chromosome TSV files from array job outputs,
concatenates them, computes overall metrics, and generates visualization plots.
"""
import logging
import argparse
import numpy as np
import pandas as pd
import session_info
import seaborn as sns
from scipy import stats
from pathlib import Path

COMPARISON_LABELS = {
    "simu_vs_rfmix": "RFMix Comparison",
    "simu_vs_flare": "FLARE Comparison",
    "rfmix_vs_flare": "LAI Comparison",
}


def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine per-chromosome phased global ancestry results."
    )
    parser.add_argument(
        "--metrics-root", type=Path, default=Path("../_m"),
        help="Root directory containing phased/unphased results.",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("combined"),
        help="Output subdirectory name for combined results.",
    )
    parser.add_argument(
        "--population", type=str, choices=["two", "three"],
        default="three",
        help="Population type to process."
    )
    return parser.parse_args()


def r2_score(y_true, y_pred):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot == 0:
        return np.nan
    return 1 - ss_res / ss_tot


def choose_corr(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    if np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return np.nan, np.nan

    nx, ny = len(x), len(y)
    if nx < 8 or ny < 8:
        corr, p_val = stats.spearmanr(x, y)
        return corr, p_val

    x_norm = stats.normaltest(x)
    y_norm = stats.normaltest(y)
    if x_norm.pvalue > 0.05 and y_norm.pvalue > 0.05:
        corr, p_val = stats.pearsonr(x, y)
        return corr, p_val

    corr, p_val = stats.spearmanr(x, y)
    return corr, p_val


def compute_comparison_metrics(df_a, df_b, label_a, label_b, comparison_label):
    merged = df_a.merge(
        df_b, on=["sample_id", "chrom", "ancestry"],
        how="inner", suffixes=(f"_{label_a}", f"_{label_b}"),
    )
    results = []
    grouped = merged.groupby(["chrom", "ancestry"], dropna=False)
    for (chrom, ancestry), group in grouped:
        x = group[f"g_anc_{label_a}"]
        y = group[f"g_anc_{label_b}"]
        r2 = r2_score(y, x)
        bias = float((x - y).mean())
        corr, p_val = choose_corr(x, y)
        results.append(
            {
                "chrom": chrom,
                "ancestry": ancestry,
                "r2": r2,
                "bias": bias,
                "corr": corr,
                "p_val": p_val,
                "lai": label_a,
                "comparison": comparison_label,
            }
        )
    return pd.DataFrame(results)


def add_identity_line(ax, data_min, data_max):
    ax.plot([data_min, data_max], [data_min, data_max],
            linestyle="--", color="gray", linewidth=1)


def plot_scatter_per_ancestry(df, output_dir):
    sns.set_theme(style="whitegrid", context="talk")
    comparisons = [
        ("simu", "rfmix", "Simulated vs RFMix"),
        ("simu", "flare", "Simulated vs FLARE"),
        ("flare", "rfmix", "FLARE vs RFMix"),
    ]
    method_labels = {
        "simu": "Simulated",
        "rfmix": "RFMix",
        "flare": "FLARE"
    }
    frames = []
    for left, right, label in comparisons:
        left_df = df[df["method"] == left]
        right_df = df[df["method"] == right]
        merged = left_df.merge(
            right_df, on=["sample_id", "chrom", "ancestry"],
            how="inner", suffixes=("_x", "_y"),
        )
        merged = merged.rename(columns={"g_anc_x": "x", "g_anc_y": "y"})
        merged.loc[:, "comparison"] = label
        merged["x_label"] = method_labels[left]
        merged["y_label"] = method_labels[right]
        frames.append(merged)
    plot_df = pd.concat(frames, ignore_index=True)

    g = sns.FacetGrid(
        plot_df, row="ancestry", col="comparison", hue="ancestry", height=3.4,
        aspect=1.1, sharex=False, sharey=False, margin_titles=True,
    )
    g.map_dataframe(sns.scatterplot, x="x", y="y", alpha=0.6, s=18)
    g.map_dataframe(sns.regplot, x="x", y="y", scatter=False,
                    ci=95, color="black")
    for (facet_keys, data), ax in zip(g.facet_data(), g.axes.flat):
        data_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
        data_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
        add_identity_line(ax, data_min, data_max)
        x_label = data["x_label"].iloc[0]
        y_label = data["y_label"].iloc[0]
        ax.set_xlabel(x_label, fontweight="bold")
        ax.set_ylabel(y_label, fontweight="bold")

    png_path = output_dir / "global_ancestry_scatter_by_ancestry.png"
    pdf_path = output_dir / "global_ancestry_scatter_by_ancestry.pdf"
    g.savefig(png_path, dpi=300, bbox_inches="tight")
    g.savefig(pdf_path, bbox_inches="tight")


def plot_scatter_by_comparison(df, output_dir):
    sns.set_theme(style="whitegrid", context="talk")
    comparisons = [
        ("simu", "rfmix", "Simulated vs RFMix"),
        ("simu", "flare", "Simulated vs FLARE"),
        ("flare", "rfmix", "FLARE vs RFMix"),
    ]
    frames = []
    for left, right, label in comparisons:
        left_df = df[df["method"] == left]
        right_df = df[df["method"] == right]
        merged = left_df.merge(
            right_df, on=["sample_id", "chrom", "ancestry"],
            how="inner", suffixes=("_x", "_y"),
        )
        merged = merged.rename(columns={"g_anc_x": "x", "g_anc_y": "y"})
        merged.loc[:, "comparison"] = label
        frames.append(merged)
    plot_df = pd.concat(frames, ignore_index=True)

    g = sns.FacetGrid(
        plot_df, col="comparison", hue="ancestry", height=4.0,
        aspect=1.1, sharex=False, sharey=False,
    )
    g.map_dataframe(sns.scatterplot, x="x", y="y", alpha=0.6, s=20)
    g.map_dataframe(sns.regplot, x="x", y="y", scatter=False,
                    ci=95, color="black")
    for ax in g.axes.flat:
        data_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
        data_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
        add_identity_line(ax, data_min, data_max)
        ax.set_xlabel("", fontweight="bold")
        ax.set_ylabel("", fontweight="bold")
    g.set_axis_labels("Global Ancestry (x)", "Global Ancestry (y)")
    g.add_legend(title="Ancestry", fontsize=10, title_fontsize=11)

    png_path = output_dir / "global_ancestry_scatter_by_comparison.png"
    pdf_path = output_dir / "global_ancestry_scatter_by_comparison.pdf"
    g.savefig(png_path, dpi=300, bbox_inches="tight")
    g.savefig(pdf_path, bbox_inches="tight")


def plot_boxplots(summary_df, output_dir):
    sns.set_theme(style="whitegrid", context="talk")
    plot_df = summary_df.copy()
    plot_df = plot_df.rename(
        columns={"comparison": "Comparison", "ancestry": "Ancestry"}
    )
    plot_df = plot_df.melt(
        id_vars=["Comparison", "Ancestry"], value_vars=["corr", "bias"],
        var_name="Metric", value_name="Value",
    )

    g = sns.catplot(
        data=plot_df, x="Comparison", y="Value", hue="Ancestry",
        col="Metric", kind="box", sharey=False, height=4.0,
        aspect=1.1,
    )

    for ax in g.axes.flat:
        sns.stripplot(
            data=plot_df[plot_df["Metric"] == ax.get_title().split(" = ")[-1]],
            x="Comparison", y="Value", hue="Ancestry", dodge=True,
            alpha=0.6, size=4, ax=ax,
        )
        ax.set_xlabel("", fontweight="bold")
        ax.set_ylabel("", fontweight="bold")
        ax.legend_.remove()
        for label in ax.get_xticklabels():
            label.set_rotation(45)
            label.set_horizontalalignment('right')

    handles, labels = g.axes.flat[0].get_legend_handles_labels()
    g.fig.legend(handles, labels, title="Ancestry", loc="upper right")

    png_path = output_dir / "global_ancestry_metrics_boxplot.png"
    pdf_path = output_dir / "global_ancestry_metrics_boxplot.pdf"
    g.savefig(png_path, dpi=300, bbox_inches="tight")
    g.savefig(pdf_path, bbox_inches="tight")


def collect_chromosome_files(pop_dir: Path):
    """Collect per-chromosome TSV files from chr* subdirectories."""
    tidy_frames = []
    metrics_frames = []

    for chrom in range(1, 23):
        chrom_dir = pop_dir / f"chr{chrom}"
        if not chrom_dir.exists():
            logging.warning("Missing chromosome directory: %s", chrom_dir)
            continue

        tidy_path = chrom_dir / "global_ancestry_per_chrom.tsv"
        metrics_path = chrom_dir / "global_ancestry_metrics_per_chrom.tsv"

        if tidy_path.exists():
            tidy_df = pd.read_csv(tidy_path, sep="\t")
            tidy_frames.append(tidy_df)
            logging.info("Loaded tidy data for chr%d: %d rows", chrom, len(tidy_df))
        else:
            logging.warning("Missing tidy file: %s", tidy_path)

        if metrics_path.exists():
            metrics_df = pd.read_csv(metrics_path, sep="\t")
            metrics_frames.append(metrics_df)
            logging.info("Loaded metrics for chr%d: %d rows", chrom, len(metrics_df))
        else:
            logging.warning("Missing metrics file: %s", metrics_path)

    if not tidy_frames:
        return None, None

    combined_tidy = pd.concat(tidy_frames, ignore_index=True)
    combined_metrics = pd.concat(metrics_frames, ignore_index=True) if metrics_frames else None

    return combined_tidy, combined_metrics


def main():
    configure_logging()
    args = parse_args()

    metrics_root = args.metrics_root.resolve()

    # Process both phased and unphased if they exist
    for phase in ["phased", "unphased"]:
        phase_dir = metrics_root / phase
        if not phase_dir.exists():
            logging.info("Skipping %s (directory not found)", phase)
            continue

        pop_dir = phase_dir / f"{args.population}_populations"
        if not pop_dir.exists():
            logging.info("Skipping %s/%s_populations (not found)", phase, args.population)
            continue

        output_dir = pop_dir / args.output
        output_dir.mkdir(parents=True, exist_ok=True)

        logging.info("Processing %s/%s_populations", phase, args.population)

        # Check if we need to combine from per-chromosome files
        combined_tidy_path = pop_dir / "global_ancestry_per_chrom.tsv"

        if combined_tidy_path.exists() and (pop_dir / "chr1").exists():
            logging.info("Found both combined file and per-chrom directories. Using per-chrom files.")
            tidy, metrics = collect_chromosome_files(pop_dir)
        elif (pop_dir / "chr1").exists():
            logging.info("Collecting per-chromosome files...")
            tidy, metrics = collect_chromosome_files(pop_dir)
        elif combined_tidy_path.exists():
            logging.info("Loading existing combined file...")
            tidy = pd.read_csv(combined_tidy_path, sep="\t")
            metrics_path = pop_dir / "global_ancestry_metrics_per_chrom.tsv"
            metrics = pd.read_csv(metrics_path, sep="\t") if metrics_path.exists() else None
        else:
            logging.warning("No data found for %s/%s_populations", phase, args.population)
            continue

        if tidy is None or tidy.empty:
            logging.warning("No tidy data found for %s/%s_populations", phase, args.population)
            continue

        # Save combined tidy data
        combined_tidy_out = output_dir / "global_ancestry_per_chrom.tsv"
        tidy.to_csv(combined_tidy_out, sep="\t", index=False)
        logging.info("Wrote combined tidy data to %s (%d rows)", combined_tidy_out, len(tidy))

        # Recompute metrics from combined data
        simu_df = tidy[tidy["method"] == "simu"].rename(columns={"g_anc": "g_anc_simu"})
        rfmix_df = tidy[tidy["method"] == "rfmix"].rename(columns={"g_anc": "g_anc_rfmix"})
        flare_df = tidy[tidy["method"] == "flare"].rename(columns={"g_anc": "g_anc_flare"})

        simu_vs_rfmix = compute_comparison_metrics(
            rfmix_df, simu_df, "rfmix", "simu",
            COMPARISON_LABELS["simu_vs_rfmix"],
        )
        simu_vs_flare = compute_comparison_metrics(
            flare_df, simu_df, "flare", "simu",
            COMPARISON_LABELS["simu_vs_flare"],
        )
        rfmix_vs_flare = compute_comparison_metrics(
            rfmix_df, flare_df, "rfmix", "flare",
            COMPARISON_LABELS["rfmix_vs_flare"],
        )

        summary_df = pd.concat(
            [simu_vs_rfmix, simu_vs_flare, rfmix_vs_flare],
            ignore_index=True,
        )
        summary_path = output_dir / "global_ancestry_metrics_per_chrom.tsv"
        summary_df.to_csv(summary_path, sep="\t", index=False)
        logging.info("Wrote metrics to %s", summary_path)

        # Generate plots
        logging.info("Generating plots for %s/%s_populations", phase, args.population)
        plot_scatter_per_ancestry(tidy, output_dir)
        plot_scatter_by_comparison(tidy, output_dir)
        plot_boxplots(summary_df, output_dir)
        logging.info("Saved plots to %s", output_dir)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
