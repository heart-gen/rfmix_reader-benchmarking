import logging
import argparse
import numpy as np
import pandas as pd
import session_info
import seaborn as sns
from scipy import stats
from pyhere import here
from pathlib import Path
from rfmix_reader import read_simu, read_rfmix, read_flare

COMPARISON_LABELS = {
    "simu_vs_rfmix": "Simulation Comparison",
    "simu_vs_flare": "Simulation Comparison",
    "rfmix_vs_flare": "LAI Comparison",
}

def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Per-chromosome global ancestry comparisons for RFMix, simulation, and FLARE."
    )
    parser.add_argument("--simu-input", type=Path, required=True)
    parser.add_argument("--rfmix-input", type=Path, required=True)
    parser.add_argument("--flare-input", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("./"))
    parser.add_argument("--population", type=str, choices=["two", "three"],
                        default="three")
    return parser.parse_args()


def sort_g_anc(df):
    return df.sort_values(["sample_id", "chrom"]).reset_index(drop=True)


def ancestry_columns(df):
    return [c for c in df.columns if c not in {"sample_id", "chrom"}]


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


def tidy_g_anc(df, method):
    anc_cols = ancestry_columns(df)
    return (
        df.melt(
            id_vars=["sample_id", "chrom"],
            value_vars=anc_cols,
            var_name="ancestry",
            value_name="g_anc",
        )
        .assign(method=method)
        .loc[:, ["sample_id", "chrom", "ancestry", "method", "g_anc"]]
    )


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
        ("sim", "rfmix", "Sim vs RFMix"),
        ("sim", "flare", "Sim vs FLARE"),
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
        plot_df, row="ancestry", col="comparison", hue="chrom",
        height=3.4, aspect=1.1, sharex=False, sharey=False,
        margin_titles=True,
    )
    g.map_dataframe(sns.scatterplot, x="x", y="y", alpha=0.6, s=18)
    g.map_dataframe(sns.regplot, x="x", y="y", scatter=False,
                    ci=95, color="black")
    for ax in g.axes.flat:
        data_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
        data_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
        add_identity_line(ax, data_min, data_max)
        ax.set_xlabel("", fontweight="bold")
        ax.set_ylabel("", fontweight="bold")
    g.set_axis_labels("Global Ancestry (x)", "Global Ancestry (y)")
    g.add_legend(title="Chromosome", fontsize=11, title_fontsize=12)

    png_path = output_dir / "global_ancestry_scatter_by_ancestry.png"
    pdf_path = output_dir / "global_ancestry_scatter_by_ancestry.pdf"
    g.savefig(png_path, dpi=300, bbox_inches="tight")
    g.savefig(pdf_path, bbox_inches="tight")


def plot_scatter_by_comparison(df, output_dir):
    sns.set_theme(style="whitegrid", context="talk")
    comparisons = [
        ("sim", "rfmix", "Sim vs RFMix"),
        ("sim", "flare", "Sim vs FLARE"),
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
    g.map_dataframe(sns.scatterplot, x="x", y="y", alpha=0.6,
                    s=20, style="ancestry")
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
        id_vars=["Comparison", "Ancestry"], value_vars=["r2", "bias"],
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

    handles, labels = g.axes.flat[0].get_legend_handles_labels()
    g.fig.legend(handles, labels, title="Ancestry", loc="upper right")

    png_path = output_dir / "global_ancestry_metrics_boxplot.png"
    pdf_path = output_dir / "global_ancestry_metrics_boxplot.pdf"
    g.savefig(png_path, dpi=300, bbox_inches="tight")
    g.savefig(pdf_path, bbox_inches="tight")


def main():
    configure_logging()
    args = parse_args()

    pop_dir = args.output / f"{args.population}_populations"
    pop_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Loading global ancestry tables.")
    binary_path = args.rfmix_input / "binary_files"
    _, g_anc_rfmix, _ = read_rfmix(here(args.rfmix_input), binary_dir=here(binary_path))
    _, g_anc_flare, _ = read_flare(here(args.flare_input))
    _, g_anc_simu, _ = read_simu(here(args.simu_input))

    g_anc_simu = sort_g_anc(g_anc_simu)
    g_anc_rfmix = sort_g_anc(g_anc_rfmix)
    g_anc_flare = sort_g_anc(g_anc_flare)
    logging.info(
        "Loaded shapes - simu: %s, rfmix: %s, flare: %s",
        g_anc_simu.shape, g_anc_rfmix.shape, g_anc_flare.shape,
    )

    simu_cols = set(g_anc_simu.columns)
    rfmix_cols = set(g_anc_rfmix.columns)
    flare_cols = set(g_anc_flare.columns)

    if simu_cols != rfmix_cols or simu_cols != flare_cols:
        raise ValueError("Input g_anc columns do not match across sources.")

    if g_anc_simu.shape != g_anc_rfmix.shape or g_anc_simu.shape != g_anc_flare.shape:
        raise ValueError("Input g_anc shapes do not match across sources.")

    logging.info("Reshaping global ancestry tables to tidy format.")
    tidy = pd.concat(
        [
            tidy_g_anc(g_anc_simu, "simu"),
            tidy_g_anc(g_anc_rfmix, "rfmix"),
            tidy_g_anc(g_anc_flare, "flare"),
        ],
        ignore_index=True,
    )

    simu_df = tidy[tidy["method"] == "simu"].rename(columns={"g_anc": "g_anc_simu"})
    rfmix_df = tidy[tidy["method"] == "rfmix"].rename(columns={"g_anc": "g_anc_rfmix"})
    flare_df = tidy[tidy["method"] == "flare"].rename(columns={"g_anc": "g_anc_flare"})

    logging.info("Computing per-chromosome metrics for simulation comparisons.")
    simu_vs_rfmix = compute_comparison_metrics(
        rfmix_df, simu_df, "rfmix", "simu",
        COMPARISON_LABELS["simu_vs_rfmix"],
    )
    simu_vs_flare = compute_comparison_metrics(
        flare_df, simu_df, "flare", "simu",
        COMPARISON_LABELS["simu_vs_flare"],
    )
    logging.info("Computing per-chromosome metrics for LAI comparison.")
    rfmix_vs_flare = compute_comparison_metrics(
        rfmix_df, flare_df, "rfmix", "flare",
        COMPARISON_LABELS["rfmix_vs_flare"],
    )

    summary_df = pd.concat(
        [simu_vs_rfmix, simu_vs_flare, rfmix_vs_flare],
        ignore_index=True,
    )

    summary_path = pop_dir / "global_ancestry_metrics_per_chrom.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    logging.info("Wrote summary metrics to %s", summary_path)

    logging.info("Generating scatter and boxplot figures.")
    plot_scatter_per_ancestry(tidy, pop_dir)
    plot_scatter_by_comparison(tidy, pop_dir)
    plot_boxplots(summary_df, pop_dir)
    logging.info("Saved plots to %s", pop_dir)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
