import logging
import argparse
import numpy as np
import pandas as pd
import session_info
from scipy import stats
from pyhere import here
from pathlib import Path
from rfmix_reader import (
    read_simu,
    read_rfmix,
    read_flare,
    phase_rfmix_chromosome_to_zarr,
)

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
        description="Per-chromosome phased global ancestry comparisons for RFMix, simulation, and FLARE."
    )
    parser.add_argument("--simu-input", type=Path, required=True)
    parser.add_argument("--rfmix-input", type=Path, required=True)
    parser.add_argument("--flare-input", type=Path, required=True)
    parser.add_argument("--sample-annot", type=Path, required=True,
                        help="Path to reference sample annotations for phasing")
    parser.add_argument("--ref-input", type=Path, required=True,
                        help="Path to reference zarr root directory (without chr suffix)")
    parser.add_argument("--output", type=Path, default=Path("./"))
    parser.add_argument("--population", type=str, choices=["two", "three"],
                        default="three")
    parser.add_argument("--chrom", type=int, default=None,
                        help="Chromosome to process (1-22). If not provided, process all.")
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


def compute_global_ancestry_from_phased(local_ancestry, labels, sample_ids, chrom):
    """
    Compute per-chromosome global ancestry from phased local ancestry array.

    Parameters
    ----------
    local_ancestry : array-like
        Phased local ancestry array of shape (n_variants, n_haplotypes)
        where n_haplotypes = 2 * n_samples. Values are ancestry labels (0, 1, 2, ...).
    labels : list
        List of ancestry label names (e.g., ["AFR", "EUR", "NAT"])
    sample_ids : array-like
        Sample IDs corresponding to the samples
    chrom : int or str
        Chromosome number

    Returns
    -------
    pd.DataFrame
        Global ancestry per sample for this chromosome with columns:
        sample_id, chrom, and one column per ancestry label
    """
    # Ensure numpy array
    if hasattr(local_ancestry, "compute"):
        local_ancestry = local_ancestry.compute()
    local_ancestry = np.asarray(local_ancestry)

    n_variants, n_haplotypes = local_ancestry.shape
    n_samples = n_haplotypes // 2
    n_ancestries = len(labels)

    # Count ancestry occurrences per haplotype, then combine per sample
    # Global ancestry = proportion of loci assigned to each ancestry
    g_anc_data = []

    for s in range(n_samples):
        # Get both haplotypes for this sample
        hap0 = local_ancestry[:, 2 * s]
        hap1 = local_ancestry[:, 2 * s + 1]

        # Count valid positions (non-negative ancestry labels)
        valid_hap0 = hap0 >= 0
        valid_hap1 = hap1 >= 0

        # Count ancestries across both haplotypes
        ancestry_counts = np.zeros(n_ancestries)
        total_valid = 0

        for a in range(n_ancestries):
            count0 = np.sum((hap0 == a) & valid_hap0)
            count1 = np.sum((hap1 == a) & valid_hap1)
            ancestry_counts[a] = count0 + count1

        total_valid = np.sum(valid_hap0) + np.sum(valid_hap1)

        # Convert to proportions
        if total_valid > 0:
            ancestry_props = ancestry_counts / total_valid
        else:
            ancestry_props = np.full(n_ancestries, np.nan)

        row = {"sample_id": sample_ids[s], "chrom": chrom}
        for i, label in enumerate(labels):
            row[label] = ancestry_props[i]
        g_anc_data.append(row)

    return pd.DataFrame(g_anc_data)


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
        # Set dynamic axis labels
        x_label = data["x_label"].iloc[0]
        y_label = data["y_label"].iloc[0]
        ax.set_xlabel(x_label, fontweight="bold")
        ax.set_ylabel(y_label, fontweight="bold")
    # Save figure
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


def main():
    configure_logging()
    args = parse_args()

    pop_dir = args.output / f"{args.population}_populations"
    pop_dir.mkdir(parents=True, exist_ok=True)

    phased_dir = pop_dir / "phased_files"
    phased_dir.mkdir(parents=True, exist_ok=True)

    binary_path = args.rfmix_input / "binary_files"

    # Determine which chromosomes to process
    chroms = [args.chrom] if args.chrom is not None else list(range(1, 23))
    single_chrom = args.chrom is not None

    all_g_anc_simu = []
    all_g_anc_rfmix = []
    all_g_anc_flare = []

    for chrom in chroms:
        chrom_label = f"chr{chrom}"
        logging.info("Processing chromosome %s", chrom_label)

        # Load simulated global ancestry (ground truth)
        logging.info("Loading simulated global ancestry for %s", chrom_label)
        _, g_anc_simu, _ = read_simu(here(args.simu_input), chrom=chrom)
        all_g_anc_simu.append(g_anc_simu)

        # Load FLARE global ancestry
        logging.info("Loading FLARE global ancestry for %s", chrom_label)
        _, g_anc_flare, _ = read_flare(here(args.flare_input), chrom=chrom)
        all_g_anc_flare.append(g_anc_flare)

        # Phase RFMix outputs and compute global ancestry
        logging.info("Phasing RFMix outputs for %s", chrom_label)

        # Get sample IDs from RFMix
        loci_rfmix, g_anc_rfmix_unphased, _ = read_rfmix(
            here(args.rfmix_input), binary_dir=here(binary_path), chrom=chrom
        )
        sample_ids = g_anc_rfmix_unphased["sample_id"].values
        labels = [c for c in g_anc_rfmix_unphased.columns if c not in {"sample_id", "chrom"}]

        # Build reference path for this chromosome
        ref_zarr_path = args.ref_input / f"1kGP.chr{chrom}.filtered.snpsOnly.afr_washington.zarr"

        output_path = phased_dir / f"phased_chr{chrom}.zarr"

        local_array = phase_rfmix_chromosome_to_zarr(
            file_prefix=here(args.rfmix_input),
            ref_zarr_root=here(ref_zarr_path),
            binary_dir=here(binary_path),
            sample_annot_path=here(args.sample_annot),
            output_path=here(output_path),
            chrom=str(chrom),
        )

        # Extract phased local ancestry
        logging.info("Computing phased global ancestry for %s", chrom_label)
        phased_local = local_array["local_ancestry"]

        # Compute global ancestry from phased local ancestry
        g_anc_rfmix = compute_global_ancestry_from_phased(
            phased_local, labels, sample_ids, chrom
        )
        all_g_anc_rfmix.append(g_anc_rfmix)

    # Concatenate all chromosomes
    logging.info("Concatenating global ancestry across processed chromosomes")
    g_anc_simu = pd.concat(all_g_anc_simu, ignore_index=True)
    g_anc_rfmix = pd.concat(all_g_anc_rfmix, ignore_index=True)
    g_anc_flare = pd.concat(all_g_anc_flare, ignore_index=True)

    logging.info(
        "Loaded shapes - simu: %s, rfmix: %s, flare: %s",
        g_anc_simu.shape, g_anc_rfmix.shape, g_anc_flare.shape,
    )

    simu_cols = set(g_anc_simu.columns)
    rfmix_cols = set(g_anc_rfmix.columns)
    flare_cols = set(g_anc_flare.columns)

    if simu_cols != rfmix_cols or simu_cols != flare_cols:
        logging.warning(
            "Column mismatch - simu: %s, rfmix: %s, flare: %s",
            simu_cols, rfmix_cols, flare_cols
        )

    logging.info("Reshaping global ancestry tables to tidy format.")

    # Save per-chromosome or combined file based on mode
    if single_chrom:
        chrom_dir = pop_dir / f"chr{args.chrom}"
        chrom_dir.mkdir(parents=True, exist_ok=True)
        tidy_path = chrom_dir / "global_ancestry_per_chrom.tsv"
    else:
        tidy_path = pop_dir / "global_ancestry_per_chrom.tsv"

    tidy = pd.concat(
        [
            tidy_g_anc(g_anc_simu, "simu"),
            tidy_g_anc(g_anc_rfmix, "rfmix"),
            tidy_g_anc(g_anc_flare, "flare"),
        ],
        ignore_index=True,
    )
    tidy.to_csv(tidy_path, sep="\t", index=False)

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

    if single_chrom:
        summary_path = chrom_dir / "global_ancestry_metrics_per_chrom.tsv"
    else:
        summary_path = pop_dir / "global_ancestry_metrics_per_chrom.tsv"

    summary_df.to_csv(summary_path, sep="\t", index=False)
    logging.info("Wrote summary metrics to %s", summary_path)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
