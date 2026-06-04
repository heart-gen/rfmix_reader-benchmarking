import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

DEFAULT_GENES = ["APOE", "ABCA7", "BIN1", "CR1", "MS4A6A", "CLU"]


def main():
    parser = argparse.ArgumentParser(description="Plot AD SNP local ancestry summaries")
    parser.add_argument("--ancestry", default="functional_results/ad_local_ancestry/_m/ad_snp_ancestry_aggregate.tsv")
    parser.add_argument("--outdir", default="functional_results/ad_local_ancestry/_m/figures")
    parser.add_argument("--genes", nargs="*", default=DEFAULT_GENES)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.ancestry, sep="\t")
    frac_cols = [c for c in df.columns if c.endswith("_fraction")]
    if not frac_cols:
        raise ValueError("No ancestry fraction columns found in aggregate table")

    plot_df = df.copy()
    if "gene" in plot_df.columns:
        plot_df = plot_df[plot_df["gene"].isin(args.genes) | plot_df["gene"].isna()]
    id_col = "variant_id" if "variant_id" in plot_df.columns else "target_id"
    long = plot_df.melt(
        id_vars=[c for c in ["gene", id_col, "chrom", "pos", "matched"] if c in plot_df.columns],
        value_vars=frac_cols,
        var_name="ancestry",
        value_name="fraction",
    ).dropna(subset=["fraction"])
    long["ancestry"] = long["ancestry"].str.replace("_fraction", "", regex=False)
    long["locus"] = long.get("gene", pd.Series(index=long.index, dtype=object)).fillna(long[id_col].astype(str))

    sns.set(style="whitegrid", context="talk")
    plt.figure(figsize=(max(10, 0.35 * long[id_col].nunique()), 6))
    ax = sns.barplot(data=long, x=id_col, y="fraction", hue="ancestry")
    ax.set_xlabel("Target SNP")
    ax.set_ylabel("Ancestry haplotype fraction")
    ax.tick_params(axis="x", rotation=90)
    plt.tight_layout()
    plt.savefig(outdir / "ad_target_ancestry_fractions.png", dpi=300)
    plt.savefig(outdir / "ad_target_ancestry_fractions.pdf")
    print(f"Wrote figures to {outdir}")


if __name__ == "__main__":
    main()
