## Script to visualize global ancestry for TWO population ground-truth data
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import session_info
import re

def plot_global_ancestry(avg_data, save_path="global_ancestry_ground_truth"):
    # Ensure correct order of columns
    avg_data = avg_data[["Sample", "EUR", "AFR"]]

    # Set up figure
    plt.figure(figsize=(12, 6))

    # Stacked barplot: AFR on bottom, EUR on top
    plt.bar(
        avg_data["Sample"],
        avg_data["AFR"],
        label="AFR",
        color="red"
    )
    plt.bar(
        avg_data["Sample"],
        avg_data["EUR"],
        bottom=avg_data["AFR"],
        label="EUR",
        color="blue"
    )

    plt.title("Global Ancestry Proportions")
    plt.ylabel("Ancestry Proportion")
    plt.xlabel("Individuals")
    plt.xticks([], [])  # remove cluttered x-ticks
    plt.legend(title="Population")
    plt.tight_layout()

    # Save in multiple formats
    for ext in ["png", "pdf", "svg"]:
        out_file = f"{save_path}.{ext}"
        plt.savefig(out_file, dpi=300)
        print(f"Saved {out_file}")
    plt.close()

def plot_ancestry_whiskers(folder_path, out_prefix="global_ancestry_whiskers"):
    import re
    from matplotlib.patches import Patch

    folder_path = Path(folder_path).resolve()
    out_prefix = Path(out_prefix).resolve()

    all_data = []
    file_pattern = re.compile(r"chr(\d+)\.tsv")

    for file in folder_path.glob("*.tsv"):
        match = file_pattern.search(file.name)
        if not match:
            continue
        chrom = int(match.group(1))

        df = pd.read_csv(file, sep="\t")
        df = df.rename(columns={"YRI": "AFR", "CEU": "EUR"})
        df["Chromosome"] = chrom

        df_melt = df.melt(
            id_vars=["Sample", "Chromosome"],
            value_vars=["AFR", "EUR"],
            var_name="Ancestry",
            value_name="Proportion"
        )
        all_data.append(df_melt)

    if not all_data:
        raise ValueError(f"No valid .tsv files found in {folder_path}")

    data = pd.concat(all_data, ignore_index=True)
    data = data[data["Chromosome"].between(1, 22)]
    data["Chromosome"] = pd.Categorical(
        data["Chromosome"], categories=list(range(1, 22)), ordered=True
    )

    colors = {"AFR": "red", "EUR": "blue"}

    plt.figure(figsize=(14, 6))
    ax = sns.boxplot(
        data=data,
        x="Chromosome",
        y="Proportion",
        hue="Ancestry",
        palette=colors,
        showcaps=True,
        fliersize=1,
        linewidth=1.2
    )

    ax.set_title("Ancestry Proportion per Chromosome", fontsize=14, weight="bold")
    ax.set_ylabel("Ancestry Proportion")
    ax.set_xlabel("Chromosome")
    ax.set_ylim(0, 1)

    # Move legend outside top-right
    ax.legend(title="Ancestry", loc="upper left", bbox_to_anchor=(1.02, 1))

    plt.tight_layout()

    for ext in ["pdf", "png", "svg"]:
        plt.savefig(out_prefix.with_suffix(f".{ext}"), bbox_inches="tight")
        print(f"Saved {out_prefix.with_suffix(f'.{ext}')}")

    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description="Visualize global ancestry proportions for TWO population data."
    )
    parser.add_argument(
        "--folder", type=str, required=True,
        help="Path to the folder containing the TSV file"
    )
    parser.add_argument(
        "--file", type=str, required=True,
        help="Name of the TSV file with ancestry averages"
    )
    parser.add_argument(
        "--output", type=str, default="global_ancestry_ground_truth",
        help="Output filename prefix (default: global_ancestry_ground_truth)"
    )
    parser.add_argument(
        "--chromosome_plots", action="store_true",  # [NEW] flag to also produce whisker plots
        help="If set, also generate per-chromosome ancestry whisker plots"
    )
    
    args = parser.parse_args()

    # Build file path
    file_path = Path(args.folder) / args.file

    # Load TSV
    avg_data = pd.read_csv(file_path, sep="\t")

    # Plot and save
    plot_global_ancestry(avg_data, save_path=args.output)

    # [NEW] Optional per-chromosome whisker plots
    if args.chromosome_plots:
        plot_ancestry_whiskers(args.folder, out_prefix="global_ancestry_whiskers")
        
    # Session info
    session_info.show()


if __name__ == "__main__":
    main()
