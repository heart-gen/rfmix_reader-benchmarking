## Script to visualize global ancestry for TWO population ground-truth data
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import session_info


def plot_global_ancestry(avg_data, save_prefix="global_ancestry_ground_truth"):
    """
    avg_data: DataFrame with columns ['Sample', 'EUR', 'AFR']
    save_prefix: output filename prefix (without extension)
    """
    df_melt = avg_data.melt(
        id_vars="Sample",
        value_vars=["EUR", "AFR"],
        var_name="Population",
        value_name="Proportion"
    )

    # Define custom color mapping
    palette = {"EUR": "blue", "AFR": "red"}
    
    plt.figure(figsize=(12,6))
    sns.barplot(
        data=df_melt,
        x="Sample",
        y="Proportion",
        hue="Population",
        palette=palette
    )
    plt.title("Global Ancestry Proportions")
    plt.ylabel("Ancestry Proportion")
    plt.xlabel("Sample")
    plt.ylim(0, 1.0)
    plt.xticks(rotation=90)
    plt.legend(title="Population")
    plt.tight_layout()

    # Save in multiple formats
    for ext in ["png", "pdf", "svg"]:
        out_file = f"{save_prefix}.{ext}"
        plt.savefig(out_file, dpi=300, bbox_inches="tight")
        print(f"Saved: {out_file}")

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
    args = parser.parse_args()

    # Build file path
    file_path = Path(args.folder) / args.file

    # Load TSV
    avg_data = pd.read_csv(file_path, sep="\t")

    # Plot and save
    plot_global_ancestry(avg_data, save_prefix=args.output)

    # Session info
    session_info.show()


if __name__ == "__main__":
    main()
