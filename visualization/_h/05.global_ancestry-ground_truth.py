## Script to visualize global ancestry for TWO population ground-truth data
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import session_info

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
    plot_global_ancestry(avg_data, save_path=args.output)

    # Session info
    session_info.show()


if __name__ == "__main__":
    main()
