# --------------------------
# Load Libraries
# --------------------------

# Standard libraries
import session_info
from pyhere import here
from pathlib import Path

# Project specific
import argparse
import pandas as pd

# --------------------------
# Compute global ancestry
# average per sample
# --------------------------

def compute_global_ancestry(folder_path: Path, file_names: list, output_file: str = "global_ancestry.tsv"):
   
    dfs = []
    for fname in file_names:
        file_path = folder_path / fname
        df = pd.read_csv(file_path, sep="\t")

        # Rename CEU -> EUR, YRI -> AFR
        df = df.rename(columns={"CEU": "EUR", "YRI": "AFR"})

        # Keep only required columns
        df = df[["Sample", "EUR", "AFR"]]
        dfs.append(df)

    # Combine all chromosomes
    all_data = pd.concat(dfs, ignore_index=True)

    # Compute average per sample
    avg_data = all_data.groupby("Sample", as_index=False).mean()

    # Save output
    output_path = folder_path / output_file
    avg_data.to_csv(output_path, sep="\t", index=False)

    print(f" Global ancestry averages saved to {output_path}")

# --------------------------
# MAIN
#--------------------------

def main():
    parser = argparse.ArgumentParser(description="Compute average ancestry proportions from TSV files.")
    parser.add_argument(
        "--folder_path",
        type=str,
        required=True,
        help="Path to folder containing TSV files (relative to project root)."
    )
    parser.add_argument(
        "--file_names",
        nargs="+",
        required=True,
        help="List of TSV file names to process."
    )
    parser.add_argument(
        "--output_file",
        type=str,
        default="global_ancestry.tsv",
        help="Output TSV filename."
    )
    args = parser.parse_args()

    # Resolve folder path relative to project root
    folder_path = Path(args.folder_path).resolve()

    compute_global_ancestry(folder_path, args.file_names, args.output_file)
    session_info.show()
    
if __name__ == '__main__':
    main()