## Script to visualize local ancestry for the ground-truth data
import session_info
from pyhere import here
from pathlib import Path
from rfmix_reader import (
    read_rfmix,
    generate_tagore_bed,
    plot_local_ancestry_tagore
)

def _load_two_pop_ground_truth():
    prefix_path = Path(
        here("input/simulations/two_populations/ground_truth/_m/")
    )
    binary_dir = prefix_path / "binary_files"

    if binary_dir.exists():
        return read_rfmix(prefix_path, binary_dir=binary_dir)
    else:
        return read_rfmix(
            prefix_path,
            binary_dir=binary_dir,
            generate_binary=True
        )

def main():
    # General configuration
    sample_num = 1
    build = "hg38"
    prefix = f"local_ancestry.{build}"

    # Load the two-population ground truth data
    loci, rf_q, admix = _load_two_pop_ground_truth()

    # Generate BED dataframe for plotting
    bed_df = generate_tagore_bed(loci, rf_q, admix, sample_num)

    # Output file prefix
    out_prefix = f"{prefix}.ground_truth_2pop"

    # Plot in both PNG and PDF formats
    for oformat in ["png", "pdf"]:
        plot_local_ancestry_tagore(bed_df, out_prefix, build, oformat)

    # Show session info
    session_info.show()

if __name__ == "__main__":
    main()
