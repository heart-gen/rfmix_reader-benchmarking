## Script to visualize global ancestry for TWO population ground-truth data
import session_info
from pyhere import here
from pathlib import Path
from rfmix_reader import (
    read_rfmix,
    plot_global_ancestry,
    plot_ancestry_by_chromosome
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
    # Load the two-population ground truth data
    _, rf_q, _ = _load_two_pop_ground_truth()

    # Plot global ancestry
    plot_global_ancestry(
        rf_q,
        dpi=300,
        bbox_inches="tight",
        palette="Set1",      # Two-pop palette
        figsize=(6, 6),
        save_path="global_ancestry.ground_truth_2pop"
    )

    # Plot ancestry by chromosome
    plot_ancestry_by_chromosome(
        rf_q,
        dpi=300,
        bbox_inches="tight",
        palette="Set1",      # Two-pop palette
        save_path="chromosome_summary.ground_truth_2pop"
    )

    # Session information
    session_info.show()

if __name__ == "__main__":
    main()
