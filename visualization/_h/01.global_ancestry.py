## Script to visualize global ancestry
import session_info
from pyhere import here
from pathlib import Path
from rfmix_reader import (
    read_rfmix,
    plot_global_ancestry,
    plot_ancestry_by_chromosome
)

def _load_real_data():
    prefix_path = Path(here("input/real_data/rfmix-version/_m/"))
    binary_dir = prefix_path / "binary_files/"
    if binary_dir.exists():
        return read_rfmix(prefix_path, binary_dir=binary_dir)
    else:
        return read_rfmix(prefix_path, binary_dir=binary_dir,
                          generate_binary=True)


def _load_simu_data(pop_num):
    pop_loc = "two_populations" if pop_num == 2 else "three_populations"
    prefix_path = Path(here("input/simulations")) / pop_loc / "_m/rfmix-out/"
    binary_dir = prefix_path / "binary_files"
    if binary_dir.exists():
        return read_rfmix(prefix_path, binary_dir=binary_dir)
    else:
        return read_rfmix(prefix_path, binary_dir=binary_dir,
                          generate_binary=True)


def main():
    # Plot real data
    _, rf_q, _ = _load_real_data()
    plot_global_ancestry(
        rf_q, dpi=300, bbox_inches="tight", palette="Set1",
        figsize=(6,6), save_path="global_ancestry.real_data_2pop"
    )
    plot_ancestry_by_chromosome(
        rf_q, dpi=300, bbox_inches="tight", palette="Set1",
        save_path="chromosome_summary.real_data_2pop"
    )
    # Plot simulated data, three populations
    _, rf_q, _ = _load_simu_data(3)
    plot_global_ancestry(
        rf_q, dpi=300, bbox_inches="tight", palette="tab10",
        figsize=(6,6), save_path="global_ancestry.simu_data_3pop"
    )
    plot_ancestry_by_chromosome(
        rf_q, dpi=300, bbox_inches="tight", palette="tab10",
        save_path="chromosome_summary.simu_data_3pop"
    )
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
