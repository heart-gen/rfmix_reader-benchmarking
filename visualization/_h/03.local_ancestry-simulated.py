## Script to visualize global ancestry
import session_info
from pyhere import here
from pathlib import Path
from rfmix_reader import (
    read_rfmix,
    generate_tagore_bed,
    plot_local_ancestry_tagore
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
    # General configuration
    sample_num = 12 ## This will be the 13th index
    build = "hg38"; prefix = f"local_ancestry.{build}"
    verbose = True; force = True ## Overwrite

    # Plot simulated data, three populations
    loci, rf_q, admix = _load_simu_data(3)
    bed_df = generate_tagore_bed(loci, rf_q, admix, sample_num)
    for oformat in ["png", "pdf"]:
        out_prefix = f"{prefix}.simu_data_3pop"
        plot_local_ancestry_tagore(bed_df, out_prefix, build,
                                   oformat, verbose, force)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
