## This script generates a BED file of local ancestry
## for use by tagore visualization.

import session_info
from pyhere import here
from rfmix_reader import read_rfmix, generate_tagore_bed

def main():
    # Load the data
    prefix_path = here("input/real_data/_m/")
    binary_dir = here("real_data/gpu_version/_m/binary_files/")
    loci, rf_q, admix = read_rfmix(prefix_path, binary_dir=binary_dir)
    # Extract BED format
    bed = generate_tagore_bed(loci, rf_q, admix, 0)
    # Save data
    bed.to_csv("sample_num_0.bed", index=False)
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
