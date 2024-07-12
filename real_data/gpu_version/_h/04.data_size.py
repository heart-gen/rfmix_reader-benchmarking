# This script will test the memory usage and executive time for
# RFMix-reader.
import argparse
from time import time
from pyhere import here
from rfmix_reader import read_rfmix

def read_data(prefix_path, BINARIES):
    loci, rf_q, admix = read_rfmix(prefix_path, generate_binary=BINARIES)
    print(f"Loci size: {loci.shape[0]}")


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description='Measure memory usage and execution time.')
    parser.add_argument('--binaries',  action="store_true",
                        help='Create binaries')
    # Parse the arguments
    args = parser.parse_args()
    prefix_path = here("input/real_data/_m/")
    start_time = time()
    _ = read_data(prefix_path, args.binaries)
    end_time = time()
    print(f"Execution time: {end_time - start_time} seconds")


if __name__ == "__main__":
    main()
