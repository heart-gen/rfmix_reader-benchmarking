# This script will test the memory usage and executive time for
# RFMix-reader.
import argparse
from time import time
from pyhere import here
from rfmix_reader import read_flare
from memory_profiler import profile

@profile
def read_data(prefix_path):
    loci, rf_q, admix = read_rfmix(prefix_path)
    admix = admix.compute()


def main():
    prefix_path = here("input/simulations/two_populations/_m",
                       "flare-out")
    start_time = time()
    _ = read_data(prefix_path)
    end_time = time()
    print(f"Execution time: {end_time - start_time} seconds")


if __name__ == "__main__":
    main()
