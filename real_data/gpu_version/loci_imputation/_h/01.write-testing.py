# This script will test the memory usage and executive time for
# RFMix-reader.
import session_info
from time import time
from pyhere import here
from rfmix_reader import (
    read_rfmix,
    write_data,
    write_imputed,
    interpolate_array
)

def _read_data(prefix_path, BIN_DIR):
    loci, rf_q, admix = read_rfmix(prefix_path, binary_dir=BIN_DIR)
    loci.rename(columns={"chromosome": "chrom",
                         "physical_position": "pos"},
                inplace=True)
    return loci, rf_q, admix


def _load_genotypes(plink_prefix_path):
    from tensorqtl import pgen
    pgr = pgen.PgenReader(plink_prefix_path)
    variant_df = pgr.variant_df
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pgr.load_genotypes(), variant_df


def _clean_data(loci):
    plink_prefix = here("input/genotypes/TOPMed_LIBD")
    _, variant_df = _load_genotypes(plink_prefix)
    variant_df = variant_df.drop_duplicates(subset=["chrom","pos"],keep='first')
    variant_loci_df = variant_df.merge(loci.to_pandas(), on=["chrom", "pos"],
                                       how="outer", indicator=True)\
                                .loc[:, ["chrom", "pos", "i", "_merge"]]
    return variant_loci_df


def main():
    # Load loci data
    prefix_path = here("input/real_data/_m/")
    binary_dir = here(prefix_path, "binary_files")
    loci, rf_q, admix = _read_data(prefix_path, binary_dir)
    # Write data
    start_time = time()
    write_data(loci, rf_q, admix, prefix="local-ancestry")
    end_time = time()
    print(f"Execution time: {end_time - start_time} seconds")

    # Interpolate data
    variant_loci = _clean_data(loci)
    start_time = time()
    z = interpolate_array(variant_loci, admix, prefix_path)
    end_time = time()
    print(f"Execution time: {end_time - start_time} seconds")

    # Write imputed data
    start_time = time()
    write_imputed(rf_q, admix, variant_loci, z, target_rows=200_000,
                  prefix="imputed-ancestry", verbose=True)
    end_time = time()
    print(f"Execution time: {end_time - start_time} seconds")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
