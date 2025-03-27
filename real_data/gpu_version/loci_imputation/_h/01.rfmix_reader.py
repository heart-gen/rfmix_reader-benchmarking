# This script will test the memory usage and executive time for
# RFMix-reader.
from time import time
from pyhere import here
from memory_profiler import profile
from rfmix_reader import read_rfmix, interpolate_array

def read_data(prefix_path, BIN_DIR):
    return read_rfmix(prefix_path, binary_dir=BIN_DIR)


def _load_genotypes(plink_prefix_path):
    from tensorqtl import pgen
    pgr = pgen.PgenReader(plink_prefix_path)
    variant_df = pgr.variant_df
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pgr.load_genotypes(), variant_df


def main():
    # Load loci data
    prefix_path = here("input/real_data/_m/")
    binary_dir = here(prefix_path, "binary_files")
    loci, rf_q, admix = read_data(prefix_path, binary_dir)
    loci.rename(columns={"chromosome": "chrom",
                         "physical_position": "pos"},
                inplace=True)
    # Load variant information
    plink_prefix = here("input/genotypes/TOPMed_LIBD")
    _, variant_df = _load_genotypes(plink_prefix)
    variant_df = variant_df.drop_duplicates(subset=["chrom","pos"],keep='first')
    variant_loci_df = variant_df.merge(loci.to_pandas(), on=["chrom", "pos"],
                                       how="outer", indicator=True)\
                                .loc[:, ["chrom", "pos", "i", "_merge"]]
    z = interpolate_array(variant_loci_df, admix, prefix_path)
