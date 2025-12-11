import logging
import argparse
import pandas as pd
import session_info
import xarray as xr
from pyhere import here
from pathlib import Path

from localqtl import PgenReader
from rfmix_reader import read_rfmix, interpolate_array

def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_parameters():
    parser = argparse.ArgumentParser(description="Locus-Level Imputation Accuracy")
    parser.add_argument("--rfmix-input", type=Path, required=True)
    parser.add_argument("--phased-zarr", type=Path, required=True)
    parser.add_argument("--population", type=str, choices=["two","three"], default="three")
    parser.add_argument("--method", type=str, choices=["linear", "nearest", "stepwise"],
                        default="linear")
    return parser.parse_args()


def load_variants(plink_path):
    pgr = PgenReader(plink_path)
    variant_df = pgr.variant_df.copy()
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return variant_df[~variant_df.duplicated(subset=["chrom", "pos"], keep="first")]


def standardize_variant_columns(df):
    return df.rename(columns={"chromosome": "chrom", "physical_position": "pos"})


def impute_data(loci_df, variants, admix, zarr_path, method):
    if hasattr(loci_df, "to_pandas"):
        loci_df = loci_df.to_pandas()
    variant_df = variants.merge(loci_df, on=["chrom", "pos"], how="outer")
    variant_df = variant_df.loc[:, ["chrom", "pos", "i"]]
    _ = interpolate_array(variant_df, admix, zarr_outdir=zarr_path,
                          interpolation=method, use_bp_positions=True)
    return variant_df


def align_variants(df_ref, df_target):
    ref_idx = df_ref.set_index(["chrom", "pos"])
    target_idx = df_target.set_index(["chrom", "pos"])
    common = target_idx.index.intersection(ref_idx.index)
    return target_idx.loc[common].reset_index()


def main():
    configure_logging()
    args = parse_parameters()
    pop_dir = Path(f"{args.population}_populations")

    # Load and impute data
    logging.info("Loading PLINK variants...")
    plink_path = here("input/simulations", pop_dir.name, "_m/plink-files/simulated")
    variants_df = load_variants(plink_path)

    logging.info("Reading RFMix outputs...")
    binary_path = args.rfmix_input / "binary_files"
    loci_df, _, _ = read_rfmix(here(args.rfmix_input), binary_dir=here(binary_path))
    loci_df = standardize_variant_columns(loci_df)

    logging.info("Reading phased data...")
    phased_path = args.phased_zarr / "phased_data.zarr"
    local_ancestry = xr.open_zarr(here(phased_path))
    admix = local_ancestry["local_ancestry"].chunk({"variant": 20_000, "sample": 100})

    logging.info("Interpolating ancestry data...")
    method_path = here(args.rfmix_input) / "phased" / args.method
    method_path.mkdir(parents=True, exist_ok=True)
    zarr_path = method_path / "imputed_local_ancestry"

    variant_loci_df = impute_data(
        loci_df, variants_df, admix, zarr_path, args.method
    )
    variant_loci_df = align_variants(variants_df, variant_loci_df)
    variant_path = f"{method_path}/imputed_variant.parquet"
    variant_loci_df.to_parquet(variant_path)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
