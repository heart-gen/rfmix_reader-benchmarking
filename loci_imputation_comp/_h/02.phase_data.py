import logging
import argparse
import session_info
from pyhere import here
from pathlib import Path
from rfmix_reader import phase_rfmix_chromosome_to_zarr


def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_parameters():
    parser = argparse.ArgumentParser(description="Locus-Level Imputation Accuracy")
    parser.add_argument("--rfmix-input", type=Path, required=True)
    parser.add_argument("--ref-input", type=Path, default=Path("input/references/_m"))
    parser.add_argument("--chrom", type=int, default=21)
    return parser.parse_args()


def main():
    configure_logging()
    args = parse_parameters()

    # Phase output per chromosome
    logging.info("Phase RFMix outputs per chromosome...")
    binary_path = args.rfmix_input / "binary_files"
    phased_path = args.rfmix_input / "phased_files"
    ref_zarr = args.ref_input / "reference_zarr"
    sample_annot_path = args.ref_input / "samples_id2"
    output_path = f"{phased_path}/phased_chr{args.chrom}.zarr"
    _ = phase_rfmix_chromosome_to_zarr(
        file_prefix=here(args.rfmix_input),
        ref_zarr_root=here(ref_zarr),
        binary_dir=here(args.ref_zarr),
        sample_annot_path=here(sample_annot_path),
        output_path=here(output_path),
        chrom=str(args.chrom),
    )

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
