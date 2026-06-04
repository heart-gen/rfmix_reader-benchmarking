import argparse
from pathlib import Path

import pandas as pd
from rfmix_reader import extract_locus_ancestry

DEFAULT_RFMIX_DIR = "input/aanri_data/rfmix-version/_m"


def read_samples(path: str | None) -> list[str] | None:
    if path is None:
        return None
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def main():
    parser = argparse.ArgumentParser(description="Extract AD SNP-level MSP local ancestry")
    parser.add_argument("--targets", default="functional_results/ad_local_ancestry/_m/ad_target_snps.tsv")
    parser.add_argument("--rfmix-dir", default=DEFAULT_RFMIX_DIR)
    parser.add_argument("--outdir", default="functional_results/ad_local_ancestry/_m")
    parser.add_argument("--samples", default=None, help="Optional one-sample-id-per-line subset")
    parser.add_argument("--sample-level", action="store_true", help="Also write sample-level ancestry calls")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    targets = pd.read_csv(args.targets, sep="\t")
    samples = read_samples(args.samples)

    aggregate_parts = []
    sample_parts = []
    for chrom, chrom_targets in targets.groupby(targets["chrom"].astype(str), sort=True):
        aggregate = extract_locus_ancestry(args.rfmix_dir, chrom_targets, samples=samples, aggregate=True)
        aggregate_parts.append(aggregate)
        aggregate.to_csv(outdir / f"ad_snp_ancestry_chr{chrom}.tsv", sep="\t", index=False)
        if args.sample_level:
            sample_df = extract_locus_ancestry(args.rfmix_dir, chrom_targets, samples=samples, aggregate=False)
            sample_parts.append(sample_df)

    aggregate_all = pd.concat(aggregate_parts, ignore_index=True) if aggregate_parts else pd.DataFrame()
    aggregate_all.to_csv(outdir / "ad_snp_ancestry_aggregate.tsv", sep="\t", index=False)
    try:
        aggregate_all.to_parquet(outdir / "ad_snp_ancestry_aggregate.parquet", index=False)
    except Exception as exc:
        print(f"Skipping aggregate parquet: {exc}")

    if args.sample_level:
        sample_all = pd.concat(sample_parts, ignore_index=True) if sample_parts else pd.DataFrame()
        try:
            sample_all.to_parquet(outdir / "ad_snp_ancestry_sample_level.parquet", index=False)
        except Exception as exc:
            fallback = outdir / "ad_snp_ancestry_sample_level.tsv.gz"
            sample_all.to_csv(fallback, sep="\t", index=False, compression="gzip")
            print(f"Wrote sample-level TSV fallback because parquet failed: {exc}")

    print(f"Wrote aggregate ancestry for {len(aggregate_all)} target rows")


if __name__ == "__main__":
    main()
