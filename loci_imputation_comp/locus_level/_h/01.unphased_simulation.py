import json
import logging
import argparse
import numpy as np
import pandas as pd
import session_info
from pyhere import here
from pathlib import Path

from sklearn.metrics import (
    f1_score,
    recall_score,
    precision_score,
    confusion_matrix,
    matthews_corrcoef
)

from localqtl import PgenReader
from rfmix_reader import read_rfmix, read_simu, interpolate_array

def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_parameters():
    parser = argparse.ArgumentParser(description="Locus-Level Imputation Accuracy")
    parser.add_argument("--rfmix-input", type=Path, required=True)
    parser.add_argument("--simu-input", type=Path, required=True)
    parser.add_argument("--population", type=str, choices=["two","three"], default="three")
    parser.add_argument("--method", type=str, choices=["linear", "nearest", "stepwise"],
                        default="linear")
    return parser.parse_args()


def load_variants(plink_path):
    pgr = PgenReader(plink_path)
    variant_df = pgr.variant_df.copy()
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return variant_df[~variant_df.duplicated(subset=["chrom", "pos"], keep="first")]


def impute_data(loci_df, variants, admix, zarr_path, method):
    if hasattr(loci_df, "to_pandas"):
        loci_df = loci_df.to_pandas()
    variant_df = variants.merge(loci_df, on=["chrom", "pos"], how="outer")
    variant_df = variant_df.loc[:, ["chrom", "pos", "i"]]
    z = interpolate_array(variant_df, admix, zarr_outdir=zarr_path,
                          interpolation=method, use_bp_positions=True)
    return variant_df, z


def admix_to_haplotypes(admix):
    V, S, A = admix.shape
    haplo = np.full((V, S * 2), -1, dtype=np.int8)
    for v in range(V):
        for s in range(S):
            counts = admix[v, s]
            labels = []
            for a, count in enumerate(counts):
                labels.extend([a] * int(count))
            if len(labels) == 2:
                haplo[v, s * 2] = labels[0]
                haplo[v, s * 2 + 1] = labels[1]
    return haplo


def _to_hard_calls(arr):
    if arr.ndim == 3:
        return admix_to_haplotypes(arr)
    return arr


def compute_metrics_one_sample(t, p, labels):
    cm = confusion_matrix(t, p, labels=labels)
    prec = precision_score(t, p, labels=labels, average=None, zero_division=0)
    rec = recall_score(t, p, labels=labels, average=None, zero_division=0)
    f1 = f1_score(t, p, labels=labels, average=None, zero_division=0)
    acc = np.mean(t == p)
    mcc = matthews_corrcoef(t, p)
    return {
        "accuracy": acc,
        "mcc": mcc,
        "precision": prec,
        "recall": rec,
        "f1": f1,
        "confusion_matrix": cm
    }


def compute_per_sample_metrics(true_anc, inferred_anc):
    true_hard = _to_hard_calls(true_anc)
    pred_hard = _to_hard_calls(inferred_anc)
    labels = np.unique(true_hard[true_hard >= 0])
    _, n_samples = true_hard.shape
    metrics = {}
    for s in range(n_samples):
        t = true_hard[:, s]
        p = pred_hard[:, s]
        mask = (t >= 0) & (p >= 0)
        metrics[s] = compute_metrics_one_sample(t[mask], p[mask], labels)
    return metrics


def write_sample_metrics_tsv(sample_metrics, outfile):
    rows = []
    for s, m in sample_metrics.items():
        row = {
            "haplotype": s,
            "accuracy": m["accuracy"],
            "mcc": m["mcc"],
            "precision": ",".join(map(str, m["precision"])),
            "recall": ",".join(map(str, m["recall"])),
            "f1": ",".join(map(str, m["f1"])),
            "confusion_matrix": ";".join(
                ",".join(map(str, row)) for row in m["confusion_matrix"]
            )
        }
        rows.append(row)
    pd.DataFrame(rows).to_csv(outfile, sep="\t", index=False)


def compute_locus_metrics_json(t, p, method, outfile=None):
    t_hard = _to_hard_calls(t)
    p_hard = _to_hard_calls(p)
    labels = list(np.unique(t_hard[t_hard >= 0]))
    accuracy = np.mean(t_hard == p_hard)
    cm = confusion_matrix(t_hard.ravel(), p_hard.ravel(), labels=labels).tolist()
    precision = precision_score(t_hard.ravel(), p_hard.ravel(), labels=labels, average=None, zero_division=0)
    recall = recall_score(t_hard.ravel(), p_hard.ravel(), labels=labels, average=None, zero_division=0)
    f1 = f1_score(t_hard.ravel(), p_hard.ravel(), labels=labels, average=None, zero_division=0)
    mcc = matthews_corrcoef(t_hard.ravel(), p_hard.ravel())
    metrics = {
        "method": method,
        "labels": labels,
        "overall_accuracy": float(accuracy),
        "confusion_matrix": cm,
        "precision": {int(k): float(v) for k, v in zip(labels, precision)},
        "recall": {int(k): float(v) for k, v in zip(labels, recall)},
        "f1": {int(k): float(v) for k, v in zip(labels, f1)},
        "mcc": float(mcc),
        "shape": list(map(int, t.shape))
    }
    if outfile:
        with open(outfile, "w") as f:
            json.dump(metrics, f, indent=4)
    return metrics


def main():
    configure_logging()
    args = parse_parameters()
    pop_dir = Path(f"{args.population}_populations")
    pop_dir.mkdir(parents=True, exist_ok=True)

    # Load and impute data
    logging.info("Loading PLINK variants...")
    plink_path = here("input/simulations", pop_dir.name, "_m/plink-files/simulated")
    variants_df = load_variants(plink_path)

    logging.info("Reading RFMix outputs...")
    binary_path = args.rfmix_input / "binary_files"
    loci_df,_,admix = read_rfmix(here(str(args.rfmix_input)), binary_dir=str(binary_path))
    method_path = here(args.rfmix_input / args.method)
    method_path.mkdir(parents=True, exist_ok=True)
    zarr_path = method_path / "imputed_local_ancestry"

    logging.info("Interpolating ancestry data...")
    loci_df = loci_df.rename(columns={"chromosome": "chrom", "physical_position": "pos"})
    variant_loci_df, inferred_anc = impute_data(
        loci_df, variants_df, admix, here(str(zarr_path)), args.method
    )

    # Filter and align
    variant_loci_df = variant_loci_df[variant_loci_df.set_index(["chrom", "pos"]).index.isin(
        variants_df.set_index(["chrom", "pos"]).index
    )]
    idx = variant_loci_df.index.values
    inferred_anc = da.from_zarr(inferred_anc)[idx, :, :]
    inferred_anc = inferred_anc.rechunk((50_000, 100, -1)).compute()

    # Load ground truth data
    logging.info("Loading ground truth data")
    loci_gt, _, admix_gt = read_simu(here(str(args.simu_input)))
    loci_gt = loci_gt.rename(columns={"chromosome": "chrom", "physical_position": "pos"})
    loci_gt = loci_gt[loci_gt.set_index(["chrom", "pos"]).index.isin(
        variant_loci_df.set_index(["chrom", "pos"]).index
    )]
    loci_gt = loci_gt.reset_index(drop=True)
    true_anc = true_anc[loci_gt.index.values, :, :]
    true_anc = admix_gt.compute()

    # Calculate metrics
    logging.info("Computing locus-level metrics...")
    json_outfile = args.output / "locus_metrics.json"
    compute_locus_metrics_json(true_anc, inferred_anc, args.method, outfile=json_outfile)

    logging.info("Computing per-haplotype metrics...")
    sample_metrics = compute_per_sample_metrics(true_anc, inferred_anc)
    tsv_outfile = args.output / "per_haplotype_metrics.tsv"
    write_sample_metrics_tsv(sample_metrics, tsv_outfile)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
