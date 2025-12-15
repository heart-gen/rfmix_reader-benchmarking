import json
import logging
import argparse
import numpy as np
import pandas as pd
import session_info
from pyhere import here
import dask.array as da
from pathlib import Path

from localqtl import PgenReader
from rfmix_reader import read_simu, read_rfmix, interpolate_array

def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_parameters():
    parser = argparse.ArgumentParser(description="Imputation Accuracy")
    parser.add_argument("--rfmix-input", type=Path, required=True)
    parser.add_argument("--simu-input", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("./"))
    parser.add_argument("--population", type=str, choices=["two","three"], default="three")
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
    admix = interpolate_array(variant_df, admix, zarr_outdir=zarr_path,
                              interpolation=method, use_bp_positions=True)
    return variant_df, admix


def admix_counts_to_haplotypes_numpy(admix):
    """
    admix: (V, S, A) integer counts summing to 2
    returns: (V, 2S) haplotype labels or -1
    """
    V, S, A = admix.shape
    out = np.full((V, 2 * S), -1, dtype=np.int8)

    admix = admix.astype(np.int16, copy=False)
    valid = admix.sum(axis=2) == 2

    csum = np.cumsum(admix, axis=2)
    hap0 = (csum >= 1).argmax(axis=2)
    hap1 = (csum >= 2).argmax(axis=2)

    out[:, 0::2] = np.where(valid, hap0, -1)
    out[:, 1::2] = np.where(valid, hap1, -1)
    return out


def to_hard_calls(arr):
    """Convert (V,S,A) -> (V, 2S)"""
    if arr.ndim == 3:
        return admix_counts_to_haplotypes_numpy(arr)
    return arr


def confusion_counts(t_hard, p_hard, n_labels):
    """Compute confusion matrix counts for hard calls."""
    t = np.asarray(t_hard).ravel()
    p = np.asarray(p_hard).ravel()

    mask = (t >= 0) & (p >= 0)
    if not np.any(mask):
        return np.zeros((n_labels, n_labels), dtype=np.int64)

    idx = t[mask].astype(np.int64) * n_labels + p[mask].astype(np.int64)
    cm = np.bincount(idx, minlength=n_labels * n_labels).reshape(n_labels, n_labels)

    return cm


def metrics_from_confusion(cm):
    K = cm.shape[0]
    total = cm.sum()
    acc = np.trace(cm) / total if total else np.nan

    tp = np.diag(cm).astype(float)
    fp = cm.sum(axis=0) - tp
    fn = cm.sum(axis=1) - tp

    precision = np.divide(tp, tp + fp, out=np.zeros_like(tp), where=(tp + fp) > 0)
    recall    = np.divide(tp, tp + fn, out=np.zeros_like(tp), where=(tp + fn) > 0)
    f1        = np.divide(2 * precision * recall,
                          precision + recall,
                          out=np.zeros_like(tp),
                          where=(precision + recall) > 0)

    # Multiclass MCC (Gorodkin)
    t_sum = cm.sum(axis=1).astype(float)
    p_sum = cm.sum(axis=0).astype(float)
    s = cm.sum().astype(float)
    c = np.trace(cm).astype(float)

    denom = np.sqrt((s**2 - (p_sum**2).sum()) * (s**2 - (t_sum**2).sum()))
    mcc = ((c * s) - (t_sum * p_sum).sum()) / denom if denom else np.nan

    return acc, mcc, precision, recall, f1


def compute_locus_metrics(true_anc, inferred_anc, labels, method, outfile):
    t_hard = to_hard_calls(true_anc)
    p_hard = to_hard_calls(inferred_anc)

    cm = confusion_counts(t_hard, p_hard, len(labels))
    acc, mcc, precision, recall, f1 = metrics_from_confusion(cm)

    per_ancestry_accuracy = {
        labels[i]: float(cm[i, i] / cm[i].sum()) if cm[i].sum() > 0 else np.nan
        for i in range(len(labels))
    }

    metrics = {
        "method": method,
        "labels": labels,
        "overall_accuracy": float(acc),
        "confusion_matrix": cm.tolist(),
        "precision": {i: float(v) for i, v in enumerate(precision)},
        "recall": {i: float(v) for i, v in enumerate(recall)},
        "f1": {i: float(v) for i, v in enumerate(f1)},
        "mcc": float(mcc),
        "shape": list(map(int, true_anc.shape)),
        "per_ancestry_accuracy": per_ancestry_accuracy,
    }
    with open(outfile, "w") as f:
        json.dump(metrics, f, indent=4)


def compute_breakpoint_distances_linear(anc, positions):
    switches = np.where(anc[1:] != anc[:-1])[0] + 1
    if len(switches) == 0:
        return np.full(len(positions), np.inf)

    bp_pos = positions[switches]

    left  = np.full(len(positions), -np.inf)
    right = np.full(len(positions),  np.inf)

    left[switches]  = bp_pos
    right[switches] = bp_pos

    left  = np.maximum.accumulate(left)
    right = np.minimum.accumulate(right[::-1])[::-1]
    return np.minimum(np.abs(positions - left), np.abs(right - positions))


def breakpoint_error_analysis(true_anc, inferred_anc, positions, bins, outfile):
    t0 = to_hard_calls(true_anc)[:, 0]
    p0 = to_hard_calls(inferred_anc)[:, 0]

    distances = compute_breakpoint_distances_linear(t0, positions)
    bin_ids = np.digitize(distances, bins) - 1
    n_bins = len(bins) - 1

    err = (t0 != p0)
    err_sum = np.array([(err[bin_ids == b]).sum() for b in range(n_bins)], dtype=np.int64)
    n_sum   = np.array([(bin_ids == b).sum() for b in range(n_bins)], dtype=np.int64)

    error_rate = np.divide(err_sum, n_sum, out=np.full(n_bins, np.nan), where=n_sum>0)

    pd.DataFrame({
        "distance_bin": [f"{bins[i]}-{bins[i+1]}" for i in range(n_bins)],
        "error_rate": error_rate,
        "n_sites": n_sum,
    }).to_csv(outfile, sep="\t", index=False)


def segment_metrics(true_anc, inferred_anc, positions, outfile):
    t0 = to_hard_calls(true_anc)[:, 0]
    p0 = to_hard_calls(inferred_anc)[:, 0]

    segments, start = [], 0
    for i in range(1, len(t0)):
        if t0[i] != t0[i - 1]:
            segments.append((start, i))
            start = i

    segments.append((start, len(t0)))

    seg_acc, seg_len = [], []
    for s, e in segments:
        seg_len.append(positions[e - 1] - positions[s])
        acc = (p0[s:e] == t0[s:e]).mean()
        seg_acc.append(acc)

    metrics = {
        "mean_segment_accuracy": float(np.mean(seg_acc)),
        "median_segment_accuracy": float(np.median(seg_acc)),
        "mean_tract_length_bp": float(np.mean(seg_len)),
        "median_tract_length_bp": float(np.median(seg_len)),
        "n_segments": len(segments),
    }
    with open(outfile, "w") as f:
        json.dump(metrics, f, indent=4)


def switch_error_rate(true_anc, inferred_anc, outfile):
    t0 = to_hard_calls(true_anc)[:, 0]
    p0 = to_hard_calls(inferred_anc)[:, 0]

    true_switches = set(np.where(t0[1:] != t0[:-1])[0] + 1)
    pred_switches = set(np.where(p0[1:] != p0[:-1])[0] + 1)

    if len(pred_switches) == 0:
        ser = np.nan
    else:
        false_switches = pred_switches - true_switches
        ser = len(false_switches) / len(pred_switches)

    metrics = {
        "switch_error_rate": float(ser),
        "n_true_switches": len(true_switches),
        "n_predicted_switches": len(pred_switches),
    }
    with open(outfile, "w") as f:
        json.dump(metrics, f, indent=4)


def main():
    configure_logging()
    args = parse_parameters()

    pop_dir = args.output / f"{args.population}_populations"
    pop_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Loading PLINK variants...")
    plink_path = here("input/simulations", pop_dir.name, "_m/plink-files/simulated")
    variants_df = load_variants(plink_path)
    variants_df = variants_df[(variants_df["chrom"] == "chr21")]
    
    logging.info("Loading ground truth data")
    # Load ground truth data
    loci_gt, g_anc_gt, admix_gt = read_simu(here(args.simu_input), chrom=21)
    labels   = list(g_anc_gt.drop(["sample_id", "chrom"], axis=1).columns)
    loci_gt  = standardize_variant_columns(loci_gt)
    dup_mask = ~loci_gt.duplicated(subset=["chrom", "pos"])
    loci_gt  = loci_gt.loc[dup_mask].reset_index(drop=True)
    loci_gt  = loci_gt.reset_index(names="_gt_pos")
    admix_gt = admix_gt[dup_mask.to_numpy()]

    logging.info("Reading RFMix outputs...")
    binary_path = args.rfmix_input / "binary_files"
    loci_df, _, admix = read_rfmix(here(args.rfmix_input),
                                   binary_dir=here(binary_path), chrom=21)
    loci_df = standardize_variant_columns(loci_df)

    for method in ["linear", "stepwise", "nearest"]:
        logging.info(f"Processing method: {method}")
        method_path =  pop_dir / method
        method_path.mkdir(parents=True, exist_ok=True)

        zarr_path = method_path / "imputed_local_ancestry" / "local-ancestry.zarr"
        variant_gt, admix_imp = impute_data(
            loci_df, variants_df, admix, zarr_path, method
        )

        variant_gt = variant_gt.reset_index(drop=True)
        variant_gt["_zarr_pos"] = np.arange(len(variant_gt))

        loci_aligned = loci_gt.merge(
            variant_gt[["chrom", "pos", "_zarr_pos"]],
            on=["chrom", "pos"], how="inner", sort=False
        )
        zarr_idx = loci_aligned["_zarr_pos"].to_numpy(dtype=np.int64)
        gt_idx   = loci_aligned["_gt_pos"].to_numpy(dtype=np.int64)

        # Lazy reads
        inferred = admix_imp[zarr_idx, :, :]
        true_anc = admix_gt[gt_idx, :, :].compute()
        positions = loci_aligned["pos"].to_numpy()

        # Calculate metrics
        logging.info(f"[{method.upper()}] Computing locus-level metrics")
        compute_locus_metrics(
            true_anc, inferred, labels, method, method_dir / "locus_metrics.json"
        )

        logging.info(f"[{method}] Computing breakpoint distance error")
        bins = [0, 1_000, 5_000, 10_000, 50_000, 100_000, np.inf]
        breakpoint_error_analysis(
            true_anc, inferred, positions, bins, method_dir / "breakpoint_error.tsv"
        )

        logging.info(f"[{method}] Computing segment-level metrics")
        segment_metrics(
            true_anc, inferred, positions, method_dir / "segment_metrics.json",
        )

        logging.info(f"[{method}] Computing switch error rate")
        switch_error_rate(
            true_anc, inferred, method_dir / "switch_error_rate.json",
        )

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
