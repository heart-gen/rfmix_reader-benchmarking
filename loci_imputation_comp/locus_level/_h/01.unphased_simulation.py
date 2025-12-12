import json
import logging
import argparse
import numpy as np
import pandas as pd
import session_info
from pyhere import here
import dask.array as da
from pathlib import Path

from rfmix_reader import read_simu

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


def discover_methods(rfmix_root):
    """Find valid method subdirectories under rfmix_root."""
    methods = [
        p.name for p in rfmix_root.iterdir()
        if p.is_dir() and (p / "imputed_local_ancestry").exists()
    ]
    return sorted(methods)


def standardize_variant_columns(df):
    return df.rename(columns={"chromosome": "chrom", "physical_position": "pos"})


def align_variants(variant_df, loci_gt):
    return loci_gt.merge(
        variant_df[["chrom", "pos", "_zarr_pos"]], on=["chrom", "pos"],
        how="inner", sort=False
    )


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


def _to_hard_calls(arr):
    """Convert (V,S,A) -> (V, 2S) lazily"""
    if arr.ndim == 3:
        return da.map_blocks(
            admix_counts_to_haplotypes_numpy,
            arr,
            dtype=np.int8,
            chunks=(arr.chunks[0], (arr.shape[1] * 2,)),
            drop_axis=2,
        )
    return arr


def _confmat_block(t, p, n_labels):
    t = t.ravel()
    p = p.ravel()
    mask = (t >= 0) & (p >= 0)

    if mask.sum() == 0:
        return np.zeros(n_labels * n_labels, dtype=np.int64)

    idx = t[mask].astype(np.int64) * n_labels + p[mask].astype(np.int64)
    return np.bincount(idx, minlength=n_labels * n_labels).astype(np.int64)


def confusion_counts_dask(t_hard, p_hard, n_labels):
    counts = da.map_blocks(
        _confmat_block,
        t_hard,
        p_hard,
        dtype=np.int64,
        chunks=(n_labels * n_labels,),
        drop_axis=list(range(t_hard.ndim)),
        new_axis=0,
        n_labels=n_labels,
    )
    return counts.sum(axis=0)


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


def compute_locus_metrics(true_anc, inferred_anc, n_labels, method, outfile):
    t_hard = to_hard_calls(true_anc)
    p_hard = to_hard_calls(inferred_anc)

    flat = confusion_counts_dask(t_hard, p_hard, n_labels).compute()
    cm = flat.reshape((n_labels, n_labels))

    acc, mcc, precision, recall, f1 = metrics_from_confusion(cm)
    
    per_ancestry_accuracy = {
        i: float(cm[i, i] / cm[i].sum()) if cm[i].sum() > 0 else np.nan
        for i in range(n_labels)
    }

    metrics = {
        "method": method,
        "labels": list(range(n_labels)),
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

    return metrics


def per_haplotype_metrics(true_anc, inferred_anc, n_labels, outfile):
    t_hard = to_hard_calls(true_anc)
    p_hard = to_hard_calls(inferred_anc)

    H = t_hard.shape[1]
    rows = []

    for h in range(H):
        flat = confusion_counts_dask(
            t_hard[:, h:h+1],
            p_hard[:, h:h+1],
            n_labels
        ).compute()

        cm = flat.reshape((n_labels, n_labels))
        acc, mcc, precision, recall, f1 = metrics_from_confusion(cm)

        rows.append({
            "haplotype": h,
            "accuracy": acc,
            "mcc": mcc,
            "precision": ",".join(map(str, precision)),
            "recall": ",".join(map(str, recall)),
            "f1": ",".join(map(str, f1)),
            "confusion_matrix": ";".join(",".join(map(str, r)) for r in cm),
        })

    pd.DataFrame(rows).to_csv(outfile, sep="\t", index=False)


def compute_breakpoint_distances(true_hard_hap0, positions):
    anc = true_hard_hap0.compute() if hasattr(true_hard_hap0, "compute") else true_hard_hap0
    switches = np.where(anc[1:] != anc[:-1])[0] + 1

    if len(switches) == 0:
        return np.full(len(positions), np.inf)

    bp_pos = positions[switches]
    return np.min(np.abs(positions[:, None] - bp_pos[None, :]), axis=1)


def breakpoint_error_analysis(true_anc, inferred_anc, positions, bins, outfile):
    t_hard = to_hard_calls(true_anc)
    p_hard = to_hard_calls(inferred_anc)

    # One haplotype is sufficient for breakpoint
    t0 = t_hard[:, 0]
    p0 = p_hard[:, 0]

    distances = compute_breakpoint_distances(t0, positions)
    bin_ids = np.digitize(distances, bins) - 1
    n_bins = len(bins) - 1

    def _block_error(t_block, p_block, bin_block):
        err = (t_block != p_block)
        out_err = np.zeros(n_bins, dtype=np.int64)
        out_n   = np.zeros(n_bins, dtype=np.int64)

        for b in range(n_bins):
            m = bin_block == b
            out_err[b] += err[m].sum()
            out_n[b]   += m.sum()

        return np.stack([out_err, out_n], axis=0)

    stats = da.map_blocks(
        _block_error, t0, p0, da.from_array(bin_ids, chunks=t0.chunks),
        dtype=np.int64, chunks=(2, n_bins), new_axis=0,
    ).sum(axis=0)

    err_sum, n_sum = stats.compute()
    error_rate = err_sum / n_sum

    df = pd.DataFrame({
        "distance_bin": [f"{bins[i]}-{bins[i+1]}" for i in range(n_bins)],
        "error_rate": error_rate,
        "n_sites": n_sum,
    })

    df.to_csv(outfile, sep="\t", index=False)


def segment_metrics(true_anc, inferred_anc, positions, outfile):
    t_hard = to_hard_calls(true_anc)
    p_hard = to_hard_calls(inferred_anc)

    # Use hap0 only for segmentation
    t0 = t_hard[:, 0].compute()
    p0 = p_hard[:, 0]

    segments, start = [], 0
    for i in range(1, len(t0)):
        if t0[i] != t0[i - 1]:
            segments.append((start, i))
            start = i
    segments.append((start, len(t0)))

    seg_acc = seg_len = []
    for s, e in segments:
        seg_len.append(positions[e - 1] - positions[s])
        acc = (p0[s:e].compute() == t0[s:e]).mean()
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
    t_hard = to_hard_calls(true_anc)
    p_hard = to_hard_calls(inferred_anc)

    t0 = t_hard[:, 0].compute()
    p0 = p_hard[:, 0].compute()

    true_switches = set(np.where(t0[1:] != t0[:-1])[0] + 1)
    pred_switches = set(np.where(p0[1:] != p0[:-1])[0] + 1)

    if len(pred_switches) == 0:
        ser = 0.0
    else:
        false_switches = pred_switches - true_switches
        ser = len(false_switches) / len(pred_switches)

    with open(outfile, "w") as f:
        json.dump(
            {
                "switch_error_rate": float(ser),
                "n_true_switches": len(true_switches),
                "n_predicted_switches": len(pred_switches),
            }, f, indent=4,
        )


def main():
    configure_logging()
    args = parse_parameters()

    n_labels = 2 if args.population == "two" else 3
    pop_dir = args.output / f"{args.population}_populations"
    pop_dir.mkdir(parents=True, exist_ok=True)


    methods = discover_methods(here(args.rfmix_input))
    logging.info(f"Evaluating methods: {methods}")
    
    logging.info("Loading ground truth data")
    # Load ground truth data
    loci_gt, _, admix_gt = read_simu(here(args.simu_input))
    loci_gt  = standardize_variant_columns(loci_gt)
    dup_mask = ~loci_gt.duplicated(subset=["chrom", "pos"])
    loci_gt  = loci_gt.loc[dup_mask].reset_index(drop=True)
    admix_gt = admix_gt[dup_mask.to_numpy()]

    for method in methods:
        logging.info(f"Processing method: {method}")
        method_dir = pop_dir / method
        method_dir.mkdir(exist_ok=True)

        method_path  = here(args.rfmix_input) / method
        zarr_path    = method_path / "imputed_local_ancestry" / "local-ancestry.zarr"

        variant_df = pd.read_parquet(method_path / "imputed_variant.parquet")
        variant_df = variant_df.reset_index(drop=True)
        variant_df["_zarr_pos"] = np.arange(len(variant_df))

        loci_aligned = align_variants(variant_df, loci_gt)
        idx = loci_aligned["_zarr_pos"].to_numpy(dtype=np.int64)

        # Lazy reads
        inferred = da.take(da.from_zarr(zarr_path), idx, axis=0)
        true_anc = admix_gt[idx, :, :]

        # Calculate metrics
        logging.info(f"[{method}] Computing locus-level metrics")
        compute_locus_metrics(
            true_anc, inferred, n_labels, method, method_dir / "locus_metrics.json"
        )

        logging.info(f"[{method}] Computing per-haplotype metrics")
        per_haplotype_metrics(
            true_anc, inferred, n_labels, method_dir / "per_haplotype_metrics.tsv"
        )

        logging.info(f"[{method}] Computing breakpoint distance error")
        bins = [0, 1_000, 5_000, 10_000, 50_000, 100_000, np.inf]
        breakpoint_error_analysis(
            true_anc, inferred, variant_df["pos"].values[idx], bins,
            method_dir / "breakpoint_error.tsv"
        )

        logging.info(f"[{method}] Computing segment-level metrics")
        segment_metrics(
            true_anc, inferred, variant_df["pos"].values[idx],
            method_dir / "segment_metrics.json",
        )

        logging.info(f"[{method}] Computing switch error rate")
        switch_error_rate(
            true_anc, inferred, method_dir / "switch_error_rate.json",
        )

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
