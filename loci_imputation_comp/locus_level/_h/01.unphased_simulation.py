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
    parser = argparse.ArgumentParser(description="Locus-Level Imputation Accuracy")
    parser.add_argument("--rfmix-input", type=Path, required=True)
    parser.add_argument("--simu-input", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("./"))
    parser.add_argument("--population", type=str, choices=["two","three"], default="three")
    parser.add_argument("--method", type=str, choices=["linear", "nearest", "stepwise"],
                        default="linear")
    return parser.parse_args()


def standardize_variant_columns(df):
    return df.rename(columns={"chromosome": "chrom", "physical_position": "pos"})


def align_variants(df_ref, df_target):
    ref_idx = df_ref.set_index(["chrom", "pos"])
    target_idx = df_target.set_index(["chrom", "pos"])
    common = target_idx.index.intersection(ref_idx.index)
    return df_target.loc[common].reset_index()


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


def main():
    configure_logging()
    args = parse_parameters()

    pop_dir = args.output / f"{args.population}_populations"
    pop_dir.mkdir(parents=True, exist_ok=True)

    n_labels = 2 if args.population == "two" else 3

    # Load inferred ancestry (lazy)
    logging.info("Loading imputed variants and ancestry data")
    method_path  = here(args.rfmix_input) / args.method
    zarr_path    = method_path / "imputed_local_ancestry" / "local-ancestry.zarr"

    variant_df = pd.read_parquet(method_path / "imputed_variant.parquet")
    inferred   = da.from_zarr(zarr_path)
    
    logging.info("Loading ground truth data")
    # Load ground truth data
    loci_gt, _, admix_gt = read_simu(here(args.simu_input))
    loci_gt  = standardize_variant_columns(loci_gt)
    loci_gt  = align_variants(variant_df, loci_gt)

    logging.info("Align variants")
    idx = loci_gt.index.to_numpy()
    inferred = da.take(inferred, idx, axis=0)
    true_anc = admix_gt[idx, :, :]

    # Calculate metrics
    logging.info("Computing locus-level metrics")
    compute_locus_metrics(
        true_anc,
        inferred,
        n_labels,
        args.method,
        pop_dir / "locus_metrics.json"
    )

    logging.info("Computing per-haplotype metrics")
    per_haplotype_metrics(
        true_anc,
        inferred,
        n_labels,
        pop_dir / "per_haplotype_metrics.tsv"
    )

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
