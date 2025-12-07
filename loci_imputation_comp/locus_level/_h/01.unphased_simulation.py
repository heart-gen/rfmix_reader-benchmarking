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
    parser.add_argument("--rfmix-input", type=Path, required=True,
                        help="Path to rfmix input data prefix.")
    parser.add_argument("--simu-input", type=Path, required=True,
                        help="Path to simulation input data prefix.")
    parser.add_argument("--output", type=str, required=True, help="Output directory")
    parser.add_argument("--population", type=str, choices=["two", "three"],
                        default="three", help="Number of populations for admixture.")
    parser.add_argument("--method", type=str, choices=["linear", "nearest", "stepwise"],
                        default="linear", help="Imputation method")
    return parser.parse_args()


def load_variants(plink_path):
    pgr = PgenReader(plink_path)
    variant_df = pgr.variant_df
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    mask = ~variant_df.duplicated(subset=["chrom", "pos"], keep="first")
    return variant_df[mask]
    

def impute_data(loci_df, variants, admix, zarr_path, method):
    variant_df = variants.merge(loci_df.to_pandas(), on=["chrom", "pos"],
                                how="outer", indicator=True)\
                         .loc[:, ["chrom", "pos", "i", "_merge"]]
    z = interpolate_array(variant_df, admix, zarr_outdir=zarr_path,
                          interpolation=method, use_bp_positions=True,
                          chunk_size=50_000)
    return z


def compute_locus_accuracy(true_anc, inferred_anc):
    """Compute locus-level accuracy."""
    total = true_anc.size
    correct = np.sum(true_anc == inferred_anc)
    return correct / total


def compute_locus_matrices(true_anc, inferred_anc):
    """Compute confusion matrix for ancestry imputation."""
    true_anc = true_anc.ravel()
    inferred_anc = inferred_anc.ravel()
    
    labels = np.unique(true_anc)
    cm = confusion_matrix(true_anc, inferred_anc, labels=labels)
    return cm


def compute_per_ancestry_precision(true_anc, inferred_anc):
    true_anc = true_anc.ravel()
    inferred_anc = inferred_anc.ravel()
    
    labels = np.unique(true_anc)
    precisions = precision_score(
        true_anc, inferred_anc, labels=labels, average=None, zero_division=0
    )
    return {label: p for label, p in zip(labels, percisions)}


def compute_per_ancestry_recall(true_anc, inferred_anc):
    """Compute recall for each ancestry."""
    true_anc = true_anc.ravel()
    inferred_anc = inferred_anc.ravel()
    
    labels = np.unique(true_anc)

    rec = recall_score(
        true_anc, inferred_anc, labels=labels, average=None, zero_division=0
    )
    return {label: r for label, r in zip(labels, rec)}


def compute_per_ancestry_f1(true_anc, inferred_anc):
    """Compute F1 score for each ancestry."""
    true_anc = true_anc.ravel()
    inferred_anc = inferred_anc.ravel()

    labels = np.unique(true_anc)

    f1 = f1_score(
        true_anc, inferred_anc, labels=labels, average=None, zero_division=0
    )
    return {label: f for label, f in zip(labels, f1)}


def compute_mcc(true_anc, inferred_anc):
    """Compute global MCC for all ancestries."""
    true_anc = true_anc.ravel()
    inferred_anc = inferred_anc.ravel()

    return matthews_corrcoef(true_anc, inferred_anc)


def compute_per_ancestry_accuracy(true_anc, inferred_anc):
    """
    Compute ancestry-specific accuracy: P(predicted == true | true == ancestry).
    """
    true_anc = true_anc.ravel()
    inferred_anc = inferred_anc.ravel()

    labels = np.unique(true_anc)

    results = {}
    for lab in labels:
        mask = (true_anc == lab)
        if mask.sum() == 0:
            results[lab] = np.nan
        else:
            results[lab] = np.mean(inferred_anc[mask] == lab)
    return results


def compute_locus_metrics_json(t, p, method, outfile=None):
    """
    Compute all locus-level metrics and return a JSON-serializable dict.
    """
    # Flatten
    t_flat = t.ravel()
    p_flat = p.ravel()

    labels = list(np.unique(t_flat))

    # Compute metrics
    accuracy = compute_locus_accuracy(t, p)
    cm = compute_confusion_matrices(t, p, labels=labels).tolist()
    precision = compute_per_ancestry_precision(t, p, labels=labels)
    recall = compute_per_ancestry_recall(t, p, labels=labels)
    f1 = compute_per_ancestry_f1(t, p, labels=labels)
    mcc_global = compute_mcc(t, p)
    mcc_per = compute_per_ancestry_mcc(t, p, labels=labels)
    acc_pop = compute_per_ancestry_accuracy(t, p, labels=labels)

    # Bundle into JSON-serializable dict
    metrics = {
        "method": method, "labels": labels,
        "overall_accuracy": float(accuracy),
        "confusion_matrix": cm,
        "precision_per_ancestry": {int(k): float(v) for k, v in precision.items()},
        "recall_per_ancestry": {int(k): float(v) for k, v in recall.items()},
        "f1_per_ancestry": {int(k): float(v) for k, v in f1.items()},
        "accuracy_per_ancestry": {int(k): float(v) for k, v in acc_pop.items()},
        "mcc_global": float(mcc_global),
        "mcc_per_ancestry": {int(k): float(v) for k, v in mcc_per.items()},
        "n_variants": int(t.shape[0]),
        "n_samples": int(t.shape[1]),
        "shape": list(map(int, t.shape))
    }

    if outfile is not None:
        with open(outfile, "w") as f:
            json.dump(metrics, f, indent=4)

    return metrics


def main():
    # Parameters
    configure_logging()
    args = parse_parameters()
    pop_levels = f"{args.population}_populations"

    outdir = args.outdir
    outdir.mkdir(exist_ok=True, parents=True)

    # Load ground truth data
    loci_gt, _, admix_gt = read_simu(str(args.simu_input))
    true_anc = admix_gt.compute()

    # Load genotypes and variants
    plink_path = here("input/simulations", pop_levels, "_m/plink-files/simulated")
    variants_df = load_variants(plink_path)

    # Load ancestry
    loci_df, _, admix = read_rfmix(str(args.rfmix_input))
    method_path = args.rfmix_input / args.method
    method_path.mkdir(exist_ok=True, parents=True)
    
    zarr_path = method_path / "imputed_local_ancestry"
    z = impute_data(loci_df, variants_df, admix, zarr_path, args.method)

    # Align ground truth and inferred ancestry
    ## TODO
    inferred_anc = z
    
    # Calculate accuracy
    outfile = outdir / "some-name.json"
    compute_locus_metrics_json(true_anc, inferred_anc, args.method, outfile=outfile)
    
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
    
