import os, re
import logging
import platform
import argparse
import numpy as np
import random, json
import psutil, time
import pandas as pd
from typing import Callable, List

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def collect_metadata(parser_version, task_id, replicate, label):
    return {
        "parser": parser_version,
        "task": task_id,
        "replicate": replicate,
        "label": label,
        "hardware": platform.platform(),
        "software_versions": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "psutil": psutil.__version__,
        }
    }


def get_prefixes(input_dir: str, task: int, verbose: bool = True) -> list[dict]:
    if task == 1:
        chroms = [21]
    elif task == 2:
        chroms = [1]
    elif task == 3:
        chroms = list(range(1, 23))
    else:
        raise ValueError(f"Unsupported task id: {task}")

    prefixes = []
    for chrom in chroms:
        prefixes.append({
            "rfmix.Q": os.path.join(input_dir, f"chr{chrom}.Q"),
            "fb.tsv":  os.path.join(input_dir, f"chr{chrom}.fb.tsv"),
        })

    if verbose:
        logging.info("Task %d -> %d chromosomes: %s", task, len(chroms), chroms)

    return prefixes


def _read_csv(fn: str) -> pd.DataFrame:
    return pd.read_csv(fn, sep="\t", engine="c")


def _read_file(fn_list: List[dict], read_func: Callable) -> List[pd.DataFrame]:
    return [read_func(f) for f in fn_list]


def concat_tables(tables: List[pd.DataFrame]) -> pd.DataFrame:
    return pd.concat(tables, ignore_index=True)


def _read_Q(fn_dict):
    fn = fn_dict["rfmix.Q"]
    df = _read_csv(fn)
    chrom_match = re.search(r'chr(\d+)', fn)
    chrom = chrom_match.group(0) if chrom_match else None
    df["chrom"] = chrom
    return df


def _read_fb(fn_dict):
    fn = fn_dict["fb.tsv"]
    return _read_csv(fn)


def table_to_numpy(df: pd.DataFrame, dtype: np.dtype = np.float32) -> np.ndarray:
    return df.to_numpy(dtype=dtype)


def _subset_populations(X: np.ndarray, npops: int) -> np.ndarray:
    pop_subset = []
    ncols = X.shape[1]
    if ncols % npops != 0:
        raise ValueError("The number of columns in X must be divisible by npops.")

    for pop_start in range(npops):
        X0 = X[:, pop_start::npops]
        if X0.shape[1] % 2 != 0:
            raise ValueError("Number of columns must be even.")
        X0_summed = X0[:, ::2] + X0[:, 1::2]
        pop_subset.append(X0_summed)
    return np.stack(pop_subset, axis=2)


def simulate_analysis(input_dir: str, task: int):
    fn_list = get_prefixes(input_dir, task)
    g_anc_tables = _read_file(fn_list, _read_Q)
    g_anc = concat_tables(g_anc_tables)

    pops = [col for col in g_anc.columns if col not in ("sample_id", "chrom")]

    X_tables = _read_file(fn_list, _read_fb)
    X = concat_tables(X_tables)

    loci = X[["chromosome", "physical_position"]]
    ancestry_df = X.drop(columns=["chromosome", "physical_position", "genetic_position", "genetic_marker_index"])

    ancestry_matrix = table_to_numpy(ancestry_df, dtype=np.float32)
    admix = _subset_populations(ancestry_matrix, len(pops))

    return loci, g_anc, admix


def get_peak_cpu_memory_mb() -> float:
    import resource
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_maxrss / 1024.0


def _dict_to_df_single_row(d: dict[str, float]) -> pd.DataFrame:
    return pd.DataFrame([d])


def run_task(input_dir: str, label: str, task: int):
    output_dir = os.path.join("output", label)
    os.makedirs(output_dir, exist_ok=True)

    for replicate in range(3, 6):
        seed = replicate
        random.seed(seed)
        np.random.seed(seed)

        logging.info(f"Replicate {replicate}: Reading files from {input_dir}")

        start = time.time()
        _, g_anc, _ = simulate_analysis(input_dir, task)
        wall_time = time.time() - start
        peak_mem = get_peak_cpu_memory_mb()

        # Compute means per column excluding non-numeric ones
        g_anc_means = g_anc.select_dtypes(include=[np.number]).mean().to_dict()

        result_path = os.path.join(output_dir, f"result_replicate_{replicate}.csv")
        df_result = _dict_to_df_single_row(g_anc_means)
        df_result.to_csv(result_path, index=False)

        meta = collect_metadata("pandas", task, replicate, label)
        meta["wall_time_sec"] = wall_time
        meta["peak_cpu_memory_MB"] = peak_mem

        meta_path = os.path.join(output_dir, f"meta_replicate_{replicate}.json")
        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=2)

        logging.info("Finished replicate %d in %.2fs, peak RSS: %.2f MB",
                     replicate, wall_time, peak_mem)
