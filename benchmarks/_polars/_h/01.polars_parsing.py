import os, re
import logging
import platform
import argparse
import numpy as np
import random, json
import psutil, time
import polars as pl
from typing import Callable, List
from collections import OrderedDict as odict

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
            "polars": pl.__version__,
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


def _types(fn: str, Q: bool = False, n_rows: int = 100) -> dict[str, pl.DataType]:
    """Infer dtypes using Polars sample read."""
    sample = pl.read_csv(
        fn, separator=r"\s+", n_rows=n_rows, comment_prefix="#", has_header=False
    )

    dtypes = sample.dtypes_dict()

    if Q:
        # First column = categorical
        ((first_name, first_type), *rest) = dtypes.items()
        return {"sample_id": pl.Categorical, **dict(rest)}

    return dtypes


def _read_csv(fn: str, dtypes=None) -> pl.LazyFrame:
    return pl.scan_csv(
        fn, separator=r"\s+", comment_prefix="#", has_header=False, dtypes=dtypes
    )


def _read_file(fn_list: List[dict], read_func: Callable) -> List[pl.DataFrame]:
    return [read_func(f) for f in fn_list]


def concat_tables(tables: List[pl.LazyFrame]) -> pl.DataFrame:
    return pl.concat(tables, how="vertical_relaxed")


def _read_Q(fn_dict):
    fn = fn_dict["rfmix.Q"]
    header = _types(fn, Q=True)
    df = _read_csv(fn, dtypes=header)

    chrom_match = re.search(r'chr(\d+)', fn)
    chrom = chrom_match.group(0) if chrom_match else None

    return df.with_columns([pl.lit(chrom).alias("chrom")])


def _read_fb(fn_dict):
    fn = fn_dict["fb.tsv"]
    header = _types(fn, Q=False)
    return _read_csv(fn, dtypes=header)


def _subset_populations(df: pl.LazyFrame, npops: int) -> np.ndarray:
    col_names = df.schema.keys()
    ncols = len(col_names)

    if ncols % npops != 0:
        raise ValueError("The number of columns in X must be divisible by npops.")

    loci_per_pop = ncols // npops
    if loci_per_pop % 2 != 0:
        raise ValueError("Number of columns per population must be even.")

    loci = loci_per_pop // 2
    pop_frames = []
    for pop_idx in range(npops):
        pop_cols = col_names[pop_idx::npops]
        pop_sum_exprs = [
            (pl.col(pop_cols[i]) + pl.col(pop_cols[i + 1])).alias(f"pop{pop_idx}_locus{j}")
            for j, i in enumerate(range(0, len(pop_cols), 2))
        ]
        # Lazy select; defer collection until end
        pop_lazy = df.select(pop_sum_exprs)
        pop_frames.append(pop_lazy)

    # Collect all at once
    pop_dfs = [lazy.collect() for lazy in pop_frames]
    pop_arrays = [df.to_numpy() for df in pop_dfs]

    return np.stack(pop_arrays, axis=2)


def simulate_analysis(input_dir: str, task: int):
    fn_list = get_prefixes(input_dir, task)

    g_anc_tables = _read_file(fn_list, _read_Q)
    g_anc = concat_tables(g_anc_tables).collect()

    pops = [col for col in g_anc.columns if col not in ("sample_id", "chrom")]

    X_tables = _read_file(fn_list, _read_fb)
    X = concat_tables(X_tables)

    loci = X.select(["chromosome", "physical_position"])

    ancestry_df = X.select(pl.col(pl.NUMERIC_DTYPES))

    admix = _subset_populations(ancestry_df, len(pops))

    return loci, g_anc, admix


def get_peak_cpu_memory_mb() -> float:
    import resource
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_maxrss / 1024.0


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

        # Numeric columns only
        g_anc_means = {
            col: g_anc[col].mean() for col in g_anc.columns
            if g_anc[col].dtype in pl.NUMERIC_DTYPES
        }

        result_path = os.path.join(output_dir, f"result_replicate_{replicate}.csv")
        pl.DataFrame(g_anc_means).write_csv(result_path)

        meta = collect_metadata("polars", task, replicate, label)
        meta["wall_time_sec"] = wall_time
        meta["peak_cpu_memory_MB"] = peak_mem

        meta_path = os.path.join(output_dir, f"meta_replicate_{replicate}.json")
        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=2)

        logging.info(
            "Finished replicate %d in %.2fs, peak RSS: %.2f MB",
            replicate, wall_time, peak_mem,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Large File Polars Processor")
    parser.add_argument("--input", type=str, required=True,
                        help="Input directory with data files")
    parser.add_argument("--label", type=str, required=True,
                        help="Label for output directory")
    parser.add_argument("--task", type=int, choices=[1, 2, 3], required=True,
                        help="Task type: 1=small chrom, 2=large chrom, 3=all chroms")
    args = parser.parse_args()

    run_task(args.input, args.label, args.task)
