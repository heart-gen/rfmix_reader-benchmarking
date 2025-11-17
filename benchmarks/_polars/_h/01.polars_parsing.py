import os, re
import logging
import platform
import argparse
import numpy as np
import random, json
import psutil, time
import polars as pl
import polars.selectors as cs
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
            "rfmix.Q": os.path.join(input_dir, f"chr{chrom}.rfmix.Q"),
            "fb.tsv":  os.path.join(input_dir, f"chr{chrom}.fb.tsv"),
        })

    if verbose:
        logging.info("Task %d -> %d chromosomes: %s", task, len(chroms), chroms)

    return prefixes


def _types(fn: str) -> list[str]:
    """Extract the Q header from a .rfmix.Q file."""
    header_line: str | None = None
    hash_lines_seen = 0
    with open(fn, "r") as f:
        for line in f:
            if not line.startswith("#"):
                continue

            hash_lines_seen += 1
            if hash_lines_seen == 1:
                continue

            header_line = line[1:].strip()
            break

    if header_line is None:
        raise ValueError(f"No suitable header line starting with '#' found in {fn}")

    if "\t" in header_line:
        header = header_line.split("\t")
    else:
        header = header_line.split()

    if header:
        header[0] = "sample_id"
    return header


def _read_csv(fn: str, header: list[str]=None, Q: bool=True) -> pl.LazyFrame:
    return pl.scan_csv(
        fn, separator="\t", comment_prefix="#", has_header=not Q, 
        new_columns=header if Q else None,
    )


def _read_file(fn_list: List[dict], read_func: Callable) -> List[pl.LazyFrame]:
    return [read_func(f) for f in fn_list]


def concat_tables(tables: List[pl.LazyFrame]) -> pl.LazyFrame:
    return pl.concat(tables, how="vertical_relaxed")


def _read_Q(fn_dict):
    fn = fn_dict["rfmix.Q"]
    header = _types(fn)
    lf = _read_csv(fn, header=header, Q=True)
    chrom_match = re.search(r'chr(\d+)', fn)
    chrom = chrom_match.group(0) if chrom_match else None
    return lf.with_columns([pl.lit(chrom).alias("chrom")])


def _read_fb(fn_dict):
    fn = fn_dict["fb.tsv"]
    return _read_csv(fn, Q=False)


def _subset_populations(df: pl.LazyFrame, npops: int) -> np.ndarray:
    col_names = list(df.collect_schema())
    ncols = len(col_names)

    if ncols % npops != 0:
        raise ValueError("The number of columns in X must be divisible by npops.")

    loci_per_pop = ncols // npops
    if loci_per_pop % 2 != 0:
        raise ValueError("Number of columns per population must be even.")

    pop_frames = []
    for pop_idx in range(npops):
        pop_cols = col_names[pop_idx::npops]
        pop_sum_exprs = [
            (pl.col(pop_cols[i]) + pl.col(pop_cols[i + 1])).alias(f"pop{pop_idx}_locus{j}")
            for j, i in enumerate(range(0, len(pop_cols), 2))
        ]
        pop_lazy = df.select(pop_sum_exprs)
        pop_frames.append(pop_lazy)

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

    loci = X.select(["chromosome", "physical_position"]).collect()
    meta_cols = ["chromosome", "physical_position", "genetic_position",
                 "genetic_marker_index"]
    meta_in = [c for c in meta_cols if c in X.collect_schema().names()]

    ancestry_lf = X.drop(meta_in).select(cs.numeric())
    admix = _subset_populations(ancestry_lf, len(pops))

    return loci, g_anc, admix


def get_peak_cpu_memory_mb() -> float:
    import resource
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_maxrss / 1024.0


def run_task(input_dir: str, label: str, task: int):
    output_dir = os.path.join("output", label)
    os.makedirs(output_dir, exist_ok=True)

    for replicate in range(1, 6):
        seed = replicate + 13
        random.seed(seed); np.random.seed(seed)

        logging.info(f"Replicate {replicate}: Reading files from {input_dir}")

        start = time.time()
        _, _, _ = simulate_analysis(input_dir, task)
        wall_time = time.time() - start
        peak_mem = get_peak_cpu_memory_mb()

        # Output meta data
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
