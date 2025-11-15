import re
import psutil
import logging
import platform
import argparse
import random, os
import json, time
import numpy as np
import pyarrow as pa
from numpy import stack
import pyarrow.csv as pacsv
import pyarrow.compute as pc
from typing import Callable, List

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Metadata logging
def collect_metadata(parser_version, task_id, replicate, label):
    return {
        "parser": parser_version,
        "task": task_id,
        "replicate": replicate,
        "label": label,
        "hardware": platform.platform(),
        "software_versions": {
            "python": platform.python_version(),
            "pyarrow": pa.__version__,
            "psutil": psutil.__version__,
        }
    }


# Simulate file discovery (replace this with your real get_prefixes function)
def get_prefixes(prefix_path, verbose=True):
    return [
        {"rfmix.Q": f"{prefix_path}/chr1.Q", "fb.tsv": f"{prefix_path}/chr1.fb.tsv"},
        {"rfmix.Q": f"{prefix_path}/chr22.Q", "fb.tsv": f"{prefix_path}/chr22.fb.tsv"},
    ]

# PyArrow file reader
def _read_arrow_table(fn: str) -> pa.Table:
    return pacsv.read_csv(
        fn,
        read_options=pacsv.ReadOptions(use_threads=True),
        parse_options=pacsv.ParseOptions(delimiter="\t")
    )


# General-purpose file reader
def _read_file(fn_list: List[dict], read_func: Callable) -> List[pa.Table]:
    return [read_func(f) for f in fn_list]


def concat_tables(tables: List[pa.Table]) -> pa.Table:
    return pa.concat_tables(tables, promote=True)

# Read Q matrix with chromosome info
def _read_Q(fn_dict):
    fn = fn_dict["rfmix.Q"]
    table = _read_arrow_table(fn)
    chrom_match = re.search(r'chr(\d+)', fn)
    chrom = chrom_match.group(0) if chrom_match else None

    # Add chrom column to arrow Table
    chrom_col = pa.array([chrom] * table.num_rows)
    return table.append_column("chrom", chrom_col)

# Read fb.tsv file
def _read_fb(fn_dict):
    fn = fn_dict["fb.tsv"]
    return _read_arrow_table(fn)


# Subset population matrix
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
    return stack(pop_subset, 1)

# 🧠 Updated simulate_analysis
def simulate_analysis(prefix_path: str):
    fn_list = get_prefixes(prefix_path)

    # Read and concatenate Q files
    g_anc_tables = _read_file(fn_list, _read_Q)
    g_anc = concat_tables(g_anc_tables)

    # Infer population labels from column names (skip 'sample_id' and 'chrom')
    col_names = g_anc.column_names
    pops = [name for name in col_names if name not in ("sample_id", "chrom")]

    # Read and concatenate fb.tsv files
    X_tables = _read_file(fn_list, _read_fb)
    X = concat_tables(X_tables)

    # Extract loci and ancestry matrix
    loci = X.select(["chromosome", "physical_position"])
    ancestry_table = X.drop(["chromosome", "physical_position", "sample_id", "marker"])

    # Convert Arrow Table to NumPy array for numerical operation
    ancestry_matrix = np.column_stack([ancestry_table[column].to_numpy() for column in ancestry_table.column_names])

    # Subset
    admix = _subset_populations(ancestry_matrix, len(pops))

    # Optional: Compute Arrow summary statistics (or convert to Pandas only if needed)
    # Example using Arrow: mean
    g_anc_means = {name: pc.mean(g_anc[name]).as_py() for name in pops}
    admix_means = np.mean(admix, axis=0)

    return g_anc_means, admix_means


# Simulated file selection based on task
def select_file(input_dir: str, task: int):
    if task == 1:
        return os.path.join(input_dir, "chr1.tsv")
    elif task == 2:
        return os.path.join(input_dir, "chr22.tsv")
    else:
        return os.path.join(input_dir, "all_chr.tsv")

# Track peak CPU memory in MB
def get_peak_cpu_memory():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss / (1024 ** 2)  # RSS: Resident Set Size in MB


# Main processing function
def run_task(input_dir: str, label: str, task: int):
    output_dir = os.path.join("output", label)
    os.makedirs(output_dir, exist_ok=True)

    for replicate in range(3, 6):  # Replicates 3 to 5
        seed = replicate
        random.seed(seed)
        np.random.seed(seed)

        file_path = select_file(input_dir, task)
        logging.info(f"Replicate {replicate}: Reading file {file_path}")

        start_mem = get_peak_cpu_memory()
        start = time.time()
        g_anc_means, _ = simulate_analysis(file_path)
        wall_time = time.time() - start
        end_mem = get_peak_cpu_memory()
        peak_mem = max(start_mem, end_mem)

        # Save result
        result_path = os.path.join(output_dir, f"result_replicate_{replicate}.csv")
        table = pa.table([g_anc_means])
        pacsv.write_csv(table, result_path)

        # Save metadata
        meta = collect_metadata("pyarrow", task, replicate, label)
        meta["wall_time_sec"] = wall_time
        meta["peak_cpu_memory_MB"] = peak_mem
        meta_path = os.path.join(output_dir, f"meta_replicate_{replicate}.json")
        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=2)

        logging.info(f"Finished replicate {replicate} in {wall_time:.2f}s, CPU Peak: {peak_mem:.2f} MB")

# CLI entrypoint
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Large File PyArrow Processor")
    parser.add_argument("--input", type=str, required=True, help="Input directory with data files")
    parser.add_argument("--label", type=str, required=True, help="Label for output directory")
    parser.add_argument("--task", type=int, choices=[1, 2, 3], required=True,
                        help="Task type: 1=small chrom, 2=large chrom, 3=all chroms")
    args = parser.parse_args()

    run_task(args.input, args.label, args.task)
