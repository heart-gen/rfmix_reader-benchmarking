import re, io
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
            "pyarrow": pa.__version__,
            "psutil": psutil.__version__,
        }
    }


def get_prefixes(input_dir: str, task: int, verbose: bool = True) -> list[dict]:
    if task == 1:
        chroms = [21]
    elif task == 2:
        chroms = [1]
    elif task == 3:
        chroms = list(range(1, 23))  # 1..22
    else:
        raise ValueError(f"Unsupported task id: {task}")

    prefixes: list[dict[str, str]] = []
    for chrom in chroms:
        prefixes.append({
            "rfmix.Q": os.path.join(input_dir, f"chr{chrom}.rfmix.Q"),
            "fb.tsv":  os.path.join(input_dir, f"chr{chrom}.fb.tsv"),
        })

    if verbose:
        logging.info("Task %d -> %d chromosome(s): %s", task, len(chroms), chroms)

    return prefixes


def _split_header_and_data(fn: str, Q: bool = False) -> tuple[list[str], list[str]]:
    """
    Return (column_names, data_lines) for an RFMix Q or fb.tsv file.

    For Q=True:
      - Treat the *second* line starting with '#' as the header.
    For Q=False (fb.tsv):
      - Treat the first non-empty, non-comment line as the header.
    """
    with open(fn, "r") as f:
        lines = f.readlines()

    if len(lines) < 1:
        raise ValueError(f"{fn} is empty")

    header_line = header_idx = None
    if Q:
        hash_seen = 0
        for i, line in enumerate(lines):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                hash_seen += 1
                if hash_seen == 2:
                    header_line = stripped[1:]  # drop leading '#'
                    header_idx = i
                    break
        if header_line is None:
            raise ValueError(f"Could not find Q header line (second '#') in {fn}")
    else:
        for i, line in enumerate(lines):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                continue
            header_line = stripped
            header_idx = i
            break
        if header_line is None:
            raise ValueError(f"Could not find header line in {fn}")

    # Split header into column names (whitespace or tabs)
    if "\t" in header_line:
        column_names = header_line.split("\t")
    else:
        column_names = header_line.split()

    if Q and column_names: # Normalize first column name for Q
        if column_names[0] != "sample_id":
            column_names[0] = "sample_id"

    # All lines after the header are data
    data_start = header_idx + 1
    data_lines = lines[data_start:]

    return column_names, data_lines


def _types(fn: str, Q: bool = False, n_rows: int = 100) -> dict[str, pa.DataType]:
    """
    Infer column types using a small sample of data, with header removed.
    """
    column_names, data_lines = _split_header_and_data(fn, Q=Q)

    sample_lines = data_lines[:n_rows]
    buffer = io.BytesIO("".join(sample_lines).encode("utf-8"))

    table = pacsv.read_csv(
        buffer,
        read_options=pacsv.ReadOptions(
            column_names=column_names,  # no header row in buffer
            use_threads=True,
        ),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(column_types=None),  # let Arrow infer
    )

    schema = table.schema
    types: dict[str, pa.DataType] = {
        field.name: field.type for field in schema
    }

    if Q and "sample_id" in types:
        types["sample_id"] = pa.string()

    return types


def _read_arrow_table(fn: str, Q: bool = False) -> pa.Table:
    column_names, data_lines = _split_header_and_data(fn, Q=Q)
    inferred_types = _types(fn, Q=Q, n_rows=100)

    buffer = io.BytesIO("".join(data_lines).encode("utf-8"))
    table = pacsv.read_csv(
        buffer,
        read_options=pacsv.ReadOptions(
            column_names=column_names,  # names for all columns
            use_threads=True,
        ),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(column_types=inferred_types),
    )
    return table


def _read_file(fn_list: List[dict], read_func: Callable) -> List[pa.Table]:
    return [read_func(f) for f in fn_list]


def concat_tables(tables: List[pa.Table]) -> pa.Table:
    return pa.concat_tables(tables, promote_options='default')


def _read_Q(fn_dict):
    fn = fn_dict["rfmix.Q"]
    table = _read_arrow_table(fn, Q=True)
    chrom_match = re.search(r'chr(\d+)', fn)
    chrom = chrom_match.group(0) if chrom_match else None
    chrom_col = pa.array([chrom] * table.num_rows)
    return table.append_column("chrom", chrom_col)


def _read_fb(fn_dict):
    fn = fn_dict["fb.tsv"]
    return _read_arrow_table(fn, Q=False)


def table_to_numpy(table: pa.Table, dtype: np.dtype = np.float32) -> np.ndarray:
    cols = [np.asarray(table[col].to_numpy(), dtype=dtype) for col in table.column_names]
    return np.column_stack(cols)


def _subset_populations(X: np.ndarray, npops: int) -> np.ndarray:
    pop_subset: list[np.ndarray] = []
    ncols = X.shape[1]
    if ncols % npops != 0:
        raise ValueError("The number of columns in X must be divisible by npops.")

    for pop_start in range(npops):
        X0 = X[:, pop_start::npops]
        if X0.shape[1] % 2 != 0:
            raise ValueError("Number of columns must be even.")
        X0_summed = X0[:, ::2] + X0[:, 1::2]
        pop_subset.append(X0_summed)
    return stack(pop_subset, axis=2)


def simulate_analysis(input_dir: str, task: int):
    fn_list = get_prefixes(input_dir, task)

    # Read and concatenate Q files
    g_anc_tables = _read_file(fn_list, _read_Q)
    g_anc = concat_tables(g_anc_tables)

    col_names = g_anc.column_names
    pops = [name for name in col_names if name not in ("sample_id", "chrom")]

    # Read and concatenate fb.tsv files
    X_tables = _read_file(fn_list, _read_fb)
    X = concat_tables(X_tables)

    loci = X.select(["chromosome", "physical_position"])
    drop_cols = [c for c in ["chromosome", "physical_position",
                             "genetic_position", "genetic_marker_index"]
                 if c in X.column_names]
    ancestry_table = X.drop(columns=drop_cols)
    ancestry_matrix = table_to_numpy(ancestry_table, dtype=np.float32)
    admix = _subset_populations(ancestry_matrix, len(pops))

    return loci, g_anc, admix


def get_peak_cpu_memory_mb() -> float:
    import resource
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_maxrss / 1024.0


def run_task(input_dir: str, output_path: str, label: str, task: int):
    output_dir = os.path.join(output_path, label)
    os.makedirs(output_dir, exist_ok=True)

    for replicate in range(1, 6):
        seed = replicate + 13
        random.seed(seed); np.random.seed(seed)

        logging.info(f"Replicate {replicate}: Reading files from {input_dir}")

        start = time.time()
        _, _, _ = simulate_analysis(input_dir, task)
        wall_time = time.time() - start
        peak_mem = get_peak_cpu_memory_mb()

        # Save metadata
        meta = collect_metadata("pyarrow", task, replicate, label)
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
    parser = argparse.ArgumentParser(description="Large File PyArrow Processor")
    parser.add_argument("--input", type=str, required=True,
                        help="Input directory with data files")
    parser.add_argument("--output", type=str, required=True, help="Output directory")
    parser.add_argument("--label", type=str, required=True,
                        help="Label for output directory")
    parser.add_argument("--task", type=int, choices=[1, 2, 3], required=True,
                        help="Task type: 1=small chrom, 2=large chrom, 3=all chroms")
    args = parser.parse_args()

    run_task(args.input, args.output, args.label, args.task)
