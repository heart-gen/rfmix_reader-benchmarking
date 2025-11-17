import rmm
import cudf
import os, re
import logging
import platform
import argparse
import cupy as cp
import json, random
import psutil, time
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
            "cudf": cudf.__version__,
            "cupy": cp.__version__,
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

    prefixes = [{"rfmix.Q": os.path.join(input_dir, f"chr{chrom}.Q"),
                 "fb.tsv": os.path.join(input_dir, f"chr{chrom}.fb.tsv")}
                for chrom in chroms]

    if verbose:
        logging.info("Task %d -> %d chromosomes: %s", task, len(chroms), chroms)

    return prefixes


def _types(fn: str, Q: bool = False, n_rows: int = 100) -> dict:
    # Rough type inference using cudf sample
    sample = cudf.read_csv(fn, sep='\t', skiprows=1, nrows=n_rows)
    if Q:
        header = {"sample_id": "category"}
        header.update({k: str(v) for k, v in zip(sample.columns[1:], sample.dtypes[1:])})
    else:
        header = {k: str(v) for k, v in zip(sample.columns, sample.dtypes)}
    return header


def _read_csv(fn: str, header=None, Q: bool = False) -> cudf.DataFrame:
    col_names = list(header.keys())
    dtypes = {i: header[name] for i, name in enumerate(col_names)}
    if Q:
        df = cudf.read_csv(fn, sep='\t', header=None, comment="#",
                           names=col_names, dtype=dtypes)
    else:
        df = cudf.read_csv(fn, sep="\t", comment="#", dtype=dtypes)
    return df


def _read_file(fn_list: List[dict], read_func: Callable) -> List[cudf.DataFrame]:
    return [read_func(f) for f in fn_list]


def concat_tables(tables: List[cudf.DataFrame]) -> cudf.DataFrame:
    return cudf.concat(tables, ignore_index=True)


def _read_Q(fn_dict):
    fn = fn_dict["rfmix.Q"]
    header = _types(fn, Q=True)
    df = _read_csv(fn, header=header, Q=True)
    chrom_match = re.search(r'chr(\d+)', fn)
    chrom = chrom_match.group(0) if chrom_match else None
    df["chrom"] = chrom
    return df


def _read_fb(fn_dict):
    fn = fn_dict["fb.tsv"]
    header = _types(fn, Q=False)
    return _read_csv(fn, header=header, Q=False)


def table_to_cupy(df: cudf.DataFrame, dtype=cp.float32) -> cp.ndarray:
    return df.to_cupy().astype(dtype, copy=False)


def _subset_populations(X: cp.ndarray, npops: int) -> cp.ndarray:
    pop_subset = []
    ncols = X.shape[1]
    if ncols % npops != 0:
        raise ValueError("Number of columns in X must be divisible by npops.")
    for pop_start in range(npops):
        X0 = X[:, pop_start::npops]
        if X0.shape[1] % 2 != 0:
            raise ValueError("Number of columns per population must be even.")
        X0_summed = X0[:, ::2] + X0[:, 1::2]
        pop_subset.append(X0_summed)
    return cp.stack(pop_subset, axis=2)


def simulate_analysis(input_dir: str, task: int):
    fn_list = get_prefixes(input_dir, task)
    g_anc_tables = _read_file(fn_list, _read_Q)
    g_anc = concat_tables(g_anc_tables)

    pops = [col for col in g_anc.columns if col not in ("sample_id", "chrom")]

    X_tables = _read_file(fn_list, _read_fb)
    X = concat_tables(X_tables)

    loci = X[["chromosome", "physical_position"]]
    drop_cols = [c for c in ["chromosome", "physical_position", "genetic_position",
                             "genetic_marker_index"] if c in X.columns]
    ancestry_df = X.drop(columns=drop_cols)
    ancestry_df = ancestry_df.select_dtypes(include=["float64", "float32",
                                                     "int32", "int64"])

    ancestry_matrix = table_to_cupy(ancestry_df)
    admix = _subset_populations(ancestry_matrix, len(pops))

    return loci, g_anc, admix


def get_peak_cpu_memory_mb() -> float:
    import resource
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_maxrss / 1024.0


def init_rmm_with_stats():
    base = rmm.mr.CudaMemoryResource()
    pool = rmm.mr.PoolMemoryResource(base)
    stats = rmm.mr.StatisticsResourceAdapter(pool)
    rmm.mr.set_current_device_resource(stats)
    return stats


def _dict_to_df_single_row(d: dict[str, float]) -> cudf.DataFrame:
    return cudf.DataFrame([d])


def run_task(input_dir: str, label: str, task: int):
    output_dir = os.path.join("output", label)
    os.makedirs(output_dir, exist_ok=True)

    global stats_resource

    for replicate in range(3, 6):
        seed = replicate
        random.seed(seed)
        cp.random.seed(seed)

        logging.info(f"Replicate {replicate}: Reading files from {input_dir}")

        try:
            start = time.time()
            _, g_anc, _ = simulate_analysis(input_dir, task)
            wall_time = time.time() - start
            peak_mem = get_peak_cpu_memory_mb()
            gpu_peak_bytes = stats_resource.get_high_watermark()
            gpu_peak_mb = gpu_peak_bytes / 1024**2

            g_anc_means = g_anc.select_dtypes(include=["float64", "float32", "int32", "int64"]).mean().to_pandas().to_dict()

            result_path = os.path.join(output_dir, f"result_replicate_{replicate}.csv")
            df_result = _dict_to_df_single_row(g_anc_means)
            df_result.to_pandas().to_csv(result_path, index=False)

            meta = collect_metadata("cudf", task, replicate, label)
            meta["wall_time_sec"] = wall_time
            meta["peak_cpu_memory_MB"] = peak_mem
            meta["peak_gpu_memory_MB"] = gpu_mem

            meta_path = os.path.join(output_dir, f"meta_replicate_{replicate}.json")
            with open(meta_path, "w") as f:
                json.dump(meta, f, indent=2)

            logging.info(
                "Finished replicate %d in %.2fs, peak RSS: %.2f MB, peak GPU: %.2f MB",
                replicate, wall_time, peak_mem, gpu_mem
            )

        except (MemoryError, RuntimeError) as e:
            if "out of memory" in str(e).lower():
                logging.error("GPU OOM on replicate %d: %s", replicate, str(e))
            else:
                logging.exception("Unexpected error on replicate %d", replicate)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Large File cuDF Processor")
    parser.add_argument("--input", type=str, required=True,
                        help="Input directory with data files")
    parser.add_argument("--label", type=str, required=True,
                        help="Label for output directory")
    parser.add_argument("--task", type=int, choices=[1, 2, 3], required=True,
                        help="Task type: 1=small chrom, 2=large chrom, 3=all chroms")
    args = parser.parse_args()

    run_task(args.input, args.label, args.task)
