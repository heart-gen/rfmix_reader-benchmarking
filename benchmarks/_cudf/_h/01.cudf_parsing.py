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
import rmm.statistics as rmm_stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def configure_rmm():
    rmm.reinitialize(pool_allocator=True)


def collect_cpu_info() -> dict:
    info: dict[str, object] = {}
    # OS / psutil view of cores
    info["os_cpu_count"] = os.cpu_count()
    try:
        info["psutil_cpu_count_logical"] = psutil.cpu_count(logical=True)
        info["psutil_cpu_count_physical"] = psutil.cpu_count(logical=False)
    except Exception:
        pass

    # SLURM environment (if present)
    for name in (
        "SLURM_CPUS_PER_TASK",
        "SLURM_JOB_CPUS_PER_NODE",
        "SLURM_NTASKS",
        "SLURM_TASKS_PER_NODE",
    ):
        val = os.environ.get(name)
        if val is not None:
            info[name] = val

    # Thread-related environment variables
    for name in (
        "POLARS_MAX_THREADS",
        "RAYON_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
    ):
        val = os.environ.get(name)
        if val is not None:
            info[f"env_{name}"] = val

    # Library-specific thread pools (if available)
    try:
        import polars as pl
        if hasattr(pl, "thread_pool_size"):
            info["polars_thread_pool_size"] = pl.thread_pool_size()
        elif hasattr(pl, "threadpool_size"):
            info["polars_thread_pool_size"] = pl.threadpool_size()
    except Exception:
        pass

    try:
        import pyarrow as pa
        info["pyarrow_cpu_count"] = pa.cpu_count()
        if hasattr(pa, "io_thread_count"):
            info["pyarrow_io_thread_count"] = pa.io_thread_count()
    except Exception:
        pass

    return info


def collect_metadata(parser_version, task_id, replicate, label):
    meta = {
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
    meta["cpu_info"] = collect_cpu_info()
    return meta


def is_oom_error(e: BaseException) -> tuple[bool, str]:
    """Heuristic classification of OOM vs other error."""
    msg = str(e).lower()
    if isinstance(e, MemoryError):
        return True, "cpu"
    # cuDF / CuPy / RMM / CUDA messages
    oom_keywords = ["out of memory", "cuda error: out of memory", "rmm_alloc"]
    if any(k in msg for k in oom_keywords):
        return True, "gpu"
    return False, "unknown"


def get_prefixes(input_dir: str, task: int, verbose: bool = True) -> list[dict]:
    if task == 1:
        chroms = [21]
    elif task == 2:
        chroms = [1]
    elif task == 3:
        chroms = list(range(1, 23))
    else:
        raise ValueError(f"Unsupported task id: {task}")

    prefixes = [{"rfmix.Q": os.path.join(input_dir, f"chr{chrom}.rfmix.Q"),
                 "fb.tsv": os.path.join(input_dir, f"chr{chrom}.fb.tsv")}
                for chrom in chroms]

    if verbose:
        logging.info("Task %d -> %d chromosome(s): %s", task, len(chroms), chroms)

    return prefixes


def _types(fn: str, Q: bool = False, n_rows: int = 100) -> dict:
    """Rough type inference using cudf sample."""
    if Q:
        sample = cudf.read_csv(fn, sep="\t", skiprows=1, nrows=n_rows)
        dtypes = {}
        dtypes["sample_id"] = "category"
        for col, dt in zip(sample.columns[1:], sample.dtypes[1:]):
            dtypes[str(col)] = dt
        return dtypes
    else:
        sample = cudf.read_csv(fn, sep="\t", nrows=n_rows, header=0)
        return {str(col): dt for col, dt in zip(sample.columns, sample.dtypes)}


def _read_csv(fn: str, dtypes: dict, Q: bool = False) -> cudf.DataFrame:
    if Q:
        col_names = list(dtypes.keys())
        df = cudf.read_csv(
            fn, sep='\t', header=None, comment="#",
            names=col_names, dtype=dtypes
        )
    else:
        df = cudf.read_csv(fn, sep="\t", comment="#", dtype=dtypes)
    return df


def _read_file(fn_list: List[dict], read_func: Callable) -> List[cudf.DataFrame]:
    return [read_func(f) for f in fn_list]


def concat_tables(tables: List[cudf.DataFrame]) -> cudf.DataFrame:
    return cudf.concat(tables, ignore_index=True)


def _read_Q(fn_dict):
    fn = fn_dict["rfmix.Q"]
    dtypes = _types(fn, Q=True)
    df = _read_csv(fn, dtypes=dtypes, Q=True)
    chrom_match = re.search(r'chr(\d+)', fn)
    chrom = chrom_match.group(0) if chrom_match else None
    df["chrom"] = chrom
    return df


def _read_fb(fn_dict):
    fn = fn_dict["fb.tsv"]
    dtypes = _types(fn, Q=False)
    return _read_csv(fn, dtypes=dtypes, Q=False)


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


def _dict_to_df_single_row(d: dict[str, float]) -> cudf.DataFrame:
    return cudf.DataFrame([d])


def run_task(input_dir: str, output_path: str, label: str, task: int):
    output_dir = os.path.join(output_path, label)
    os.makedirs(output_dir, exist_ok=True)
    configure_rmm()

    for replicate in range(1, 6):
        seed = replicate + 13
        random.seed(seed); cp.random.seed(seed)
        meta_path = os.path.join(output_dir, f"meta_replicate_{replicate}.json")

        logging.info(f"Replicate {replicate}: Reading files from {input_dir}")
        logging.info("CPU info: %s", collect_cpu_info())

        # Initial meta
        meta = collect_metadata("cudf", task, replicate, label)
        meta["status"] = "running"
        meta["oom_type"] = None
        meta["wall_time_sec"] = None
        meta["peak_cpu_memory_MB"] = None
        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=2)

        with rmm_stats.statistics():
            status = "success"
            oom_kind = error_msg = None
            start = time.time()
            try:
                _, _, _ = simulate_analysis(input_dir, task)
            except Exception as e:
                is_oom, kind = is_oom_error(e)
                status = "oom" if is_oom else "error"
                oom_kind = kind if is_oom else None
                error_msg = str(e)[:500]  # truncate
                logging.exception("Error on replicate %d: %s", replicate, e)
            finally:
                wall_time = time.time() - start
                stats = rmm_stats.get_statistics()
                peak_gpu = float(stats.peak_bytes) / (1024 ** 2)            
                peak_cpu = get_peak_cpu_memory_mb()

                # Update meta (overwrite file)
                meta = collect_metadata("cudf", task, replicate, label)
                meta["status"] = status
                meta["oom_type"] = oom_kind
                meta["wall_time_sec"] = wall_time
                meta["peak_cpu_memory_MB"] = peak_cpu
                meta["peak_gpu_memory_MB"] = peak_gpu
                if error_msg:
                    meta["error"] = error_msg
            
                with open(meta_path, "w") as f:
                    json.dump(meta, f, indent=2)

        logging.info(
            "Finished replicate %d in %.2fs, peak RSS: %.2f MB, peak GPU: %.2f MB",
            replicate, wall_time, peak_cpu, peak_gpu
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Large File cuDF Processor")
    parser.add_argument("--input", type=str, required=True,
                        help="Input directory with data files")
    parser.add_argument("--output", type=str, required=True,
                        help="Output directory")
    parser.add_argument("--label", type=str, required=True,
                        help="Label for output directory")
    parser.add_argument("--task", type=int, choices=[1, 2, 3], required=True,
                        help="Task type: 1=small chrom, 2=large chrom, 3=all chroms")
    args = parser.parse_args()

    run_task(args.input, args.output, args.label, args.task)
