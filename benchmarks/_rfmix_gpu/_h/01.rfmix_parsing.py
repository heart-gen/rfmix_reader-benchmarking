import cupy as cp
import pandas as pd
import os, logging, rmm
import argparse, platform
import psutil, time, json, random
import rmm.statistics as rmm_stats
from rfmix_reader import read_rfmix

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


def collect_metadata(parser_version, task_id, replicate, label, GPU):
    versions = {
        "python": platform.python_version(),
        "psutil": psutil.__version__,
    }
    if GPU:
        import cudf
        import cupy as cp
        versions["cudf"] = cudf.__version__
        versions["cupy"] = cp.__version__
        backend = "GPU"
    else:
        versions["pandas"] = pd.__version__
        backend = "CPU"
    meta = {
        "parser": parser_version,
        "task": task_id,
        "replicate": replicate,
        "label": label,
        "backend": backend,
        "hardware": platform.platform(),
        "software_versions": versions,
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


def load_data(prefix_path, binary_dir, BINARIES):
    loci, g_anc, admix = read_rfmix(
        prefix_path, binary_dir=binary_dir, generate_binary=BINARIES
    )
    admix = admix.compute()
    return g_anc


def get_peak_cpu_memory_mb() -> float:
    import resource
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_maxrss / 1024.0


def run_task(input_dir: str, output_path: str, label: str, task: int,
             BINARIES: bool, GPU: bool):
    output_dir = os.path.join(output_path, label)
    os.makedirs(output_dir, exist_ok=True)
    if GPU: configure_rmm()

    for replicate in range(1, 6):
        seed = replicate + 13
        random.seed(seed); cp.random.seed(seed)
        meta_path = os.path.join(output_dir, f"meta_replicate_{replicate}.json")

        logging.info(f"Replicate {replicate}: Reading files from {input_dir}")
        logging.info("CPU info: %s", collect_cpu_info())

        # Initial meta
        meta = collect_metadata("rfmix_reader", task, replicate, label, GPU)
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
                binaries_path = output_dir.replace("no_binaries", "binaries")
                _ = load_data(input_dir, binaries_path, BINARIES)
            except Exception as e:
                is_oom, kind = is_oom_error(e)
                status = "oom" if is_oom else "error"
                oom_kind = kind if is_oom else None
                error_msg = str(e)[:500]  # truncate
                logging.exception("Error on replicate %d: %s", replicate, e)
            finally:
                wall_time = time.time() - start
                stats = rmm_stats.get_statistics() if GPU else 0.0
                peak_gpu = float(stats.peak_bytes) / (1024 ** 2) if GPU else 0.0
                peak_cpu = get_peak_cpu_memory_mb()
            
                # Update meta (overwrite file)
                meta = collect_metadata("rfmix_reader", task, replicate, label, GPU)
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
    parser = argparse.ArgumentParser(description="Large File RFMix-Reader Processor")
    parser.add_argument("--input", type=str, required=True,
                        help="Input directory with data files")
    parser.add_argument("--output", type=str, required=True, help="Output directory")
    parser.add_argument("--label", type=str, required=True,
                        help="Label for output directory")
    parser.add_argument("--task", type=int, choices=[1, 2, 3], required=True,
                        help="Task type: 1=small chrom, 2=large chrom, 3=all chroms")
    parser.add_argument('--binaries',  action="store_true", help='Create binaries')
    parser.add_argument('--gpu', action="store_true", help="Run with GPU")
    args = parser.parse_args()

    run_task(args.input, args.output, args.label, args.task, args.binaries, args.gpu)
