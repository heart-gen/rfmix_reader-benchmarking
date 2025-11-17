import pandas as pd
import os, logging, rmm
import argparse, platform
import psutil, time, json, random
from rfmix_reader import read_flare

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

stats_resource = None # Global variable

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
        import numpy as cp
        versions["pandas"] = pd.__version__
        backend = "CPU"
    return {
        "parser": parser_version,
        "task": task_id,
        "replicate": replicate,
        "label": label,
        "backend": backend,
        "hardware": platform.platform(),
        "software_versions": versions,
    }


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


def load_data(prefix_path):
    loci, g_anc, admix = read_flare(prefix_path)
    admix = admix.compute()
    return g_anc


def get_peak_cpu_memory_mb() -> float:
    import resource
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_maxrss / 1024.0


def init_gpu_memory_tracking():
    global stats_resource
    rmm.reinitialize(pool_allocator=True)
    base = rmm.mr.get_current_device_resource()
    stats_resource = rmm.mr.StatisticsResourceAdapter(base)
    rmm.mr.set_current_device_resource(stats_resource)


def get_peak_gpu_memory_mb() -> float:
    if stats_resource is None:
        return 0.0
    return stats_resource.get_high_watermark() / 1024**2


def run_task(input_dir: str, label: str, task: int, GPU: bool):
    output_dir = os.path.join("output", label)
    os.makedirs(output_dir, exist_ok=True)

    init_gpu_memory_tracking()

    for replicate in range(3, 6):
        if stats_resource is not None:
            stats_resource.reset()

        seed = replicate
        random.seed(seed)
        cp.random.seed(seed)

        logging.info(f"Replicate {replicate}: Reading files from {input_dir}")

        start = time.time()
        status = "success"
        oom_kind = None
        error_msg = None

        try:
            g_anc = load_data(input_dir)
        except Exception as e:
            wall_time = time.time() - start
            is_oom, kind = is_oom_error(e)
            status = "oom" if is_oom else "error"
            oom_kind = kind if is_oom else None
            error_msg = str(e)[:500]  # truncate
            logging.exception("Error on replicate %d: %s", replicate, e)
        else:
            wall_time = time.time() - start
            
        peak_cpu = get_peak_cpu_memory_mb() # Always measure
        try:
            peak_gpu = get_peak_gpu_memory_mb()
        except Exception:
            peak_gpu = 0.0

        meta = collect_metadata("rfmix_reader", task, replicate, label, GPU)
        meta["status"] = status
        meta["oom_type"] = oom_kind
        meta["wall_time_sec"] = wall_time
        meta["peak_cpu_memory_MB"] = peak_cpu
        meta["peak_gpu_memory_MB"] = peak_gpu

        if status == "success":
            numeric = g_anc.select_dtypes(include=["float64", "float32", "int32", "int64"])
            if hasattr(numeric, "to_pandas"): # normalize to pandas
                numeric = numeric.to_pandas()
            g_anc_means = numeric.mean().to_dict()
            df_result = pd.DataFrame([g_anc_means])
            df_result.to_csv(
                os.path.join(output_dir, f"result_replicate_{replicate}.csv"),
                index=False,
            )
            
        meta_path = os.path.join(output_dir, f"meta_replicate_{replicate}.json")
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
    parser.add_argument("--label", type=str, required=True,
                        help="Label for output directory")
    parser.add_argument("--task", type=int, choices=[1, 2, 3], required=True,
                        help="Task type: 1=small chrom, 2=large chrom, 3=all chroms")
    parser.add_argument('--gpu', action="store_true", help="Run with GPU")
    args = parser.parse_args()

    run_task(args.input, args.label, args.task, args.gpu)
