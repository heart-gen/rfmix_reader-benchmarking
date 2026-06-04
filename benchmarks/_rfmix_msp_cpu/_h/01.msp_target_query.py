import argparse
import json
import logging
import os
import platform
import random
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
from rfmix_reader import extract_locus_ancestry, read_rfmix

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def collect_cpu_info() -> dict:
    info: dict[str, object] = {"os_cpu_count": os.cpu_count()}
    try:
        info["psutil_cpu_count_logical"] = psutil.cpu_count(logical=True)
        info["psutil_cpu_count_physical"] = psutil.cpu_count(logical=False)
    except Exception:
        pass
    for name in (
        "SLURM_CPUS_PER_TASK",
        "SLURM_JOB_CPUS_PER_NODE",
        "SLURM_NTASKS",
        "SLURM_TASKS_PER_NODE",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
    ):
        val = os.environ.get(name)
        if val is not None:
            info[name] = val
    return info


def peak_cpu_memory_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def collect_metadata(task_id: int, replicate: int, label: str, operation: str) -> dict:
    return {
        "parser": "rfmix_reader",
        "task": task_id,
        "replicate": replicate,
        "label": label,
        "backend": "CPU",
        "input_format": "msp",
        "operation": operation,
        "hardware": platform.platform(),
        "software_versions": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "psutil": psutil.__version__,
        },
        "cpu_info": collect_cpu_info(),
    }


def is_oom_error(exc: BaseException) -> bool:
    return isinstance(exc, MemoryError) or "out of memory" in str(exc).lower()


def msp_files(input_dir: str) -> list[Path]:
    root = Path(input_dir)
    files = sorted(root.glob("*.msp.tsv")) + sorted(root.glob("*.msp.tsv.gz"))
    if not files:
        raise FileNotFoundError(f"No .msp.tsv or .msp.tsv.gz files found under {input_dir}")
    return files


def read_msp_intervals(fn: Path) -> pd.DataFrame:
    return pd.read_csv(
        fn,
        sep="\t",
        skiprows=1,
        usecols=["#chm", "spos", "epos"],
        compression="infer",
    ).rename(columns={"#chm": "chrom"})


def derive_targets(input_dir: str, target_count: int, seed: int) -> pd.DataFrame:
    frames = []
    for fn in msp_files(input_dir):
        df = read_msp_intervals(fn)
        df["pos"] = ((df["spos"].astype(np.int64) + df["epos"].astype(np.int64)) // 2).astype(np.int64)
        frames.append(df[["chrom", "pos"]])
    targets = pd.concat(frames, ignore_index=True).drop_duplicates()
    if len(targets) > target_count:
        targets = targets.sample(n=target_count, random_state=seed).sort_values(["chrom", "pos"])
    targets = targets.reset_index(drop=True)
    targets.insert(0, "target_id", [f"target_{i + 1}" for i in range(len(targets))])
    return targets


def load_targets(path: str | None, input_dir: str, target_count: int, seed: int) -> pd.DataFrame:
    if path is None:
        return derive_targets(input_dir, target_count, seed)
    targets = pd.read_csv(path, sep="\t")
    required = {"chrom", "pos"}
    missing = required - set(targets.columns)
    if missing:
        raise ValueError(f"Target file is missing required columns: {sorted(missing)}")
    return targets


def write_meta(output_dir: Path, replicate: int, operation: str, meta: dict) -> None:
    with open(output_dir / f"meta_replicate_{replicate}_{operation}.json", "w") as fh:
        json.dump(meta, fh, indent=2)


def run_operation(output_dir: Path, replicate: int, operation: str, task: int, label: str, func):
    meta = collect_metadata(task, replicate, label, operation)
    meta.update({"status": "running", "oom_type": None, "wall_time_sec": None, "peak_cpu_memory_MB": None, "peak_gpu_memory_MB": 0.0})
    write_meta(output_dir, replicate, operation, meta)

    start = time.time()
    status = "success"
    oom_type = None
    error_msg = None
    extra = {}
    try:
        extra = func()
    except Exception as exc:
        status = "oom" if is_oom_error(exc) else "error"
        oom_type = "cpu" if status == "oom" else None
        error_msg = str(exc)[:500]
        logging.exception("%s failed on replicate %d", operation, replicate)
    finally:
        meta = collect_metadata(task, replicate, label, operation)
        meta.update({
            "status": status,
            "oom_type": oom_type,
            "wall_time_sec": time.time() - start,
            "peak_cpu_memory_MB": peak_cpu_memory_mb(),
            "peak_gpu_memory_MB": 0.0,
            **extra,
        })
        if error_msg:
            meta["error"] = error_msg
        write_meta(output_dir, replicate, operation, meta)


def run_task(input_dir: str, output_path: str, label: str, task: int, targets_path: str | None, target_count: int):
    output_dir = Path(output_path) / label
    output_dir.mkdir(parents=True, exist_ok=True)

    for replicate in range(1, 6):
        seed = replicate + 13
        random.seed(seed)
        np.random.seed(seed)
        logging.info("Replicate %d: MSP benchmark from %s", replicate, input_dir)

        run_operation(
            output_dir,
            replicate,
            "read_msp",
            task,
            label,
            lambda: _read_msp_operation(input_dir),
        )

        targets = load_targets(targets_path, input_dir, target_count, seed)
        targets_out = output_dir / f"targets_replicate_{replicate}.tsv"
        targets.to_csv(targets_out, sep="\t", index=False)

        run_operation(
            output_dir,
            replicate,
            "target_query",
            task,
            label,
            lambda: _target_query_operation(input_dir, targets, output_dir, replicate),
        )


def _read_msp_operation(input_dir: str) -> dict:
    loci, g_anc, admix = read_rfmix(input_dir, verbose=False, read_q=True)
    return {
        "n_loci_segments": int(len(loci)),
        "n_samples": int(admix.shape[1]),
        "n_ancestries": int(admix.shape[2]),
        "has_q": bool(g_anc is not None),
    }


def _target_query_operation(input_dir: str, targets: pd.DataFrame, output_dir: Path, replicate: int) -> dict:
    ancestry = extract_locus_ancestry(input_dir, targets, aggregate=True)
    out = output_dir / f"target_query_replicate_{replicate}.tsv"
    ancestry.to_csv(out, sep="\t", index=False)
    return {
        "n_targets": int(len(targets)),
        "n_matched_targets": int(ancestry["matched"].sum()) if "matched" in ancestry else 0,
        "target_query_output": str(out),
    }


def main():
    parser = argparse.ArgumentParser(description="CPU MSP RFMix-reader benchmark")
    parser.add_argument("--input", required=True, help="Directory with .msp.tsv(.gz) files")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--label", required=True, help="Label for output directory")
    parser.add_argument("--task", type=int, choices=[1, 2, 3], required=True)
    parser.add_argument("--targets", default=None, help="Optional TSV with chrom and pos columns")
    parser.add_argument("--target-count", type=int, default=1000, help="Number of interval-derived targets when --targets is absent")
    args = parser.parse_args()
    run_task(args.input, args.output, args.label, args.task, args.targets, args.target_count)


if __name__ == "__main__":
    main()
