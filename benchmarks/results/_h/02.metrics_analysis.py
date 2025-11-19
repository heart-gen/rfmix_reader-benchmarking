"""
Analyze benchmark metadata and generate summary table + figure.
"""
import numpy as np
import session_info
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt

def infer_backend(row: pd.Series) -> str:
    """Infer CPU vs GPU backend."""
    parser_raw = str(row.get("parser", "")).lower()
    backend_raw = str(row.get("backend", "")).strip()
    extras_raw = str(row.get("extras", "")).lower()

    # Keep clean CPU/GPU if already present
    if backend_raw.upper() in {"CPU", "GPU"}:
        return backend_raw.upper()

    # Heuristics based on parser
    if parser_raw in {"pandas", "pyarrow", "polars"}:
        return "CPU"
    if parser_raw == "cudf":
        return "GPU"

    # Fallback to extras tagging
    if "gpu" in extras_raw:
        return "GPU"
    if "cpu" in extras_raw:
        return "CPU"

    # Default
    return "CPU"


def clean_parser(row: pd.Series) -> str:
    """Normalize parser names."""
    p_raw = str(row.get("parser", "")).lower()
    extras_raw = str(row.get("extras", "")).lower()

    if "rfmix_reader" in p_raw:
        if extras_raw == "binaries":
            return "rfmix-reader [binaries]"
        else:
            return "rfmix-reader [no-binaries]"
    return p_raw


def add_missing_replicates(df: pd.DataFrame) -> pd.DataFrame:
    # Ensure replicate is integer
    df["replicate"] = df["replicate"].astype(int)

    # Define the expected combinations
    all_reps = pd.DataFrame({"replicate": [1, 2, 3, 4, 5]})
    group_keys = ["parser", "task", "population_model"]
    
    # Create a full index of all group x replicate combinations
    full_index = (
        df[group_keys].drop_duplicates().merge(all_reps, how="cross")
    )

    # Merge full index with data to include missing rows
    df_full = full_index.merge(df, on=group_keys + ["replicate"], how="left")

    # Update backend and success
    df_full["backend"] = df_full.apply(infer_backend, axis=1)
    df_full["status"] = df_full["status"].fillna("running")
    return df_full


def load_and_clean_metadata(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")

    # Ensure expected columns exist
    required = [
        "parser", "task", "replicate", "backend", "status", "wall_time_sec",
        "peak_cpu_memory_MB", "peak_gpu_memory_MB", "population_model",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in metadata: {missing}")

    # Infer backend (CPU/GPU)
    df["backend"] = df.apply(infer_backend, axis=1)

    # Clean parser names
    df["parser_clean"] = df.apply(clean_parser, axis=1)

    # Add missing replicates that failed
    df = add_missing_replicates(df)

    # Coerce to convenient dtypes
    df["task"] = df["task"].astype(str)
    df["population_model"] = df["population_model"].astype(str)
    return df


def summarize(df: pd.DataFrame, expected_reps: int = 5) -> pd.DataFrame:
    """Summarize by parser_clean, backend, task, and population_model."""
    group_cols = ["parser_clean", "backend", "task", "population_model"]

    def summary_stats(series: pd.Series):
        """Return (median, q1, q3, iqr) after dropping NaNs; all NaN if no data."""
        vals = pd.to_numeric(series, errors="coerce").dropna()
        if len(vals) == 0:
            return np.nan, np.nan, np.nan, np.nan
        q1 = vals.quantile(0.25)
        q3 = vals.quantile(0.75)
        median = vals.median()
        iqr = q3 - q1
        return median, q1, q3, iqr

    records = []
    for key, gdf in df.groupby(group_cols):
        parser_clean, backend, task, pop_model = key

        # Logged counts
        n_logged = len(gdf)
        n_success_logged = (gdf["status"] == "success").sum()
        status_lower = gdf["status"].astype(str).str.lower()
        is_oom = gdf["oom_type"].fillna("").ne("") | (status_lower == "running")
        n_oom_logged = is_oom.sum()
        n_other_fail_logged = n_logged - n_success_logged - n_oom_logged

        # Reliability counts
        if n_success_logged == 0 and n_oom_logged > 0:
            # Early-stop: first replicate OOM => assume all expected_reps are OOM
            n_total_reliability = expected_reps
            n_oom_reliability = expected_reps
            n_success_reliability = 0
            n_other_fail_reliability = 0
        else:
            # Actually run
            n_total_reliability = n_logged
            n_oom_reliability = n_oom_logged
            n_success_reliability = n_success_logged
            n_other_fail_reliability = n_other_fail_logged

        # Stats from successful runs only
        g_success = gdf[gdf["status"] == "success"]
        if g_success.empty:
            wall_time_median_sec = wall_time_q1_sec = wall_time_q3_sec = wall_time_iqr_sec = np.nan
            cpu_mem_median_mb = cpu_mem_q1_mb = cpu_mem_q3_mb = cpu_mem_iqr_mb = np.nan
            gpu_mem_median_mb = gpu_mem_q1_mb = gpu_mem_q3_mb = gpu_mem_iqr_mb = np.nan
        else:
            (
                wall_time_median_sec, wall_time_q1_sec, wall_time_q3_sec, wall_time_iqr_sec,
            ) = summary_stats(g_success["wall_time_sec"])

            (
                cpu_mem_median_mb, cpu_mem_q1_mb, cpu_mem_q3_mb, cpu_mem_iqr_mb,
            ) = summary_stats(g_success["peak_cpu_memory_MB"])

            (
                gpu_mem_median_mb, gpu_mem_q1_mb, gpu_mem_q3_mb, gpu_mem_iqr_mb,
            ) = summary_stats(g_success["peak_gpu_memory_MB"])

        rec = {
            "parser_clean": parser_clean, "backend": backend,
            "task": task, "population_model": pop_model,

            # logged counts
            "n_logged": n_logged, "n_success_logged": n_success_logged,
            "n_oom_logged": n_oom_logged, "n_other_fail_logged": n_other_fail_logged,

            # reliability counts
            "n_total_reliability": n_total_reliability,
            "n_success_reliability": n_success_reliability,
            "n_oom_reliability": n_oom_reliability,
            "n_other_fail_reliability": n_other_fail_reliability,

            # performance stats
            "wall_time_median_sec": wall_time_median_sec,
            "wall_time_q1_sec": wall_time_q1_sec,
            "wall_time_q3_sec": wall_time_q3_sec,
            "wall_time_iqr_sec": wall_time_iqr_sec,
            "cpu_mem_median_mb": cpu_mem_median_mb,
            "cpu_mem_q1_mb": cpu_mem_q1_mb,
            "cpu_mem_q3_mb": cpu_mem_q3_mb,
            "cpu_mem_iqr_mb": cpu_mem_iqr_mb,
            "gpu_mem_median_mb": gpu_mem_median_mb,
            "gpu_mem_q1_mb": gpu_mem_q1_mb,
            "gpu_mem_q3_mb": gpu_mem_q3_mb,
            "gpu_mem_iqr_mb": gpu_mem_iqr_mb,
        }
        records.append(rec)

    summary = pd.DataFrame.from_records(records)
    summary = summary.sort_values(
        ["parser_clean", "population_model", "task", "backend"]
    )
    return summary


def make_figure(summary: pd.DataFrame, out_prefix: str | Path) -> None:
    """
    Create a figure with:
      * Panel A: runtime median (+ IQR) by parser/backend/task (faceted by population_model)
      * Panel B: CPU memory median (+ IQR)
      * Panel C: GPU memory median (+ IQR, GPU only)
    """
    sns.set(style="whitegrid", context="talk")

    summary = summary.copy()
    summary["task_pop"] = summary["task"].astype(str) + " (" + summary["population_model"].astype(str) + ")"

    fig, axes = plt.subplots(
        nrows=3, ncols=1, figsize=(12, 16), sharex=True, constrained_layout=True,
    )

    def add_iqr_errorbars(ax, data, median_col, q1_col, q3_col):
        for patch, (_, row) in zip(ax.patches, data.iterrows()):
            x = patch.get_x() + patch.get_width() / 2
            y = row[median_col]
            lower = y - row[q1_col]
            upper = row[q3_col] - y
            ax.errorbar(x, y, yerr=[[lower], [upper]], fmt="none", capsize=4)

    # Panel A: Runtime
    ax = axes[0]
    runtime_df = summary
    sns.barplot(
        data=runtime_df, x="task_pop", y="wall_time_median_sec", hue="backend", ax=ax,
    )
    add_iqr_errorbars(
        ax, runtime_df, "wall_time_median_sec", "wall_time_q1_sec", "wall_time_q3_sec",
    )
    ax.set_ylabel("Runtime (s, median [IQR])")
    ax.set_xlabel("")
    ax.set_title("A. Runtime by task, population model, parser, and backend")
    ax.legend(title="Backend")

    # Panel B: CPU memory
    ax = axes[1]
    cpu_df = summary
    sns.barplot(
        data=cpu_df, x="task_pop", y="cpu_mem_median_mb", hue="backend", ax=ax,
    )
    add_iqr_errorbars(
        ax, cpu_df, "cpu_mem_median_mb", "cpu_mem_q1_mb", "cpu_mem_q3_mb",
    )
    ax.set_ylabel("Peak CPU memory (MB)")
    ax.set_xlabel("")
    ax.set_title("B. Peak CPU memory")
    ax.legend_.remove()

    # Panel C: GPU memory (GPU only)
    ax = axes[2]
    gpu_df = summary[summary["backend"] == "GPU"].copy()
    sns.barplot(
        data=gpu_df, x="task_pop", y="gpu_mem_median_mb", hue="parser_clean", ax=ax,
    )
    add_iqr_errorbars(
        ax, gpu_df, "gpu_mem_median_mb", "gpu_mem_q1_mb", "gpu_mem_q3_mb",
    )
    ax.set_ylabel("Peak GPU memory (MB)")
    ax.set_xlabel("Task (population model)")
    ax.set_title("C. Peak GPU memory (GPU runs only)")
    ax.legend(title="Parser")

    png_path = out_prefix.parent / f"{out_prefix.name}_summary.png"
    pdf_path = out_prefix.parent / f"{out_prefix.name}_summary.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"[INFO] Wrote performance figures to: {png_path} and {pdf_path}")


def make_reliability_figure(summary: pd.DataFrame, out_prefix: str | Path) -> None:
    """Fraction of replicates in each outcome category."""
    df = summary.copy()
    for col in ["n_success_reliability", "n_oom_reliability", "n_other_fail_reliability"]:
        df[col] = df[col].astype(float)

    df["success_frac"] = df["n_success_reliability"] / df["n_total_reliability"]
    df["oom_frac"] = df["n_oom_reliability"] / df["n_total_reliability"]
    df["other_fail_frac"] = df["n_other_fail_reliability"] / df["n_total_reliability"]

    long = df.melt(
        id_vars=["parser_clean", "backend", "task", "population_model"],
        value_vars=["success_frac", "oom_frac", "other_fail_frac"],
        var_name="outcome", value_name="fraction",
    )
    mapping = {
        "success_frac": "Success", "oom_frac": "OOM",
        "other_fail_frac": "Other failure",
    }
    long["outcome"] = long["outcome"].map(mapping)

    sns.set(style="whitegrid", context="talk")
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(
        data=long, x="parser_clean", y="fraction", hue="outcome",
    )
    ax.set_ylabel("Fraction of replicates")
    ax.set_xlabel("Parser")
    ax.set_title("Reliability across replicates (per parser)")
    plt.legend(title="Outcome")
    plt.tight_layout()

    png_path = out_prefix.parent / f"{out_prefix.name}_reliability.png"
    pdf_path = out_prefix.parent / f"{out_prefix.name}_reliability.pdf"
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()
    print(f"[INFO] Wrote reliability figures to: {png_path} and {pdf_path}")


def main():
    metadata_path = Path("../_m/metadata.tsv")
    out_prefix = Path("benchmark")

    df = load_and_clean_metadata(metadata_path)

    # Write cleaned per-run table (for R)
    clean_tsv = out_prefix.parent / f"{out_prefix.name}_clean.tsv"
    cols_for_r = [
        "parser", "parser_clean", "backend", "task", "population_model",
        "replicate", "label", "status", "oom_type", "wall_time_sec",
        "peak_cpu_memory_MB", "peak_gpu_memory_MB", "extras",
    ]
    cols_for_r = [c for c in cols_for_r if c in df.columns]
    df[cols_for_r].to_csv(clean_tsv, sep="\t", index=False)
    print(f"[INFO] Wrote cleaned per-run table to: {clean_tsv}")

    # Summary table
    summary = summarize(df)
    summary_tsv = out_prefix.parent / f"{out_prefix.name}_summary.tsv"
    summary.to_csv(summary_tsv, sep="\t", index=False)
    print(f"[INFO] Wrote summary table to: {summary_tsv}")

    # Figure
    make_figure(summary, out_prefix)
    make_reliability_figure(summary, out_prefix)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
