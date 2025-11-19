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

    if "flare" in p_raw:
        return "rfmix-reader [flare]"
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


def format_task_pop(task, population_model) -> str:
    # Map raw population_model to pretty text
    pop_map = {
        "two_pop": "Two Ancestries",
        "three_pop": "Three Ancestries",
    }
    pretty_pop = pop_map.get(
        str(population_model),
        str(population_model).replace("_", " ").title()
    )
    return f"Task {task}\n({pretty_pop})"


def make_figure(summary: pd.DataFrame, out_prefix: str | Path) -> None:
    """
    Performance figure faceted by backend, colored by parser_clean.
    """
    sns.set(style="whitegrid", context="talk")

    df = summary.copy()
    df["task_pop"] = df.apply(
        lambda r: format_task_pop(r["task"], r["population_model"]),
        axis=1,
    )

    backends = ["CPU", "GPU"]
    metric_specs = [
        {
            "name": "Runtime",
            "y": "wall_time_median_sec",
            "q1": "wall_time_q1_sec",
            "q3": "wall_time_q3_sec",
            "title": "A. Runtime (median [IQR])",
            "ylabel": "Runtime (s)",
        },
        {
            "name": "CPU memory",
            "y": "cpu_mem_median_mb",
            "q1": "cpu_mem_q1_mb",
            "q3": "cpu_mem_q3_mb",
            "title": "B. Peak CPU memory (median [IQR])",
            "ylabel": "Peak CPU memory (MB)",
        },
        {
            "name": "GPU memory",
            "y": "gpu_mem_median_mb",
            "q1": "gpu_mem_q1_mb",
            "q3": "gpu_mem_q3_mb",
            "title": "C. Peak GPU memory (median [IQR])",
            "ylabel": "Peak GPU memory (MB)",
        },
    ]
    fig, axes = plt.subplots(
        nrows=3, ncols=2, figsize=(14, 16), sharex='col', constrained_layout=True,
    )

    def add_iqr_errorbars(ax, data, median_col, q1_col, q3_col):
        if data.empty:
            return
        for patch, (_, row) in zip(ax.patches, data.iterrows()):
            x = patch.get_x() + patch.get_width() / 2
            y = row[median_col]
            lower = y - row[q1_col]
            upper = row[q3_col] - y
            ax.errorbar(x, y, yerr=[[lower], [upper]], fmt="none", capsize=4)

    for row_idx, spec in enumerate(metric_specs):
        ycol = spec["y"]
        q1col = spec["q1"]
        q3col = spec["q3"]

        for col_idx, backend in enumerate(backends):
            ax = axes[row_idx, col_idx]
            sub = df[df["backend"] == backend].copy()

            if spec["name"] == "GPU memory" and backend == "CPU":
                ax.text(
                    0.5, 0.5, "No GPU memory for CPU runs",
                    ha="center", va="center",
                    transform=ax.transAxes,
                )
                ax.set_axis_off()
                continue

            # If no data for this backend/metric, annotate and skip
            if sub.empty or sub[ycol].isna().all():
                ax.text(
                    0.5,
                    0.5,
                    f"No data for {backend}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_axis_off()
                continue

            sns.barplot(
                data=sub,
                x="task_pop",
                y=ycol,
                hue="parser_clean",
                ax=ax,
            )
            add_iqr_errorbars(ax, sub, ycol, q1col, q3col)

            if col_idx == 0:
                ax.set_ylabel(spec["ylabel"])
            else:
                ax.set_ylabel("")

            if row_idx == len(metric_specs) - 1:
                ax.set_xlabel("Task / population model")
            else:
                ax.set_xlabel("")

            ax.set_title(f"{spec['title']} - {backend}")
            if col_idx == 1:
                ax.legend(title="Parser", fontsize=10)
            else:
                ax.legend_.remove()
            ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    png_path = out_prefix.parent / f"{out_prefix.name}_summary.png"
    pdf_path = out_prefix.parent / f"{out_prefix.name}_summary.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"[INFO] Wrote performance figures to: {png_path} and {pdf_path}")


def make_reliability_figure(summary: pd.DataFrame, out_prefix: str | Path) -> None:
    """Stacked barplot of reliability per parser/backend"""
    sns.set(style="whitegrid", context="talk")
    df = summary.copy()
    for col in ["n_success_reliability", "n_oom_reliability", "n_other_fail_reliability"]:
        df[col] = df[col].astype(float)

    # Aggregate reliability counts over all tasks/pop models
    agg = (
        df.groupby(["parser_clean", "backend"], as_index=False)[
            [
                "n_success_reliability",
                "n_oom_reliability",
                "n_other_fail_reliability",
                "n_total_reliability",
            ]
        ]
        .sum()
    )
    
    # Avoid division-by-zero
    agg = agg[agg["n_total_reliability"] > 0].copy()

    # Fractions
    agg["success_frac"] = agg["n_success_reliability"] / agg["n_total_reliability"]
    agg["oom_frac"] = agg["n_oom_reliability"] / agg["n_total_reliability"]
    agg["other_fail_frac"] = (
        agg["n_other_fail_reliability"] / agg["n_total_reliability"]
    )

    # Label for x-axis: "rfmix-reader [binaries] (GPU)" etc.
    agg["label"] = agg["parser_clean"] + " (" + agg["backend"] + ")"

    # Stacked barplot using matplotlib
    x = np.arange(len(agg))
    width = 0.8

    fig, ax = plt.subplots(figsize=(10, 6))

    bottom = np.zeros(len(agg))

    for outcome, col in [("Success", "success_frac"),
                         ("OOM", "oom_frac"),
                         ("Other failure", "other_fail_frac")]:
        ax.bar(x, agg[col].values, width, bottom=bottom, label=outcome)
        bottom += agg[col].values

    ax.set_ylabel("Fraction of replicates")
    ax.set_xlabel("Parser / backend")
    ax.set_title("Reliability across replicates (per parser/backend)")
    ax.set_xticks(x)
    ax.set_xticklabels(agg["label"], rotation=45, ha="right")
    ax.set_ylim(0, 1.0)

    ax.legend(title="Outcome")
    plt.tight_layout()

    png_path = out_prefix.parent / f"{out_prefix.name}_reliability.png"
    pdf_path = out_prefix.parent / f"{out_prefix.name}_reliability.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
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
