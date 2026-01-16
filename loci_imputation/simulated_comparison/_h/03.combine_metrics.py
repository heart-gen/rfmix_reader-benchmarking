#!/usr/bin/env python3
import argparse
import pandas as pd
import session_info
import logging, json
from pathlib import Path

def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Combine locus-level metrics by phase, population count, and "
            "interpolation method."
        )
    )
    parser.add_argument(
        "--metrics-root", type=Path, default=Path("../_m"),
        help="Root directory that contains phased/unphased results.",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("./combined_metrics"),
        help="Output directory for combined CSV files.",
    )
    return parser.parse_args()


def load_json(path: Path):
    if not path.exists():
        logging.warning("Missing metrics file: %s", path)
        return None
    with path.open() as handle:
        return json.load(handle)


def normalize_method(method_name: str) -> str:
    METHOD_LABELS = {
        "stepwise": "stepwise", "linear": "linear", "nearest": "nearest",
    }
    return METHOD_LABELS.get(method_name, method_name)


def collect_metrics(metrics_root: Path):
    summary_rows, breakpoint_rows = [], []
    for phase_dir in sorted(metrics_root.iterdir()):
        if not phase_dir.is_dir():
            continue

        phase = phase_dir.name
        for population_dir in sorted(phase_dir.glob("*_populations")):
            if not population_dir.is_dir():
                continue

            population = population_dir.name
            for method_dir in sorted(population_dir.iterdir()):
                if not method_dir.is_dir():
                    continue

                method = method_dir.name
                interpolation_method = normalize_method(method)
                locus_metrics = load_json(method_dir / "locus_metrics.json")
                segment_metrics = load_json(method_dir / "segment_metrics.json")
                switch_metrics = load_json(method_dir / "switch_error_rate.json")

                row = {
                    "phase": phase, "population": population, "method": method,
                    "interpolation_method": interpolation_method,
                }

                if locus_metrics:
                    row.update(
                        {
                            "overall_accuracy": locus_metrics.get("overall_accuracy"),
                            "mcc": locus_metrics.get("mcc"),
                        }
                    )
                    shape = locus_metrics.get("shape")
                    if shape and len(shape) == 3:
                        row.update(
                            {"n_variants": shape[0],
                             "n_samples": shape[1],
                             "n_labels": shape[2]}
                        )

                    per_ancestry = locus_metrics.get("per_ancestry_accuracy", {})
                    for label, value in per_ancestry.items():
                        row[f"per_ancestry_accuracy_{label}"] = value

                    precision = locus_metrics.get("precision", {})
                    recall    = locus_metrics.get("recall", {})
                    f1 = locus_metrics.get("f1", {})
                    labels = locus_metrics.get("labels", [])
                    for idx, label in enumerate(labels):
                        row[f"precision_{label}"] = precision.get(str(idx))
                        row[f"recall_{label}"] = recall.get(str(idx))
                        row[f"f1_{label}"] = f1.get(str(idx))

                if segment_metrics:
                    row.update(segment_metrics)

                if switch_metrics:
                    row.update(switch_metrics)

                summary_rows.append(row)
                breakpoint_path = method_dir / "breakpoint_error.tsv"
                if breakpoint_path.exists():
                    breakpoint_df = pd.read_csv(breakpoint_path, sep="\t")
                    breakpoint_df.insert(0, "phase", phase)
                    breakpoint_df.insert(1, "population", population)
                    breakpoint_df.insert(2, "method", method)
                    breakpoint_df.insert(3, "interpolation_method", interpolation_method)
                    breakpoint_rows.append(breakpoint_df)
                else:
                    logging.warning("Missing breakpoint error file: %s", breakpoint_path)

    return summary_rows, breakpoint_rows


def main():
    configure_logging()
    args = parse_args()

    # Parameter setup
    metrics_root = args.metrics_root.resolve()
    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Collect data
    summary_rows, breakpoint_rows = collect_metrics(metrics_root)
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_df.sort_values(
            ["phase", "population", "interpolation_method"], inplace=True
        )
        summary_df.to_csv(output_dir / "summary_metrics.csv", index=False)
        logging.info("Wrote summary metrics to %s", output_dir / "summary_metrics.csv")
    else:
        logging.warning("No summary metrics found under %s", metrics_root)

    if breakpoint_rows:
        breakpoint_df = pd.concat(breakpoint_rows, ignore_index=True)
        breakpoint_df.to_csv(output_dir / "breakpoint_metrics.csv", index=False)
        logging.info(
            "Wrote breakpoint metrics to %s", output_dir / "breakpoint_metrics.csv"
        )
    else:
        logging.warning("No breakpoint metrics found under %s", metrics_root)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
