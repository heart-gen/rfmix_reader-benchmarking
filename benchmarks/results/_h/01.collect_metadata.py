"""
Collect all per-parser meta_replicate_[1-5].json into a single metadata table.
"""
import json
import argparse
import session_info
import pandas as pd
from pathlib import Path

def collect_general_metrics(rootdir: str, outfile: str):
    rows = []
    for f in Path(rootdir).rglob("meta_replicate_*.json"):
        try:
            meta = json.load(open(f))
            path_parts = f.parts # Extract components
            try:
                m_index = path_parts.index("_m")
                framework = path_parts[m_index - 1].lstrip("_")
                pop_model = path_parts[m_index + 1] # two or three population
                # Detect if backend (cpu, gpu, binaries, no_binaries)
                maybe_backend = path_parts[m_index + 2]
                if maybe_backend.startswith("task_"):
                    backend = ""  # No intermediate backend dir
                    task = maybe_backend
                else:
                    backend = maybe_backend
                    task = path_parts[m_index + 3]
            except (ValueError, IndexError):
                framework, pop_model, backend, task = "", "", "", ""

            # Attach new fields
            meta["framework"] = framework # should match parser
            meta["population_model"] = pop_model
            meta["extras"] = backend

            rows.append(meta)

        except Exception as e:
            print(f"Skipping {f}: {e}")

    if not rows:
        raise RuntimeError(f"No meta_replicate_*.json files found under {rootdir}")

    df = pd.DataFrame(rows)
    df.to_csv(outfile, sep="\t", index=False)
    print(f"Collected {len(df)} runs → {outfile}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rootdir", required=True,
                        help="Root directory containing run subdirs")
    parser.add_argument("--out", default="metadata.tsv", help="Output file (TSV)")
    args = parser.parse_args()

    collect_general_metrics(args.rootdir, args.out)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
