# This script will test the memory usage and execution time for
# RFMix-reader using Polars.
import numpy as np
import polars as pl
from time import time
from pyhere import here
from re import search as rsearch
from typing import Callable, List
from memory_profiler import profile
from rfmix_reader import get_prefixes
from collections import OrderedDict as odict

@profile
def read_data(prefix_path, verbose=True):
    fn = get_prefixes(prefix_path, verbose)

    # Global ancestry files
    g_anc = _read_file(fn, lambda f: _read_Q(f["rfmix.Q"]))
    pops = g_anc[0].drop(["sample_id", "chrom"]).columns
    g_anc = pl.concat(g_anc, how="vertical")

    # Local ancestry files
    X = _read_file(fn, lambda f: _read_fb(f["fb.tsv"]))
    X = pl.concat(X, how="vertical")

    loci = X.select(["chromosome", "physical_position"])
    X = X.select(pl.all().exclude(["chromosome", "physical_position", "sample_id"]))
    admix = _subset_populations(X.to_numpy(), len(pops))


def _read_file(fn: List[str], read_func: Callable) -> List[pl.DataFrame]:
    """
    Read data from multiple files using a provided read function.
    """
    return [read_func(file_name) for file_name in fn]


def _read_fb(fn: str) -> pl.DataFrame:
    header = odict(_types(fn, False))
    return _read_tsv(fn, header)


def _read_g_anc(fn: str) -> pl.DataFrame:
    """
    Read the global ancestry matrix from a file and add chromosome information.
    """
    header = _types(fn)
    df = _read_csv(fn, header)
    m = rsearch(r"chr(\d+)", fn)
    chrom = m.group(0) if m else None
    if chrom:
        df = df.with_columns(pl.lit(chrom).alias("chrom"))
    else:
        print(f"Warning: Could not extract chromosome information from '{fn}'")
    return df


def _read_csv(fn: str, header: dict) -> pl.DataFrame:
    """
    Read a CSV/TSV file into a Polars DataFrame with specified data types.
    """
    return pl.read_csv(
        fn,
        has_header=False,
        separator="\t",
        new_columns=list(header.keys()),
        dtypes=header,
        comment_char="#",
    )


def _read_tsv(fn: str, header: dict) -> pl.DataFrame:
    """
    Read a TSV file into a Polars DataFrame.
    """
    return pl.read_csv(
        fn,
        has_header=True,
        separator="\t",
        new_columns=list(header.keys()),
        dtypes=header,
        comment_char="#",
    )


def _types(fn: str, Q: bool = True) -> dict:
    """
    Infer the data types of columns in a TSV file using Polars.
    """
    df = pl.read_csv(fn, has_header=True, separator="\t", n_rows=2, skip_rows=1)
    if Q:
        header = {"sample_id": pl.Utf8}
        header.update({col: df[col].dtype for col in df.columns[1:]})
    else:
        header = {col: df[col].dtype for col in df.columns}
    return header


def _subset_populations(X: np.ndarray, npops: int) -> np.ndarray:
    """
    Subset and process the input array X based on populations.
    """
    pop_subset = []
    ncols = X.shape[1]
    if ncols % npops != 0:
        raise ValueError("The number of columns in X must be divisible by npops.")
    for pop_start in range(npops):
        X0 = X[:, pop_start::npops]
        if X0.shape[1] % 2 != 0:
            raise ValueError("Number of columns must be even.")
        X0_summed = X0[:, ::2] + X0[:, 1::2]
        pop_subset.append(X0_summed)
    return np.stack(pop_subset, 1)


def main():
    prefix_path = here("input/simulations/two_populations/_m/rfmix-out/")
    start_time = time()
    _ = read_data(prefix_path)
    end_time = time()
    print(f"Execution time: {end_time - start_time:.2f} seconds")


if __name__ == "__main__":
    main()
