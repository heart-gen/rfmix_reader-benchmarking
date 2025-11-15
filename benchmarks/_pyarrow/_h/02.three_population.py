# This script will test the memory usage and executive time for
# RFMix-reader.
from time import time
from pyhere import here
from typing import Callable, List
from memory_profiler import profile
from rfmix_reader import get_prefixes
from collections import OrderedDict as odict
from pandas import DataFrame, read_csv, concat

@profile
def read_data(prefix_path, verbose=True):
    fn = get_prefixes(prefix_path, verbose)
    rf_q = _read_file(fn, lambda f: _read_Q(f["rfmix.Q"]))
    pops = rf_q[0].drop(["sample_id", "chrom"], axis=1).columns.values
    rf_q = concat(rf_q, axis=0, ignore_index=True)
    X = _read_file(fn, lambda f: _read_fb(f["fb.tsv"]))
    X = concat(X, axis=0, ignore_index=True)
    loci = X.loc[:, ["chromosome", "physical_position"]]
    X = X.iloc[:, 4:].values
    admix = _subset_populations(X, len(pops))


def is_available():
    return False


def _read_file(fn: List[str], read_func: Callable) -> List:
    """
    Read data from multiple files using a provided read function.

    Parameters:
    ----------
    fn (List[str]): A list of file paths to read.
    read_func (Callable): A function to read data from each file.

    Returns:
    -------
    List: A list containing the data read from each file.
    """
    data = [];
    for file_name in fn:
        data.append(read_func(file_name))
    return data


def _read_fb(fn: str) -> DataFrame:
    header = odict(_types(fn, False))
    return _read_tsv(fn, header)


def _read_Q(fn: str) -> DataFrame:
    """
    Read the Q matrix from a file and add the chromosome information.

    Parameters:
    ----------
    fn (str): The file path of the Q matrix file.

    Returns:
    -------
    DataFrame: The Q matrix with the chromosome information added.
    """
    from re import search
    header = odict(_types(fn))
    df = _read_csv(fn, header)
    match = search(r'chr(\d+)', fn)
    if match:
        chrom = match.group(0)
        df["chrom"] = chrom
    else:
        print(f"Warning: Could not extract chromosome information from '{fn}'")
    return df


def _read_csv(fn: str, header: dict) -> DataFrame:
    """
    Read a CSV file into a pandas DataFrame with specified data types.

    Parameters:
    ----------
    fn (str): The file path of the CSV file.
    header (dict): A dictionary mapping column names to data types.

    Returns:
    -------
    DataFrame: The data read from the CSV file as a pandas DataFrame.
    """
    if is_available():
        df = read_csv(fn,sep="\t",header=None,names=list(header.keys()),
                      dtype=header,comment="#")
    else:
        df = read_csv(fn,delim_whitespace=True,header=None,
                      names=list(header.keys()),dtype=header,comment="#",
                      compression=None, engine="c",iterator=False)
    return df


def _read_tsv(fn: str, header: dict) -> DataFrame:
    """
    Read a TSV file into a pandas DataFrame.

    Parameters:
    ----------
    fn (str): File name of the TSV file.

    Returns:
    -------
    DataFrame: DataFrame containing specified columns from the TSV file.
    """
    if is_available():
        df = read_csv(fn,sep="\t",header=0, names=list(header.keys()),
                      dtype=header,comment="#")
    else:
        chunks = read_csv(fn,delim_whitespace=True, header=0,
                          names=list(header.keys()),dtype=header,comment="#",
                          chunksize=100000)
        # Concatenate chunks into single DataFrame
        df = concat(chunks, ignore_index=True)
    return df


def _types(fn: str, Q: bool = True) -> dict:
    """
    Infer the data types of columns in a TSV file.

    Parameters:
    ----------
    fn (str) : File name of the TSV file.

    Returns:
    -------
    dict : Dictionary mapping column names to their inferred data types.
    """
    from pandas import StringDtype
    if is_available():
        df = read_csv(fn,sep="\t",nrows=2,skiprows=1)
    else:
        df = read_csv(fn,delim_whitespace=True,nrows=2,skiprows=1)
    if Q:
        # Initialize the header dictionary with the sample_id column
        header = {"sample_id": StringDtype()}
        # Update the header dictionary with the data types of the remaining columns
        header.update(df.dtypes[1:].to_dict())
    else:
        header = df.dtypes.to_dict()
    return header


def _subset_populations(X, npops):
    """
    Subset and process the input array X based on populations.

    Parameters:
    X (ndarray): Input array where columns represent data for different populations.
    npops (int): Number of populations for column processing.

    Returns:
    ndarray: Processed array with adjacent columns summed for each population subset.
    """
    from numpy import concatenate
    pop_subset = []
    pop_start = 0
    ncols = X.shape[1]
    if ncols % npops != 0:
        raise ValueError("The number of columns in X must be divisible by npops.")
    while pop_start < npops:
        X0 = X[:, pop_start::npops] # Subset based on populations
        if X0.shape[1] % 2 != 0:
            raise ValueError("Number of columns must be even.")
        X0_summed = X0[:, ::2] + X0[:, 1::2] # Sum adjacent columns
        pop_subset.append(X0_summed)
        pop_start += 1
    return concatenate(pop_subset, 1)


def main():
    prefix_path = here("input/simulations/three_populations/_m/")
    start_time = time()
    _ = read_data(prefix_path)
    end_time = time()
    print(f"Execution time: {end_time - start_time} seconds")


if __name__ == "__main__":
    main()
