import numpy as np
import session_info
import pandas as pd
from json import dumps
import sys, re, gzip, csv
from io import TextIOWrapper
from numba import njit, prange
from collections import defaultdict

def group_haplotypes(columns):
    """
    Group haplotype columns into pairs per individual.
    Returns a dict: {Sample_1: ('Sample_1_1', 'Sample_1_2'), ...}
    """
    sample_pairs = defaultdict(list)
    for col in columns:
        m = re.match(r"(Sample_\d+)_\d+", col)
        if m:
            sample_name = m.group(1)
            sample_pairs[sample_name].append(col)
    return {k: v for k, v in sample_pairs.items() if len(v) == 2}


def get_ref_ancestry_and_map(filename, nrows=1000):
    sample_df = pd.read_csv(filename, sep="\t", nrows=nrows, index_col=0)
    ref_ancestry = sample_df.iloc[0].value_counts().idxmax()
    all_ancestries = sorted(set(sample_df.values.ravel()))
    return ref_ancestry, {anc: i for i, anc in enumerate(all_ancestries)}


@njit(parallel=True)
def batch_encode_pairwise(a1, a2, n_ancestries):
    n_rows, n_samples = a1.shape
    result = np.zeros((n_rows, n_samples, n_ancestries), dtype=np.uint8)
    for i in prange(n_rows):
        for j in range(n_samples):
            counts = np.zeros(n_ancestries, dtype=np.uint8)
            counts[a1[i, j]] += 1
            counts[a2[i, j]] += 1
            result[i, j] = counts
    return result


def stream_encode(input_file, output_file, chunksize=10000, chrom="chr1"):
    open_fn = gzip.open if output_file.endswith('.gz') else open
    with open_fn(output_file, 'wt', newline='') as raw_out:
        out_f = TextIOWrapper(raw_out) if isinstance(raw_out, gzip.GzipFile) else raw_out
        writer = csv.writer(out_f, delimiter='\t')

        # Read only column headers and determine haplotype groups
        with (gzip.open if input_file.endswith(".gz") else open)(input_file, 'rt') as f:
            header = f.readline().strip().split('\t')[1:]
        grouped = group_haplotypes(header)
        if not grouped:
            raise ValueError("No haplotype pairs found in input.")

        # Efficient ancestry map estimation
        _, ancestry_map = get_ref_ancestry_and_map(input_file)
        ancestries = list(sorted(ancestry_map.keys()))
        n_ancestries = len(ancestries)

        # Write metadata and header
        out_f.write(f"##Ancestries: {dumps(ancestries)}\n")
        out_f.write(f"##AncestryMap: {dumps(ancestry_map)}\n")
        writer.writerow(['#CHROM', 'POS'] + list(grouped.keys()))

        # Define dtypes
        hap_cols = [col for pair in grouped.values() for col in pair]
        dtype_dict = {col: 'category' for col in hap_cols}
        usecols = ["Position"] + hap_cols

        # Stream data
        for chunk in pd.read_csv(input_file, sep="\t", index_col=0,
                                 usecols=usecols, chunksize=chunksize,
                                 dtype=dtype_dict):

            # Ensure category ordering matches ancestry_map
            for col in hap_cols:
                chunk[col] = chunk[col].astype(pd.CategoricalDtype(categories=ancestries,
                                                                   ordered=False))
                chunk[col] = chunk[col].cat.codes.astype(np.uint8)

            # Create arrays for both haplotype sets
            a1_array = chunk[[pair[0] for pair in grouped.values()]].to_numpy()
            a2_array = chunk[[pair[1] for pair in grouped.values()]].to_numpy()

            encoded = batch_encode_pairwise(a1_array, a2_array, n_ancestries)
            pos_array = chunk.index.to_numpy()

            rows_out = []
            for i in range(encoded.shape[0]):
                encoded_row = [chrom, pos_array[i]]
                encoded_row += ["|".join(map(str, encoded[i, j]))
                                for j in range(encoded.shape[1])]
                rows_out.append(encoded_row)

            writer.writerows(rows_out)


def main(input_file, output_file, chunksize, chrom):
    stream_encode(input_file, output_file, chunksize, chrom)


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python 01.convert_hap_to_geno.py <chrom> <input_tsv[gz]> <output_tsv[gz]> [chunksize]")
        sys.exit(1)

    chrom = f"chr{sys.argv[1]}"
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    chunksize = int(sys.argv[4]) if len(sys.argv) >= 5 else 10000

    main(input_file, output_file, chunksize, chrom)
