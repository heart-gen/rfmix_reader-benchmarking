import numpy as np
import session_info
import pandas as pd
from json import dumps
from numba import njit
import sys, re, gzip, csv
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


@njit
def encode_pairwise_numba(a1, a2, n_ancestries):
    counts = np.zeros(n_ancestries, dtype=np.uint8)
    counts[a1] += 1
    counts[a2] += 1
    return counts


def stream_encode(input_file, output_file, chunksize=10000, chrom="chr1"):
    open_fn = gzip.open if output_file.endswith('.gz') else open
    with open_fn(output_file, 'wt', newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')

        # Read only column headers and determine haplotype groups
        with (gzip.open if input_file.endswith(".gz") else open)(input_file, 'rt') as f:
            header = f.readline().strip().split('\t')[1:]
        grouped = group_haplotypes(header)
        if not grouped:
            raise ValueError("No haplotype pairs found in input.")

        # Efficient ancestry map estimation
        _, ancestry_map = get_ref_ancestry_and_map(input_file, nrows=1000)
        ancestries = list(sorted(ancestry_map.keys()))
        n_ancestries = len(ancestries)

        print(f"Ancestries: {ancestry_map.keys()}")
        print(f"Encoding: {n_ancestries} ancestries ({'|'.join(ancestries)})")

        # Write metadata and header
        out_f.write(f"##Ancestries: {dumps(ancestries)}\n")
        out_f.write(f"##Encoding locations: {dumps('|'.join(ancestries))}\n")
        out_f.write("##Legend: 0 - absent; 1 - one allele; 2 - both alleles\n")
        writer.writerow(['#CHROM', 'POS'] + list(grouped.keys()))

        # Define dtypes
        hap_cols = [col for pair in grouped.values() for col in pair]
        dtype_dict = {col: 'category' for col in hap_cols}

        # Stream data
        for chunk in pd.read_csv(input_file, sep="\t", index_col=0,
                                 usecols=["Position"] + hap_cols,
                                 chunksize=chunksize, dtype=dtype_dict):
            # Convert ancestries to integer codes
            chunk_int = chunk.apply(lambda col: col.astype('category').cat.codes)
            for pos, row in chunk_int.iterrows():
                encoded_row = [chrom, pos]
                for sample, (h1, h2) in grouped.items():
                    a1 = row[h1]
                    a2 = row[h2]
                    counts = encode_pairwise_numba(a1, a2, n_ancestries)
                    encoded = "|".join(str(c) for c in counts)
                    encoded_row.append(encoded)
                writer.writerow(encoded_row)


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
