import session_info
import pandas as pd
from json import dumps
import sys, re, gzip, csv
from collections import defaultdict

def group_haplotypes(columns):
    """
    Group haplotype columns into pairs per individual.
    Returns a dict: {Sample_1: ('Sample_1_1', 'Sample_1_2'), ...}
    """
    sample_pairs = defaultdict(list)
    for col in columns:
        match = re.match(r"(Sample_\d+)_\d+", col)
        if match:
            sample_name = match.group(1)
            sample_pairs[sample_name].append(col)
    # Ensure we only include samples with exactly 2 haplotypes
    return {k: v for k, v in sample_pairs.items() if len(v) == 2}


def get_ref_ancestry_and_map(filename, nrows=1000):
    sample_df = pd.read_csv(filename, sep="\t", nrows=nrows, index_col=0)
    ref_ancestry = sample_df.iloc[0].value_counts().idxmax()
    all_ancestries = sorted(set(sample_df.values.ravel()))
    return ref_ancestry, {anc: i for i, anc in enumerate(all_ancestries)}


def stream_encode(input_file, output_file, mode='auto', chunksize=1000,
                  chrom="chr1"):
    # Open output handle
    open_fn = gzip.open if output_file.endswith('.gz') else open
    with open_fn(output_file, 'wt', newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')

        # Grab column names and ancestry map
        sample_df = pd.read_csv(input_file, sep="\t", nrows=1)
        hap_cols = sample_df.columns[1:]  # skip index
        grouped = group_haplotypes(hap_cols)
        if not grouped:
            raise ValueError("No haplotype pairs found in input.")

        # Determine reference ancestry and mapping
        _, ancestry_map = get_ref_ancestry_and_map(input_file)
        ancestries = sorted(ancestry_map.keys())

        print(f"Ancestry map: {ancestry_map}")
        print(f"Encoding format: presence/absence per ancestry: {'|'.join(ancestries)}")

        # Write metadata comments
        out_f.write(f"##ancestry_mapping={dumps(ancestries)}\n")

        # Write header
        writer.writerow(['#CHROM', 'POS'] + list(grouped.keys()))

        # Stream SNP rows
        for chunk in pd.read_csv(input_file, sep="\t", index_col=0,
                                 chunksize=chunksize):
            for pos, row in chunk.iterrows():
                encoded_row = [chrom, pos]
                for sample, (h1, h2) in grouped.items():
                    a1 = row[h1]
                    a2 = row[h2]
                    present = {a1, a2}
                    encoded_vec = [str(int(a in present)) for a in ancestries]
                    encoded_row.append("|".join(encoded_vec))
                writer.writerow(encoded_row)


def main(input_file, output_file, mode, chunksize, chrom):
    stream_encode(input_file, output_file, mode, chunksize, chrom)


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python 01.convert_hap_to_geno.py <chrom> <input_tsv[gz]> <output_tsv[gz]> [biallelic|multiallelic|auto] [chunksize]")
        sys.exit(1)

    chrom = f"chr{sys.argv[1]}"
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    mode = sys.argv[4] if len(sys.argv) >= 5 else 'auto'
    chunksize = int(sys.argv[5]) if len(sys.argv) >= 6 else 1000

    main(input_file, output_file, mode, chunksize, chrom)
