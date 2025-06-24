import sys
import pandas as pd
import re
from collections import defaultdict

def group_haplotypes(columns):
    """
    Group haplotype columns into pairs per individual.
    Assumes format: Sample_1_H1, Sample_1_H2, etc.
    Returns a dict: {Sample_1: ('Sample_1_H1', 'Sample_1_H2'), ...}
    """
    sample_pairs = defaultdict(list)
    for col in columns:
        match = re.match(r"(Sample_\d+)_H\d+", col)
        if match:
            sample_name = match.group(1)
            sample_pairs[sample_name].append(col)
    # Ensure we only include samples with exactly 2 haplotypes
    return {k: v for k, v in sample_pairs.items() if len(v) == 2}


def encode_biallelic_genotypes(df, ref_ancestry=None):
    """
    Convert haplotype ancestry calls to 0/1/2 per individual.
    ref_ancestry defines which ancestry is coded as 0.
    """
    grouped = group_haplotypes(df.columns)
    out_df = pd.DataFrame(index=df.index)
    if not ref_ancestry:
        # Automatically pick most frequent ancestry at the first position
        ref_ancestry = df.iloc[0].value_counts().idxmax()
    print(f"Using reference ancestry: {ref_ancestry}")
    for sample, (hap1, hap2) in grouped.items():
        def encode_pair(a1, a2):
            if a1 == a2:
                return 0 if a1 == ref_ancestry else 2
            else:
                return 1

        encoded = [
            encode_pair(h1, h2)
            for h1, h2 in zip(df[hap1], df[hap2])
        ]
        out_df[sample] = encoded
    return out_df


def encode_multiallelic(df):
    """
    Encodes multiallelic ancestries as categorical integers.
    Returns: DataFrame where each ancestry is encoded as an int sum per individual.
    """
    grouped = group_haplotypes(df.columns)
    out_df = pd.DataFrame(index=df.index)
    # Build consistent ancestry-to-int map
    all_ancestries = sorted(set(df.values.ravel()))
    ancestry_map = {anc: i for i, anc in enumerate(all_ancestries)}
    print(f"Ancestry map: {ancestry_map}")
    for sample, (hap1, hap2) in grouped.items():
        encoded = [
            ancestry_map[a1] + ancestry_map[a2]
            for a1, a2 in zip(df[hap1], df[hap2])
        ]
        out_df[sample] = encoded
    return out_df, ancestry_map


def main(input_file, output_file, mode='auto'):
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    ancestries = sorted(set(df.values.ravel()))
    print(f"Detected ancestries: {ancestries}")
    if mode == 'biallelic' or (mode == 'auto' and len(ancestries


if __name__ == "__main__":
    main()
