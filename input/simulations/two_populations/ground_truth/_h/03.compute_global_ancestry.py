# --------------------------
# Load Libraries
# --------------------------

# Standard libraries
import os
import session_info
from pyhere import here
from pathlib import Path

# Project specific
import gzip
import pandas as pd
from typing import List, Tuple
import ast
import argparse
import sys


# --------------------------
# Parsing the file data
# --------------------------

def open_vcf_gz(path: str):
    """
    Open .vcf.gz file types.
    """
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'r')

def parse_header_and_ancestries(path) -> Tuple[dict, List[str], List[str]]:
    """
    Extract ancestries and sample names from the VCF header.

    Returns:
        ancestries: list of ancestry names (if available)
        sample_names: list of sample names from #CHROM header
    """
    ancestries = None
    sample_names = None
    with open_vcf_gz(path) as f:
        for line in f:
            if line.startswith('##'):
                if 'Ancestries' in line:
                    try:
                        ancestries = ast.literal_eval(line.split(':', 1)[1].strip())
                    except Exception:
                        ancestries = None             
                continue
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                try:
                    start_idx = parts.index('Sample_1')
                except ValueError:
                    # fallback: use previous default slicing
                    start_idx = 2
                sample_names = parts[start_idx:]
                break
    return ancestries, sample_names

def parse_ancestry_token(token: str) -> List[int]:
    """
    Parse a sample token like '0|2' or '0|2:otherinfo' -> [0, 2]
    Returns a list of ints (one per ancestry).
    """
    if token is None or token == '.':
        return []
    token = token.split(':', 1)[0]
    parts = token.split('|')
    try:
        return [int(x) for x in parts]
    except ValueError:
        # if parsing fails, return empty to indicate missing data
        return []

# --------------------------
# Processing files and
# Computing global ancestry
# --------------------------

def process_file(path: Path,
                 ancestries: List[str],
                 weight_intervals: bool,
                 sample_sums: dict,
                 sample_total_weight: dict):
    """Process one input file and update sample_sums and sample_total_weight in-place.
    
    Also prints example numeric counts for the first variant row (for visual inspection).
    """
    num_ances = len(ancestries)
    
    with open_vcf_gz(path) as fh:
        # === PARSE HEADER ===
        sample_names = None
        sample_start_idx = 2
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                try:
                    sample_start_idx = parts.index('Sample_1')
                except ValueError:
                    sample_start_idx = 2
                sample_names = parts[sample_start_idx:]
                break
        if sample_names is None:
            raise RuntimeError(f"No #CHROM header found in {path}")

        # === INITIALIZE SAMPLE SUMS ===
        for s in sample_names:
            if s not in sample_sums:
                sample_sums[s] = [0.0] * num_ances
                sample_total_weight[s] = 0.0

        prev_pos = None
        prev_counts = {s: None for s in sample_names}

        # === ITERATE VARIANT LINES ===
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            pos = int(parts[1])
            sample_fields = parts[sample_start_idx : sample_start_idx + len(sample_names)]

            # Parse ancestry tokens
            parsed = [parse_ancestry_token(tok) for tok in sample_fields]

            # Skip rows for samples with missing calls (do not add denominator)
            # We assume correct length (num_ances) for non-missing tokens
            for i, arr in enumerate(parsed):
                if not arr:  # Missing or invalid token
                    parsed[i] = None

            # === ADD COUNTS ===
            if not weight_intervals:
                for sname, counts in zip(sample_names, parsed):
                    if counts is None:
                        continue  # skip missing data
                    for a in range(num_ances):
                        sample_sums[sname][a] += counts[a]
                    sample_total_weight[sname] += 2.0
                continue

            # Interval weighting
            if prev_pos is None:
                prev_pos = pos
                for sname, counts in zip(sample_names, parsed):
                    prev_counts[sname] = counts
                continue

            interval_len = pos - prev_pos
            if interval_len < 0:
                print(f"Warning: non-positive interval (length={interval_len}) at pos {pos} in {path.name}", file=sys.stderr)
                interval_len = 0
            
            for sname in sample_names:
                counts = prev_counts[sname]
                if counts is None:
                    continue
                for a in range(num_ances):
                    sample_sums[sname][a] += counts[a] * interval_len
                sample_total_weight[sname] += 2.0 * interval_len

            prev_pos = pos
            for sname, counts in zip(sample_names, parsed):
                prev_counts[sname] = counts

def compute_proportions(
    sample_sums: Dict[str, List[float]],
    sample_total_weight: Dict[str, float],
    ancestries: List[str]
    ) -> Dict[str, Dict[str, float]]:

    # Convert sample_sums into DataFrame: rows=samples, columns=ancestries
    df_sums = pd.DataFrame.from_dict(sample_sums, orient='index', columns=ancestries)

    # Convert total weights into Series
    total_weight_series = pd.Series(sample_total_weight)

    # Divide each ancestry count by total weight per sample
    df_props = df_sums.div(total_weight_series, axis=0)

    # Convert back to nested dict: sample -> {ancestry: proportion}
    return df_props.to_dict(orient='index')

def write_tsv(output_path: str, proportions: dict, ancestries: List[str]):
    with open(output_path, 'w') as out:
        header = ['Sample'] + ancestries
        out.write('\t'.join(header) + '\n')
        for s in sorted(proportions.keys()):
            row = [s] + [f"{proportions[s][a]:.6f}" for a in ancestries]
            out.write('\t'.join(row) + '\n')

            
# --------------------------
# MAIN
# --------------------------
    
    
def main():
    p = argparse.ArgumentParser(
        description='Compute global ancestry proportions from local ancestry VCF-like files'
    )
    p.add_argument(
        'inputs',
        nargs='+',
        help='Input per-chromosome files (vcf or vcf.gz).'
    )
    p.add_argument(
        '--ancestries',
        nargs='+',
        help='Ordered list of ancestry labels (e.g. CEU YRI). If omitted, tries to read from first file header.'
    )
    p.add_argument(
        '--weight',
        action='store_true',
        help='Weight calls by interval (pos[i+1]-pos[i]) instead of simple site counts.'
    )
    p.add_argument(
        '--out',
        default='global_ancestry.tsv',
        help='Output TSV file'
    )
    args = p.parse_args()

    # -------------------------
    # Prepare containers
    # -------------------------
    sample_sums = {}          # sample -> [sum_per_ancestry]
    sample_total_weight = {}  # sample -> total_weight (allele-units)

    # -------------------------
    # Determine ancestries
    # -------------------------
    ancestries = args.ancestries
    if ancestries is None:
        # try to read from first file
        first_file = args.inputs[0]
        file_ancestries, sample_names = parse_header_and_ancestries(first_file)
        if file_ancestries is None:
            p.error('No ancestries provided and none found in file header. Pass --ancestries.')
        ancestries = file_ancestries

    print(f"Using ancestries: {ancestries}", file=sys.stderr)

    # -------------------------
    # Process each per-chromosome file
    # -------------------------
    for path in args.inputs:
        print(f"Processing {path} ...", file=sys.stderr)
        process_file(
            path=path,
            ancestries=ancestries,
            weight_intervals=args.weight,
            sample_sums=sample_sums,
            sample_total_weight=sample_total_weight
        )

    # -------------------------
    # Compute proportions
    # -------------------------
    proportions = compute_proportions(sample_sums, sample_total_weight, ancestries)

    # -------------------------
    # Write output
    # -------------------------
    write_tsv(args.out, proportions, ancestries)
    print(f"Wrote global ancestry table to {args.out}", file=sys.stderr)

session_info.show()
    
if __name__ == '__main__':
    main()

if __name__ == '__main__':
    main()