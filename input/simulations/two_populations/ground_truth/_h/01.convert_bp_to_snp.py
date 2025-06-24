import sys
import pysam
import pandas as pd

def parse_bp_file(bp_file):
    ancestry_blocks = {}
    current_sample = None
    with open(bp_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("Sample"):
                current_sample = line
                ancestry_blocks[current_sample] = []
            else:
                ancestry, chrom, bp, cm = line.split()
                ancestry_blocks[current_sample].append((int(bp), ancestry))
    return ancestry_blocks


def get_snp_positions_from_vcf(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    positions = [records.pos for records in vcf.fetch()]
    return sorted(positions)


def assign_ancestry_to_snps(snp_positions, ancestry_blocks):
    """
    Assigns ancestry to each SNP based on the block breakpoint structure.
    """
    ancestry_matrix = {}
    for sample_hap, blocks in ancestry_blocks.items():
        ancestry_per_snp = []
        block_idx = 0
        for snp in snp_positions:
            while block_idx < len(blocks) and snp > blocks[block_idx][0]:
                block_idx += 1
            if block_idx == 0:
                ancestry = blocks[0][1]
            else:
                ancestry = blocks[block_idx - 1][1]
            ancestry_per_snp.append(ancestry)
        ancestry_matrix[sample_hap] = ancestry_per_snp
    results_df = pd.DataFrame.from_dict(ancestry_matrix, orient="columns")
    results_df.index = snp_positions
    results_df.index.name = "Position"
    return results_df


def main(bp_file, snp_file, output_file):
    # Read SNP positions file (1 column, no header)
    snp_positions = get_snp_positions_from_vcf(snp_file)
    ancestry_blocks = parse_bp_file(bp_file)
    per_snp_ancestry = assign_ancestry_to_snps(snp_positions, ancestry_blocks)
    per_snp_ancestry.to_csv(output_file, sep='\t')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python 01.convert_bp_to_snp.py <bp_file> <snp_file> <output_file>")
        sys.exit(1)
    bp_file = sys.argv[1]
    snp_file = sys.argv[2]
    output_file = sys.argv[3]
    main(bp_file, snp_file, output_file)
