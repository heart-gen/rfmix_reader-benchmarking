import sys
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


def assign_ancestry_to_snps(snp_positions, ancestry_blocks):
    """
    Assigns ancestry to each SNP based on the block breakpoint structure.
    """
    results = pd.DataFrame(index=snp_positions)
    for sample, blocks in ancestry_blocks.items():
        ancestry_per_snp = []
        block_idx = 0
        for snp in snp_positions:
            # Move to the correct block
            while block_idx < len(blocks) and snp > blocks[block_idx][0]:
                block_idx += 1
            # Assign ancestry of the *previous* block
            if block_idx == 0:
                ancestry = blocks[0][1]
            else:
                ancestry = blocks[block_idx - 1][1]
            ancestry_per_snp.append(ancestry)
        results[sample] = ancestry_per_snp
    return results


def main(bp_file, snp_file, output_file):
    # Read SNP positions file (1 column, no header)
    snp_positions = pd.read_csv(snp_file, header=None)[0].tolist()

    ancestry_blocks = parse_bp_file(bp_file)
    per_snp_ancestry = assign_ancestry_to_snps(snp_positions, ancestry_blocks)
    per_snp_ancestry.index.name = 'Position'
    per_snp_ancestry.to_csv(output_file, sep='\t')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python 01.convert_bp_to_snp.py <bp_file> <snp_file> <output_file>")
        sys.exit(1)
    bp_file = sys.argv[1]
    snp_file = sys.argv[2]
    output_file = sys.argv[3]
    main(bp_file, snp_file, output_file)
