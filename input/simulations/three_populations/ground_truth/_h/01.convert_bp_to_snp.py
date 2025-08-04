import sys
import gzip
import pysam

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


def iter_snp_positions_vcf(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    for record in vcf.fetch():
        yield record.pos



def get_ancestry_positions(snp_pos, blocks):
    block_idx = 0
    while block_idx < len(blocks) and snp_pos > blocks[block_idx][0]:
        block_idx += 1
    if block_idx == 0:
        return blocks[0][1]
    else:
        return blocks[block_idx - 1][1]


def generate_ancestry_table(bp_file, vcf_file, output_file):
    """
    Assigns ancestry to each SNP based on the block breakpoint structure.
    """
    ancestry_blocks = parse_bp_file(bp_file)
    haplotypes = list(ancestry_blocks.keys())
    # Open output file
    open_fn = gzip.open if output_file.endswith(".gz") else open
    with open_fn(output_file, "wt") as out_f:
        out_f.write("Position\t" + "\t".join(haplotypes) + "\n")
        for snp_pos in iter_snp_positions_vcf(vcf_file):
            row = [str(snp_pos)]
            for hap in haplotypes:
                ancestry = get_ancestry_positions(snp_pos, ancestry_blocks[hap])
                row.append(ancestry)
            out_f.write("\t".join(row) + "\n")


def main(bp_file, vcf_file, output_file):
    generate_ancestry_table(bp_file, vcf_file, output_file)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python 01.convert_bp_to_snp.py <bp_file> <vcf_file> <output_file>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])
