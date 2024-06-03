import argparse
import os
import shutil
import gzip
import collections

from utilis import *


def main():
    parser = argparse.ArgumentParser(description="isoforms finder")
    # parser.add_argument("-l", "--long_read", required=True, help="Use long read sequencing")
    parser.add_argument("-1", "--read1", required=True, help="Path to the first input FastQ file")
    parser.add_argument("-2", "--read2", required=True, help="Path to the second input FastQ file")
    parser.add_argument("-m", "--merge", required=True, help="Path to the second input FastQ file")

    parser.add_argument("-u", "--umi", required=True, help="Path to the reference genome file")
    parser.add_argument("-c", "--count_table", required=True, help="medaka module")
    # parser.add_argument("-t", "--num_threads", type=int, default=8, help="Number of threads for processing")
    parser.add_argument("-o", "--output", required=True, help="Output directory for results")

    args = parser.parse_args()

    read_1 = args.read1
    read_2 = args.read2
    merge_fasta = args.merge
    umi_fasta = args.umi
    output_fasta = args.output
    count_table = args.count_table
    # files = os.listdir()

    read1 = decompress_files(read_1)
    read2 = decompress_files(read_2)

    read1_fasta = fq2fa(read1)
    read2_fasta = fq2fa(read2)

    merge_file = merge_seqs(read1_fasta, read2_fasta, merge_fasta)

    umi_fasta, output_fasta = find_UMI(merge_file, umi_fasta, output_fasta)

    umi_cluster(umi_fasta, output_fasta, count_table)


if __name__ == "__main__":
    main()
