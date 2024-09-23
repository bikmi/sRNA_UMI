import argparse
import time
import os

from utilis import *


def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description="UMI finder", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-1", "--read1", required=True, help="Read1 FastQ file", type=str)
    parser.add_argument("-2", "--read2", required=True, help="Read2 FastQ file", type=str)
    parser.add_argument("-l", "--umi_length", required=False, help="Length of UMI sequences",
                        type=int, default=12)
    parser.add_argument("-c", "--count_table", required=True, help="Table of UMI count of each miRNA", type=str)
    parser.add_argument("-o", "--output_fasta", required=True, help="Final output miRNA fasta file", type=str)
    parser.add_argument("--min_length", required=False, help="Minimum length of miRNA sequences", type=int, default=15)
    parser.add_argument("--max_length", required=False, help="Maximum length of miRNA sequences", type=int, default=55)
    parser.add_argument("-n", "--mismatch_base", required=False,
                        help="Number of the mismatch base in common sequences", type=int, default=0)
    parser.add_argument("-s", "--sum_mismatch", required=False,
                        help="Sum of the mismatch base in common sequences", type=int, default=0)
    parser.add_argument("-fs", "--front_window_size", required=False, help="Size of front slide window", type=int,
                        default=1)
    parser.add_argument("-fq", "--front_mean_quality", required=False,
                        help="Average base quality of front slide window", type=int, default=20)
    parser.add_argument("-ts", "--tail_window_size", required=False, help="Size of tail slide window", type=int,
                        default=1)
    parser.add_argument("-tq", "--tail_mean_quality", required=False,
                        help="Average base quality of tail slide window", type=int, default=20)
    parser.add_argument("-aq", "--avg_quality", required=False, help="Average base quality of slide window", type=int,
                        default=15)
    parser.add_argument("-N", "--max_N", required=False, help="Maximum N base in read", type=int, default=0)
    parser.add_argument("-q", "--base_quality", required=False, help="Minimum base quality", type=int, default=20)
    parser.add_argument("-u", "--unqualified_percent_limit", required=False,
                        help="Percentage of unqualified base in read", type=float, default=2)



    args = parser.parse_args()

    read_1 = args.read1
    read_2 = args.read2
    count_table = args.count_table
    umi_length = args.umi_length
    min_len = args.min_length
    max_len = args.max_length
    average_quality = args.avg_quality
    front_window_size = args.front_window_size
    front_mean_quality = args.front_mean_quality
    tail_window_size = args.tail_window_size
    tail_mean_quality = args.tail_mean_quality
    ambiguous_base = args.max_N
    base_quality = args.base_quality
    low_base_percentage = args.unqualified_percent_limit
    mismatch_base = args.mismatch_base
    sum_mismatch = args.sum_mismatch
    output_fasta = args.output_fasta

    read_1 = os.path.abspath(read_1)
    read_2 = os.path.abspath(read_2)

    sample_name = os.path.basename(read_1)
    file_prefix = os.path.splitext(sample_name)[0].split('_')[0]
    directory_path = os.path.join(os.path.dirname(read_1), file_prefix)

    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

    os.chdir(directory_path)

    read1_decompress = decompress_files(read_1, directory_path)
    read2_decompress = decompress_files(read_2, directory_path)

    read1_clean, read2_clean = paired_end_quality_control(read1_decompress, read2_decompress,
                                                          file_prefix + '.R1.clean.fq', file_prefix + '.R2.clean.fq',
                                                          average_quality, front_window_size, front_mean_quality,
                                                          tail_window_size, tail_mean_quality, ambiguous_base,
                                                          base_quality,
                                                          low_base_percentage)
    # os.remove(read1_decompress)
    # os.remove(read2_decompress)

    read1_fasta = fq2fa(read1_clean)
    read2_fasta = fq2fa(read2_clean)

    # os.remove(read1_clean)
    # os.remove(read2_clean)

    merge_file = merge_seqs(read1_fasta, read2_fasta)

    # os.remove(read1_fasta)
    # os.remove(read2_fasta)

    umi_fasta1, output_fasta1 = find_umi(merge_file, umi_length, min_len, max_len, mismatch_base, sum_mismatch)

    umi_table = umi_cluster(umi_fasta1, output_fasta1)

    dereplicated_table = remove_error_seqs(umi_table)

    count_umi(dereplicated_table, count_table, output_fasta)

    end_time = time.time()

    execution_time = end_time - start_time

    minutes = int(execution_time // 60)
    seconds = int(execution_time % 60)

    print(f"Processing completed in {minutes}m{seconds}s")


if __name__ == "__main__":
    main()
