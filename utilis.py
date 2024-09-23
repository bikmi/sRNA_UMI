import collections
import shutil
import gzip
import os


def decompress_files(gzip_file, subfolder_path):
    if gzip_file.endswith('.gz'):
        # decompress gzip files
        output_file = os.path.join(subfolder_path, os.path.basename(gzip_file).replace('.gz', ''))
        with gzip.open(gzip_file, 'rb') as gf, open(output_file, 'wb') as ugf:
            shutil.copyfileobj(gf, ugf)
        return output_file


def fq2fa(fastq):
    # convert fastq to fasta file
    fasta_file = fastq.replace('.fq', '.fa')
    with open(fastq, 'r') as fq, open(fasta_file, 'w') as fa:
        for lineID, line in enumerate(fq, 1):
            if lineID % 4 == 1:
                sequence_id = line
                fa.write(sequence_id.replace('@', '>'))
            elif lineID % 4 == 2:
                sequence = line
                fa.write(sequence)

    return fasta_file


def reverse_complement(dna_sequence):  # reverse complement
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(dna_sequence))


def merge_seqs(reads1, reads2, merge_reads='merge.fasta'):
    """Merge sequences from two FastA files based on matching criteria."""
    with open(reads1, 'r') as r1, open(reads2, 'r') as r2, open(merge_reads, 'w+') as r3:
        while True:
            read1_line1 = r1.readline().rstrip()
            read1_line2 = r1.readline().rstrip()
            read2_line1 = r2.readline().rstrip()
            read2_line2 = r2.readline().rstrip()
            read2_line2_reverse_complete = reverse_complement(read2_line2)
            read_name = read1_line1.split(' ')[0]
            if not (read1_line1 and read1_line2) or not (read2_line1 and read2_line2):
                break

            max_overlap_len = min(len(read1_line2), len(read2_line2_reverse_complete))
            for i in range(max_overlap_len, 0, -1):
                if read2_line2_reverse_complete[-i:] == read1_line2[:i]:
                    if len(read1_line2[-i:]) >= 10:
                        merged_read = read2_line2_reverse_complete + read1_line2[i:]
                        r3.write(read_name + '\n')
                        r3.write(merged_read + '\n')

    return merge_reads


# def find_umi(fasta, umi_length, min_len, max_len, umi_fasta='UMI.fasta', out_put='UMI_miRNA.fasta'):
#     str1 = 'AACTGTAGGCACCATCAAT'
#     str2 = 'AGATCGGAAGAGCACACGTCT'
#     str3 = 'GTTCAGAGTTCTACAGTCCGACGATC'
#
#     with open(fasta, 'r') as fa, open(umi_fasta, 'w+') as umi, open(out_put, 'w+') as o:
#         while True:
#             line1 = fa.readline().rstrip()
#             line2 = fa.readline().rstrip()
#             if not line1:
#                 break
#
#             if str1 in line2 and str2 in line2 and str3 in line2:
#                 str1_index = line2.find(str1)
#                 str2_index = line2.find(str2)
#                 str3_index = line2.find(str3)
#                 str1_end_index = str1_index + len(str1) - 1
#                 str3_end_index = str3_index + len(str3) - 1
#
#                 umi_seq = line2[str1_end_index + 1: str2_index]
#
#                 if abs(len(umi_seq) - 12) <= abs(umi_length - 12):
#                     umi.write(line1 + '\n')
#                     umi.write(umi_seq + '\n')
#
#                 insert_seq = line2[str3_end_index + 1: str1_index]
#
#                 if max_len >= len(insert_seq) >= min_len:
#                     o.write(line1 + ' ' + umi_seq + '\n')
#                     o.write(insert_seq + '\n')
#
#     return umi_fasta, out_put

def hamming_distance(s1, s2):
    """Calculate the Hamming distance between two strings of equal length."""
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


def find_with_mismatch(line, pattern, max_mismatch):
    """Find the first occurrence of a pattern in a line allowing for a maximum number of mismatches."""
    for i in range(len(line) - len(pattern) + 1):
        if hamming_distance(line[i:i + len(pattern)], pattern) <= max_mismatch:
            return i
    return -1


def find_umi(fasta, umi_length, min_len, max_len, mismatch, sum_mismatch, umi_fasta='UMI.fasta',
             out_put='UMI_miRNA.fasta'):
    str1 = 'AACTGTAGGCACCATCAAT'
    str2 = 'AGATCGGAAGAGCACACGTCT'
    str3 = 'GTTCAGAGTTCTACAGTCCGACGATC'

    len_str1 = len(str1)
    len_str2 = len(str2)
    len_str3 = len(str3)

    if mismatch == 0:
        with open(fasta, 'r') as fa, open(umi_fasta, 'w+') as umi, open(out_put, 'w+') as o:
            while True:
                line1 = fa.readline().rstrip()
                line2 = fa.readline().rstrip()
                if not line1:
                    break

                if str1 in line2 and str2 in line2 and str3 in line2:
                    str1_index = line2.find(str1)
                    str2_index = line2.find(str2)
                    str3_index = line2.find(str3)
                    str1_end_index = str1_index + len(str1) - 1
                    str3_end_index = str3_index + len(str3) - 1

                    umi_seq = line2[str1_end_index + 1: str2_index]

                    if abs(len(umi_seq) - 12) <= abs(umi_length - 12):
                        umi.write(line1 + '\n')
                        umi.write(umi_seq + '\n')

                    insert_seq = line2[str3_end_index + 1: str1_index]

                    if max_len >= len(insert_seq) >= min_len:
                        o.write(line1 + ' ' + umi_seq + '\n')
                        o.write(insert_seq + '\n')

        return umi_fasta, out_put
    else:
        with open(fasta, 'r') as fa, open(umi_fasta, 'w+') as umi, open(out_put, 'w+') as o:
            while True:
                line1 = fa.readline().rstrip()
                line2 = fa.readline().rstrip()
                if not line1:
                    break

                # Calculate mismatches for each pattern
                str1_index = find_with_mismatch(line2, str1, mismatch)
                if str1_index != -1:
                    str1_mismatches = hamming_distance(line2[str1_index:str1_index + len_str1], str1)
                else:
                    continue

                str2_index = find_with_mismatch(line2, str2, mismatch)
                if str2_index != -1:
                    str2_mismatches = hamming_distance(line2[str2_index:str2_index + len_str2], str2)
                else:
                    continue

                str3_index = find_with_mismatch(line2, str3, mismatch)
                if str3_index != -1:
                    str3_mismatches = hamming_distance(line2[str3_index:str3_index + len_str3], str3)
                else:
                    continue

                # Ensure the total number of mismatches is within the allowed limit
                if (str1_mismatches + str2_mismatches + str3_mismatches) <= sum_mismatch:
                    str1_end_index = str1_index + len_str1 - 1
                    str3_end_index = str3_index + len_str3 - 1

                    umi_seq = line2[str1_end_index + 1: str2_index]

                    if abs(len(umi_seq) - 12) <= abs(umi_length - 12):
                        umi.write(line1 + '\n')
                        umi.write(umi_seq + '\n')

                    insert_seq = line2[str3_end_index + 1: str1_index]

                    if max_len >= len(insert_seq) >= min_len:
                        o.write(line1 + ' ' + umi_seq + '\n')
                        o.write(insert_seq + '\n')

        return umi_fasta, out_put


def umi_cluster(umi_fasta, seqs_fasta, output_file='cluster.tsv', seq_number=1):
    umi_count = collections.defaultdict(int)
    seqs_count = collections.defaultdict(int)

    with open(umi_fasta, 'r') as umi, open(seqs_fasta, 'r') as seqs:

        while True:  # 统计UMI的个数
            umi_name = umi.readline().rstrip()
            umi_read = umi.readline().rstrip()
            if len(umi_name) == 0 or len(umi_read) == 0:
                break
            umi_count[umi_read] += 1

        while True:  # 统计UMI下各个序列的个数
            seq_name = seqs.readline().rstrip()
            umi_read1 = seq_name.split(' ')[-1]
            seq_read = seqs.readline().rstrip()
            if len(seq_name) == 0 or len(seq_read) == 0:
                break

            if umi_read1 in umi_count and umi_count[umi_read1] >= seq_number:
                '''筛选出UMI序列的统计数量大于一定阈值的序列，当前的序列数量的阈值设置为1'''
                seqs_count[(umi_read1, seq_read)] += 1

    sorted_seqs = sorted(seqs_count.items(), key=lambda item: (item[0][0], item[1]))
    '''按照seqs_count中的umi_read1和数量进行排序'''

    with open(output_file, 'w') as out:
        for (umi_read, seq_read), count in sorted_seqs:
            out.write(f'{umi_read}\t{seq_read}\t{count}\n')

    return output_file


def remove_error_seqs(input_table, output_table='clean_cluster.tsv'):
    data_dict = {}

    with open(input_table, 'r') as ib, open(output_table, 'w') as ob:
        for line in ib:
            # 分割每一行并获取col1, col2, col3
            col1, col2, col3 = line.strip().split('\t')
            col3 = int(col3)

            # 如果col1不在字典中，添加新的词典条目
            if col1 not in data_dict:
                data_dict[col1] = [{'col2': col2, 'col3': col3}]
            else:
                '''如果col1在字典中，比较并更新col3值；如果新col3值大于字典中的最大值，则替换；如果col3值等于字典中最大值，则追加'''
                max_col3 = max(value['col3'] for value in data_dict[col1])
                if col3 > max_col3:
                    data_dict[col1] = [{'col2': col2, 'col3': col3}]
                elif col3 == max_col3:
                    data_dict[col1].append({'col2': col2, 'col3': col3})

        ob.write('UMI\tSequence\tCount\n')

        for col1, values in data_dict.items():
            for value in values:
                if value["col3"] != 1:
                    ob.write(f'{col1}\t{value["col2"]}\t{value["col3"]}\n')

    return output_table


def count_umi(input_table, output_table, output_fasta):
    # count the sequences with same umi

    umi_count = collections.defaultdict(int)
    with open(input_table, 'r') as ib, open(output_table, 'w') as ob, open(output_fasta, 'w') as of:
        for i, line in enumerate(ib):
            if i == 0:
                continue
            # 分割每一行并获取col1, col2, col3
            line = line.strip().split('\t')
            col1, col2, col3 = line
            umi_count[col2] += 1
            if len(line) != 3:
                continue
        ob.write('Sequences\tCount\n')

        i = 0
        for seqs, count in umi_count.items():
            i += 1
            ob.write(f'{seqs}\t{count}\n')
            of.write(f'>miRNA_{i}\n{seqs}\n')


# def format_fastq_record(record):
#     """Format a single FASTQ record with newline characters."""
#     return f"{record[0]}\n{record[1]}\n{record[2]}\n{record[3]}\n"


def paired_end_quality_control(read1_filename, read2_filename, output_read1_filename, output_read2_filename,
                               avg_quality_threshold, cut_front_window_size, cut_front_mean_quality,
                               cut_tail_window_size, cut_tail_mean_quality, max_n_bases,
                               qualified_quality_phred, unqualified_percent_limit):
    def average_quality(quality_scores):
        return sum(quality_scores) / len(quality_scores)

    def count_n_bases(sequence):
        return sequence.upper().count('N')

    def trim_window_front(sequence, quality_scores, front_window_size, front_quality_threshold):
        """Trim from the 5' end."""
        for i in range(0, len(quality_scores) - front_window_size + 1):
            window_avg = average_quality(quality_scores[i:i + front_window_size])
            if window_avg >= front_quality_threshold:
                return sequence[i:], quality_scores[i:]
        return '', []

    def trim_window_tail(sequence, quality_scores, tail_window_size, tail_quality_threshold):
        """Trim from the 3' end."""
        for i in range(len(quality_scores) - tail_window_size, -1, -1):
            window_avg = average_quality(quality_scores[i:i + tail_window_size])
            if window_avg >= tail_quality_threshold:
                return sequence[:i + tail_window_size], quality_scores[:i + tail_window_size]
        return '', []

    def percent_unqualified_bases(quality_scores, qualified_quality_value):
        unqualified_count = sum(1 for q in quality_scores if q < qualified_quality_value)
        return (unqualified_count / len(quality_scores)) * 100

    # Read and process read1 and read2
    with open(read1_filename, 'r') as f1, \
            open(read2_filename, 'r') as f2, \
            open(output_read1_filename, 'w') as out_file1, \
            open(output_read2_filename, 'w') as out_file2:
        while True:
            try:
                identifier_read1 = next(f1).strip()
                sequence_read1 = next(f1).strip()
                plus_line_read1 = next(f1).strip()
                quality_read1 = next(f1).strip()

                identifier_read2 = next(f2).strip()
                sequence_read2 = next(f2).strip()
                plus_line_read2 = next(f2).strip()
                quality_read2 = next(f2).strip()

                # Decode quality scores from Phred33 encoding
                quality_scores_read1 = [ord(char) - 33 for char in quality_read1]
                quality_scores_read2 = [ord(char) - 33 for char in quality_read2]

                # Trim the sequences from the 5' end
                sequence_trimmed_read1, quality_scores_trimmed_read1 = trim_window_front(
                    sequence_read1, quality_scores_read1, cut_front_window_size, cut_front_mean_quality)
                sequence_trimmed_read2, quality_scores_trimmed_read2 = trim_window_front(
                    sequence_read2, quality_scores_read2, cut_front_window_size, cut_front_mean_quality)

                # Trim the sequences from the 3' end
                sequence_trimmed_read1, quality_scores_trimmed_read1 = trim_window_tail(
                    sequence_trimmed_read1, quality_scores_trimmed_read1, cut_tail_window_size, cut_tail_mean_quality)
                sequence_trimmed_read2, quality_scores_trimmed_read2 = trim_window_tail(
                    sequence_trimmed_read2, quality_scores_trimmed_read2, cut_tail_window_size, cut_tail_mean_quality)

                # After trimming, check if the sequences are not empty
                if len(sequence_trimmed_read1) == 0 or len(sequence_trimmed_read2) == 0:
                    continue  # Skip if any sequence is empty after trimming

                # Check the number of 'N' bases
                n_bases_read1 = count_n_bases(sequence_trimmed_read1)
                n_bases_read2 = count_n_bases(sequence_trimmed_read2)
                if n_bases_read1 > max_n_bases or n_bases_read2 > max_n_bases:
                    continue  # Skip if the number of 'N' bases exceeds the threshold

                # Check the percentage of unqualified bases
                unqualified_percent_read1 = percent_unqualified_bases(quality_scores_trimmed_read1,
                                                                      qualified_quality_phred)
                unqualified_percent_read2 = percent_unqualified_bases(quality_scores_trimmed_read2,
                                                                      qualified_quality_phred)
                if unqualified_percent_read1 > unqualified_percent_limit \
                        or unqualified_percent_read2 > unqualified_percent_limit:
                    continue  # Skip if the percentage of unqualified bases exceeds the limit

                # Check if the average quality meets the threshold
                avg_quality_read1 = average_quality(quality_scores_trimmed_read1)
                avg_quality_read2 = average_quality(quality_scores_trimmed_read2)
                if avg_quality_read1 < avg_quality_threshold or avg_quality_read2 < avg_quality_threshold:
                    continue  # Skip if average quality is below threshold

                # Append the processed FASTQ records to the output lists
                out_file1.write(
                    f"{identifier_read1}\n{sequence_trimmed_read1}\n{plus_line_read1}\n"
                    f"{''.join([chr(q + 33) for q in quality_scores_trimmed_read1])}\n")
                out_file2.write(
                    f"{identifier_read2}\n{sequence_trimmed_read2}\n{plus_line_read2}\n"
                    f"{''.join([chr(q + 33) for q in quality_scores_trimmed_read2])}\n")
            except StopIteration:
                break

    return output_read1_filename, output_read2_filename
