import collections
import os
import shutil
import gzip


def decompress_files(gzip_file):
    # decompress gzip files
    output_file = gzip_file.replace('.gz', '')
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


def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(dna_sequence))


def merge_seqs(reads1, reads2, merge_reads):
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


def find_UMI(fasta, umi_fasta, out_put):
    str1 = 'AACTGTAGGCACCATCAAT'
    str2 = 'AGATCGGAAGAGCACACGTCT'
    str3 = 'GTTCAGAGTTCTACAGTCCGACGATC'

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
                str3_end_index = str3_index + len(str1) - 1

                umi_seq = line2[str1_end_index + 1: str2_index]
                umi.write(line1 + '\n')
                umi.write(umi_seq + '\n')

                insert_seq = line2[str3_end_index + 1: str1_index]
                o.write(line1 + ' ' + umi_seq + '\n')
                o.write(insert_seq + '\n')

    return umi_fasta, out_put


def umi_cluster(umi_fasta, seqs_fasta, output_file, seq_number=1):
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
                # print(umi_read1)
                seqs_count[(umi_read1, seq_read)] += 1

    sorted_seqs = sorted(seqs_count.items(), key=lambda item: item[0][0] and item[1])

    with open(output_file, 'w') as out:
        for (umi_read, seq_read), count in sorted_seqs:
            out.write(f'{umi_read}\t{seq_read}\t{count}\n')