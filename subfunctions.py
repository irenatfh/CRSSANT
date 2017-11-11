import numpy as np
import re


################################################################################
def process_cigar(cigar_str):
    # parse the cigar string into operations (ops) and operation lengths (lens)
    ops_raw = re.findall('\D+', cigar_str)
    lens_strs = re.findall('\d+', cigar_str)
    lens_raw = [int(i) for i in lens_strs]
    # merge duplicate consecutive operations
    ops = []
    lens = []
    lens.append(lens_raw[0])
    ops.append(ops_raw[0])
    for i in range(1, len(ops_raw)):
        if ops_raw[i] == ops_raw[i-1]:
            lens[-1] += lens_raw[i]
        else:
            lens.append(lens_raw[i])
            ops.append(ops_raw[i])
    # Standardize all cigar strings to start and end with soft-clipped regions
    if ops[0] != 'S':
        lens = [0] + lens
        ops = ['S'] + ops
    if ops[-1] != 'S':
        lens = lens + [0]
        ops = ops + ['S']
    
    return ops, lens


################################################################################
def get_reads(sam_file):
    traits = ['L start pos', 'L stop pos', 'R start pos', 'R stop pos', 'DG']
    reads_info = []
    reads_ids = []
    line_counter = 1
    cigar_ops_set = set()
    with open(sam_file, 'r') as f:
        for line in f:
            if line[0] != '@':
                data = line.split('\t')
                read_id = data[0]
                pos_align = int(data[3])
                cigar_str = data[5]
                xg = data[19]
                dg = int(data[20].split(':')[-1])
                if (xg == 'XG:i:0') or (xg == 'XG:i:1'):
                    cigar_ops, cigar_lens = process_cigar(cigar_str)
                    cigar_ops_set.update({tuple(cigar_ops)})
                    if cigar_ops == ['S', 'M', 'N', 'M', 'S']:
                        l_pos_start = pos_align + cigar_lens[0]
                        l_pos_stop = l_pos_start + cigar_lens[1] - 1
                        r_pos_start = l_pos_stop + cigar_lens[2] + 1
                        r_pos_stop = r_pos_start + cigar_lens[3] - 1
                        reads_ids.append(read_id)
                        reads_info.append((l_pos_start, l_pos_stop, r_pos_start, r_pos_stop, dg))
            line_counter += 1
    reads_info = np.asarray(reads_info)
    reads_ids = np.asarray(reads_ids)
    reads_dict = dict(zip(reads_ids, reads_info))
    return reads_info, reads_ids, reads_dict