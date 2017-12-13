import numpy as np
import re
import networkx as nx


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
                        l_pos_start = pos_align
                        l_pos_stop = l_pos_start + cigar_lens[1] - 1
                        r_pos_start = l_pos_stop + cigar_lens[2] + 1
                        r_pos_stop = r_pos_start + cigar_lens[3] - 1
                        reads_ids.append(read_id)
                        reads_info.append((l_pos_start, l_pos_stop, 
                                           r_pos_start, r_pos_stop, dg))
            line_counter += 1
    reads_info = np.asarray(reads_info)
    reads_ids = np.asarray(reads_ids)
    reads_dict = dict(zip(reads_ids, reads_info))
    return reads_info, reads_ids, reads_dict


################################################################################
def get_overlaps(read_1, read_2):
    overlap_l = min(read_1[1], read_2[1]) - max(read_1[0], read_2[0]) + 1
    overlap_r = min(read_1[3], read_2[3]) - max(read_1[2], read_2[2]) + 1
    overlap_g = min(read_1[2], read_2[2]) - max(read_1[1], read_2[1]) + 1
    span_l = max(read_1[1], read_2[1]) - min(read_1[0], read_2[0]) + 1
    span_r = max(read_1[3], read_2[3]) - min(read_1[2], read_2[2]) + 1
    span_g = max(read_1[2], read_2[2]) - min(read_1[1], read_2[1]) + 1
    return overlap_l, overlap_r, overlap_g, span_l, span_r, span_g


################################################################################
def count_crosslinks(seq, fc, mfe):
    l_bp_inds = [i for i in range(len(seq)) if fc[i] == '(' ]
    r_bp_inds = [i for i in range(len(seq)) if fc[i] == ')' ][::-1]
    cl_counter = 0
    for (l_ind, r_ind) in zip(l_bp_inds, r_bp_inds):
        if (l_ind != 0) or (r_ind != len(seq)):
            if (seq[l_ind:l_ind+2] == 'TA') and (seq[r_ind-1:r_ind+1] == 'TA'):
                cl_counter += 1  # "after" condition
        if (seq[l_ind-1:l_ind+1] == 'AT') and (seq[r_ind:r_ind+2] == 'AT'):
            cl_counter += 1  # "before" condition
    print("%s uridine cross-linking sites found in %s-bp stem" %(cl_counter, len(l_bp_inds)))