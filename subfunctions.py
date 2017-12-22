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
    uu_cl_counter = 0
    uc_cl_counter = 0
    for (l_ind, r_ind) in zip(l_bp_inds, r_bp_inds):
        if (l_ind != 0) or (r_ind != len(seq)):
            if (set(seq[l_ind]+seq[r_ind-1]) == set('TC')) or \
                (set(seq[l_ind]+seq[r_ind-1]) == set('T')):
                # "after" condition
                if (seq[l_ind] == 'T') and (seq[r_ind-1] == 'T'):
                    uu_cl_counter += 1
                else:
                    uc_cl_counter += 1
            if (set(seq[l_ind]+seq[r_ind+1]) == set('TC')) or \
                (set(seq[l_ind]+seq[r_ind+1]) == set('T')):
                # "before" condition
                if (seq[l_ind] == 'T') and (seq[r_ind+1] == 'T'):
                    uu_cl_counter += 1
                else:
                    uc_cl_counter += 1
    return [uu_cl_counter, uc_cl_counter], len(l_bp_inds)