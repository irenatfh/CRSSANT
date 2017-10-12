import numpy as np
import re


################################################################################
def process_cigar(pos_align, cigar_str):
    lens_strs = re.findall('\d+', cigar_str)
    lens = [int(i) for i in lens_strs]
    ops = re.findall('\D+', cigar_str)
    M_last = ''.join(ops).rfind('M')
    # Standardize all cigar strings to be of format S,M,N,M,S
    if ops[0] != 'S':
        lens = [0] + lens
        ops = ['S'] + ops
    if ops[-1] != 'S':
        lens = lens + [0]
        ops = ops + ['S']
    l_pos_start = pos_align + lens[0]
    l_pos_stop = l_pos_start + lens[1]
    r_pos_start = l_pos_stop + lens[2]
    r_pos_stop = r_pos_start + lens[3]
    
    return l_pos_start, l_pos_stop, r_pos_start, r_pos_stop