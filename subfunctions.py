import numpy as np
import re


################################################################################
def process_cigar(cigar_str):
    """
    Process the CIGAR string accompanying a read.

    Parameters
    ----------
    cigar_str : str
        CIGAR string
        
    Returns
    -------
    list, list
        List of operations, list of lengths associated with the operations

    """
    # Parse the cigar string into operations (ops) and operation lengths (lens)
    ops_raw = re.findall('\D+', cigar_str)
    lens_strs = re.findall('\d+', cigar_str)
    lens_raw = [int(i) for i in lens_strs]
    # Merge duplicate consecutive operations
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
def get_overlaps(read_1_inds, read_2_inds):
    """
    Calculate the overlaps between two reads and the distance spanned by them.
    
    The reads are each assumed to comprise two arms, a left and a right arm.

    Parameters
    ----------
    read_1_inds : np array
        np array containing information of read.
        [read left start, read left stop, read right start, read right stop]
    read_2_inds : np array
        np array containing information of read.
        [read left start, read left stop, read right start, read right stop]

    Returns
    -------
    int, int, int, int, int, int
        overlap of left arms, overlap of right arms, overlap of the gaps,
        span of the left arm, span of the right arm, span of the gaps

    """
    overlap_l = min(read_1[1], read_2[1]) - max(read_1[0], read_2[0]) + 1
    overlap_r = min(read_1[3], read_2[3]) - max(read_1[2], read_2[2]) + 1
    overlap_g = min(read_1[2], read_2[2]) - max(read_1[1], read_2[1]) + 1
    span_l = max(read_1[1], read_2[1]) - min(read_1[0], read_2[0]) + 1
    span_r = max(read_1[3], read_2[3]) - min(read_1[2], read_2[2]) + 1
    span_g = max(read_1[2], read_2[2]) - min(read_1[1], read_2[1]) + 1
    
    return overlap_l, overlap_r, overlap_g, span_l, span_r, span_g


################################################################################
def get_intersection(read_1_inds, read_2_inds):
    """
    Calculate the intersections between two reads, defined as the total number
    of overlapping bases over any of the arms (left arm and left arm, left arm 
    and right arm, etc.)
    
    The reads are each assumed to comprise two arms, a left and a right arm.

    Parameters
    ----------
    read_1_inds : np array
        np array containing information of read.
        [read left start, read left stop, read right start, read right stop]
    read_2_inds : np array
        np array containing information of read.
        [read left start, read left stop, read right start, read right stop]

    Returns
    -------
    int
        overlap of left arms, overlap of right arms, overlap of the gaps,
        span of the left arm, span of the right arm, span of the gaps

    """
    range_1 = np.concatenate((np.arange(read_1_inds[0], read_1_inds[1] + 1),
                              np.arange(read_1_inds[2], read_1_inds[3] + 1)), 
                             axis=0)
    range_2 = np.concatenate((np.arange(read_2_inds[0], read_2_inds[1] + 1),
                              np.arange(read_2_inds[2], read_2_inds[3] + 1)), 
                             axis=0)
    intersection = set(range_1).intersection(set(range_2))
    
    return len(intersection)


################################################################################
def read_inds_to_vec(read_inds, gene_length, gene_start):
    """
    Converts a read index array into a read vector

    Parameters
    ----------
    read_inds : np array
        np array containing information of read.
        [read left start, read left stop, read right start, read right stop]

    Returns
    -------
    np array
        vector of length gene_length with read arm positions indicated by 1s 

    """
    read_vec = np.zeros(gene_length)
    l_arm_range = np.arange(read_inds[0], read_inds[1] + 1) - (gene_start - 1)
    r_arm_range = np.arange(read_inds[2], read_inds[3] + 1) - (gene_start - 1)
    total_range = np.concatenate((l_arm_range, r_arm_range), axis=0)
    for i in total_range:
        read_vec[i] = 1
    
    return read_vec


################################################################################
def count_crosslinks(seq, fc):
    """
    Count uridine crosslinking sites in the DG sequence.

    Valid crosslinks are U-U and U-C.

    Parameters
    ----------
    seq : str
        DG sequence, left and right arms concatenated
    fc : str
        ViennaRNA folding structure

    Returns
    -------
    np array
        [uu crosslink count, uc crosslink count, helix length]

    """
    l_bp_inds = [i for i in range(len(seq)) if fc[i] == '(' ]
    r_bp_inds = [i for i in range(len(seq)) if fc[i] == ')' ][::-1]
    uu_cl_counter = 0
    uc_cl_counter = 0
    for (l_ind, r_ind) in zip(l_bp_inds, r_bp_inds):
        if l_ind < max(l_bp_inds): # Check previous base on right arm
            if (set(seq[l_ind] + seq[r_ind-1]) == set('TC')) or \
                (set(seq[l_ind] + seq[r_ind-1]) == set('T')):
                if (seq[l_ind] == 'T') and (seq[r_ind-1] == 'T'):
                    uu_cl_counter += 1
                else:
                    uc_cl_counter += 1
                    
        if l_ind > min(l_bp_inds):# Check following base on right arm
            if (set(seq[l_ind] + seq[r_ind+1]) == set('TC')) or \
                (set(seq[l_ind] + seq[r_ind+1]) == set('T')):
                if (seq[l_ind] == 'T') and (seq[r_ind+1] == 'T'):
                    uu_cl_counter += 1
                else:
                    uc_cl_counter += 1
                    
    return np.array([uu_cl_counter, uc_cl_counter, len(l_bp_inds)], dtype=np.int)