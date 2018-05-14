import numpy as np
import re
import sys
sys.path.append('/home/ihwang/software/ViennaRNA-2.4.3/interfaces/Python3')
import RNA
import ushuffle


################################################################################
def process_cigar(cigar_str):
    """
    Function to process CIGAR string in sequencing reads

    Parameters
    ----------
    cigar_str : str
        CIGAR string
        
    Returns
    -------
    ops, lens : list, list
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


def get_overlap_ratios(inds_1, inds_2):
    """
    Function to calculate the overlap ratio between two reads

    Parameters
    ----------
    inds_1 : np array
        np array containing read indices
        [read left start, read left stop, read right start, read right stop]
    inds_2 : np array
        np array containing read indices
        [read left start, read left stop, read right start, read right stop]

    Returns
    -------
    float, float

    """
    overlap_l = min(inds_1[1], inds_2[1]) - \
                max(inds_1[0], inds_2[0]) + 1
    overlap_r = min(inds_1[3], inds_2[3]) - \
                max(inds_1[2], inds_2[2]) + 1
    span_l = max(inds_1[1], inds_2[1]) - \
             min(inds_1[0], inds_2[0]) + 1
    span_r = max(inds_1[3], inds_2[3]) - \
             min(inds_1[2], inds_2[2]) + 1
    ratio_l = overlap_l / span_l
    ratio_r = overlap_r / span_r
    return ratio_l, ratio_r


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


################################################################################
def generate_shifted_inds(gene_inds, struct_len, struct_inds, n_shifts):
    gene_range = range(gene_inds[0], gene_inds[1] + 1 - struct_len)
    valid_inds = list(set(gene_range).difference(struct_inds[0:1]))
    inds_shift = np.random.choice(valid_inds, min(len(valid_inds), n_shifts),
                                  replace=False)
    return inds_shift

################################################################################
def shift_dg(dg_inds, inds_shift, ref_seq, ARM_FLAG=0):
    mfes_shifted = []
    for ind_shift in inds_shift:
        inds_shifted = np.copy(dg_inds)
        if (ARM_FLAG) == 0 or (ARM_FLAG == 1):
            shift = ind_shift - dg_inds[0]
            inds_shifted[:2] += shift
            if ARM_FLAG == 0:
                inds_shifted[2:] += shift
        else:
            shift = ind_shift - dg_inds[2]
            inds_shifted[2:] += shift
        inds_shifted_folded, [fc, mfe] = calculate_stem_mfe(inds_shifted,
                                                            ref_seq,
                                                            TRUNC_FLAG=0)
        mfes_shifted.append(mfe)
    return mfes_shifted


################################################################################
def shuffle_dg(l_arm, r_arm, n_shuffles):
    cut_point = len(l_arm) + 1
    shuffled_seqs = set()
    shuffled_mfes = []
    while len(shuffled_seqs) < n_shuffles:
        seq = l_arm + r_arm
        if seq not in shuffled_seqs:
            shuffled_seqs.add(seq)
            [shuffled_fc, shuffled_mfe] = RNA.fold_compound(
                l_arm + '&' + r_arm).mfe_dimer()
            l_symbols = [i for i in shuffled_fc[ : cut_point] if i != '.']
            r_symbols = [i for i in shuffled_fc[cut_point : ] if i != '.']
            if (len(l_symbols) < 2 or len(r_symbols) < 2) or \
            (set(l_symbols) != set('(') or set(r_symbols) != set(')')):
                shuffled_fc = '.' * len(seq)
                shuffled_mfe = 0.0
            shuffled_mfes.append(shuffled_mfe)
        ushuffle.shuffle1(l_arm, len(l_arm), 2)
        l_arm = ushuffle.shuffle2()
        ushuffle.shuffle1(r_arm, len(r_arm), 2)
        r_arm = ushuffle.shuffle2()
    return shuffled_mfes