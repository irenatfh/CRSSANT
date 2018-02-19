import numpy as np
import re
import sys
sys.path.append('/home/ihwang/software/ViennaRNA-2.4.3/interfaces/Python3')
import RNA


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
    overlap_l = min(read_1_inds[1], read_2_inds[1]) - \
                max(read_1_inds[0], read_2_inds[0]) + 1
    overlap_r = min(read_1_inds[3], read_2_inds[3]) - \
                max(read_1_inds[2], read_2_inds[2]) + 1
    overlap_g = min(read_1_inds[2], read_2_inds[2]) - \
                max(read_1_inds[1], read_2_inds[1]) + 1
    span_l = max(read_1_inds[1], read_2_inds[1]) - \
             min(read_1_inds[0], read_2_inds[0]) + 1
    span_r = max(read_1_inds[3], read_2_inds[3]) - \
             min(read_1_inds[2], read_2_inds[2]) + 1
    span_g = max(read_1_inds[2], read_2_inds[2]) - \
             min(read_1_inds[1], read_2_inds[1]) + 1
    
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


################################################################################
def calculate_stem_mfe(stem_inds, ref_seq, TRUNC_FLAG=1):
    """
    Calculate the minimum free energy structure and energy of a stem. Attempt
    truncation if possible (and if desired, as indicated by the TRUNC_FLAG).

    Parameters
    ----------
    stem_inds : np array
        RNA stem arm indices [left_start, left_stop, right_start, right_start]
    ref_seq : str
        Reference sequence

    Returns
    -------
    str, float
        [structure in parenthetical notation, minimum free energy of structure]

    """
    
    folded_stem_inds = np.copy(stem_inds)
    l_arm = ref_seq[stem_inds[0] : stem_inds[1] + 1]
    r_arm = ref_seq[stem_inds[2] : stem_inds[3] + 1]
    read_start = stem_inds[0]
    cut_point = len(l_arm)
    seq = l_arm + r_arm
    # Attempt folding the sequence
    res = RNA.fold_compound(l_arm + '&' + r_arm)
    [fc, mfe] = res.mfe_dimer()
    l_symbols = [i for i in fc[ : cut_point] if i != '.']
    r_symbols = [i for i in fc[cut_point : ] if i != '.']
    # Check for fatal helix folding results
    if (len(l_symbols) < 2 or len(r_symbols) < 2) or \
       (set(l_symbols) != set('(') or set(r_symbols) != set(')')):
        folded_stem_inds = np.zeros(4)
        fc = '.' * len(seq)
        mfe = 0.0
    else:
        if TRUNC_FLAG == 1:
            # Check for folding results that might be fixable by truncation
            if (set(fc[:cut_point]) != set('.(')):
                off_inds = [i for i in range(cut_point) if fc[i] == ')']
                l_arm = l_arm[off_inds[-1] + 1 : ]
                folded_stem_inds[0] += off_inds[-1] + 1
                read_start += off_inds[-1] + 1
            if (set(fc[cut_point:]) != set('.)')):
                off_inds = [i for i in range(len(r_arm)) if 
                            fc[cut_point : ][i] == '(']
                r_arm = r_arm[ : off_inds[0]]
                folded_stem_inds[3] = folded_stem_inds[2] + off_inds[0] - 1
            # Re-attempt helix folding/fold helix again if no truncation
            cut_point = len(l_arm)
            seq = l_arm + r_arm
            res = RNA.fold_compound(l_arm + '&' + r_arm)
            [fc, mfe] = res.mfe_dimer()
            # Check if truncation was successful
            l_symbols = [i for i in fc[ : cut_point] if i != '.']
            r_symbols = [i for i in fc[cut_point :] if i != '.']
            # If truncation was unsuccessful output null fc
            if (len(l_symbols) < 2 or len(r_symbols) < 2) or \
               (set(l_symbols) != set('(') or set(r_symbols) != set(')')):
                folded_stem_inds = np.zeros(4)
                fc = '.' * len(seq)
                mfe = 0.0
            else:
                pass
    
    return folded_stem_inds, [fc, mfe]