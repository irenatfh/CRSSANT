import numpy as np
import re
import math
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


def fold_optimize_stem(stem_inds, ref_seq):
    """
    Function to fold and optimize an RNA stem using RNAfold
    
    Results in the minimum free energy structure and energy of the folded stem, 
    and attempts folding optimization by truncating the stem, if possible.

    Parameters
    ----------
    stem_inds : np array
        Stem arm indices
    ref_seq : str
        Reference sequence
        
    Returns
    -------
    folded_stem_inds, fc, mfe : list, str, float
    """
    folded_stem_inds = np.copy(stem_inds)
    seq_l = ref_seq[stem_inds[0] : stem_inds[1] + 1]
    seq_r = ref_seq[stem_inds[2] : stem_inds[3] + 1]
    cut_point = len(seq_l)
    seq = seq_l + seq_r
    # Attempt folding
    res = RNA.fold_compound(seq_l + '&' + seq_r)
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
        # Check for folding results that might be fixable by truncation
        if set(fc[:cut_point]).issubset(set('.(')) is False:
            off_inds = [i for i in range(cut_point) if fc[i] == ')']
            seq_l = seq_l[off_inds[-1] + 1 : ]
            folded_stem_inds[0] += off_inds[-1] + 1
        if set(fc[cut_point:]).issubset(set('.)')) is False:
            off_inds = [i for i in range(len(seq_r)) if 
                        fc[cut_point : ][i] == '(']
            seq_r = seq_r[ : off_inds[0]]
            folded_stem_inds[3] = folded_stem_inds[2] + off_inds[0] - 1
        # Re-attempt helix folding/fold helix again if no truncation
        cut_point = len(seq_l)
        seq = seq_l + seq_r
        res = RNA.fold_compound(seq_l + '&' + seq_r)
        [fc, mfe] = res.mfe_dimer()
        # Check if truncation was successful
        l_symbols = [i for i in fc[ : cut_point] if i != '.']
        r_symbols = [i for i in fc[cut_point :] if i != '.']
        # If truncation was unsuccessful output null fc and mfe
        if (len(l_symbols) < 2 or len(r_symbols) < 2) or \
           (set(l_symbols) != set('(') or set(r_symbols) != set(')')):
            folded_stem_inds = np.zeros(4)
            fc = '.' * len(seq)
            mfe = 0.0
        else:
            pass
    return folded_stem_inds, fc, mfe


def fold_stem(seq_l, seq_r):
    """
    Function to fold an RNA stem using RNAfold
    
    This function simply attempts folding given left and right stem sequences.

    Parameters
    ----------
    seq_l : np array
        Stem left arm sequence
    seq_r : str
        Stem right arm sequence

    Returns
    -------
    fc, mfe : str, float
    """
    cut_point = len(seq_l)
    fc, mfe = RNA.fold_compound(seq_l + '&' + seq_r).mfe_dimer()
    fc_l = [i for i in fc[ : cut_point] if i != '.']
    fc_r = [i for i in fc[cut_point : ] if i != '.']
    if (len(fc_l) < 2 or len(fc_r) < 2) or \
    (set(fc_l) != set('(') or set(fc_r) != set(')')):
        fc = '.' * len(fc)
        mfe = 0.0
    return fc, mfe


def shuffle_stem(seq_l, seq_r, n):
    """
    Function to shuffle a stem up to n times

    Parameters
    ----------
    seq_l : np array
        Stem left arm sequence
    seq_r : str
        Stem right arm sequence
    n : int
        Number of shifts

    Returns
    -------
    mfes_shuffled : dict
    """
    seqs_shuffled = set()
    mfes_shuffled = []
    loop = 0
    while (len(seqs_shuffled) < n) and (loop < math.factorial(len(seq_l))**2):
        seq = seq_l + seq_r
        if seq not in seqs_shuffled:
            seqs_shuffled.add(seq)
            fc, mfe = fold_stem(seq_l, seq_r)
            mfes_shuffled.append(mfe)
        seq_l = ushuffle.shuffle(seq_l, len(seq_l), 2)
        seq_r = ushuffle.shuffle(seq_r, len(seq_r), 2)
        loop += 1
    mfes_shuffled = sorted(mfes_shuffled)[::-1]  # sort from worst to best MFE
    return mfes_shuffled


def shift_stem(stem_inds, ref_seq, gene_inds, n):
    """
    Function to shift  a stem up to n times

    Parameters
    ----------
    stem_inds : list
        Stem index list
    ref_seq : str
        Reference sequence
    gene_inds : list
        Gene index list
    n : int
        Number of shifts

    Returns
    -------
    mfes_shifted : dict
    """
    stem_len = stem_inds[3] - stem_inds[0] + 1
    gene_range = range(gene_inds[0], gene_inds[1] + 1 - stem_len)
    valid_inds = list(set(gene_range).difference([stem_inds[0]]))
    inds_shift = np.random.choice(
        valid_inds, min(len(valid_inds), n), replace=False)
    mfes_shifted = []
    for ind_shift in inds_shift:
        inds_shifted = np.copy(stem_inds)
        shift = ind_shift - stem_inds[0]
        inds_shifted += shift
        inds, fc, mfe = fold_optimize_stem(inds_shifted, ref_seq)
        mfes_shifted.append(mfe)
    mfes_shifted = sorted(mfes_shifted)[::-1]  # sort from worst to best MFE
    return mfes_shifted


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
