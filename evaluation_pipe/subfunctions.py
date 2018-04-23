import numpy as np
import sys
sys.path.append('/home/ihwang/software/ViennaRNA-2.4.3/interfaces/Python3')
import RNA
import ushuffle


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
    TRUNC_FLAG : int
        

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
            if set(fc[:cut_point]).issubset(set('.(')) is False:
                off_inds = [i for i in range(cut_point) if fc[i] == ')']
                l_arm = l_arm[off_inds[-1] + 1 : ]
                folded_stem_inds[0] += off_inds[-1] + 1
                read_start += off_inds[-1] + 1
            if set(fc[cut_point:]).issubset(set('.)')) is False:
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