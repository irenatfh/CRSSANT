# This file is part of CRSSANT:
# Computational RNA Secondary Structure Analysis using Network Techniques
#
###############################################################################
"""
This module is a collection of functions that perform DG analysis tasks
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu


import numpy as np
import sys
import subfunctions as sf
sys.path.append('/home/ihwang/software/ViennaRNA-2.4.3/interfaces/Python3')
import RNA


def get_struct_dict(dg_dict, ref_seq):
    """
    Function...

    Parameters
    ----------
    dg_dict : type
    ref_seq : str

    Returns
    -------
    struct_dict : dict
    """
    struct_dict = {}
    for (dg, dg_info) in dg_dict.items():
        struct_inds, fc, mfe = sf.fold_stem(dg_info['arm_inds'], ref_seq)
        if mfe != 0:
            struct_dict[dg] = {}
            struct_dict[dg]['struct_inds'] = struct_inds
            struct_dict[dg]['fc'] = fc
            struct_dict[dg]['mfe'] = mfe
    return struct_dict


def shuffle_dg(l_arm, r_arm, n):
    cut_point = len(l_arm) + 1
    shuffled_seqs = set()
    shuffled_mfes = []
    while len(shuffled_seqs) < n:
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