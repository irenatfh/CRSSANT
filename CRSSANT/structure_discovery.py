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


def get_struct_dict(dg_dict, ref_seq):
    """
    Function to extract valid structures to a new dictionary

    Parameters
    ----------
    dg_dict : dict
        Dictionary of DGs and info        
    ref_seq : str

    Returns
    -------
    struct_dict : dict
    """
    struct_dict = {}
    for (dg, dg_info) in dg_dict.items():
        struct_inds, fc, mfe = sf.fold_optimize_stem(
            dg_info['arm_inds'], ref_seq
        )
        if mfe != 0:
            struct_dict[dg] = {}
            struct_dict[dg]['struct_inds'] = struct_inds
            struct_dict[dg]['fc'] = fc
            struct_dict[dg]['mfe'] = mfe
    return struct_dict


def test_structs(struct_dict, ref_seq):
    """
    Function...

    Parameters
    ----------
    P1 : type
        Descr.     

    Returns
    -------
    tests_dict : dict
    """
    tests_dict = {}
    for (struct, struct_info) in struct_dict.items():
        tests_dict[struct] = {}
        struct_inds = struct_info['struct_inds']
        seq_l = region_seq[struct_inds[0] : struct_inds[1] + 1]
        seq_r = region_seq[struct_inds[2] : struct_inds[3] + 1]
        mfes_shuffled = shuffle_dg(seq_l, seq_r, n)

        
def shuffle_dg(seq_l, seq_r, n):
    """
    Function...

    Parameters
    ----------
    P1 : type
        Descr.     

    Returns
    -------
    tests_dict : dict
    """
    shuffled_seqs = set()
    shuffled_mfes = []
    while len(shuffled_seqs) < n:
        seq = seq_l + seq_r
        if seq not in shuffled_seqs:
            shuffled_seqs.add(seq)
            [shuffled_fc, shuffled_mfe] = sf.fold_stem(seq_l, seq_r)
            shuffled_mfes.append(shuffled_mfe)
        ushuffle.shuffle1(seq_l, len(seq_l), 2)
        seq_l = ushuffle.shuffle2()
        ushuffle.shuffle1(seq_r, len(seq_r), 2)
        seq_r = ushuffle.shuffle2()
    return shuffled_mfes
###############################################################################