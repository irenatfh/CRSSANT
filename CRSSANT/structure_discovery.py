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
import bisect
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


def test_structs(struct_dict, ref_seq, gene_inds, n):
    """
    Function to test structures using shifts and shuffles
    
    Shifts are random along the gene reference sequence, but shuffles
    are those that maintain dinucleotide content.

    Parameters
    ----------
    struct_dict : dict
        Dictionary of structurally-valid DGs
    ref_seq : str
        Reference sequence
    gene_inds : list
        Gene index list
    n : int
        Number of tests to perform

    Returns
    -------
    test_dict : dict
    """
    test_dict = {}
    for (struct, struct_info) in struct_dict.items():
        test_dict[struct] = {}
        struct_inds = struct_info['struct_inds']
        struct_mfe = struct_info['mfe']
        
        
        # Shuffle arm contents
        seq_l = ref_seq[struct_inds[0] : struct_inds[1] + 1]
        seq_r = ref_seq[struct_inds[2] : struct_inds[3] + 1]
        mfes_shuffled = sf.shuffle_stem(seq_l, seq_r, n)
        ranking = bisect.bisect_right(
            [abs(i) for i in mfes_shuffled], abs(struct_mfe)
        )
        test_dict[struct]['shuffles'] = ranking / len(mfes_shuffled)

        
        # Shift arms
        mfes_shifted = sf.shift_stem(struct_inds, ref_seq, gene_inds, n)
        ranking = bisect.bisect_right(
            [abs(i) for i in mfes_shifted], abs(struct_mfe)
        )
        test_dict[struct]['shifts'] = ranking / len(mfes_shifted)
    return test_dict


def filter_structs(test_dict, t):
    """
    Function that checks which structures exceed a validity threshold

    Parameters
    ----------
    test_dict : dict
        Dictionary of structure test results
    t : float
        Validity threshold

    Returns
    -------
    struct_list : dict
    """
    struct_list = []
    for (struct, tests) in test_dict.items():
        if (tests['shuffles'] > t) and (tests['shifts'] > t):
            struct_list.append(struct)
    return struct_list
        
