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


def get_stem_dict(dg_dict, ref_seq):
    """
    Function to extract valid stems to a new dictionary

    Parameters
    ----------
    dg_dict : dict
        Dictionary of DGs and info        
    ref_seq : str

    Returns
    -------
    stem_dict : dict
    """
    stem_dict = {}
    for (dg, dg_info) in dg_dict.items():
        inds, fc, mfe = sf.fold_optimize_stem(
            dg_info['arm_inds'], ref_seq
        )
        if mfe != 0:
            stem_dict[dg] = {}
            stem_dict[dg]['stem_inds'] = inds
            stem_dict[dg]['fc'] = fc
            stem_dict[dg]['mfe'] = mfe
            crosslinks, basepairs = sf.get_stem_info(inds, fc, ref_seq)
            stem_dict[dg]['crosslinks'] = crosslinks
            stem_dict[dg]['basepairs'] = basepairs
    return stem_dict


def test_stems(stem_dict, ref_seq, gene_inds, n):
    """
    Function to test stems using shifts and shuffles
    
    Shifts are random along the gene reference sequence, but shuffles
    are those that maintain dinucleotide content.

    Parameters
    ----------
    stem_dict : dict
        Dictionary of stem-forming DGs
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
    for (stem, stem_info) in stem_dict.items():
        test_dict[stem] = {}
        stem_inds = stem_info['stem_inds']
        stem_mfe = stem_info['mfe']
        
        
        # Shuffle arm contents
        seq_l = ref_seq[stem_inds[0] : stem_inds[1] + 1]
        seq_r = ref_seq[stem_inds[2] : stem_inds[3] + 1]
        mfes_shuffled = sf.shuffle_stem(seq_l, seq_r, n)
        ranking = bisect.bisect_right(
            [abs(i) for i in mfes_shuffled], abs(stem_mfe)
        )
        test_dict[stem]['shuffles'] = ranking / len(mfes_shuffled)

        
        # Shift arms
        mfes_shifted = sf.shift_stem(stem_inds, ref_seq, gene_inds, n)
        ranking = bisect.bisect_right(
            [abs(i) for i in mfes_shifted], abs(stem_mfe)
        )
        test_dict[stem]['shifts'] = ranking / len(mfes_shifted)
    return test_dict


def filter_stems(test_dict, t):
    """
    Function to check which stems pass a validity threshold

    Parameters
    ----------
    test_dict : dict
        Dictionary of structure test results
    t : float
        Validity threshold

    Returns
    -------
    stem_list : dict
    """
    stem_list = []
    for (stem, tests) in test_dict.items():
        if (tests['shuffles'] >= t) and (tests['shifts'] >= t):
            stem_list.append(stem)
    return stem_list