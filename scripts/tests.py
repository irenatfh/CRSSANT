# This file is part of CRSSANT:
# Crosslinked RNA Secondary Structure Analysis using Network Techniques
#
###############################################################################
"""
This module is a collection of functions that perform an assortment of tests
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu

import numpy as np
import subfunctions as sf


def get_sg_dict(bp_bed, ref_dict):
    """
    Function to create dictionary of gold standard structures
    
    Parameters
    ----------
    bp_bed : str
        Path to gold standard structures file (BED)
    ref_dict : dict
        Reference dictionary
    
    Returns
    -------
    sg_dict : dict
    """
    sg_dict = {}
    l_inds = []
    r_inds = []
    with open(bp_bed, 'r') as f:
        for line in f:
            if 'graphType=' not in line:
                data = line.split('\t')
                region = data[0]
                l_ind = int(data[1])
                r_ind = int(data[2])
                l_inds.append(l_ind)
                r_inds.append(r_ind)
    l_inds = np.asarray(l_inds)
    r_inds = np.asarray(r_inds)
    stem_limits = 1 + np.array(
        [
            i for i in range(len(l_inds) - 1) 
            if (np.abs(np.diff(l_inds)[i]) > 3)
            and ((np.abs(np.diff(r_inds)[i]) > 3))
        ]
    )
    stem_limits = np.concatenate(
        (np.concatenate((np.array([0]), stem_limits)), np.array([len(l_inds)]))
    )
    for i in range(len(stem_limits) - 1):
        stem_start = stem_limits[i]
        stem_stop = stem_limits[i+1]
        l_stem = -1 + np.array(
            [l_inds[j] for j in range(stem_start, stem_stop)]
        )  # python is 0-indexed
        r_stem = -1 + np.array(
            sorted([r_inds[j] for j in range(stem_start, stem_stop)])
        )
        stem_inds = np.array([l_stem[0], l_stem[-1], r_stem[0], r_stem[-1]])
        sg_dict[i] = stem_inds
    return sg_dict


def shuffle_stem(l_seq, r_seq, n):
    """
    Function to shuffle a stem up to n times

    Parameters
    ----------
    l_seq : np array
        Stem left arm sequence
    r_seq : str
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
    while (len(seqs_shuffled) < n) and (loop < math.factorial(len(l_seq))**2):
        seq = l_seq + r_seq
        if seq not in seqs_shuffled:
            seqs_shuffled.add(seq)
            fc, mfe = fold_stem(l_seq, r_seq)
            mfes_shuffled.append(mfe)
        l_seq = ushuffle.shuffle(l_seq, len(l_seq), 2)
        r_seq = ushuffle.shuffle(r_seq, len(r_seq), 2)
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