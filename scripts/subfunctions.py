#!/usr/bin/env python
#
# This file is part of CRSSANT:
# Crosslinked RNA Secondary Structure Analysis using Network Techniques
#
###############################################################################

"""
Main script for running CRSSANT analysis and discovery pipelines
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu

import numpy as np
import re
import math
import RNA
import ushuffle
from itertools import chain


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


def get_overlaps(inds_1, inds_2):
    """
    Calculate the overlaps between two reads and the distance spanned by them.
    
    The reads are each assumed to comprise two arms, a left and a right arm.
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
    overlap_l, overlap_r, span_l, span_r: int, int, int, int
    """
    overlap_l = min(inds_1[1], inds_2[1]) - \
                max(inds_1[0], inds_2[0]) + 1
    overlap_r = min(inds_1[3], inds_2[3]) - \
                max(inds_1[2], inds_2[2]) + 1
    span_l = max(inds_1[1], inds_2[1]) - \
             min(inds_1[0], inds_2[0]) + 1
    span_r = max(inds_1[3], inds_2[3]) - \
             min(inds_1[2], inds_2[2]) + 1
    return overlap_l, overlap_r, span_l, span_r


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
    ratio_l, ratio_r : float, float
    """
    overlap_l, overlap_r, span_l, span_r = get_overlaps(inds_1, inds_2)
    ratio_l = overlap_l / span_l
    ratio_r = overlap_r / span_r
    return ratio_l, ratio_r


def fold_stem(l_seq, r_seq):
    """
    Function to fold an RNA stem using RNAfold with constraints

    Parameters
    ----------
    l_seq : np array
        Stem left arm sequence
    r_seq : str
        Stem right arm sequence

    Returns
    -------
    fc, mfe : str, float
    """
    fc = RNA.fold_compound(l_seq + r_seq)
    fc.hc_add_from_db('<'*len(l_seq) + '>'*len(r_seq))
    fc, mfe = fc.pf()
    return fc, mfe


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


def get_stem_info(inds, fc, ref_seq):
    """
    Function to get crosslinking and basepairing information

    Valid crosslinks are U-U and U-C

    Parameters
    ----------
    inds : np array
        Stem indices
    fc : str
        Stem basepairing structure
    ref_seq : str
        Referene sequence

    Returns
    -------
    crosslinks, basepairs : list, list
    """
    inds_l = [
        (i + inds[0]) for i in range(len(fc)) 
        if fc[i] == '(' or fc[i] == '}'
    ]
    inds_r = [
        (inds[3] - i) for i in range(len(fc)) 
        if fc[::-1][i] == ')' or fc[::-1][i] == '}'
    ]
    inds_r = inds_r[::-1]
    
    
    # Check crosslinking sites
    uu_count = 0
    uc_count = 0
    for (ind_l, r_ind) in zip(inds_l, inds_r):
        if ind_l < max(inds_l):  # Check next base on right arm
            if (ref_seq[ind_l] == 'T') and (ref_seq[r_ind - 1] == 'T'):
                uu_count += 1
            elif set(ref_seq[ind_l] + ref_seq[r_ind - 1]) == set('T'):
                uc_count += 1
                    
        if ind_l > min(inds_l):  # Check previous base on right arm
            if (ref_seq[ind_l] == 'T') and (ref_seq[r_ind + 1] == 'T'):
                uu_count += 1
            elif (set(ref_seq[ind_l] + ref_seq[r_ind + 1]) == set('T')):
                uc_count += 1           
    crosslinks = [uu_count, uc_count]
    basepairs = [inds_l, inds_r]
    return crosslinks, basepairs


def get_gene(inds, region, genes, ref_dict):
    """
    Function to get gene

    This function assumes that the left and right arms of the stem
    are located in the same gene

    Parameters
    ----------
    inds : np array
        Stem indices
    region : str
        Genomic region
    genes : list
        List of genes
    ref_dict:
        Reference dictionary

    Returns
    -------
    gene : str
    """
    for gene in genes:
        gene_start = ref_dict[region]['genes'][gene][0]
        gene_stop = ref_dict[region]['genes'][gene][1] + 1
        if (inds[0] in range(gene_start, gene_stop)) and \
            (inds[2] in range(gene_start, gene_stop)):
            return gene