# This file is part of CRSSANT:
# Crosslinked RNA Secondary Structure Analysis using Network Techniques
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
import itertools


def dg_to_sg_dict(dg_filtered_dict, reads_dict, p=90):
    """
    Function to create stem group (SG) dict from DG dict
    
    To create SG dict, read arms from all filtered DGs in the gene are 
    aggregated, and then pth percentile arm length is calculated. In each DG,
    reads that are below that arm length threshold are taken to be the reads
    for the SG corresponding to that DG.

    Parameters
    ----------
    dg_filtered_dict : dict
        Dictionary of filtered DGs and their reads
    reads_dict : dict
        Dictionary of reads
    p : int
        Percentile for arm length cutoff

    Returns
    -------
    sg_reads_dict : dict
    """
    # Get read length cutoff based on percentile 
    all_arm_lengths = []
    for (dg, dg_reads) in dg_filtered_dict.items():
        all_arm_lengths += [
            (
                reads_dict[read_id][1] - reads_dict[read_id][0] + 1,
                reads_dict[read_id][1] - reads_dict[read_id][0] + 1
            )
         for read_id in dg_reads
        ]
    all_arm_lengths = np.array(list(itertools.chain(*all_arm_lengths)))
    cutoff_length = np.percentile(all_arm_lengths, p)
    
    
    # Filter out reads that have either arm longer than the cutoff
    sg_reads_dict = {}
    for (dg, dg_reads) in dg_filtered_dict.items():
        filtered_reads = [
            read_id for read_id in dg_reads if
            (
                reads_dict[read_id][1] - reads_dict[read_id][0] + 1 
                <= cutoff_length
            )
            and
            (
                reads_dict[read_id][3] - reads_dict[read_id][2] + 1
                <= cutoff_length
            )
                ]
        if len(filtered_reads) > 0:
            sg_reads_dict[dg] = filtered_reads
    return sg_reads_dict


def create_sg_dict(sg_reads_dict, reads_dict, ref_seq):
    """
    Function to create SG dict

    Parameters
    ----------
    sg_reads_dict : dict
        Dictionary of DGs and info  
    reads_dict : dict
        Dictionary of reads
    ref_seq : str
        Reference sequence

    Returns
    -------
    stem_dict : dict
    """
    sg_dict = {}
    for (sg, sg_reads) in sg_reads_dict.items():
        sg_reads_inds = np.array(
            [reads_dict[read_id][:4] for read_id in sg_reads]
        )
        inds = [
            np.min(sg_reads_inds[:, 0]), np.max(sg_reads_inds[:, 1]),
            np.min(sg_reads_inds[:, 2]), np.max(sg_reads_inds[:, 3])
        ]
        l_seq = ref_seq[inds[0] : inds[1] + 1]
        r_seq = ref_seq[inds[2] : inds[3] + 1]
        fc, mfe = sf.fold_stem(l_seq, r_seq)
        cut_point = len(l_seq)
        fc_l = [i for i in fc[ : cut_point] if i == '(' or i == '{']
        fc_r = [i for i in fc[cut_point : ] if i == ')' or i == '}']
        if mfe != 0:
            sg_dict[sg] = {}
            sg_dict[sg]['num_reads'] = len(sg_reads_inds)
            sg_dict[sg]['fc'] = fc
            sg_dict[sg]['mfe'] = mfe
            crosslinks, basepairs = sf.get_stem_info(inds, fc, ref_seq)
            sg_dict[sg]['crosslinks'] = crosslinks
            sg_dict[sg]['basepairs'] = basepairs
    return sg_dict


def test_stems(stem_dict, ref_seq, gene_inds, n):
    """
    Function to test stems using shifts and shuffles
    
    Shifts are random along the gene reference sequence, but shuffles
    maintain dinucleotide content.

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
        l_seq = ref_seq[stem_inds[0] : stem_inds[1] + 1]
        r_seq = ref_seq[stem_inds[2] : stem_inds[3] + 1]
        mfes_shuffled = sf.shuffle_stem(l_seq, r_seq, n)
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