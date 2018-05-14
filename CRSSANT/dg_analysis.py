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
import itertools


def get_preliminary_dgs(reads_dict, reads_dg_dict):
    """
    Function that creates inverse dictionary of reads_dg_dict

    Parameters
    ----------
    reads_dict : dict
        Dictionary of reads
    reads_dg_dict : dict
        Dictionary of reads and their spectral clustering DG assignments

    Returns
    -------
    dg_reads_dict : dict
    """
    dgs_list = set(reads_dg_dict.values())
    dg_reads_dict = {}
    for dg in dgs_list:
        dg_reads_list = [
            read_id for (read_id, read_dg) in reads_dg_dict.items() 
            if read_dg == dg
        ]
        dg_reads_dict[dg] = dg_reads_list
    return dg_reads_dict


def add_reads_to_dg(dg_reads_dict, reads_ids, reads_dict, dg_ind, t=0.3):
    """
    Function to adjust duplex groups

    The goal is to assign reads that were not sampled to DGs. Each read is 
    checked for overlap with an existing DGs.
    + If the read has no overlap with any DGs, create a new DG.
    + If the read overlaps at least two DGs equally, again create a new DG.
    + Otherwise, add the read to the DG with which it has the largest overlap
    Note that this process updates the dg_reads_dict throughout the loop.

    Parameters
    ----------
    dg_reads_dict : dict
        Dictionary of DGs and their reads
    reads_ids : list
        List of reads indices
    reads_dict : dict
        Dictionary of reads
    dg_index : int
    t : float
        Overlap threshold

    Returns
    -------
    dg_reads_dict, dg_ind : dict, int
    """
    new_dg = 0
    for read_id in reads_ids:
        read_inds = reads_dict[read_id][:4]
        dg_overlaps = {}
        for (dg, dg_reads) in dg_reads_dict.items():
            dg_inds = np.median(
                np.array(
                    [reads_dict[dg_read][:4] for dg_read in dg_reads]
                ), axis=0
            )
            ratio_l, ratio_r = sf.get_overlap_ratios(read_inds, dg_inds)
            if (ratio_l > t) and (ratio_r > t):
                dg_overlaps[dg] = ratio_l + ratio_r
        if len(dg_overlaps) == 0:
            new_dg = 1
        else:
            max_overlap = max(dg_overlaps.values())
            max_dgs_list = [
                i for i in dg_overlaps.keys() if dg_overlaps[i] == max_overlap
            ]
            if len(max_dgs_list) > 1:
                new_dg = 1
        if new_dg == 0:
            dg_reads_dict[max_dgs_list[0]].append(read_id)
        else:
            dg_reads_dict[dg_ind] = [read_id]
            dg_ind += 1
        new_dg = 0 
    return dg_reads_dict, dg_ind


def filter_dgs(dg_reads_dict):
    """
    Function to filter out invalid DGs
    
    DGs with fewer than three reads are eliminated.

    Parameters
    ----------
    dg_reads_dict : dict
        Dictionary of DGs and their reads

    Returns
    -------
    filtered_dict : dict
    """
    filtered_dict = {}
    for (dg, dg_reads) in dg_reads_dict.items():
        num_reads = len(dg_reads)
        if num_reads > 2:
            filtered_dict[dg] = dg_reads
    return filtered_dict


def create_dg_dict(dg_reads_dict, reads_dict, ng_ind):
    """
    Function that creates the DG dictionary
    
    The DG dictionary no longer has individual read IDs and instead contains
    metadata about each DG, including number of reads, coverage, and
    non-overlapping group (NG).
    
    Coverage is defined as c / sqrt(a*b) where
        + c = number of reads in a given DG
        + a = number of reads overlapping the left arm of the DG
        + b = number of reads overlapping the right arm of the DG

    Parameters
    ----------
    dg_reads_dict : dict
        Dictionary of DGs and their reads
    reads_dict : dict
        Dictionary of reads

    Returns
    -------
    dg_dict : dict
    """
    gene_reads = list(itertools.chain.from_iterable(dg_reads_dict.values()))
    dg_dict = {}
    ng_dict = {}
    for (dg, dg_reads) in dg_reads_dict.items():
        dg_min = min([reads_dict[i][0] for i in dg_reads])
        dg_max = max([reads_dict[i][3] for i in dg_reads])
        dg_inds = [int(i) for i in 
                   np.median(
                       np.array(
                           [reads_dict[i][0:4] for i in dg_reads]
                       ), axis=0
                   )
                  ]
        overlapping_l = 0
        overlapping_r = 0
        # Calculate coverage
        for read_id in gene_reads:
            read_inds = reads_dict[read_id][0:4]
            ratio_l, ratio_r = sf.get_overlap_ratios(dg_inds, read_inds)
            if ratio_l > 0:
                overlapping_l += 1
            if ratio_r > 0:
                overlapping_r += 1
        cov = len(dg_reads) / np.sqrt(overlapping_l * overlapping_r)
        # Assign NG
        if len(ng_dict) == 0:
            ng_dict[ng_ind] = [dg]
            dg_ng = ng_ind
            ng_ind += 1
        else:
            ng_assigned = 0
            for (ng, ng_dgs) in ng_dict.items():
                ng_overlaps = np.zeros(len(ng_dgs))
                for (i, ng_dg) in enumerate(ng_dgs):
                    ng_dg_reads = dg_reads_dict[ng_dg]
                    ng_dg_min = min([reads_dict[i][0] for i in ng_dg_reads])
                    ng_dg_max = max([reads_dict[i][3] for i in ng_dg_reads])
                    overlap = min(dg_max, ng_dg_max) - max(dg_min, ng_dg_min)
                    if overlap > 0:
                        ng_overlaps[i] = 1
                if sum(ng_overlaps) == 0:
                    ng_dict[ng].append(dg)
                    dg_ng = ng
                    ng_assigned += 1
                    break
            if ng_assigned == 0:
                ng_dict[ng_ind] = [dg]
                dg_ng = ng_ind
                ng_ind += 1
        dg_dict[dg] = {
            'arm_indices' : dg_inds, 'coverage': cov, 
            'num_reads' : len(dg_reads), 'NG' : dg_ng
        }
    return dg_dict, ng_ind