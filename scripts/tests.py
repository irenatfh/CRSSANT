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


def get_dg_dict(dg_bed):
    """
    Function to create DG dictionary based on DG file
    
    Parameters
    ----------
    dg_bed : str
        Path to DG file (BED)
        
    Returns
    -------
    dg_dict : dict
    """
    # Get DG dict
    dg_dict = {}
    dg_valid_structs = []
    with open(dg_bed, 'r') as f:
        for line in f:
            data = line.split('\t')
            region = data[0]
            dg = int(data[3].split('_')[1])
            dg_dict[dg] = {}
            coverage = float(data[3].split('_')[2])
            num_reads = int(data[4])
            l_start = int(data[1]) - 1
            lens = [int(i) for i in data[10].split(',')]
            r_start = [int(i) for i in data[11].split(',')][1] + l_start
            dg_inds = np.array([l_start, l_start + lens[0] - 1, 
                                r_start, r_start + lens[1] - 1])
            dg_dict[dg]['region'] = region
            dg_dict[dg]['arm_inds'] = dg_inds
            dg_dict[dg]['coverage'] = coverage
            dg_dict[dg]['num_reads'] = num_reads
    return dg_dict


def get_gs_dict(gs_bed, ref_dict, regions_genes_dict):
    """
    Function to create dictionary of gold standard structures
    
    Parameters
    ----------
    gs_bed : str
        Path to gold standard structures file (BED)
    ref_dict : dict
        Reference dictionary
    regions_genes_dict : dict
        Analysis dictionary
    
    Returns
    -------
    gs_dict : dict
    """
    gs_dict = {}
    regions = []
    l_inds = []
    r_inds = []
    with open(gs_bed, 'r') as f:
        for line in f:
            if 'graphType=' in line:
                pass
            else:
                data = line.split('\t')
                region = data[0]
                l_ind = int(data[1])
                r_ind = int(data[2])
                regions.append(region)
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
        region = regions[stem_start]
        if region in regions_genes_dict:
            genes = regions_genes_dict[region]
            gene = sf.get_gene(stem_inds, region, genes, ref_dict)
            if gene is not None:
                gs_dict[i] = {}
                gs_dict[i]['stem_inds'] = stem_inds
                gs_dict[i]['region'] = region
                gs_dict[i]['gene'] = gene
    return gs_dict


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