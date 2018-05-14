# This file is part of CRSSANT:
# Computational RNA Secondary Structure Analysis using Network Techniques
#
###############################################################################
"""
This module is a collection of functions that perform pre-processing tasks
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu

import numpy as np
import subfunctions as sf


def get_reference_dict(seq_file, gene_file):
    """
    Function to parse reference sequence and gene files into a dictionary

    Parameters
    ----------
    seq_file : str
        Path to reference sequence file (FASTA)
    gene_file : 
        Path to reference gene file (BED)

    Returns
    -------
    ref_dict : dict
        Nested dictionary of region sequence and genes
        {region:{seq: str, genes: {str: list}}
    """
    ref_dict = {}
    with open(seq_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                region = line.split('>')[-1].rstrip()
            else:
                ref_dict[region] = {}
                ref_dict[region]['sequence'] = line.rstrip()
                ref_dict[region]['genes'] = {}

    with open(gene_file, 'r') as f:
        for line in f:
            data = line.split('\t')
            region = data[0]
            gene_start = int(data[1]) - 1
            gene_stop = int(data[2]) - 1
            gene = data[3]
            if gene not in ref_dict[region]['genes'].keys():
                ref_dict[region]['genes'][gene] = [gene_start, gene_stop]
    return ref_dict


def parse_reads(reads_file, ref_dict):
    """
    Function to parse a reads file into a dictionary

    Parameters
    ----------
    reads_file : str
        Path to reads file (SAM)
    ref_dict : str
        Reference dictionary

    Returns
    -------
    reads_dict : dict
        Dictionary of read information
        {read ID:[L start, L stop, R start, R stop, 
                  region, L arm gene, R arm gene]}
    """
    reads_dict = {}
    with open(reads_file, 'r') as f:
        for line in f:
            if line[0] != '@':
                data = line.rstrip().split('\t')
                read_id = data[0]
                region = data[2]
                pos_align = int(data[3])
                cigar_str = data[5]
                if len(data) > 19:
                    xg = data[19]
                else:
                    xg = 'XG:i:0'
                if (xg == 'XG:i:0') or (xg == 'XG:i:1'):
                    cigar_ops, cigar_lens = sf.process_cigar(cigar_str)
                    if cigar_ops == ['S', 'M', 'N', 'M', 'S']:
                        l_pos_start = pos_align - 1  # biology is 0-indexed
                        l_pos_stop = l_pos_start + cigar_lens[1] - 1
                        r_pos_start = l_pos_stop + cigar_lens[2] + 1
                        r_pos_stop = r_pos_start + cigar_lens[3] - 1
                        for (gene, inds) in ref_dict[region]['genes'].items():
                            if inds[0] <= l_pos_start <= inds[1]:
                                l_gene = gene
                            if inds[0] <= r_pos_start <= inds[1]:
                                r_gene = gene
                        try:
                            l_gene, r_gene
                        except NameError:
                            pass
                        else:
                            reads_dict[read_id] = [
                                l_pos_start, l_pos_stop, r_pos_start, 
                                r_pos_stop, region, l_gene, r_gene
                            ]
    return reads_dict


def init_outputs(in_sam, out_sam, out_info, out_bp, out_aux):
    """
    Function to initialize output files

    Parameters
    ----------
    in_sam : str
        Path to input reads file (SAM)
    out_sam : str
        Path to output reads file (SAM)
    out_info : str
        Path to output info file (BED)
    out_bp : str
        Path to output basepairing file (BED)
    out_aux : str
        Path to output auxiliary file (AUX)

    Returns
    -------
    """
    with open(in_sam, 'r') as f_r, open(out_sam, 'w') as f_w:
        for line in f_r:
            if line[0] == '@':
                f_w.write(line)
            else:
                pass
            
    with open(out_info, 'w') as f_w:
        pass
    
    with open(out_bp, 'w') as f_w:
        pass
    with open(out_aux, 'w') as f_w:
        header = [
            'DG_coverage', 'UU_cl,UC_cl,helix_length', 
            'L_start_start,L_start_stop,L_start_std', 
            'L_stop_start,L_stop_stop,L_stop_std',
            'R_start_start,R_start_stop,R_start_std', 
            'R_stop_start,R_stop_stop,R_stop_std']
        f_w.write('\t'.join(header) + '\n')
        
        
def get_genes(ref_dict, regions):
    """
    Function to get list of genes when regions are specified

    Parameters
    ----------
    ref_dict : dict
        Reference dictionary
    regions : list
        Specified list of regions

    Returns
    -------
    genes : list
    """
    genes = []
    for region in regions:
        genes += list(ref_dict[region]['genes'].keys())
    return genes


def get_regions(ref_dict, genes):
    """
    Function to get list of regions when genes are specified

    Parameters
    ----------
    ref_dict : dict
        Reference dictionary
    genes : list
        Specified list of genes

    Returns
    -------
    regions : list
    """
    regions = []
    for gene in genes:
        for region in ref_dict.keys():
            if gene in ref_dict[region]['genes'].keys():
                regions.append(region)
                continue
    regions = list(set(regions))
    return regions
        

def get_analysis_dict(ref_dict, regions, genes):
    """
    Function to compose analysis regions and genes into a dict

    Parameters
    ----------
    ref_dict : dict
        Reference dictionary
    regions : list
        Specified list of regions
    genes : list
        Specified list of genes

    Returns
    -------
    analysis_dict : dict
    """
    analysis_dict = {}
    for gene in genes:
        for region in ref_dict.keys():
            if gene in ref_dict[region]['genes'].keys():
                if region not in analysis_dict.keys():
                    analysis_dict[region] = []
                analysis_dict[region].append(gene)
                continue
    for region in regions:
        if region not in analysis_dict.keys():
            analysis_dict[region] = list(ref_dict[region]['genes'].keys())
    return analysis_dict


################################################################################
def get_gs_dict(gs_bp_bed, ref_dict, region):
    """
    From a gold standard (GS) basepairing BED file, create a dictionary of GS
    structures.

    Parameters
    ----------
    gs_bp_bed : str
        file name of GS basepairing BED file
    region : str
        genomic region of interest
    ref_dict : dict
        dictionary of reference in format {gene: [start_ind, stop_ind]}

    Returns
    -------
    dict
        {struct: np.array([left_start_ind, left_stop_ind,
                           right_start_ind, right_stop_ind]
                         ), 
                 [left_arm_gene, right_arm_gene]
                 }

    """
    gs_structs_dict = {}
    left_arms = []
    right_arms = []
    with open(gs_bp_bed, 'r') as f:
        for line in f:
            if 'graphType=' in line:
                pass
            else:
                data = line.split('\t')
                bp_region = data[0]
                if bp_region == region:
                    bp_left = int(data[1])
                    bp_right = int(data[2])
                    left_arms.append(bp_left)
                    right_arms.append(bp_right)
    left_arms = np.asarray(left_arms);
    right_arms = np.asarray(right_arms)
    arm_limits = 1 + np.array(
        [i for i in range(len(left_arms) - 1) if (np.abs(np.diff(left_arms)[i]) > 3)
         and ((np.abs(np.diff(right_arms)[i]) > 3))]
    )
    arm_limits = np.concatenate(
        (np.concatenate((np.array([0]), 
                         arm_limits)
                       ),
         np.array([len(left_arms)])
        )
    )
    for i in range(len(arm_limits) - 1):
        arm_start = arm_limits[i]
        arm_stop = arm_limits[i+1]
        left_arm = -1 + np.array([left_arms[j] for j in range(arm_start, arm_stop)]
                                )  # python is 0-indexed
        right_arm = -1 + np.array(
            sorted([right_arms[j] for j in range(arm_start, arm_stop)])
        )
        struct_inds = np.array(
            [left_arm[0], left_arm[-1], right_arm[0], right_arm[-1]]
        )
        struct_ranges = []
        for j in range(2):
            for (gene, gene_inds) in ref_dict.items():
                gene_range = range(gene_inds[0], gene_inds[1] + 1)
                if (struct_inds[j*2] in gene_range):
                    struct_ranges.append(gene)
                    break
        if len(struct_ranges) == 2:
            gs_structs_dict[i] = [struct_inds, struct_ranges]
    return gs_structs_dict


def get_dg_dict(dg_info_bed, ref_dict, ref_seq):
    """
    From a duplex group (DG) info BED file, create a dictionary of DG
    structures.

    Parameters
    ----------
    dg_info_bed : str
        file name of DG info BED file
    ref_dict : dict
        dictionary of reference in format {gene: [start_ind, stop_ind]}
    ref_seq : str
        reference sequence

    Returns
    -------
    dict
        {dg: [np.array([left_start_ind, left_stop_ind, 
                        right_start_ind, right_stop_ind]
                      ), 
              [left_arm_gene, right_arm_gene],
              coverage,
              number of reads,
              RNAfold basepairing structure,
              RNAfold minimum free energy
             ]}
    list
    
    """
    # Get DG dict
    dg_dict = {}
    dg_valid_structs = []
    with open(dg_info_bed, 'r') as f:
        for line in f:
            data = line.split('\t')
            dg = int(data[3].split('_')[1])
            covg = float(data[3].split('_')[2])
            num_reads = int(data[4])
            left_start = int(data[1]) - 1
            lens = [int(i) for i in data[10].split(',')]
            right_start = [int(i) for i in data[11].split(',')][1] + left_start
            dg_inds = np.array([left_start, left_start + lens[0] - 1, 
                                right_start, right_start + lens[1] - 1])
            dg_ranges = []
            for i in range(2):
                for (gene, gene_inds) in ref_dict.items():
                    gene_range = range(gene_inds[0], gene_inds[1] + 1)
                    if (dg_inds[i*2] in gene_range):
                        dg_ranges.append(gene)
                        break
            folded_dg_inds, [fc, mfe] = sf.calculate_stem_mfe(dg_inds, ref_seq)
            if mfe != 0:
                dg_inds = folded_dg_inds
                dg_valid_structs.append(dg)
            else:
                fc = '.'
                mfe = 0
            dg_dict[dg] = [dg_inds, dg_ranges, covg, num_reads, fc, mfe]
    return dg_dict, dg_valid_structs
