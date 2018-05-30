# This file is part of CRSSANT:
# Crosslinked RNA Secondary Structure Analysis using Network Techniques
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
        region_flag = 0
        for line in f:
            if line[0] == '>':
                # Join and save sequence for previous region (if exists)
                try:
                    region
                except NameError:
                    pass
                else:
                    ref_dict[region]['sequence'] = ''.join(
                        ref_dict[region]['sequence']
                    )
                region = line.split('>')[-1].rstrip()
                region_flag = 1
            else:
                if region_flag == 1:
                    ref_dict[region] = {}
                    ref_dict[region]['genes'] = {}
                    ref_dict[region]['sequence'] = [line.rstrip()]
                    region_flag = 0
                else:
                    ref_dict[region]['sequence'].append(line.rstrip())
    # Join sequence string for final region
    ref_dict[region]['sequence'] = ''.join(ref_dict[region]['sequence'])
                  
                
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
            
            
    with open(out_info, 'w') as f:
        pass
    
    
    with open(out_bp, 'w') as f:
        f.write('track graphType=arc itemRgb=on\n')
        
        
    with open(out_aux, 'w') as f:
        header = [
            'DG_coverage', 'UU_cl,UC_cl,stem_length', 'PASS',
            'L_start_min,L_start_max,L_start_std', 
            'L_stop_min,L_stop_max,L_stop_std',
            'R_start_min,R_start_max,R_start_std', 
            'R_stop_min,R_stop_max,R_stop_std']
        f.write('\t'.join(header) + '\n')
    return

        
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
