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
                region = line.split('>')[-1].rstrip().split(' ')[0]
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
            gene = data[3].rstrip()
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
    regions_readgenes_dict : dict
        Dictionary of set of (l_gene, r_gene) tuples for each region
    """
    reads_dict = {}
    regions_readgenes_dict = {}
    with open(reads_file, 'r') as f:
        for line in f:
            if line[0] != '@':
                data = line.rstrip().split('\t')
                read_id = data[0]
                region = data[2]
                pos_align = int(data[3])
                cigar_str = data[5]
                read_seq = data[9]
                if len(data) > 19:
                    xg = data[19]
                else:
                    xg = 'XG:i:0'
                if (xg == 'XG:i:0') or (xg == 'XG:i:1'):
                    l_gene = None
                    r_gene = None
                    cigar_ops, cigar_lens = sf.process_cigar(cigar_str)
                    if cigar_ops == ['S', 'M', 'N', 'M', 'S']:                     
                        # Get read position in reference sequence
                        l_pos_start = pos_align - 1  # biology is 0-indexed
                        l_pos_stop = l_pos_start + cigar_lens[1] - 1
                        r_pos_start = l_pos_stop + cigar_lens[2] + 1
                        r_pos_stop = r_pos_start + cigar_lens[3] - 1
                        for (gene, inds) in ref_dict[region]['genes'].items():
                            if inds[0] <= l_pos_start <= inds[1]:
                                l_gene = gene
                            if inds[0] <= r_pos_start <= inds[1]:
                                r_gene = gene
                        if l_gene and r_gene:
                            # Get read sequence
                            l_ind_start = cigar_lens[0]
                            l_ind_stop = l_ind_start + cigar_lens[1]
                            r_ind_start = l_ind_stop
                            r_ind_stop = r_ind_start + cigar_lens[3]
                            l_seq = read_seq[l_ind_start : l_ind_stop]
                            r_seq = read_seq[r_ind_start : r_ind_stop]
                            
                            
                            reads_dict[read_id] = [
                                l_pos_start, l_pos_stop, r_pos_start, 
                                r_pos_stop, region, l_gene, r_gene,
                                l_seq, r_seq
                            ]
                            if region not in regions_readgenes_dict:
                                regions_readgenes_dict[region] = {(l_gene, r_gene)}
                            else:
                                regions_readgenes_dict[region].add((l_gene, r_gene))
    return reads_dict, regions_readgenes_dict


def combine_aligned_and_chimeric_reads(
    reads, reads_chimeric
):
    """
    Function to combine files of aligned and chimeric reads
    
    This function parses paired chimeric reads and adds them to the
    combined file with an XG:i:1 designation. Aligned reads are also assigned
    an XG:i:0 designation.
    Parameters
    ----------
    reads : str
        Path to aligned reads file (SAM)
    reads_chimeric : str
        Path to chimeric reads file (SAM)
    Returns
    -------
    reads_all : str
        Path to file containing both aligned and chimeric reads
    """
    reads_all = reads.split('.sam')[0] + '_combinedReads.sam'


    with open(reads, 'r') as f_aligned, \
         open(reads_all, 'w') as f:
        for line in f_aligned:
            if line[0] != '@':
                f.write(line.strip() + '\tXG:i:0\n')
            else:
                f.write(line)


    chimeric_dict = {}
    with open(reads_chimeric, 'r') as f_chimeric:
        for line in f_chimeric:
            if line[0] != '@':
                data = line.rstrip().split('\t')
                read_id = data[0]
                read_flag = int(data[1])
                if read_id not in chimeric_dict:
                    chimeric_dict[read_id] = {}
                chimeric_dict[read_id][read_flag] = data[2:]
    paired_ids = []
    for (read_id, read_dict) in chimeric_dict.items():
        if len(read_dict) == 2:
            paired_ids.append(read_id)


    with open(reads_all, 'a') as f:
        for read_id in paired_ids:
            read_dict = chimeric_dict[read_id]
            [r_flag, l_flag] = sorted(read_dict.keys())

            
            if (r_flag == 0) and (l_flag == 256):
                l_pos = int(read_dict[l_flag][1])
                l_cigar_str = read_dict[l_flag][3]
                l_cigar_ops, l_cigar_lens = sf.process_cigar(l_cigar_str)


                r_pos = int(read_dict[r_flag][1])
                r_cigar_str = read_dict[r_flag][3]
                r_cigar_ops, r_cigar_lens = sf.process_cigar(r_cigar_str)


                # Shared info
                s_region = read_dict[r_flag][0]
                s_seq = read_dict[r_flag][7]
                s_qual = read_dict[r_flag][8]
                s_as = int(
                    round(
                        np.mean(
                    [int(read_dict[r_flag][11].split(':')[-1]), 
                     int(read_dict[l_flag][11].split(':')[-1])]
                )
                    )
                )
                s_jm = read_dict[r_flag][15]
                s_ji = read_dict[r_flag][16]

                
                if (r_cigar_ops == ['S','M','S']) and \
                   (l_cigar_ops == ['S','M','S']):
                    m_l = l_cigar_lens[1]
                    x_l = l_cigar_lens[2]
                    n = r_pos - (l_pos + m_l + x_l - 1) - 1
                    m_r = r_cigar_lens[1]

                    if x_l == 0:
                        if n < 0:
                            cigar_lens = [m_l, -n, m_r]
                            cigar_ops = ['M','I','M']
                        else:
                            cigar_lens = [m_l, n, m_r]
                            cigar_ops = ['M','N','M']
                    else:
                        if n < 0:
                            cigar_lens = [m_l, x_l, -n, m_r]
                            cigar_ops = ['M','X','I','M']
                        else:
                            cigar_lens = [m_l, x_l, n, m_r]
                            cigar_ops = ['M','X','I','M']
                    cigar = ''.join(
                        str(i)+j for (i,j) in zip(cigar_lens, cigar_ops)
                    )
                    l_seq = s_seq[l_cigar_lens[0] :]
                    r_seq = s_seq[r_cigar_lens[0] : 
                                  r_cigar_lens[0] + r_cigar_lens[1]]
                    seq = l_seq + r_seq
                    l_qual = s_qual[l_cigar_lens[0] :]
                    r_qual = s_qual[r_cigar_lens[0] : 
                                    r_cigar_lens[0] + r_cigar_lens[1]]
                    qual = l_qual + r_qual
                    line_list = [
                        read_id, l_flag, s_region, l_pos, 255, cigar, '*', 0, 
                        0, seq, qual, 'NH:i:1', 'HI:i:1', 'AS:i:%d' %s_as, 
                        'nM:i:0', 'NM:i:0', 'MD:Z:%d' %(m_l + m_r), s_jm, s_ji, 
                        'XG:i:1\n'
                    ]
                    line = '\t'.join(
                        [
                            str(tag) if isinstance(tag, str) is False 
                            else tag for tag in line_list
                        ]
                    )
                    f.write(line)
    return reads_all

            
def init_dg_outputs(
    in_sam, out_sam, out_dg
):
    """
    Function to initialize output files produced after successful DG clustering

    Parameters
    ----------
    in_sam : str
        Path to input reads file (SAM)
    out_sam : str
        Path to output reads file (SAM)
    out_dg : str

    Returns
    -------
    """
    with open(in_sam, 'r') as f_r, open(out_sam, 'w') as f_w:
        for line in f_r:
            if line[0] == '@':
                f_w.write(line)
            
            
    with open(out_dg, 'w') as f:
        pass
    return


def init_sg_outputs(
    out_sg_arc, out_sg_bp, out_sg
):
    """
    Function to initialize output files produced after successful SG assembly

    Parameters
    ----------
    out_sg_arc : str
    out_sg_bp : str
    out_sg : str

    Returns
    -------
    """
    with open(out_sg_arc, 'w') as f:
        f.write('track graphType=arc itemRgb=on\n')
        
        
    with open(out_sg_bp, 'w') as f:
        f.write('track graphType=arc itemRgb=on\n')
        
        
    with open(out_sg, 'w') as f:
        header = [
            'Group_ID', 'num_reads', 'UU_cl,UC_cl,num_basepairs',
            'L_start_min,L_start_max,L_start_std', 
            'L_stop_min,L_stop_max,L_stop_std',
            'R_start_min,R_start_max,R_start_std', 
            'R_stop_min,R_stop_max,R_stop_std']
        f.write('\t'.join(header) + '\n')
    return