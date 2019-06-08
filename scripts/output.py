# This file is part of CRSSANT:
# Crosslinked RNA Secondary Structure Analysis using Network Techniques
#
###############################################################################
"""
This module is a collection of functions that deal with writing output files
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu


import numpy as np


def write_dg_ng_sam(reads_file, rna_file, dg_reads_dict, dg_dict):
    """
    Function that writes the ouput SAM file
    
    SAM file includes only reads that contributed to DG groups, i.e. if DG had
    only one read, that DG does not pass the DG filter and the read is not
    written to the SAM file.
    
    Parameters
    ----------
    reads_file : str
        Path to original SAM file
    rna_file : str
        Path to new SAM file
    dg_reads_dict : dict
        Dictionary of DG and reads
    dg_dict : dict
        Dictionary of DG metadata, including DG and NG

    
    Returns
    -------
    appends to output file
    """
    with open(reads_file, 'r') as f_read, \
         open(rna_file, 'a') as f:
        for line in f_read:
            if line[0] != '@':
                data = line.split('\n')[0].split('\t')
                read_id = data[0]
                for (dg, dg_info) in dg_dict.items():
                    ng = dg_info['NG']
                    dg_reads_list = dg_reads_dict[dg]
                    if read_id in dg_reads_list:
                        data.append('DG:i:' + str(dg))
                        data.append('NG:i:' + str(ng))
                        f.write('\t'.join(data) + '\n')
                        break          
    return


def write_dg(dg_file, dg_dict, region):
    with open(dg_file, 'a') as f:
        for (dg, dg_info) in dg_dict.items():
            dg_inds = dg_info['arm_inds']
            coverage = dg_info['coverage']
            num_reads = dg_info['num_reads']
            l_start = dg_inds[0] + 1  # biology is 1-indexed
            r_start = dg_inds[2] + 1
            r_stop = dg_inds[3] + 1
            l_len = dg_inds[1] - dg_inds[0] + 1
            r_len = dg_inds[3] - dg_inds[2] + 1
            dg_str = 'Group_%d_%.16f' %(dg, coverage)
            line = [region, str(l_start), str(r_stop), dg_str, 
                    str(num_reads), '-', str(l_start), str(l_start), 
                    '0,0,0', '2', '%d,%d' %(l_len, r_len), 
                    '0,%d' %(r_start - l_start)]
            f.write('\t'.join(line) + '\n')     
    return


def write_sg_bp(sg_bp_file, sg_dict, region):
    with open(sg_bp_file, 'a') as f:
        for (sg, sg_info) in sg_dict.items():
            l_bp, r_bp = sg_info['basepairs']
            for (l_ind, r_ind) in zip(l_bp, r_bp):
                line = [
                    region, 
                    str(l_ind + 1),  # biology is 1-indexed
                    str(r_ind + 1), 
                    str(sg)
                ]
                f.write('\t'.join(line) + '\n')          
    return


def write_sg_arc(sg_arc_file, sg_dict, region):
    with open(sg_arc_file, 'a') as f:
        for (sg, sg_info) in sg_dict.items():
            sg_inds = sg_info['arm_inds']
            line = [
                region, 
                str(int(np.mean(sg_inds[:2])) + 1),  # biology is 1-indexed
                str(int(np.mean(sg_inds[2:])) + 1),
                str(sg)
            ]
            f.write('\t'.join(line) + '\n')      
    return


def write_sg(
    sg_file, sg_dict, sg_reads_dict, dg_dict, reads_dict
):
    with open(sg_file, 'a') as f:
        for (dg, dg_info) in dg_dict.items():
            dg_str = 'Group_%d' %(dg)
            line = [dg_str]
            
            
            # Add number of reads in SG
            if dg in sg_dict.keys():
                num_reads_str = str(sg_dict[dg]['num_reads'])
            else:
                num_reads_str = '0'
            line.append(num_reads_str)
                
            
            # Add crosslinking, stem length and structure pass
            if dg in sg_dict.keys():
                sg_info = sg_dict[dg]
                crosslinks = sg_info['crosslinks']
                basepairs = sg_info['basepairs']
                crosslinks_str = [str(i) for i in crosslinks]
                basepairs_str = [str(len(basepairs[0]))]
                crosslinks_basepairs_str = ','.join(
                    crosslinks_str + basepairs_str
                )
            else:
                crosslinks_basepairs_str = '0,0,0'
            line.append(crosslinks_basepairs_str)
            
            
            # Add read edge statistics
            if dg in sg_dict.keys():
                sg_reads_list = sg_reads_dict[dg]
                sg_reads_inds = np.array(
                    [reads_dict[i][0:4] for i in sg_reads_list]
                )
                for i in range(4):
                    edge_inds = sg_reads_inds[:,i] + 1  # biology is 1-indexed
                    edge_min = str(np.min(edge_inds))
                    edge_max = str(np.max(edge_inds))
                    edge_sd = str(np.std(edge_inds))
                    edge_str = ','.join([edge_min, edge_max, edge_sd])
                    line.append(edge_str)
            else:
                line.append('\t'.join(['0,0,0']*4))
            f.write('\t'.join(line) + '\n')
            
    return