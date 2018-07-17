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


def write_info(bed_file, dg_dict, region):
    with open(bed_file, 'a') as f:
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


def write_dg_ng_sam(reads_file, rna_file, dg_reads_dict, dg_dict):
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


def write_aux(
    aux_file, dg_dict, dg_reads_dict, reads_dict, stem_dict):
    with open(aux_file, 'a') as f:
        for (dg, dg_info) in dg_dict.items():
            coverage = dg_info['coverage']
            dg_str = 'Group_%d_%.16f' %(dg, coverage)
            
            
            # Add crosslinking, stem length and structure pass
            if dg in stem_dict:
                crosslinks = stem_dict[dg]['crosslinks']
                basepairs = stem_dict[dg]['basepairs']
                crosslinks_str = [str(i) for i in crosslinks]
                basepairs_str = [str(len(basepairs[0]))]
                crosslinks_basepairs_str = ','.join(
                    crosslinks_str + basepairs_str
                )
            else:
                crosslinks_basepairs_str = '0,0,0'
            line = [dg_str, crosslinks_basepairs_str]
            
            
            # Add read edge statistics
            dg_reads_list = dg_reads_dict[dg]
            dg_reads_inds = np.array(
                [reads_dict[i][0:4] for i in dg_reads_list]
            )
            for i in range(4):
                edge_inds = dg_reads_inds[:,i] + 1  # biology is 1-indexed
                edge_min = str(np.min(edge_inds))
                edge_max = str(np.max(edge_inds))
                edge_sd = str(np.std(edge_inds))
                edge_str = ','.join([edge_min, edge_max, edge_sd])
                line.append(edge_str)
            f.write('\t'.join(line) + '\n')
    return


def write_dg_arcs(bed_file, dg_dict, region):
    with open(bed_file, 'a') as f:
        for (dg, dg_info) in dg_dict.items():
            dg_inds = dg_info['arm_inds']
            line = [
                region, 
                str(int(np.median(dg_inds[:2])) + 1),  # biology is 1-indexed
                str(int(np.median(dg_inds[2:])) + 1),
                str(dg), '1','+', 
                str(int(np.median(dg_inds[:2])) + 1),
                str(int(np.median(dg_inds[2:])) + 1),
                '0,0,0'
            ]
            f.write('\t'.join(line) + '\n')      
    return


def write_dg_bps(bed_file, stem_dict, region):
    with open(bed_file, 'a') as f:
        for (stem, stem_info) in stem_dict.items():
            l_bp, r_bp = stem_info['basepairs']
            for (l_ind, r_ind) in zip(l_bp, r_bp):
                line = [
                    region, 
                    str(l_ind + 1),  # biology is 1-indexed
                    str(r_ind + 1), 
                    str(stem), '1','+', 
                    str(l_ind), 
                    str(r_ind), 
                    '0,0,0'
                ]
                f.write('\t'.join(line) + '\n')          
    return