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
            left_start = dg_inds[0] + 1  # biology is 1-indexed
            right_start = dg_inds[2] + 1
            right_stop = dg_inds[3] + 1
            left_len = dg_inds[1] - dg_inds[0] + 1
            right_len = dg_inds[3] - dg_inds[2] + 1
            dg_str = 'Group_%d_%.16f' %(dg, coverage)
            line = [region, str(left_start), str(right_stop), dg_str, 
                    str(num_reads), '-', str(left_start), str(left_start), 
                    '0,0,0', '2', '%d,%d' %(left_len, right_len), 
                    '0,%d' %(right_start - left_start)]
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
    aux_file, dg_dict, dg_reads_dict, reads_dict, stem_dict, struct_list
):
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
                if dg in struct_list:
                    pass_str = '1'
                else:
                    pass_str = '0'
            else:
                crosslinks_basepairs_str = '0,0,0'
                pass_str = '0'
            line = [dg_str, crosslinks_basepairs_str, pass_str]
            
            
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


def write_bp(bed_file, stem_dict, struct_list, region, gene):
    with open(bed_file, 'a') as f:
        for dg in struct_list:
            bp_l, bp_r = stem_dict[dg]['basepairs']
            for (ind_l, ind_r) in zip(bp_l, bp_r):
                line = [
                    region, str(ind_l + 1), str(ind_r + 1), gene, '1','+', 
                    str(ind_l), str(ind_r), '0,0,0'
                ]  # biology is 1-indexed
                f.write('\t'.join(line) + '\n')          
    return