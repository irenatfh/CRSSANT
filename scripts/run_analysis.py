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
import datetime
import preprocessing as pp, graphing as gp, dg_analysis as da, \
stem_discovery as sd, output as op


def run_analysis(args):
    gene_start = datetime.datetime.now()
    dg_ind = 0
    dg_reads_dict = None
    dg_dict = None
    sg_reads_dict = None
    sg_dict = None
    l_gene, r_gene, region, reads_dict, ref_dict, log, max_reads, min_overlap, bin_width = args
    
    
    region_seq = ref_dict[region]['sequence']
    ng_ind = 0
    gene_ids = [
        read_id for (read_id, read_info) in reads_dict.items() 
        if (read_info[4] == region) & (read_info[5] == l_gene) 
        & (read_info[6] == r_gene)
    ]
    gene_ids = pp.filter_gene_reads(gene_ids, reads_dict)
    if l_gene == r_gene:
        print(
            'Analyzing %s reads spanning gene %s\n' 
            %(len(gene_ids), l_gene)
        )
    else:
        print(
            'Analyzing %s reads spanning genes %s and %s\n' 
            %(len(gene_ids), l_gene, r_gene)
        )
    if len(gene_ids) > 1:


        # Analyze reads to obtain DGs
        weights = gp.get_weights(
                gene_ids, reads_dict, b_w=bin_width
            )
        if len(gene_ids) > max_reads:
            inds_samp = np.random.choice(
                len(gene_ids), max_reads, replace=False, p=weights
            )
        else:
            inds_samp = range(len(gene_ids))
        graph = gp.graph_reads(
                    [gene_ids[ind] for ind in inds_samp], 
                    reads_dict, t=min_overlap
                )
        reads_dg_dict, dg_ind = gp.cluster_graph(graph, dg_ind)
        dg_reads_dict = da.get_preliminary_dgs(
            reads_dict, reads_dg_dict
        )
        if len(gene_ids) > max_reads:


            # Add back in unsampled reads
            inds_unsamp = np.setdiff1d(
                range(len(gene_ids)), inds_samp
            )
            dg_reads_dict, dg_index = da.add_reads_to_dg(
                dg_reads_dict, [
                    gene_ids[ind] for ind in inds_unsamp
                ], 
                reads_dict, dg_ind, t=min_overlap
            )
        dg_filtered_dict = da.filter_dgs(dg_reads_dict)
        dg_dict, ng_ind = da.create_dg_dict(
                dg_filtered_dict, reads_dict, ng_ind
            )
        if len(dg_filtered_dict) >= 1:


            # If there is a valid number of DGs, do SG discovery
            # Refine reads to get stem group (SG) dict
            sg_reads_dict = sd.dg_to_sg_dict(dg_filtered_dict, reads_dict)
            sg_dict = sd.create_sg_dict(sg_reads_dict, reads_dict, region_seq)
    gene_stop = datetime.datetime.now()
    gene_time = gene_stop - gene_start
    return dg_reads_dict, dg_dict, sg_reads_dict, sg_dict, region, gene_time, len(gene_ids)