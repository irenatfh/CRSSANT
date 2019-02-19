# This file is part of CRSSANT:
# Crosslinked RNA Secondary Structure Analysis using Network Techniques
#
###############################################################################
"""
This module runs the main analysis.
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu


import numpy as np
import datetime
import preprocessing as pp, graphing as gp, dg_analysis as da, \
sg_analysis as sa, output as op


def run_analysis(args):
    gene_start = datetime.datetime.now()
    dg_reads_dict = None
    dg_stats_dict = None
    sg_reads_dict = None
    sg_dict = None
    l_gene, r_gene, reads_dict, region, region_seq, \
    c, t_o, t_eig, path, reads = args

    
    # Initalize outputs
    file_name = reads.split('/')[-1].split('.sam')[0]
    file_base = '/%s_CRSSANT.g%s,%s' %(file_name, l_gene, r_gene)
    if c == 'cliques':
        clustering_str = '%s.t_o%s' %(c, t_o)
    else:
        clustering_str = '%s.t_o%s.t_eig%s' %(c, t_o, t_eig)
    file_base += '.%s' %clustering_str
    out_sam = path + file_base + '.sam'
    out_dg = path + file_base +  '_dg.bed'
    out_sg_arc = path + file_base + '_sg_arc.bed'
    out_sg_bp = path + file_base + '_sg_bp.bed'
    out_sg = path + file_base + '_sg.aux'
    
    
    # Run analysis
    gene_ids = [
        read_id for (read_id, read_info) in reads_dict.items() 
        if (read_info[5] == l_gene) & (read_info[6] == r_gene)
    ]
    print(
        'Analyzing %s reads with left arm mapped to gene %s and right arm '
        'mapped to gene %s' %(len(gene_ids), l_gene, r_gene)
    )
    graph = gp.graph_reads(gene_ids, reads_dict, t=t_o)
    if c == 'cliques':
        reads_dg_dict = gp.get_cliques(graph)
    else:
        reads_dg_dict = gp.cluster_graph(graph, t=t_eig)
    dg_reads_dict = da.get_preliminary_dgs(
        reads_dict, reads_dg_dict
    )
    dg_filtered_dict = da.filter_dgs(dg_reads_dict, reads_dict)
    dg_stats_dict = da.create_dg_dict(
            dg_filtered_dict, reads_dict
        )
    if dg_stats_dict:
        # If there exist DGs that pass the filter, check whether any of them
        # contain non-overlapping reads
        n_nonoverlap = da.check_nonoverlapping_reads(
            dg_filtered_dict, reads_dict
        )
        if n_nonoverlap > 0:
            print(
                'ERROR: Pipeline aborted for reads with left arm mapped to gene '
                '%s and right arm mapped to gene %s. %s of %s total reads were '
                'grouped into DGs containing non-overlapping reads. If using '
                'spectral clustering method, re-run the pipeline for only this '
                'gene pair and try incrementing t_o by 0.1. If incrementing t_o '
                'doesn\'t fix it, try using cliques-finding method for '
                'clustering.\n' %(l_gene, r_gene, n_nonoverlap, len(gene_ids))
            )
            return
        else:
            # Output SAM file and dg file
            pp.init_outputs(reads, out_sam, out_dg, out_sg_arc, out_sg_bp, out_sg)
            op.write_dg_ng_sam(
                reads, out_sam, dg_reads_dict, dg_stats_dict
            )
            op.write_dg(out_dg, dg_stats_dict, region)
            sg_reads_dict = sa.dg_to_sg_dict(dg_filtered_dict, reads_dict)
            sg_stats_dict = sa.create_sg_dict(
                sg_reads_dict, reads_dict, region_seq
            )
            if sg_stats_dict:
                # Write remaining output files
                op.write_sg_bp(out_sg_bp, sg_stats_dict, region)
                op.write_sg_arc(out_sg_arc, sg_stats_dict, region)
                op.write_sg(
                    out_sg, sg_stats_dict, sg_reads_dict, dg_stats_dict, 
                    reads_dict
                )
    gene_stop = datetime.datetime.now()
    gene_time = gene_stop - gene_start
    print(
        'Total analysis time for reads with left arm mapped to gene %s and '
        'right arm mapped to gene %s: %s' %(l_gene, r_gene, gene_time)
    )