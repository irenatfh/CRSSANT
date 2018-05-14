#!/usr/bin/env python
#
# This file is part of CRSSANT:
# Computational RNA Secondary Structure Analysis using Network Techniques
#
###############################################################################

"""
Main script for running CRSSANT analysis and discovery pipelines
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu

import sys
import argparse
import time
import numpy as np
from itertools import chain
import preprocessing as pp, graphing as gp, dg_analysis as da, \
structure_discovery as sd, output as op


# Global variables
max_reads = 2000
min_overlap = 0.3


class ReadsFiles(object):
    """
    """
    def __init__(self, reads, out):
        self.name = reads.split('/')[-1].split('.sam')[0] + '_CRSSANT'
        self.out_sam = out + self.name + '.sam'
        self.out_info = out + self.name + '_info.bed'
        self.out_bp = out + self.name + '_bp.bed'
        self.out_aux = out + self.name + '.aux'
        self.log = out + self.name + '.log'
        
    
def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser=argparse.ArgumentParser(description='CRSSANT analyzes PARIS reads '
                                   'for RNA structures',
                                   epilog='python ./CRSSANT.py reads.sam '
                                   'ref.fa ref.bed out region1,region2 '
                                   'gene1,gene2')
    parser.add_argument('reads', 
                        help='Path to aligned PARIS reads file (SAM)')
    parser.add_argument('ref_seq', 
                        help='Path to reference sequence file (FASTA)')
    parser.add_argument('ref_genes', 
                        help='Path to file listing genes in reference '
                        'sequence (BED)')
    parser.add_argument('-r', '--regions', 
                        help='Genomic regions of interest (separated only by '
                        'commas, should match naming system in reference)')
    parser.add_argument('-g', '--genes', 
                        help='Genes of interest (separated only by commas)')
    parser.add_argument('out', 
                        help='Path of output')
    args = parser.parse_args(sys.argv[1:])
    return args


def main():
    
    
    # Read in arguments, parse reads, and initialize files
    args = parse_args()
    files = ReadsFiles(args.reads, args.out)
    pp.init_outputs(args.reads, files.out_sam, files.out_info, files.out_bp, 
                    files.out_aux)
    
    # Initialize reference dict analysis dict and parse reads
    ref_dict = pp.get_reference_dict(args.ref_seq, args.ref_genes)
    if args.regions and not args.genes:
        args.regions = args.regions.split(',')
        args.genes = pp.get_genes(ref_dict, args.regions)
    elif args.genes and not args.regions:
        args.genes = args.genes.split(',')
        args.regions = pp.get_regions(ref_dict, args.genes)
    elif not args.regions and not args.genes:
        args.regions = list(ref_dict.keys())
        args.genes = list(
            chain.from_iterable(
                [list(ref_dict[region]['genes'].keys()) 
                 for region in ref_dict.keys()]
            )
        )
    else:
        args.regions = args.regions.split(',')
        args.genes = args.genes.split(',')
    analysis_dict = pp.get_analysis_dict(ref_dict, args.regions, args.genes)
    reads_dict = pp.parse_reads(args.reads, ref_dict)
    

    # Analyze reads to obtain DGs
    with open(files.log, 'w') as log:
        dg_ind = 0
        ng_ind = 0
        for region in analysis_dict:
            for gene in analysis_dict[region]:
                gene_ids = [
                    read_id for (read_id, read_info) in reads_dict.items() 
                    if (read_info[4] == region) & (read_info[5] == gene) 
                    & (read_info[6] == gene)
                ]
                if len(gene_ids) > 1:
                    if len(gene_ids) > max_reads:
                        inds_samp = np.random.choice(
                            len(gene_ids), max_reads, replace=False
                        )
                        graph = gp.graph_reads(
                            [gene_ids[ind] for ind in inds_samp], 
                            reads_dict, t=min_overlap
                        )
                    else:
                        graph = gp.graph_reads(
                            gene_ids, reads_dict, t=min_overlap
                        )
                    reads_dg_dict, dg_ind = gp.cluster_graph(graph, dg_ind)
                    dg_reads_dict = da.get_preliminary_dgs(
                        reads_dict, reads_dg_dict
                    )
                    if len(gene_ids) > max_reads:
                        inds_unsamp = np.setdiff1d(
                            range(len(gene_ids)), inds_samp
                        )
                        dg_reads_dict, dg_index = da.adjust_dgs(
                            dg_reads_dict, [
                                gene_ids[ind] for ind in inds_unsamp
                            ], 
                            reads_dict, dg_ind, t=min_overlap
                        )
                    dg_filtered_dict = da.filter_dgs(dg_reads_dict)
                    dg_dict, ng_ind = da.create_dg_dict(dg_filtered_dict, reads_dict, ng_ind)
    
if __name__ == '__main__':
    main()
###############################################################################