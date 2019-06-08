#!/usr/bin/env python
#
# This file is part of CRSSANT:
# Crosslinked RNA Secondary Structure Analysis using Network Techniques
#
###############################################################################

"""
Main script for running CRSSANT analysis and discovery pipelines
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu

import sys
import argparse
import datetime
import numpy as np
from itertools import chain
import multiprocessing
import preprocessing as pp, run_analysis as ra, output as op
        
    
def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser=argparse.ArgumentParser(
        description='CRSSANT groups PARIS reads into DGs and SGs',
        epilog='python ./CRSSANT.py reads.sam ref.fa ref.bed out gene1,gene2'
    )
    parser.add_argument(
        'reads', 
        help='Path to aligned PARIS reads file (SAM)'
    )
    parser.add_argument(
        'ref_seq',
        help='Path to reference sequence file (FASTA)'
    )
    parser.add_argument(
        'ref_genes',
        help='Path to file listing genes in reference sequence (BED)'
    )
    parser.add_argument(
        '-out',
        help='Path of output'
    )
    parser.add_argument(
        '-chimeric',
        help='Path to chimeric PARIS reads file (SAM)'
    )
    parser.add_argument(
        '-genes',
        help='Genes of interest (separated only by commas)'
    )
    parser.add_argument(
        '-cluster',
        help='Clustering method, valid options are "cliques" or "spectral". '
        'Default is spectral'
    )
    parser.add_argument(
        '-t_o',
        help='Overlap threshold (any float between 0 and 1, inclusive). '
        'Default is 0.5 for spectral clustering, and 0.1 for cliques-finding'
    )
    parser.add_argument(
        '-t_eig',
        help='Eigenratio threshold (any positive number). Default is 5 for '
        'spectral clustering. Not needed for cliques-finding'
    )
    parser.add_argument(
        '-n', 
        help='Number of threads. Default is 8')
    args = parser.parse_args(sys.argv[1:])
    return args
###############################################################################

def main():
    args = parse_args()
    ##### Preprocessing
    # 1a) Read in reference dictionary and parse reads
    ref_dict = pp.get_reference_dict(args.ref_seq, args.ref_genes)
    # 1b) Check for chimeric reads
    if args.chimeric:
        reads = pp.combine_aligned_and_chimeric_reads(args.reads, args.chimeric)
    reads_dict, regions_readgenes_dict = pp.parse_reads(args.reads, ref_dict)
    # 1c) Check that SAM file contains only one genomic region
    if len(regions_readgenes_dict.keys()) == 1:
        region = list(regions_readgenes_dict.keys())[0]
        analysis_genes = regions_readgenes_dict[region]
    else:
        sys.exit(
            'PIPELINE ABORTED: CRSSANT assumes that the input SAM file '
            'contains reads mapped to only one genomic region'
        )
    # 2) Check genes arguments
    if args.genes:
        l_gene, r_gene = args.genes.split(',')
        if (l_gene, r_gene) in regions_readgenes_dict.get(region):
            analysis_genes = set([(l_gene, r_gene)])
        else:
            sys.exit(
                'PIPELINE ABORTED: SAM file does not contain any reads with '
                'left arm mapped to gene %s and right arm mapped to gene %s'
                %(l_gene, r_gene)
            )
    # 3) Check clustering arguments
    if args.cluster == 'cliques':
        if not args.t_o:
            args.t_o = 0.1
        else:
            args.t_o = float(args.t_o)
        args.t_eig = None
    elif (not args.cluster) or args.cluster == 'spectral':
        args.cluster = 'spectral'
        if not args.t_o:
            args.t_o = 0.5
        else:
            args.t_o = float(args.t_o)
        if not args.t_eig:
            args.t_eig = 5
        else:
            args.t_eig = float(args.t_eig)
    else:
        sys.exit(
            'PIPELINE ABORTED: Valid clustering arguments: "cliques" or '
            '"spectral"'
        )
    # 4) Check number of threads
    if not args.n:
        args.n = 8
    else:
        args.n = int(args.n)

        
    # Run pipeline
    arg_instances = [
        (
            l_gene, r_gene, reads_dict, region, ref_dict[region]['sequence'], 
            args.cluster, args.t_o, args.t_eig,
            args.out, args.reads
        ) for (l_gene, r_gene) in analysis_genes
    ]
    pool = multiprocessing.Pool(processes=args.n)
    results = pool.map(ra.run_analysis, arg_instances)

    
if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()