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
from itertools import chain
from analysis_pipe import subfunctions as sf, preprocessing as pp, \
graphing as gp, dg_analysis as da, output as op


class ReadsFiles(object):
    """
    """
    def __init__(self, reads, out):
        self.name = reads.split('/')[-1].split('.sam')[0]
        self.out_sam = out + self.name + '_DGNG.sam'
        self.out_info = self.out_sam.split('.sam')[0] + '_info.bed'
        self.out_bp = self.out_sam.split('.sam')[0] + '_bp.bed'
        self.out_aux = self.out_sam.split('.sam')[0] + '.aux'
        self.out_log = out + self.name + '_PARIS_analysis.log'
        
    
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
    # Read in arguments and initialize files
    args = parse_args()
    files = ReadsFiles(args.reads, args.out)
    pp.init_outputs(args.reads, files.out_sam, files.out_info, files.out_bp, 
                    files.out_aux)
    # Initialize reference dictionary analysis regions/genes
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
    print(analysis_dict)

if __name__ == '__main__':
    main()