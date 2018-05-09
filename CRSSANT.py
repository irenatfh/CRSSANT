#!/usr/bin/env python

"""
Main script for running CRSSANT analysis and discovery pipelines
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu

import sys
import time
import argparse
from analysis_pipe import subfunctions as sf, preprocessing as pp, graphing as gp, dg_analysis as da, output as op


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
    parser=argparse.ArgumentParser(description='CRSSANT analyzes PARIS reads for RNA structures')
    parser.add_argument('reads', 
                        help='Path to aligned PARIS reads file (SAM)')
    parser.add_argument('ref_seq', 
                        help='Path to reference sequence file (FASTA)')
    parser.add_argument('ref_genes', 
                        help='Path to file listing genes in reference sequence (BED)')
    parser.add_argument('-r', '--region', 
                        help='Genomic region of interest (should match naming system in reference sequence)')
    parser.add_argument('-g', '--genes', 
                        help='Genes of interest (separated only by commas)')
    parser.add_argument('out', 
                        help='Path to where results should be output')
    args = parser.parse_args(sys.argv[1:])
    return args


def main():
    args = parse_args()
    files = ReadsFiles(args.reads, args.out)
    pp.get_reference_dict(args.ref_seq, args.ref_genes)
    pp.init_outputs(args.reads, files.out_sam, files.out_info, files.out_bp, files.out_aux)
    
    # if not args.region:
    #     print('no region')
    # if args.genes:
    #     args.genes = args.genes.split(',')
    # else:
    #     print('no genes')
        

if __name__ == '__main__':
    main()