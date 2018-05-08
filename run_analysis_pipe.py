import sys


################################################################################
import argparse
parser=argparse.ArgumentParser(
    description='CRSSANT is...')
parser.add_argument('reads_file', type=str, help="""SAM file of aligned
                    reads generated by PARIS""")
parser.add_argument('ref_seq', type=str, help="""FASTA file of reference
                    sequence (.fa file)""")
parser.add_argument('ref_bed', type=str, help="""BED file of all genes in
                    reference sequence""")
parser.add_argument('region', type=str, help="""Genomic region of interest
                    (should match naming system in reference sequence)""")
parser.add_argument('genes', type=str, help="""Genes of interest (separated 
                    only by commas)""")
parser.add_argument('results_path', type=str, help="""File path for results""")
args = parser.parse_args()

reads_sam = sys.argv[1]
ref_seq = sys.argv[2]
ref_bed = sys.argv[3]
region = sys.argv[4]
genes = list(sys.argv[5].split(','))
results_path = sys.argv[6]

import analysis_pipe
from analysis_pipe import preprocessing as pp, subfunctions as sf, \
    graphing as gp, dg_analysis as da, output as op
    
        
################################################################################
