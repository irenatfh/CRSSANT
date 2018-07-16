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
import preprocessing as pp, graphing as gp, dg_analysis as da, \
structure_discovery as sd, output as op


# Global variables
max_reads = 2000
min_overlap = 0.4
tests = 100
min_rank = 0.75


class ReadsFiles(object):
    """
    """
    def __init__(self, reads, out, region_gene_str):
        self.name = reads.split('/')[-1].split('.sam')[0] + '_CRSSANT_' + \
                    region_gene_str
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
    parser.add_argument('-c', '--chimeric', 
                        help='Path to chimeric PARIS reads file (SAM)')
    parser.add_argument('-j', '--junction', 
                        help='Path to junction PARIS reads file (JUNCTION)')
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
    
    
    # Read in arguments
    args = parse_args()
    ref_dict = pp.get_reference_dict(args.ref_seq, args.ref_genes)
    if args.regions and not args.genes:
        region_gene_str = 'r%s_g' %args.regions
        args.regions = args.regions.split(',')
        args.genes = pp.get_genes(ref_dict, args.regions)
    elif args.genes and not args.regions:
        region_gene_str = 'r_g%s' %args.genes
        args.genes = args.genes.split(',')
        args.regions = pp.get_regions(ref_dict, args.genes)
    elif not args.regions and not args.genes:
        region_gene_str = 'r_g'
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
        region_gene_str = 'r%s_g%s' %(
            ''.join(args.regions), ''.join(args.genes)
        )
    if args.junction and not args.chimeric:
        print('.JUNCTION file must be accompanied by chimeric reads file')
        sys.exit()
    if args.chimeric:
        if args.junction:
            args.reads = pp.combine_aligned_and_chimeric_reads(
                args.reads, args.chimeric, args.junction
            )
        else:
            args.reads = pp.combine_aligned_and_chimeric_reads(
                args.reads, args.chimeric
            )


    # Parse analysis dictionary and reads
    analysis_dict = pp.get_analysis_dict(ref_dict, args.regions, args.genes)
    reads_dict = pp.parse_reads(args.reads, ref_dict)
    
    
    # Initialize output files
    files = ReadsFiles(args.reads, args.out, region_gene_str)
    pp.init_outputs(args.reads, files.out_sam, files.out_info, files.out_bp, 
                    files.out_aux)
    
    
    # Run analysis pipeline
    with open(files.log, 'w') as log:
        start = datetime.datetime.now()
        dg_ind = 0
        for region in analysis_dict:
            region_seq = ref_dict[region]['sequence']
            for gene in analysis_dict[region]:
                gene_start = datetime.datetime.now()
                gene_inds = ref_dict[region]['genes'][gene]
                ng_ind = 0
                gene_ids = [
                    read_id for (read_id, read_info) in reads_dict.items() 
                    if (read_info[4] == region) & (read_info[5] == gene) 
                    & (read_info[6] == gene)
                ]
                log.write(
                    'Analyzing %s reads spanning gene %s\n' 
                    %(len(gene_ids), gene)
                )
                if len(gene_ids) > 1:
                    
                    
                    # Analyze reads to obtain DGs
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
                    
                    
                    # Output SAM and _info.bed files after DG analysis
                    op.write_info(files.out_info, dg_dict, region)
                    op.write_dg_ng_sam(
                        args.reads, files.out_sam, dg_reads_dict, dg_dict
                    )
    
    
                    # Analyze DGs to discover new secondary structures
                    stem_dict = sd.get_stem_dict(dg_dict, region_seq)
                    op.write_aux(
                        files.out_aux, dg_dict, dg_reads_dict, reads_dict, 
                        stem_dict
                    )
                    op.write_bp(files.out_bp, stem_dict, region)
                gene_stop = datetime.datetime.now()
                log.write('Gene analysis time: %s\n' %(gene_stop - gene_start))
        stop = datetime.datetime.now()
        log.write('\nTotal analysis time: %s\n' %(stop - start))

if __name__ == '__main__':
    main()
###############################################################################