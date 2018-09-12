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
from joblib import Parallel, delayed
import preprocessing as pp, run_analysis as ra, output as op


# Global variables
max_reads = 10000
min_overlap = 0.4
bin_width = 40
threads = 4


class ReadsFiles(object):
    """
    """
    def __init__(self, reads, out, region_gene_str):
        self.name = reads.split('/')[-1].split('.sam')[0] + '_CRSSANT_' + \
                    region_gene_str
        self.out_sam = out + self.name + '.sam'
        self.out_dg = out + self.name + '_dg.bed'
        self.out_sg_arc = out + self.name + '_sg_arc.bed'
        self.out_sg_bp = out + self.name + '_sg_bp.bed'
        self.out_sg = out + self.name + '_sg.aux'
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
                        help='Path to file listing genes in reference'
                        ' sequence (BED)')
    parser.add_argument('-c', '--chimeric',  
                        help='Path to chimeric PARIS reads file (SAM)')
    parser.add_argument('-r', '--regions', 
                        help='Genomic regions of interest (separated only by'
                        ' commas, should match naming system in reference)')
    parser.add_argument('-g', '--genes', 
                        help='Genes of interest (separated only by commas)')
    parser.add_argument('-t', '--threads', 
                        help='Number of threads')
    parser.add_argument('out', help='Path of output')
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
        region_gene_str = 'r%s_g%s' %(
            ''.join(args.regions), ''.join(args.genes)
        )
        args.regions = args.regions.split(',')
        args.genes = args.genes.split(',')
    
    
    # Initialize output files
    files = ReadsFiles(args.reads, args.out, region_gene_str)
    pp.init_outputs(
        args.reads, files.out_sam, 
        files.out_dg, files.out_sg_arc, files.out_sg_bp, files.out_sg
    )
    
    with open(files.log, 'w') as log:
        if args.chimeric:
            print('Combining aligned and chimeric reads\n')
            log.write('Combining aligned and chimeric reads\n')
            chimeric_start = datetime.datetime.now()
            args.reads = pp.combine_aligned_and_chimeric_reads(
                args.reads, args.chimeric
            )
            chimeric_stop = datetime.datetime.now()
            log.write(
                'Chimeric read combination time: %s\n\n' \
                %(chimeric_stop - chimeric_start)
            )
    
    
        # Parse analysis dictionary and reads
        print('Preparing reads for analysis\n')
        log.write('Preparing reads for analysis\n')
        prepare_start = datetime.datetime.now()
        regions_genes_dict = pp.get_regions_genes_dict(
            ref_dict, args.regions, args.genes
        )
        reads_dict, regions_readgenes_dict = pp.parse_reads(
            args.reads, ref_dict
        )
        analysis_dict = pp.get_analysis_dict(
            regions_genes_dict, regions_readgenes_dict
        )
        prepare_stop = datetime.datetime.now()
        log.write(
            'Read preparation time: %s\n\n' %(prepare_stop - prepare_start)
        )

        
        # Run analysis pipeline
        start = datetime.datetime.now()
        arg_instances = [
            (
                l_gene, r_gene, region, reads_dict, ref_dict, 
                log, max_reads, min_overlap, bin_width
            ) for (l_gene, r_gene), region in analysis_dict.items()
        ]
        results = Parallel(n_jobs=threads,  prefer="threads")(
            map(delayed(ra.run_analysis), arg_instances)
        )
        dg_offset = 0
        for (iter_args, iter_results) in zip(arg_instances, results):
            l_gene, r_gene = iter_args[:2]
            [
                dg_reads_dict, dg_dict, sg_reads_dict, sg_dict, 
                region, gene_time, num_reads
            ] = iter_results
            if l_gene == r_gene:
                log.write(
                    'Analyzed %s reads spanning gene %s\n' \
                    'Analysis time: %s\n' \
                    %(num_reads, l_gene, gene_time)
                )
            else:
                log.write(
                    'Analyzed %s reads spanning genes %s and %s\n' \
                    'Analysis time: %s\n' \
                    %(num_reads, l_gene, r_gene, gene_time)
                )
            if dg_dict:
                # Output SAM file and dg file
                op.write_dg_ng_sam(
                    args.reads, files.out_sam, 
                    dg_reads_dict, dg_dict, dg_offset
                )
                op.write_dg(files.out_dg, dg_dict, region, dg_offset)
                if sg_dict:
                    # Write remaining output files
                    op.write_sg_arc(
                        files.out_sg_arc, sg_dict, region, dg_offset
                    )
                    op.write_sg_bp(
                        files.out_sg_bp, sg_dict, region, dg_offset
                    )
                    op.write_sg(
                        files.out_sg, sg_dict, sg_reads_dict, dg_dict, 
                        reads_dict, dg_offset
                    )
                dg_offset += max(dg_dict.keys()) + 1
        stop = datetime.datetime.now()
        log.write('\nTotal analysis time: %s\n' %(stop - start))
        print('\nTotal analysis time: %s\n' %(stop - start))
    print('Analysis complete')

    
if __name__ == '__main__':
    main()