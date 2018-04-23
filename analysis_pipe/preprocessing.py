from . import subfunctions as sf


################################################################################
def parse_reference_bed(ref_file, region, genes):
    """
    Parse a reference BED file into a dictionary.

    Parameters
    ----------
    ref_file : str
        Reference file path
    region : str
        Genomic region of interest
    genes : list
        Genes of interest

    Returns
    -------
    dict
        {gene:[start position, stop position]}

    """
    ref_dict = {}
    with open(ref_file, 'r') as f:
        for line in f:
            data = line.split('\t')
            data_region = data[0]
            gene = data[3]
            if (data_region == region) and (gene in genes):
                pos_start = int(data[1]) - 1  # Biology is 1-indexed
                pos_stop = int(data[2]) - 1
                ref_dict[gene] = [pos_start, pos_stop]
            
    return ref_dict


################################################################################
def get_reference_seq(ref_file, region):
    """
    Parse a reference sequence file into a dictionary

    Parameters
    ----------
    ref_file : str
        Reference sequence file path
    region : str
        Genomic region of interest

    Returns
    -------
    dict
        {genomic region:sequence}

    """
    ref_seq = ''
    region_flag = 0
    with open(ref_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                ref_key = line.split('>')[-1][:-1]
                if ref_key == region:
                    region_flag = 1
                else:
                    region_flag = 0
            else:
                if region_flag == 1:
                    ref_seq += line[:-1]
                
    return ref_seq


################################################################################
def parse_reads(reads_file, ref_dict, output_sam='test.sam'):
    """
    Parse a reads SAM file into a dictionary.

    Parameters
    ----------
    reads_file : str
        Reads file path
    ref_dict : str
        Reference file dictionary
    output_sam : str
        Output file name

    Returns
    -------
    dict
        {read ID:[L start, L stop, R start, R stop, L arm RNA, R arm RNA]}

    """
    reads_dict = {}
    with open(reads_file, 'r') as f_read:
        for line in f_read:
            if line[0] != '@':
                data = line.split('\n')[0].split('\t')
                read_id = data[0]
                pos_align = int(data[3])
                cigar_str = data[5]
                if len(data) > 19:
                    xg = data[19]
                else:
                    xg = 'XG:i:0'
                if (xg == 'XG:i:0') or (xg == 'XG:i:1'):
                    cigar_ops, cigar_lens = sf.process_cigar(cigar_str)
                    if cigar_ops == ['S', 'M', 'N', 'M', 'S']:
                        l_pos_start = pos_align - 1  # biology is 0-indexed
                        l_pos_stop = l_pos_start + cigar_lens[1] - 1
                        r_pos_start = l_pos_stop + cigar_lens[2] + 1
                        r_pos_stop = r_pos_start + cigar_lens[3] - 1
                        l_arm_rna = [rna for (rna, [rna_start, rna_stop]) 
                                     in ref_dict.items() if l_pos_start 
                                     in range(rna_start, rna_stop)]
                        r_arm_rna = [rna for (rna, [rna_start, rna_stop]) 
                                     in ref_dict.items() if r_pos_start 
                                     in range(rna_start, rna_stop)]
                        if (len(l_arm_rna) > 0) and (len(r_arm_rna) > 0):
                            reads_dict[read_id] = [
                                l_pos_start, l_pos_stop, 
                                r_pos_start, r_pos_stop,
                                l_arm_rna[0], r_arm_rna[0]]
                            
    return reads_dict


################################################################################
def init_output_files(reads_sam,
                      output_sam, output_info_bed, output_bp_bed, output_aux):
    with open(reads_sam, 'r') as f_r, open(output_sam, 'w') as f_w:
        for line in f_r:
            if line[0] == '@':
                f_w.write(line)
            else:
                pass
            
    with open(output_info_bed, 'w') as f_w:
        pass
    
    with open(output_bp_bed, 'w') as f_w:
        pass
    
    with open(output_aux, 'w') as f_w:
        header = ['DG_coverage',
                  'UU_cl,UC_cl,helix_length', 
                  'L_start_start,L_start_stop,L_start_std', 
                  'L_stop_start,L_stop_stop,L_stop_std',
                  'R_start_start,R_start_stop,R_start_std', 
                  'R_stop_start,R_stop_stop,R_stop_std']
        f_w.write('\t'.join(header) + '\n')
        
        