import subfunctions as sf


################################################################################
def parse_reference_bed(ref_file):
    """
    Parse a reference BED file into a dictionary.

    Parameters
    ----------
    ref_file : str
        Reference file path

    Returns
    -------
    dict
        {genomic region:{rna:[start position, stop position]}}

    """
    ref_dict = {}
    with open(ref_file, 'r') as f:
        for line in f:
            data = line.split('\t')
            region = data[0]
            if region not in ref_dict.keys():
                ref_dict[region] = {}
            pos_start = int(data[1]) - 1  # Biology is 1-indexed
            pos_stop = int(data[2]) - 1
            rna = data[3]
            ref_dict[region][rna] = [pos_start, pos_stop]
            
    return ref_dict


################################################################################
def get_reference_seq(ref_file):
    """
    Parse a reference sequence file into a dictionary

    Parameters
    ----------
    ref_file : str
        Reference sequence file path

    Returns
    -------
    dict
        {genomic region:sequence}

    """
    ref_seq_dict = {}
    with open(ref_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                ref_key = line.split('>')[-1][:-1]
            else:
                ref_seq_dict[ref_key] = line[:-1]
                
    return ref_seq_dict


################################################################################
def parse_reads(reads_file, ref_dict, output_sam):
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
    with open(reads_file, 'r') as f_read, open(output_sam, 'w') as f_write:
        for line in f_read:
            if line[0] != '@':
                data = line.split('\n')[0].split('\t')
                read_id = data[0]
                region = data[2]
                if region not in reads_dict.keys():
                    reads_dict[region] = {}
                pos_align = int(data[3])
                cigar_str = data[5]
                xg = data[19]
                if (xg == 'XG:i:0') or (xg == 'XG:i:1'):
                    cigar_ops, cigar_lens = sf.process_cigar(cigar_str)
                    if cigar_ops == ['S', 'M', 'N', 'M', 'S']:
                        l_pos_start = pos_align - 1  # biology is 0-indexed
                        l_pos_stop = l_pos_start + cigar_lens[1] - 1
                        r_pos_start = l_pos_stop + cigar_lens[2] + 1
                        r_pos_stop = r_pos_start + cigar_lens[3] - 1
                        l_arm_rna = [rna for (rna, [rna_start, rna_stop]) 
                                     in ref_dict[region].items() if l_pos_start 
                                     in range(rna_start, rna_stop)]
                        r_arm_rna = [rna for (rna, [rna_start, rna_stop]) 
                                     in ref_dict[region].items() if r_pos_start 
                                     in range(rna_start, rna_stop)]
                        if (len(l_arm_rna) > 0) and (len(r_arm_rna) > 0):
                            reads_dict[region][read_id] = [
                                l_pos_start, l_pos_stop, 
                                r_pos_start, r_pos_stop,
                                l_arm_rna[0], r_arm_rna[0]]
            else:
                f_write.write(line)
                            
    return reads_dict