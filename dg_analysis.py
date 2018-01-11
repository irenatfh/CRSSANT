import numpy as np
import sys
sys.path.append('/home/ihwang/software/ViennaRNA-2.4.3/interfaces/Python3')
import RNA
import subfunctions as sf


################################################################################
def get_preliminary_dgs(reads_dict, reads_dg_dict):
    """
    Organize reads into preliminary duplex groups (DGs) based on spectral 
    clustering.

    Parameters
    ----------
    reads_dict : dict
        Dictionary of reads and reads information
    reads_dg_dict : dict
        Dictionary of reads and their spectral clustering DG assignments

    Returns
    -------
    dict
        {dg:[read ids]}

    """
    dgs_list = set(reads_dg_dict.values())
    dg_reads_dict = {}
    for dg in dgs_list:
        dg_reads_list = [i[0] for i in reads_dg_dict.items() if i[1] == dg]
        dg_reads_dict[dg] = dg_reads_list
    return dg_reads_dict


################################################################################
def adjust_dgs(dg_reads_dict, reads_ids, reads_dict, dg_index):
    """
    Adjust duplex groups by adding in the reads that were not sampled.
    
    Each read is checked for overlap with existing DGs.
    + If the read has no overlap with any DGs, add as DG.
    + If the read is equidistant between at least two DGs, add as new DG.
    Note that this process updates the dg_reads_dict throughout the loop.

    Parameters
    ----------
    dg_reads_dict : dict
        Dictionary of reads and their spectral clustering DG assignments
    reads_ids : np array
        Numpy array of reads indices
    reads_dict : dict
        Dictionary of reads and reads information
    dg_index : int

    Returns
    -------
    dict, int
        Updated dg_reads_dict, DG index to start from

    """
    for (index, read_id) in enumerate(reads_ids):
        read_inds = reads_dict[read_id][:-2]
        dg_overlaps = {}
        for (dg, dg_reads_list) in dg_reads_dict.items():
            dg_reads_info = np.median(
                np.array([reads_dict[i][:-2] for i in dg_reads_list]), axis=0)
            dg_inds = np.array([int(dg_reads_info[0]), int(dg_reads_info[1]),
                                int(dg_reads_info[2]), int(dg_reads_info[3])])
            overlaps = sf.get_overlaps(read_inds, dg_inds)
            sum_ratio = overlaps[0]/overlaps[3] + overlaps[1]/overlaps[4]
            if sum_ratio > 0:
                dg_overlaps[dg] = sum_ratio
        if len(dg_overlaps) == 0:
            dg_reads_dict[dg_index] = [read_id]
            dg_index += 1
        else:
            max_overlap = max(dg_overlaps.values())
            max_dgs_list = [i for i in dg_overlaps.keys() 
                            if dg_overlaps[i] == max_overlap]
            if len(max_dgs_list) == 1:
                dg_reads_dict[list(dg_overlaps.keys())[0]].append(read_id)
            else:
                dg_reads_dict[dg_index] = [read_id]
                dg_index += 1
                
    return dg_reads_dict, dg_index


################################################################################
def fold_dgs(dg_reads_dict, reads_dict, ref_seq):
    """
    Fold DGs and calculate other DG information
    
    + If there exists valid ViennaRNA folding structure, add it to the dict
      as the 'basepairs' nested array. Else, add it as 'NA.'
    + If there exists valid ViennaRNA folding structure, add to dict the number 
      of uridine crosslinking sites as an array of [UU, UC]. Else, add as 'NA.'

    Parameters
    ----------
    dg_reads_dict : dict
        Dictionary of reads and their spectral clustering DG assignments
    reads_dict : dict
        Dictionary of reads and reads information
    ref_seq : str
        Reference sequence

    Returns
    -------
    dict
        {dg:{'arm_indices': np array, 'basepairs': nested np array, 
             'num_reads': int}}

    """
    dg_folded_dict = {}
    for (dg, dg_reads_list) in dg_reads_dict.items():
        dg_reads_info = np.median(
            np.array([reads_dict[i][:-2] for i in dg_reads_list]), axis=0)
        dg_inds = np.array([int(dg_reads_info[0]), int(dg_reads_info[1]),
                            int(dg_reads_info[2]), int(dg_reads_info[3])])
        l_arm = ref_seq[dg_inds[0] : dg_inds[1] + 1]
        r_arm = ref_seq[dg_inds[2] : dg_inds[3] + 1]
        read_start = dg_inds[0]
        cut_point = len(l_arm)
        seq = l_arm + r_arm
        # Attempt folding the sequence
        res = RNA.fold_compound(l_arm + '&' + r_arm)
        [fc, mfe] = res.mfe_dimer()
        l_symbols = [i for i in fc[:cut_point] if i != '.']
        r_symbols = [i for i in fc[cut_point:] if i != '.']
        # Check for fatal helix folding results
        if (len(l_symbols) < 1 or len(r_symbols) < 1) or \
           (l_symbols[-1]+r_symbols[0] != '()'):
            folded_struct = np.zeros((2,1))
            cl_arr = np.zeros(2)
        else:
            # Truncate helix arms to fix folding, if necessary
            if (set(fc[:cut_point]) != set('.(')):
                off_inds = [i for i in range(cut_point) if 
                            fc[:cut_point][i] == ')']
                l_arm = l_arm[off_inds[-1]+1 : ]
                read_start += off_inds[-1] + 1
            if (set(fc[cut_point:]) != set('.)')):
                off_inds = [i for i in range(len(r_arm)) if 
                            fc[cut_point:][i] == '(']
                r_arm = r_arm[ : off_inds[0]]
            # Re-attempt helix folding
            seq = l_arm + r_arm
            cut_point = len(l_arm)
            res = RNA.fold_compound(l_arm + '&' + r_arm)
            [fc, mfe] = res.mfe_dimer()
            # Fixing the helix folding by truncation was unsuccessful
            if (set(fc[:cut_point]) != set('.(')) or \
               (set(fc[cut_point:]) != set('.)')):
                folded_struct = np.zeros((2,1))
                cl_arr = np.zeros(2)
            else:
                uu, stem_len = sf.count_crosslinks(seq, fc)
                cl_arr = uu
                l_bps = [i for i in range(cut_point) if fc[i] == '(']
                l_bps = read_start + np.array(l_bps)
                r_bps = [i - cut_point for i in range(cut_point, len(fc)) 
                         if fc[i] == ')']
                r_bps = dg_inds[2] + np.array(r_bps)
                r_bps = r_bps[::-1]
                folded_struct = np.array([l_bps, r_bps])
        dg_folded_dict[dg] = {'arm_indices': dg_inds, 
                              'basepairs': folded_struct,
                              'num_reads': len(dg_reads_list),
                              'cl_sites': cl_arr}
                        
    return dg_folded_dict


################################################################################
def filter_dgs(dg_dict):
    """
    Filter out DGs with fewer than 2 reads.

    Parameters
    ----------
    dg_dict : dict
        DG dictionary {DG:{DG atttributes}}, where one of the attributes
        should be 'num_reads'

    Returns
    -------
    dict
        Filtered DG dictionary

    """
    filtered_dict = {}
    for (dg, dg_reads_list) in dg_dict.items():
        num_reads = dg_reads_list['num_reads']
        if num_reads > 2:
            filtered_dict[dg] = dg_reads_list
            
    return filtered_dict


################################################################################
def add_coverage_to_dgs(dg_dict, reads_dict):
    """
    Add coverage threshold to DG dictionary.
    
    Coverage is defined as c / sqrt(a*b) where
        + c = number of reads in a given DG
        + a = number of reads overlapping the left arm of the DG
        + b = number of reads overlapping the right arm of the DG

    Parameters
    ----------
    dg_dict : dict
        DG dictionary {DG:{DG atttributes}}, where the attributes
        should include valid folding structures and crosslinking sites
    reads_dict : dict
        Dictionary of reads and reads information

    Returns
    -------
    dict
        {dg:{'arm_indices': np array, 'basepairs': nested np array, 
             'num_reads': int, 'coverage': float}}

    """
    dg_cov_dict = {}
    for (dg, dg_info) in dg_dict.items():
        dg_inds = dg_info['arm_indices']
        num_reads_dg = dg_info['num_reads']
        num_reads_l = 0
        num_reads_r = 0
        for read in reads_dict:
            read_inds = reads_dict[read][:-2]
            overlaps = sf.get_overlaps(dg_inds, read_inds)
            overlap_l = overlaps[0]
            overlap_r = overlaps[1]
            if overlap_l > 0:
                num_reads_l += 1
            if overlap_r > 0:
                num_reads_r += 1
        cov = num_reads_dg / np.sqrt(num_reads_l * num_reads_r)
        dg_cov_dict[dg] = {**dg_info, **{'coverage': cov}}
            
    return dg_cov_dict


################################################################################
def add_ngs_to_dgs(dg_dict, dg_reads_dict, reads_dict, ng_index):
    """
    Add non-overlapping groups (NGs) to DG dict.

    Parameters
    ----------
    dg_dict : dict
        DG dictionary {DG:{DG atttributes}}
    dg_reads_dict : dict
        Dictionary of reads and their spectral clustering DG assignments
    reads_dict : dict
        Dictionary of reads and reads information 

    Returns
    -------
    dict
        dg_dict updated with NG attribute

    """
    ng_dict = {}
    for dg in dg_dict.keys():
        dg_reads_list = dg_reads_dict[dg]
        dg_start = min(np.array([reads_dict[i][0] for i in dg_reads_list]))
        dg_stop = max(np.array([reads_dict[i][3] for i in dg_reads_list]))
        if len(ng_dict) == 0:
            ng_dict[ng_index] = [dg]
            dg_dict[dg] = {**dg_dict[dg], **{'NG': ng_index}}
            ng_index += 1
        else:
            ng_flag = 0
            for (ng, ng_dgs) in ng_dict.items():
                overlaps_arr = np.zeros(len(ng_dgs))
                for (i, dg_other) in enumerate(ng_dgs):
                    dg_other_start = min(np.array(
                        [reads_dict[i][0] for i in dg_reads_dict[dg_other]]))
                    dg_other_stop = max(np.array(
                        [reads_dict[i][3] for i in dg_reads_dict[dg_other]]))
                    overlap = min(dg_stop, dg_other_stop) - max(dg_start, 
                                                                dg_other_start)
                    if overlap >= 0:
                        overlaps_arr[i] = 1
                if sum(overlaps_arr) == 0:
                    ng_dgs.append(dg)
                    dg_dict[dg] = {**dg_dict[dg], **{'NG': ng}}
                    ng_flag += 1
                    break
            if ng_flag == 0:
                ng_dict[ng_index] = [dg]
                dg_dict[dg] = {**dg_dict[dg], **{'NG': ng_index}}
                ng_index += 1
                
    return dg_dict, ng_index
        
            