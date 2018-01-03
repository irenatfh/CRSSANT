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
def fold_and_filter_dgs(reads_dict, dg_reads_dict, ref_seq):
    """
    Fold and filter DGs, and calculate other DG information
    
    + Filter out DGs with fewer than 2 reads
    + If there exists valid ViennaRNA folding structure, add it to the dict
      as the 'basepairs' nested array. Else, add it as 'NA.'
    + If there exists valid ViennaRNA folding structure, add to dict the number 
      of uridine crosslinking sites as an array of [UU, UC]. Else, add as 'NA.'

    Parameters
    ----------
    reads_dict : dict
        Dictionary of reads and reads information
    reads_dg_dict : dict
        Dictionary of reads and their spectral clustering DG assignments
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
        if len(dg_reads_list) >= 2:
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
                folded_struct = 'NA'
                cl_arr = 'NA'
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
                    folded_struct = 'NA'
                    cl_arr = 'NA'
                else:
                    uu, stem_len = sf.count_crosslinks(seq, fc)
                    cl_arr = uu
                    l_bps = [i for i in range(cut_point) if fc[i] == '(']
                    l_bps = read_start + np.array(l_bps)
                    r_bps = [i - cut_point 
                             for i in range(cut_point, len(fc)) 
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
def add_coverage_to_dgs(dg_folded_dict, reads_dict,
                        COV_THRESH=0):
    """
    Filter folded duplex groups by coverage threshold.
    
    Coverage is defined as c / sqrt(a*b) where
        + c = number of reads in a given DG
        + a = number of reads overlapping the left arm of the DG
        + b = number of reads overlapping the right arm of the DG
    All DG with coverage < COV_THRESH are discarded.

    Parameters
    ----------
    dg_folded_dict : int
        Dictionary of DG with valid folding structures and crosslinking sites
    reads_dict : dict
        Dictionary of reads and reads information
    COV_THRESH : coverage threshold

    Returns
    -------
    dict
        {dg:{'arm_indices': np array, 'basepairs': nested np array, 
             'num_reads': int, 'coverage': float}}

    """
    dg_cov_dict = {}
    for (dg, dg_info) in dg_folded_dict.items():
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
        if cov >= COV_THRESH:
            dg_cov_dict[dg] = {**dg_info, **{'coverage': cov}}
            
    return dg_cov_dict
        