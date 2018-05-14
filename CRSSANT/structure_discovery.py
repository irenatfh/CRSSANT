# This file is part of CRSSANT:
# Computational RNA Secondary Structure Analysis using Network Techniques
#
###############################################################################
"""
This module is a collection of functions that perform DG analysis tasks
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu


import numpy as np
import sys
import subfunctions as sf
sys.path.append('/home/ihwang/software/ViennaRNA-2.4.3/interfaces/Python3')
import RNA


def fold_dgs(dg_reads_dict, reads_dict, ref_seq):
    """
    Fold DGs and calculate other DG information
    
    + If there exists valid ViennaRNA folding structure (i.e. stem of at least
      3 bp), add it to the dictionary with arm indices, basepairs, number of
      reads and uridine crosslinking sites.
    + If the initial ViennaRNA folding structure is invalid, try truncating the
      DG.

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
             'num_reads': int, 'cl_sites': np array}}

    """
    dg_folded_dict = {}
    for (dg, dg_reads_list) in dg_reads_dict.items():
        dg_reads_info = np.median(
            np.array([reads_dict[i][ : -2] for i in dg_reads_list]), axis=0)
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
        l_symbols = [i for i in fc[ : cut_point] if i != '.']
        r_symbols = [i for i in fc[cut_point : ] if i != '.']
        # Check for fatal helix folding results
        if (len(l_symbols) < 2 or len(r_symbols) < 2) or \
           (set(l_symbols) != set('(') or set(r_symbols) != set(')')):
            folded_struct = np.zeros((2,1), dtype=np.int)
            cl_sites = np.zeros(3, dtype=np.int)
        else:
            # Check for folding results that might be fixable by truncation
            if set(fc[:cut_point]).issubset(set('.(')) is False:
                off_inds = [i for i in range(cut_point) if fc[i] == ')']
                l_arm = l_arm[off_inds[-1] + 1 : ]
                read_start += off_inds[-1] + 1
            if set(fc[cut_point:]).issubset(set('.)')) is False:
                off_inds = [i for i in range(len(r_arm)) if 
                            fc[cut_point : ][i] == '(']
                r_arm = r_arm[ : off_inds[0]]
            # Re-attempt helix folding/fold helix again if no truncation
            cut_point = len(l_arm)
            seq = l_arm + r_arm
            res = RNA.fold_compound(l_arm + '&' + r_arm)
            [fc, mfe] = res.mfe_dimer()
            # Check if truncation was successful
            l_symbols = [i for i in fc[ : cut_point] if i != '.']
            r_symbols = [i for i in fc[cut_point :] if i != '.']
            # If truncation was unsuccessful output null crosslinking sites
            if (len(l_symbols) < 2 or len(r_symbols) < 2) or \
               (set(l_symbols) != set('(') or set(r_symbols) != set(')')):
                folded_struct = np.zeros((2,1), dtype=np.int)
                cl_sites = np.zeros(3, dtype=np.int)
            else:
                cl_sites = sf.count_crosslinks(seq, fc)
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
                              'cl_sites': cl_sites}
                        
    return dg_folded_dict