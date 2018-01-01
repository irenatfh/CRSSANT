import numpy as np
import sys
sys.path.append('/home/ihwang/software/ViennaRNA-2.4.3/interfaces/Python3')
import RNA
import subfunctions as sf


################################################################################
def get_preliminary_dgs(reads_dict, reads_dg_dict):
    dgs_list = set(reads_dg_dict.values())
    dg_reads_dict = {}
    for dg in dgs_list:
        dg_reads_list = [i[0] for i in reads_dg_dict.items() if i[1] == dg]
        dg_reads_dict[dg] = dg_reads_list
    return dg_reads_dict


################################################################################
def filter_dgs(reads_dict, dg_reads_dict, ref_seq):
    dg_filtered_dict = {}
    for (dg, dg_reads_list) in dg_reads_dict.items():
        dg_reads_info = np.median(
            np.array([reads_dict[i][:-2] for i in dg_reads_list]), axis=0)
        if len(dg_reads_list) > 1:
            dg_inds = np.array([int(dg_reads_info[0]), int(dg_reads_info[1]),
                                int(dg_reads_info[2]), int(dg_reads_info[3])])
            l_arm = ref_seq[dg_inds[0] : dg_inds[1] + 1]
            r_arm = ref_seq[dg_inds[2] : dg_inds[3] + 1]
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
                continue
            else:
                # Truncate helix arms to fix folding, if necessary
                if (set(fc[:cut_point]) != set('.(')):
                    print('yodel! left!')
                    off_inds = [i for i in range(cut_point) if 
                                fc[:cut_point][i] == ')']
                    l_arm = l_arm[off_inds[-1]+1 : ]
                if (set(fc[cut_point:]) != set('.)')):
                    print('yodel! right!')
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
                    continue
                else:
                    uu, stem_len = sf.count_crosslinks(seq, fc, mfe)
                    if sum(uu) > 0:
                        print('ViennaRNA folding results for DG %s, %s reads' %(dg, len(dg_reads_list)))
                        print(l_arm + '&' + r_arm + '\n' + \
                              fc[:cut_point] + '&' + fc[cut_point:], mfe)
                        print("%s U-U and %s U-C cross-linking sites found in %s-bp stem" \
                              %(uu[0], uu[1], stem_len))
            # dg_filtered_dict[dg] = [l_arm, r_arm]
    return dg_filtered_dict