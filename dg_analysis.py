import numpy as np
import sys
sys.path.append('/home/ihwang/software/ViennaRNA-2.4.3/interfaces/Python3')
import RNA

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
            dg_filtered_dict[dg] = [l_arm, r_arm]
    return dg_filtered_dict