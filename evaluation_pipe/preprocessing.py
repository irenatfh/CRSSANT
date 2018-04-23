import numpy as np
from . import subfunctions as sf


################################################################################
def get_gs_dict(gs_bp_bed, ref_dict, region):
    """
    From a gold standard (GS) basepairing BED file, create a dictionary of GS
    structures.

    Parameters
    ----------
    gs_bp_bed : str
        file name of GS basepairing BED file
    region : str
        genomic region of interest
    ref_dict : dict
        dictionary of reference in format {gene: [start_ind, stop_ind]}

    Returns
    -------
    dict
        {struct: np.array([left_start_ind, left_stop_ind,
                           right_start_ind, right_stop_ind]
                         ), 
                 [left_arm_gene, right_arm_gene]
                 }

    """
    gs_structs_dict = {}
    left_arms = []
    right_arms = []
    with open(gs_bp_bed, 'r') as f:
        for line in f:
            if 'graphType=' in line:
                pass
            else:
                data = line.split('\t')
                bp_region = data[0]
                if bp_region == region:
                    bp_left = int(data[1])
                    bp_right = int(data[2])
                    left_arms.append(bp_left)
                    right_arms.append(bp_right)
    left_arms = np.asarray(left_arms);
    right_arms = np.asarray(right_arms)
    arm_limits = 1 + np.array(
        [i for i in range(len(left_arms) - 1) if (np.abs(np.diff(left_arms)[i]) > 3)
         and ((np.abs(np.diff(right_arms)[i]) > 3))]
    )
    arm_limits = np.concatenate(
        (np.concatenate((np.array([0]), 
                         arm_limits)
                       ),
         np.array([len(left_arms)])
        )
    )
    for i in range(len(arm_limits) - 1):
        arm_start = arm_limits[i]
        arm_stop = arm_limits[i+1]
        left_arm = -1 + np.array([left_arms[j] for j in range(arm_start, arm_stop)]
                                )  # python is 0-indexed
        right_arm = -1 + np.array(
            sorted([right_arms[j] for j in range(arm_start, arm_stop)])
        )
        struct_inds = np.array(
            [left_arm[0], left_arm[-1], right_arm[0], right_arm[-1]]
        )
        struct_ranges = []
        for j in range(2):
            for (gene, gene_inds) in ref_dict.items():
                gene_range = range(gene_inds[0], gene_inds[1] + 1)
                if (struct_inds[j*2] in gene_range):
                    struct_ranges.append(gene)
                    break
        if len(struct_ranges) == 2:
            gs_structs_dict[i] = [struct_inds, struct_ranges]
    return gs_structs_dict


################################################################################
def get_dg_dict(dg_info_bed, ref_dict, ref_seq):
    """
    From a duplex group (DG) info BED file, create a dictionary of DG
    structures.

    Parameters
    ----------
    dg_info_bed : str
        file name of DG info BED file
    ref_dict : dict
        dictionary of reference in format {gene: [start_ind, stop_ind]}
    ref_seq : str
        reference sequence

    Returns
    -------
    dict
        {dg: [np.array([left_start_ind, left_stop_ind, 
                        right_start_ind, right_stop_ind]
                      ), 
              [left_arm_gene, right_arm_gene],
              coverage,
              number of reads,
              RNAfold basepairing structure,
              RNAfold minimum free energy
             ]}
    list
    
    """
    # Get DG dict
    dg_dict = {}
    dg_valid_structs = []
    with open(dg_info_bed, 'r') as f:
        for line in f:
            data = line.split('\t')
            dg = int(data[3].split('_')[1])
            covg = float(data[3].split('_')[2])
            num_reads = int(data[4])
            left_start = int(data[1]) - 1
            lens = [int(i) for i in data[10].split(',')]
            right_start = [int(i) for i in data[11].split(',')][1] + left_start
            dg_inds = np.array([left_start, left_start + lens[0] - 1, 
                                right_start, right_start + lens[1] - 1])
            dg_ranges = []
            for i in range(2):
                for (gene, gene_inds) in ref_dict.items():
                    gene_range = range(gene_inds[0], gene_inds[1] + 1)
                    if (dg_inds[i*2] in gene_range):
                        dg_ranges.append(gene)
                        break
            folded_dg_inds, [fc, mfe] = sf.calculate_stem_mfe(dg_inds, ref_seq)
            if mfe != 0:
                dg_inds = folded_dg_inds
                dg_valid_structs.append(dg)
            else:
                fc = '.'
                mfe = 0
            dg_dict[dg] = [dg_inds, dg_ranges, covg, num_reads, fc, mfe]
    return dg_dict, dg_valid_structs