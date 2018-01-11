import numpy as np


################################################################################
def write_dg_ng_sam(reads_file, rna_file, dg_reads_dict, dg_dict):
    with open(reads_file, 'r') as f_read, \
         open(rna_file, 'w') as f_write:
        for line in f_read:
            if line[0] == '@':
                f_write.write(line)
            else:
                data = line.split('\n')[0].split('\t')
                read_id = data[0]
                for (dg, dg_info) in dg_dict.items():
                    ng = dg_info['NG']
                    dg_reads_list = dg_reads_dict[dg]
                    if read_id in dg_reads_list:
                        data.append('DG:i:' + str(dg))
                        data.append('NG:i:' + str(ng))
                        f_write.write('\t'.join(data) + '\n')
                        break
                        
    return


################################################################################
def write_info_bed(bed_file, dg_dict, region):
    with open(bed_file, 'w') as f_write:
        f_write.write('track graphType=arc\n')
        for (dg, dg_info) in dg_dict.items():
            dg_inds = dg_info['arm_indices']
            coverage = dg_info['coverage']
            num_reads = dg_info['num_reads']
            left_start = dg_inds[0] + 1  # biology is 1-indexed
            right_start = dg_inds[2] + 1
            right_stop = dg_inds[3] + 1
            left_len = dg_inds[1] - dg_inds[0]
            right_len = dg_inds[3] - dg_inds[2]
            dg_str = 'Group_%d_%.16f' %(dg, coverage)
            line = [region, str(left_start), str(right_stop), dg_str, 
                    str(num_reads), '-', str(left_start), str(left_start), 
                    '0,0,0', '2', '%d,%d' %(left_len, right_len), 
                    '0,%d' %(right_start - left_start)]
            f_write.write('\t'.join(line) + '\n')
        
    return


################################################################################
def write_helix_bed(bed_file, dg_dict, region, rna):
    with open(bed_file, 'w') as f_write:
        f_write.write('track graphType=arc\n')
        for (dg, dg_info) in dg_dict.items():
            helix_inds = dg_info['basepairs']
            len_helix = np.shape(helix_inds)[1]
            for [left_ind, right_ind] in helix_inds.T:
                line = [region, str(left_ind), str(right_ind), rna, '1', '+',
                        str(left_ind), str(left_ind), '0,0,0']
                f_write.write('\t'.join(line) + '\n')
             
    return