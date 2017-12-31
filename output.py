import numpy as np


################################################################################
def rna_dg_sam(reads_file, rna_file, kmeans_dict):
    with open(reads_file, 'r') as f_read, \
         open(rna_file, 'w') as f_write:
        for line in f_read:
            if line[0] == '@':
                f_write.write(line)
            else:
                data = line.split('\n')[0].split('\t')
                read_id = data[0]
                if read_id in kmeans_dict.keys():
                    dg = kmeans_dict[read_id]
                    data.append('DG:i:' + str(dg))
                    f_write.write('\t'.join(data) + '\n')
    return