import numpy as np
import networkx as nx
from sklearn.cluster import KMeans
import scipy as sp
import subfunctions as sf


################################################################################
def graph_reads(rna_ids, rna_info, reads_dict, t=0.25):
    graph = nx.Graph()
    graph.add_nodes_from(rna_ids)
    sorted_l_start = np.array([read_id for (read_loc, read_id) 
                               in sorted(zip(rna_info[:, 0], rna_ids))])
    sorted_l_stop = np.array([read_id for (read_loc, read_id) 
                              in sorted(zip(rna_info[:, 1], rna_ids))])
    sorted_r_start = np.array([read_id for (read_loc, read_id) 
                               in sorted(zip(rna_info[:, 2], rna_ids))])
    sorted_r_stop = np.array([read_id for (read_loc, read_id) 
                              in sorted(zip(rna_info[:, 3], rna_ids))])
    sorted_reads = np.array([sorted_l_start, sorted_l_stop, 
                             sorted_r_start, sorted_r_stop])
    for read_1_id in rna_ids:
        read_1 = reads_dict[read_1_id]
        for i in range(4):
            read_list = sorted_reads[i]
            j = np.where(read_list == read_1_id)[0][0]
            loop_flag = 1
            if i%2 == 0:
                index_condition = (j < len(rna_info) - 1)
            else:
                index_condition = (j > 0)
            while (loop_flag == 1) and index_condition:
                if i%2 == 0:
                    j += 1
                    index_condition = (j < len(rna_info) - 1)
                else:
                    j -= 1
                    index_condition = (j > 0)
                read_2_id = read_list[j]
                read_2 = reads_dict[read_2_id]
                if i == 0:
                    read_condition = (read_2[0] <= read_1[1])
                elif i == 1:
                    read_condition = (read_2[1] >= read_1[0])
                elif i == 2:
                    read_condition = (read_2[2] <= read_1[3])
                else:
                    read_condition = (read_2[3] >= read_1[2])
                if read_condition:
                    overlap_l, overlap_r, overlap_g, \
                    span_l, span_r, span_g = sf.get_overlaps(read_1, read_2)
                    if (overlap_l/span_l > t) and (overlap_r/span_r > t):
                        graph.add_edge(read_1_id, read_2_id,
                                   weight=(overlap_l/span_l + \
                                           overlap_r/span_r))
                else:
                    loop_flag = 0
    return graph


################################################################################
def spectral_clustering(graph, dg_index):
    kmeans_dict = {}
    subgraphs = list(nx.connected_component_subgraphs(graph))
    for subgraph in subgraphs:
        if len(subgraph.nodes()) > 1:
            L = nx.laplacian_matrix(subgraph).todense()
            D = np.diag(list(dict(subgraph.degree()).values()))
            w, v = sp.linalg.eigh(L, D)
            eigengaps = np.diff(w)
            k = np.argmax(eigengaps) + 1
            Y = v[:,:k]
            kmeans = KMeans(n_clusters=k, random_state=0).fit(Y)
            kmeans.labels_ += dg_index
            subgraph_dict = dict(zip(
                set(kmeans.labels_),
                [[] for i in range(len(set(kmeans.labels_)))]))
            for i in range(len(subgraph.nodes())):
                subgraph_dict[kmeans.labels_[i]].append(
                    list(subgraph.nodes())[i])
            kmeans_dict = {**kmeans_dict, 
                           **subgraph_dict}
            dg_index += k
        else:
            read_id = list(subgraph.nodes())[0]
            kmeans_dict[dg_index] = read_id
            dg_index +=1
    return kmeans_dict, dg_index