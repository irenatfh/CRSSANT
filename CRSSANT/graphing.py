# This file is part of CRSSANT:
# Computational RNA Secondary Structure Analysis using Network Techniques
#
###############################################################################
"""
This module is a collection of functions that perform graphing-related tasks
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu


import numpy as np
import networkx as nx
from sklearn.cluster import KMeans
import scipy as sp
import subfunctions as sf


def graph_reads(gene_ids, reads_dict, t=0.3):
    """
    Function that creates a weighted graph representation of the reads

    Reads are represented by nodes. Two reads that have left and right overlap
    ratios > threshold t are connected by an edge of weight overlap / span.

    Parameters
    ----------
    gene_ids : np array
        Read IDs
    reads_dict : dict
        Dictionary of reads
    t : float
        Overlap threshold

    Returns
    -------
    graph : NetworkX graph
    """
    graph = nx.Graph()
    graph.add_nodes_from(gene_ids)
    gene_inds = np.array([reads_dict[read_id][:4] for read_id in gene_ids])
    sorted_ids = []
    for i in range(4):
        sorted_ids.append(
            [read_id for (ind, read_id)
             in sorted(zip(gene_inds[:, i], gene_ids))]
        )
    for (id_1, inds_1) in zip(gene_ids, gene_inds):
        for i in range(4):
            id_list = sorted_ids[i]
            j = id_list.index(id_1)
            loop_flag = 1
            if i%2 == 0:
                check_loop = (j < len(gene_inds) - 1)
            else:
                check_loop = (j > 0)
            while (loop_flag == 1) and check_loop:
                if i%2 == 0:
                    j += 1
                    check_loop = (j < len(gene_inds) - 1)
                else:
                    j -= 1
                    check_loop = (j > 0)
                id_2 = id_list[j]
                inds_2 = reads_dict[id_2]
                if i == 0:
                    check_inds = (inds_2[0] <= inds_1[1])
                elif i == 1:
                    check_inds = (inds_2[1] >= inds_1[0])
                elif i == 2:
                    check_inds = (inds_2[2] <= inds_1[3])
                else:
                    check_inds = (inds_2[3] >= inds_1[2])
                if check_inds:
                    ratio_l, ratio_r = sf.get_overlap_ratios(inds_1, inds_2)
                    if (ratio_l > t) and (ratio_r > t):
                        graph.add_edge(
                            id_1, id_2, weight=(ratio_l + ratio_r)
                        )
                else:
                    loop_flag = 0
    return graph


def cluster_graph(graph, dg_ind):
    """
    Function that performs spectral clustering on the weighted graph

    Parameters
    ----------
    graph : NetworkX graph
        Weighted graph representation of all reads
    dg_ind : int
        DG index to start from

    Returns
    -------
    kmeans_dict, dg_ind : dict, int
    """
    kmeans_dict = {}
    subgraphs = list(nx.connected_component_subgraphs(graph))
    for subgraph in subgraphs:
        if len(subgraph.nodes()) > 1:
            L = nx.laplacian_matrix(subgraph).todense()
            D = np.diag(list(dict(subgraph.degree()).values()))
            w, v = sp.linalg.eigh(L, D)
            eigengaps = np.diff(w)
            k = np.argmax(eigengaps) + 1
            Y = v[:,:k]  # first k eigenvectors for clustering
            kmeans = KMeans(n_clusters=k, random_state=0).fit(Y)
            kmeans.labels_ += dg_ind
            subgraph_dict = dict(zip(subgraph.nodes(), kmeans.labels_))
            kmeans_dict = {**kmeans_dict, **subgraph_dict}
            dg_ind += k
        else:
            read_id = list(subgraph.nodes())[0]
            kmeans_dict[read_id] = dg_ind
            dg_ind +=1 
    return kmeans_dict, dg_ind