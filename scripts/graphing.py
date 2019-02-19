# This file is part of CRSSANT:
# Crosslinked RNA Secondary Structure Analysis using Network Techniques
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


def graph_reads(gene_ids, reads_dict, t):
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


def cluster_graph(graph, n=10, t=5):
    """
    Function that performs spectral clustering on the weighted graph.
    Cluster number, k, is determined by finding the first eigengap that is
    some amount t larger than preceding eigengaps.

    Parameters
    ----------
    graph : NetworkX graph
        Weighted graph representation of all reads
    n : int
        Number of eigenvalues (DG splits) to consider threshold
    t : int
        Multiplicity of median eigengap threshold

    Returns
    -------
    reads_dg_dict, dg_ind : dict, int
    """
    dg_ind = 0
    reads_dg_dict = {}
    subgraphs = list(nx.connected_component_subgraphs(graph))
    for subgraph in subgraphs:
        k = 1
        if len(subgraph) > 1:
            L = nx.laplacian_matrix(
                subgraph, nodelist=sorted(subgraph.nodes())
            ).todense()
            D = np.diag(
                [subgraph.degree[node] for node in sorted(subgraph.nodes())]
            )
            w, v = sp.linalg.eigh(L, D, type=1)  # Since L always symmetric
            eigengaps = np.diff(w[:(n + 1)])
            if len(eigengaps) > 2:
                if (w[1] > 1) and (w[1] >= 10*np.median(eigengaps[1:])):
                    k = 2
                else:
                    # ignore divide by 0 warning if eigengaps median is 0
                    np.seterr(divide='ignore', invalid='ignore')
                    eigenratios = np.copy(eigengaps)
                    eigenratios[1:] = np.array([
                        eigengaps[i] / np.median(eigengaps[:i])
                        for i in range(1, len(eigengaps))
                    ])
                    if max(eigenratios) >= t:
                        k = np.argmax(
                            eigenratios >= t
                        ) + 2
            Y = np.transpose(v[:k])
            kmeans = KMeans(n_clusters=k, random_state=0).fit(Y)
            kmeans.labels_ += dg_ind
            subgraph_dict = dict(zip(sorted(subgraph.nodes()), kmeans.labels_))
        else:
            subgraph_dict = {list(subgraph)[0]: dg_ind}       
        reads_dg_dict = {**reads_dg_dict, **subgraph_dict}
        dg_ind += k
    return reads_dg_dict


def get_cliques(graph):
    """
    Function that performs gets cliques from subgraphs of connected components.

    Parameters
    ----------
    graph : NetworkX graph
        Weighted graph representation of all reads

    Returns
    -------
    reads_dg_dict : dict, int
    """
    dg_ind = 0
    reads_dg_dict = {}
    subgraphs = [graph.subgraph(c) for c in nx.connected_components(graph)]
    for subgraph in subgraphs:
        if len(subgraph) > 1:
            cliques_kept = []
            cliques_all = list(nx.find_cliques(subgraph))
            cliques_all.sort(key=len)
            while len(cliques_all) > 0:
                cliques_nodes = set(
                    [node for clique in cliques_kept for node in clique]
                )
                clique_test = cliques_all.pop()
                if not set(list(clique_test)).intersection(cliques_nodes):
                    cliques_kept.append(clique_test)
            dg_inds = [
                [dg_ind + i]*len(clique) 
                for i, clique in enumerate(cliques_kept)
            ]
            subgraph_dict = dict(
                zip(
                    [node for clique in cliques_kept for node in clique],
                    [ind for inds_list in dg_inds for ind in inds_list]
                )
            )
            reads_dg_dict = {**reads_dg_dict, **subgraph_dict}
            dg_ind += len(cliques_kept)
        else:
            read_id = list(subgraph)[0]
            reads_dg_dict[read_id] = dg_ind
            dg_ind += 1
    return reads_dg_dict