# This file is part of CRSSANT:
# Crosslinked RNA Secondary Structure Analysis using Network Techniques
#
###############################################################################
"""
This module is a collection of functions that produce outputs that are useful 
for understanding aspects of the pipeline that are typically returned.
"""

# Author: Irena Fischer-Hwang
# Contact: ihwang@stanford.edu

import numpy as np
import subfunctions as sf


def get_subgraphs(graph, dg_ind):
    """
    Description

    Parameters
    ----------
    p1 : type
        Description

    Returns
    -------
    r1, r2, r3
    """
    subgraph_eig_dict = {}
    reads_subgraph_dict = {}
    subgraphs = [graph.subgraph(c) for c in nx.connected_components(graph)]
    for subgraph in subgraphs:
        if len(subgraph) > 1:
            L = nx.laplacian_matrix(subgraph).todense()
            D = np.diag([subgraph.degree[node] for node in subgraph])
            w, v = sp.linalg.eigh(L, D)
            subgraph_eig_dict[dg_ind] = w
            subgraph_dict = dict(
                zip(subgraph, [dg_ind]*len(subgraph))
            )
            reads_subgraph_dict = {**reads_subgraph_dict, **subgraph_dict}
            dg_ind += 1
        else:
            read_id = list(subgraph)[0]
            reads_subgraph_dict[read_id] = dg_ind
            dg_ind +=1 
    return reads_subgraph_dict, dg_ind, subgraph_eig_dict