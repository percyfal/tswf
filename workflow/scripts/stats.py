#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import scipy.stats
import scipy.cluster.hierarchy

def cluster_gnn_map(X, by="population"):
    '''Cluster GNN values.

    Args:
        X (:class:`~pd.DataFrame`) :
            A DataFrame consisting of gnn proportions. Each row
            corresponds to a chromosome from an individual in a
            population


    '''
    # dfg = X.groupby(by).mean()
    # Zscore normalise
    for col in list(X):
        X[col] = scipy.stats.zscore(X[col])

    row_linkage = scipy.cluster.hierarchy.linkage(X, method="average")
    order = scipy.cluster.hierarchy.leaves_list(row_linkage)
    x_pop = X.index.values[order]

    return X[x_pop], row_linkage
