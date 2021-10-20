#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import scipy.stats
import scipy.cluster.hierarchy


def cluster_gnn_map(X, by="population"):
    """Cluster GNN values.

    Args:
        X (:class:`~pd.DataFrame`) :
            A DataFrame consisting of gnn proportions. Each row
            corresponds to a chromosome from an individual in a
            population


    """
    # Zscore normalise
    for col in list(X):
        X[col] = scipy.stats.zscore(X[col])

    row_linkage = scipy.cluster.hierarchy.linkage(
        X, method="average", optimal_ordering=True
    )
    order = scipy.cluster.hierarchy.leaves_list(row_linkage)
    x_pop = X.index.values[order]
    X = X.reindex(x_pop)

    return X, row_linkage
