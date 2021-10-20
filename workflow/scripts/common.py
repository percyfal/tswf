#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import sys
import json
import collections
import pandas as pd
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
import tskit
import matplotlib

matplotlib.use("Agg")  # NOQA
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition, inset_axes
import seaborn as sns
import scipy

# Could potentially be parallelized if we want all individuals to be focal
def local_gnn(ts, focal, reference_sets):
    """Calculate the local gnn for a chromosome"""
    reference_set_map = np.zeros(ts.num_nodes, dtype=int) - 1
    for k, reference_set in enumerate(reference_sets):
        for u in reference_set:
            if reference_set_map[u] != -1:
                raise ValueError("Duplicate value in reference sets")
            reference_set_map[u] = k

    K = len(reference_sets)
    A = np.zeros((len(focal), ts.num_trees, K))
    lefts = np.zeros(ts.num_trees, dtype=float)
    rights = np.zeros(ts.num_trees, dtype=float)
    parent = np.zeros(ts.num_nodes, dtype=int) - 1
    sample_count = np.zeros((ts.num_nodes, K), dtype=int)

    # Set the intitial conditions.
    for j in range(K):
        sample_count[reference_sets[j], j] = 1

    for t, ((left, right), edges_out, edges_in) in tqdm(enumerate(ts.edge_diffs())):
        for edge in edges_out:
            parent[edge.child] = -1
            v = edge.parent
            while v != -1:
                sample_count[v] -= sample_count[edge.child]
                v = parent[v]
        for edge in edges_in:
            parent[edge.child] = edge.parent
            v = edge.parent
            while v != -1:
                sample_count[v] += sample_count[edge.child]
                v = parent[v]

        # Process this tree.
        for j, u in enumerate(focal):
            focal_reference_set = reference_set_map[u]
            p = parent[u]
            lefts[t] = left
            rights[t] = right
            while p != tskit.NULL:
                total = np.sum(sample_count[p])
                if total > 1:
                    break
                p = parent[p]
            if p != tskit.NULL:
                scale = 1 / (total - int(focal_reference_set != -1))
                for k, reference_set in enumerate(reference_sets):
                    n = sample_count[p, k] - int(focal_reference_set == k)
                    A[j, t, k] = n * scale
    return (A, lefts, rights)


def group_samples_ts(ts, by="population"):
    """Group samples by tree sequence metadata 'population' property

    NB: we could want to do grouping by some other property.
    Preferably reserve a key for groups and if none is given on tree
    sequence generation, use species (population) as proxy

    See treeseq-inference/src/analyse_human_data.py

    """
    sample_group_set_map = collections.defaultdict(list)
    for population in ts.populations():
        md = json.loads(population.metadata.decode())
        key = md[by]
        sample_group_set_map[key].extend(list(ts.samples(population=population.id)))
    groups = list(sample_group_set_map.keys())
    sample_group_sets = [sample_group_set_map[k] for k in groups]
    return groups, sample_group_sets


def get_samples(sample_data, population=None):
    if population is None:
        pop_index = np.array(range(len(sample_data.individuals_population)))
    else:
        pop_index = np.where(
            np.array(sample_data.individuals_population) == population
        )[0]
    individuals = np.array(list(sample_data.individuals()))[pop_index]
    samples = list()
    for ind in individuals:
        samples.extend(ind.samples)
    return samples


def group_samples(sample_data, by="population"):
    sample_group_set_map = collections.defaultdict(list)
    for pop in sample_data.populations():
        key = pop.metadata[by]
        sample_group_set_map[key].extend(get_samples(sample_data, pop.id))
    groups = list(sample_group_set_map.keys())
    sample_group_sets = [sample_group_set_map[k] for k in groups]
    return groups, sample_group_sets


def chromosome_gnn(ts, sample, groups, sample_group_sets):
    dflist = []
    for ind in ts.individuals():
        md = json.loads(ind.metadata.decode())
        if md["name"] != sample:
            continue

        for j, node in enumerate(ind.nodes):
            A, left, right = local_gnn(ts, [node], sample_group_sets)
            tmp = pd.DataFrame(data=A[0], columns=groups)
            tmp["left"] = left
            tmp["right"] = right
            tmp["haplotype"] = j
            # Remove rows with no difference in GNN to next row
            keep_rows = ~(tmp.iloc[:, 0:5].diff(axis=0) == 0).all(axis=1)
            tmp = tmp[keep_rows]
            dflist.append(tmp)
    df = pd.concat(dflist)
    return df


def make_colour_map(groups):
    colours = {}
    # cmap = sns.color_palette("coolwarm", n_colors=15)
    cmap = sns.diverging_palette(150, 275, s=80, l=55, n=len(groups))
    for i, g in enumerate(groups):
        colours[g] = cmap[i]
    return colours


def plot_sample_gnn(df_all, groups, outfile, sample, title=None):
    gs = matplotlib.gridspec.GridSpec(
        2, 2, height_ratios=[1, 1], hspace=0.6, width_ratios=[5, 1]
    )
    fig = plt.figure(figsize=(17, 4))
    _, ax_main = plt.subplots()
    ax_upper = plt.subplot(gs[0, 0])
    ax_lower = plt.subplot(gs[1, 0])
    colours = make_colour_map(groups)

    for j, ax in enumerate([ax_upper, ax_lower]):
        df = df_all[df_all.haplotype == j]
        left = df.left
        width = np.diff(np.append(left, df.right.max()))
        total = np.zeros_like(width)
        for group in tqdm(groups):
            # Stack bar at position left, height = df[group],
            # width = width, bottom = cumulative total
            ax.bar(
                left,
                df[group],
                bottom=total,
                width=width,
                align="edge",
                label=group,
                color=colours[group],
            )
            total += df[group]
        ax.set_title("{} haplotype ({})".format(sample, j + 1))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(0, df.right.max())
        ax.set_ylim(0, 1)
        ax.axis("off")

    # Add legende
    ax_right = plt.subplot(gs[0, 1])
    ax_right.legend()
    # Save figure
    print("Saving figure")
    plt.savefig(outfile, bbox_inches="tight", dpi=400)


def _mean_cluster(df, zscore=True):
    if zscore:
        for col in list(df):
            df[col] = scipy.stats.zscore(df[col])

    row_linkage = scipy.cluster.hierarchy.linkage(
        df, method="average", optimal_ordering=True
    )

    order = scipy.cluster.hierarchy.leaves_list(row_linkage)
    x_pop = df.index.values[order]
    return row_linkage, df[x_pop]


def plot_mean_cluster(df, outfile=None, title=None, zscore=True):
    row_linkage, df = _mean_cluster(df, zscore)
    cg = sns.clustermap(df, row_linkage=row_linkage, col_cluster=False, rasterized=True)
    cg.ax_heatmap.set_ylabel("")
    cg.fig.suptitle(title)
    for tick in cg.ax_heatmap.get_xticklabels():
        tick.set_rotation(-45)
        tick.set_ha("left")
        tick.set_rotation_mode("anchor")
    if outfile is not None:
        plt.savefig(outfile, bbox_inches="tight", dpi=400)


def plot_gnn(df, outfile, groups, sample_group_sets, title=None):
    # Track gnn proportion for each individual
    A = np.zeros((len(sample_group_sets), len(df)))
    colours = make_colour_map(groups)

    for j, group in enumerate(groups):
        # For now we assume group is a population; in tsinfer paper
        # the groupings consist of populations; we want the
        # corresponding population submatrices of df
        A[j, :] = np.sum([df[pop].values for pop in [group]], axis=0)

    index = np.argsort(A[0])[::-1]
    inds = df.Individual.values

    inds = list(inds[index])
    A = A[:, index]

    margin = 0.05
    width = (1.0 - 2.0 * margin) / len(groups)

    fig, ax = plt.subplots(1, 1, figsize=(17, 3))
    # Insert extra x coordinate between groups; also order according to group
    i = 0
    x = []
    xticks = []
    for name, group in df.groupby(level=0):
        x.extend(np.arange(i, i + len(group)))
        xticks.append(float(i + int(len(group) / 2)))
        i = i + len(group) + 1

    for j, group in enumerate(groups):
        ax.bar(
            x,
            A[j],
            bottom=np.sum(A[:j, :], axis=0),
            label=group,
            width=1,
            color=colours[group],
            align="edge",
        )
    ax.set_xlim(0, max(x))
    ax.set_ylim(0, 1)
    ax.set_xticks(np.array(xticks))
    ax.set_xticklabels(groups, rotation="vertical")
    ax.set_yticks([])
    ax.set_title(title)
    # ax.legend(groups)
    # ax.axis('off')
    plt.savefig(outfile, bbox_inches="tight", dpi=400)
