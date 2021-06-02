#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import tskit


from common import plot_sample_gnn, make_colour_map, group_samples_by_populations, plot_mean_cluster, plot_gnn

# Load tree sequence
ts = tskit.load(snakemake.input.trees)

# Load sample group sets
groups, sample_group_sets = group_samples_by_populations(ts)

# Read input data frame
df = pd.read_csv(snakemake.input.csv)


# Plot type
sample = snakemake.wildcards.sample
datatype = snakemake.wildcards.datatype

if sample == "":
    df.set_index("Species", drop=True, inplace=True)
    # mean cluster
    if datatype == "mean":
        plot_mean_cluster(df, snakemake.output.png)
    elif datatype == "gnn":
        plot_gnn(df, snakemake.output.png, groups, sample_group_sets)
else:
    plot_sample_gnn(df, groups, snakemake.output.png, sample)
