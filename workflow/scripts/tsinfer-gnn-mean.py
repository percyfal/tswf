#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

index_col = 0
if snakemake.wildcards.type == "gnn":
    index_col = [0, 1, 2]

def load_files(infiles):
    for fn in infiles:
        yield pd.read_csv(fn, index_col=index_col)

df = pd.concat(load_files(snakemake.input.csv))
df_mean = df.groupby(df.index).mean()

if snakemake.wildcards.type == "gnn":
    df_mean.set_index(pd.MultiIndex.from_tuples(
        df_mean.index, names=df.index.names), inplace=True)

df_mean.to_csv(snakemake.output.csv, header=True)
