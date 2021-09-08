#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

def load_files(infiles):
    for fn in infiles:
        yield pd.read_csv(fn)

df = pd.concat(load_files(snakemake.input.csv))
df_mean = df.groupby(df.index).mean()

df_mean.to_csv(snakemake.output.csv, header = True)
