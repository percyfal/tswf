#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import os
import sys
import pandas as pd
import tskit
import itertools
import json
import numpy as np
from collections import defaultdict
from bokeh.resources import CDN
from bokeh.embed import file_html, components
from bokeh.models import Paragraph, Div, Panel, Tabs, Dropdown, CustomJS, Select
from bokeh.document import Document
from bokeh.layouts import row, column
from bokeh.io.doc import curdoc
import bokehutils
from common import group_samples_ts
from snakemake.io import regex

# Main width
docwidth=1200
doc = Document()
# Add some general information

d = Div(text="""

<h2>About</h2>

<p>Collection of plots based on results from running <a
href="https://github.com/tskit-dev/tsinfer">tsinfer</a> on phased
bi-allelic snps from samples originating from multiple populations.
tsinfer generates a tree sequence along a reference genome, which by
algorithmic design guarantees that every position in the genome has a
complete genealogy. Every position can therefore be represented as a
genealogical tree, albeit with polytomies.</p>

<p>The plots have been created using <a
href="https://docs.bokeh.org/en/latest/index.html">bokeh</a>. To the
right of each plot, there is a toolbar with tools to pan, box zoom,
wheel zoom, save graphics, reset graphic, help, and hover tool, each
of which can be turned on and off as required. The hover tool provides
the mouse pointer with extra information when hovering the plot. The
included plots are described in the following sections.</p>



<h3>GNN clustering plot</h3>

<p>Z-score normalised GNN proportions. Following <a
href="https://www.nature.com/articles/s41588-019-0483-y#Fig4">Kelleher
et al 2019, Fig 4</a> columns have first been Z-score normalised,
followed by hierarchical clustering on rows. For each genomic
position, genealogical nearest neighbours are assigned by looking at a
focal node (where each node is a haplotype sample). The sister nodes
(i.e. nodes sharing same parent) are the genealogical nearest
neighbours. By labelling the nodes with the corresponding population
of these sister samples, we can calculate the "population" proportion
for the focal node. The rows correspond to the focal population, which
is the average of all individuals in a population, and the columns are
the corresponding GNN proportion assigned to the focal population.
Note therefore that the matrix need not necessarily be symmetric; a
focal node may predominantly be found together with a certain
composition of neighbours, but a neighbour may predominantly be found
in a completly different environment.</p>

<h3>Mean chromosome Fst plot</h3>

<p>Mean chromosome Fst values calculated from tree sequences.<p>

<h3>GNN proportions plot</h3>

<p>Plot showing average GNN proportions for all <b>individuals</b>.
</p>



""", width=1000)
doc.add_root(d)


popkey = "Species"

gnn = defaultdict(dict)
treefiles = defaultdict(dict)

rgx = re.compile(re.sub("\$$", "", regex(snakemake.params.fmt)))

for csvfile, treefile in zip(snakemake.input.csv, snakemake.input.trees):
    bn = os.path.basename(csvfile)
    try:
        k = rgx.search(os.path.basename(treefile)).group(0)
    except AttributeError:
        raise
    df = pd.read_csv(csvfile)
    df["Sample"] = list(map(lambda x: f"{x[0]}/{x[1]}", zip(df.Individual, list(map(int, df.Individual.duplicated())))))
    df = df.set_index(["Species", "Sample", "Sample node", "Individual"]).droplevel(["Sample node", "Individual"])
    gnn[k] = df
    treefiles[k] = treefile


gnn_all = pd.concat(gnn.values(), keys=['gnn-{}'.format(i+1) for i in
                                        range(len(gnn.keys()))]).mean(level=[1, 2])

# Make list of heatmaps
def _heatmap(infile, gnn, title="GNN clustering: {infile}",
             cbar_title="GNN proportion (Z-score)", plot_height=500,
             plot_width=700, visible=True):
    dfg = bokehutils.Matrix(gnn.groupby(popkey).mean())
    dfg.rescale()
    f = bokehutils.figure(dfg, plot_height=plot_height,
                          plot_width=plot_width,
                          y_axis_location="right",
                          toolbar_location="right",
                          title=title.format(infile=infile),
                          name=infile, visible=visible)

    hm = f.heatmap(row_cluster=True, row_colors=True, cbar_title=cbar_title)
    return hm

def _fst_heatmap(key, data, plot_width=700, plot_height=500, visible=True):
    f = bokehutils.figure(data, plot_width=plot_width,
                          plot_height=plot_height, title=f"Mean chromosome fst: {key}",
                          visible=visible, name=key)
    return f.heatmap(cbar_title="Fst", row_colors=True, row_cluster=True, col_cluster=True, dendrogram=False)


# GNN prop plot
def _gnnprop(infile, gnn, plot_width=1800, plot_height=400, visible=True):
    f = bokehutils.figure(bokehutils.Matrix(gnn),
                          plot_width=plot_width,
                          plot_height=plot_height,
                          title=f"GNN proportion: {infile}",
                          y_axis_label="GNN proportion", name=infile, visible=visible)
    return f.vbar_stack()

# Need to parallelize if ts inference
def _fst(key, plot_width=700, plot_height=500, visible=True):
    treefile = treefiles[key]
    print("loading treefile for key ", key)
    ts = tskit.load(treefile)
    tsdata = bokehutils.TSData(ts)
    tsdata.fst()
    return _fst_heatmap(key, tsdata.data, plot_width=plot_width, plot_height=plot_height, visible=visible), tsdata.data

hm = {}
gnnprop = {}
fst = {}
fst_data = {}
for k in list(gnn.keys()):
    print("Analysing ", k)
    v = gnn[k]
    hm[k] = _heatmap(k, v)
    gnnprop[k] = _gnnprop(k, v)
    fst[k], fst_data[k] = _fst(k)


# Add ALL chromosome results
fst_all = bokehutils.Matrix(pd.concat([v.data for v in fst_data.values()],
                    keys=['fst-{}'.format(i+1) for i in range(len(fst_data.keys()))]).mean(level=1))

fig_hm_all = _heatmap("All chromosomes", gnn_all, visible=True)
fig_gnnprop_all = _gnnprop("All chromosomes", gnn_all, visible=True)
fig_fst_all = _fst_heatmap("All chromosomes", fst_all, visible=True)

doc.add_root(Div(text="""<h3>Average over all chromosomes</h3><p><br></p>"""))
doc.add_root(row(fig_hm_all, fig_fst_all))
doc.add_root(row(fig_gnnprop_all))
doc.add_root(Div(text="""<p><br></p>"""))

##############################
# Single chromosomes
#
# NB: Currently selector does not work
#
##############################
doc.add_root(Div(text="""<h3>Single chromosomes</h3>"""))

def create_selector():
    menu = list(hm.keys())
    selector = Select(title="Select chromosome", value=menu[0], options=menu, width=400)
    figures = []
    visible = True
    for k in list(hm.keys()):
        figs = [hm[k], fst[k], gnnprop[k]]
        if visible:
            gnnprop[k].visible = visible
        figs = [gnnprop[k]]
        figures.extend(figs)
        if visible:
            visible = False

    callback = CustomJS(
        args=dict(figs=figs),
        code="""
        let selected = cb_obj.value;
        for(let fig of figs){
        fig.visible = fig.name == selected;
        }
        """,
    )
    selector.js_on_change("value", callback)

    return [selector] + figures

# selector = column(create_selector())
# doc.add_root(selector)

for k in list(hm.keys()):
    doc.add_root(Div(text=f"<h4>{k}</h4>"))
    doc.add_root(row(hm[k], fst[k]))
    doc.add_root(row(gnnprop[k]))

html = file_html(doc, CDN, "Tsinfer EDA")
with open(snakemake.output.html, "w") as fh:
    fh.write(html)
