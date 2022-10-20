#!/usr/bin/env python3
import json
import os
import re
from collections import defaultdict

import bokehutils
import pandas as pd
import tskit
from bokeh.document import Document
from bokeh.embed import file_html
from bokeh.io import export_png
from bokeh.layouts import row
from bokeh.models import CustomJS
from bokeh.models import Div
from bokeh.models import Select
from bokeh.resources import CDN
from snakemake.io import regex

# Main width
docwidth = 1200
doc = Document()
# Add some general information

d = Div(
    text="""

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
box select, wheel zoom, save graphics, reset graphic, help, and hover
tool, each of which can be turned on and off as required. The hover
tool provides the mouse pointer with extra information when hovering
the plot. The included plots are described in the following
sections.</p>



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



""",
    width=docwidth,
)
dmap = Div(
    text="""
<h3>Choropleth world map plot</h3>

<p> The world map plot shows a choropleth overlayed with sampling
sites. The plot is linked to the mean GNN clustering plot such that
selections in one plot translate to selections in the other. </p> """,
    width=docwidth,
)
doc.add_root(d)

# FIXME! One should be able to use different key in plot.
popkey = "sample_node_population"

gnn = defaultdict(dict)
treefiles = defaultdict(dict)

rgx = re.compile(re.sub(r"\$$", "", regex(snakemake.params.fmt)))  # noqa: F821

individuals = None
first = True
has_lng_lat = False
for csvfile, treefile in zip(snakemake.input.csv, snakemake.input.trees):  # noqa: F821
    bn = os.path.basename(csvfile)
    try:
        k = rgx.search(os.path.basename(treefile)).group(0)
    except AttributeError:
        raise
    df = pd.read_csv(csvfile)
    df = df.set_index(
        [
            "sample_node_population",
            "sample_node_population_id",
            "sample_name_unique",
            "sample_node_id",
            "sample_name",
        ]
    ).droplevel(["sample_node_id", "sample_node_population_id"])
    gnn[k] = df
    treefiles[k] = treefile
    # Parse first tree file for metadata
    if first:
        ts = tskit.load(treefile)
        individuals = list(ts.individuals())
        if any(
            ("longitude" in json.loads(ind.metadata.decode()).keys())
            or ("latitude" in json.loads(ind.metadata.decode()).keys())
            for ind in individuals
        ):
            has_lng_lat = True

gnn_all = (
    pd.concat(gnn.values(), keys=[f"gnn-{i + 1}" for i in range(len(gnn.keys()))])
    .groupby(level=[1, 2, 3])
    .mean()
)

gnn_all = gnn_all.sort_values(by=gnn_all.columns.to_list(), ascending=False)

# Add longitude, latitude to gnn_all
if has_lng_lat:
    doc.add_root(dmap)
    gnn_all["longitude"] = None
    gnn_all["latitude"] = None
    for ind in individuals:
        md = json.loads(ind.metadata.decode())
        sample = md["SM"]
        lng = md["longitude"]
        lat = md["latitude"]
        gnn_all.longitude.loc[:, :, sample] = lng
        gnn_all.latitude.loc[:, :, sample] = lat


# Make list of heatmaps
def _heatmap(
    infile,
    gnn,
    title="GNN clustering: {infile}",
    cbar_title="GNN proportion (Z-score)",
    plot_height=500,
    plot_width=700,
    visible=True,
):
    dfg = bokehutils.Matrix(gnn.groupby(popkey).mean())
    dfg.rescale()
    f = bokehutils.figure(
        dfg,
        plot_height=plot_height,
        plot_width=plot_width,
        y_axis_location="right",
        toolbar_location="right",
        title=title.format(infile=infile),
        name=infile,
        visible=visible,
    )

    hm = f.heatmap(row_cluster=True, row_colors=True, cbar_title=cbar_title)
    return hm


def _fst_heatmap(key, data, plot_width=700, plot_height=500, visible=True):
    f = bokehutils.figure(
        data,
        plot_width=plot_width,
        plot_height=plot_height,
        title=f"Mean chromosome fst: {key}",
        visible=visible,
        name=key,
    )
    return f.heatmap(
        cbar_title="Fst",
        row_colors=True,
        row_cluster=True,
        col_cluster=True,
        dendrogram=False,
    )


# GNN prop plot
def _gnnprop(infile, gnn, plot_width=1800, plot_height=400, visible=True):
    f = bokehutils.figure(
        bokehutils.Matrix(gnn),
        plot_width=plot_width,
        plot_height=plot_height,
        title=f"GNN proportion: {infile}",
        y_axis_label="GNN proportion",
        name=infile,
        visible=visible,
    )
    groups = [x for x in sorted(gnn.columns) if x not in ["longitude", "latitude"]]
    p = f.vbar_stack(factor_levels=[0, 1], groups=groups)
    f._fig.title.text_font_size = "18pt"
    return p, f


# Need to parallelize if ts inference
def _fst(key, plot_width=700, plot_height=500, visible=True):
    treefile = treefiles[key]
    print("loading treefile for key ", key)
    ts = tskit.load(treefile)
    tsdata = bokehutils.TSData(ts)
    tsdata.fst()
    return (
        _fst_heatmap(
            key,
            tsdata.data,
            plot_width=plot_width,
            plot_height=plot_height,
            visible=visible,
        ),
        tsdata.data,
    )


hm = {}
gnnprop = {}
fst = {}
fst_data = {}
for k in list(gnn.keys()):
    print("Analysing ", k)
    v = gnn[k]
    hm[k] = _heatmap(k, v)
    gnnprop[k], _ = _gnnprop(k, v)
    fst[k], fst_data[k] = _fst(k)


# Add ALL chromosome results
fst_all = bokehutils.Matrix(
    pd.concat(
        [v.data for v in fst_data.values()],
        keys=[f"fst-{i + 1}" for i in range(len(fst_data.keys()))],
    )
    .groupby(level=1)
    .mean()
)

fig_hm_all = _heatmap("All chromosomes", gnn_all, visible=True)
fig_gnnprop_all, gnnprop_all = _gnnprop("All chromosomes", gnn_all, visible=True)
fig_fst_all = _fst_heatmap("All chromosomes", fst_all, visible=True)
fig_worldmap = gnnprop_all.world_map()

doc.add_root(Div(text="""<h3>Average over all chromosomes</h3><p><br></p>"""))
doc.add_root(row(fig_hm_all, fig_fst_all))
if fig_worldmap is not None:
    doc.add_root(row(fig_worldmap))
doc.add_root(row(fig_gnnprop_all))
doc.add_root(Div(text="""<p><br></p>"""))

##############################
# Save high-res versions of fst, gnnclust, map, and gnnprop
##############################
try:
    export_png(fig_hm_all, filename="gnnprop.png")
    export_png(fig_fst_all, filename="gnnfst.png")
    export_png(fig_worldmap, filename="worldmap.png")
    export_png(fig_gnnprop_all, filename="gnnpropall.png")
except Exception as e:
    print(e)


##############################
# Single chromosomes
#
# NB: Currently selector does not work
#
##############################
doc.add_root(
    Div(
        text=(
            f"""<br></br><hr style="width:{docwidth}px;">"""
            """</hr><br></br><h3>Single chromosomes</h3>"""
        )
    )
)


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
with open(snakemake.output.html, "w") as fh:  # noqa: F821
    fh.write(html)
