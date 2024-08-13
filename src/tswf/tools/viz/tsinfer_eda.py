"""tsinfer-eda - exploratory data analysis

Given gnn CSV files and tree sequence files TS generated for the same
chromosomes, summarize and cluster GNN components, calculate Fst, and
optionally display samples on a map and output results to interactive
html file.

Note that CSV and TS files are paired according the order in which
they are passed to the options, so make sure they follow the order
they were generated.

"""

import logging
import os
from collections import defaultdict

import click
import pandas as pd
import tskit
from bokeh.document import Document
from bokeh.embed import file_html
from bokeh.io import export_png
from bokeh.layouts import row
from bokeh.models import Div  # type: ignore
from bokeh.resources import CDN

import tswf.viz.bokehutils as bokehutils
from tswf.cli import pass_environment


def make_doc(docwidth):
    # Main width
    doc = Document(title="Tsinfer exploratory data analysis")
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

    return doc, dmap


# Make list of heatmaps
def _heatmap(
    infile,
    gnn,
    population_key="sample_node_population",
    title="GNN clustering: {infile}",
    cbar_title="GNN proportion (Z-score)",
    height=500,
    width=700,
    visible=True,
):
    dfg = bokehutils.Matrix(gnn.groupby(population_key).mean())
    dfg.rescale()
    f = bokehutils.figure(
        dfg,
        height=height,
        width=width,
        y_axis_location="right",
        toolbar_location="right",
        title=title.format(infile=infile),
        name=infile,
        visible=visible,
    )

    hm = f.heatmap(row_cluster=True, row_colors=True, cbar_title=cbar_title)
    return hm


def _fst_heatmap(key, data, population_key, width=700, height=500, visible=True):
    f = bokehutils.figure(
        data,
        width=width,
        height=height,
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
def _gnnprop(infile, gnn, width=1800, height=400, visible=True, metadata=None):
    f = bokehutils.figure(
        bokehutils.Matrix(gnn),
        width=width,
        height=height,
        title=f"GNN proportion: {infile}",
        y_axis_label="GNN proportion",
        name=infile,
        visible=visible,
        metadata=metadata,
    )
    groups = [x for x in sorted(gnn.columns)]
    p = f.vbar_stack(factor_levels=[0, 1], groups=groups)
    f._fig.title.text_font_size = "18pt"
    return p, f


# Need to parallelize if ts inference
def _fst(key, treefile, population_key, width=700, height=500, visible=True):
    ts = tskit.load(treefile)
    tsdata = bokehutils.TSData(ts)
    tsdata.fst()
    return (
        _fst_heatmap(
            key,
            tsdata.data,
            population_key,
            width=width,
            height=height,
            visible=visible,
        ),
        tsdata.data,
    )


def _ts_individuals(ts, sample_id):
    """Get tree sequence metadata for individuals"""
    individuals = []
    for ind in list(ts.individuals()):
        md = ind.metadata
        individuals.append(md)
    return pd.DataFrame(individuals).set_index(sample_id)


@click.command(help=__doc__)
@click.option(
    "--gnn",
    type=click.Path(exists=True),
    multiple=True,
    help="GNN csv file for a chromosome",
)
@click.option(
    "--ts",
    type=click.Path(exists=True),
    multiple=True,
    help="tree sequence (TS) file for a chromosome",
)
@click.option(
    "--gnn-ts",
    type=(click.Path(exists=True), click.Path(exists=True)),
    multiple=True,
    help="<GNN csv, TS file> pairs for a chromosome",
)
@click.option("--output-file", type=click.Path(), help="output (html) file name")
@click.option(
    "--population-key",
    type=str,
    default="sample_node_population",
    help="tree sequence metadata key that defines population name",
)
@click.option(
    "--sample-id",
    type=str,
    default="variant_data_sample_id",
    help="Sample ID variable name in individual metadata",
)
@pass_environment
def main(env, gnn, ts, gnn_ts, output_file, population_key, sample_id):  # noqa: C901
    docwidth = 1200

    assert len(gnn) == len(ts), "must supply same number of GNN csv and TS files"
    assert (len(gnn_ts) == 0 and len(gnn) > 0) or (len(gnn_ts) > 0 and len(gnn) == 0), (
        "either supply GNN, TS file pairs with option "
        "--gnn-ts or separately via --gnn and --ts"
    )

    first = True
    has_lng_lat = False
    if len(gnn) > 0:
        pairs = zip(gnn, ts)
    else:
        pairs = gnn_ts

    gnn = defaultdict(dict)
    treefiles = defaultdict(dict)
    individuals = pd.DataFrame()
    for csvfile, treefile in pairs:
        # Chromosome is the directory of the file name
        k = os.path.basename(os.path.dirname(treefile))
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
        if first:
            individuals = _ts_individuals(tskit.load(treefile), sample_id)
            first = False

    if "longitude" in individuals.columns and "latitude" in individuals.columns:
        has_lng_lat = True

    gnn_all = (
        pd.concat(gnn.values(), keys=[f"gnn-{i + 1}" for i in range(len(gnn.keys()))])
        .groupby(level=[1, 2, 3])
        .mean()
    )

    gnn_all = gnn_all.sort_values(by=gnn_all.columns.to_list(), ascending=False)

    doc, dmap = make_doc(docwidth)
    if has_lng_lat:
        doc.add_root(dmap)

    hm = {}
    gnnprop = {}
    fst = {}
    fst_data = {}
    for k in list(gnn.keys()):
        v = gnn[k]
        hm[k] = _heatmap(k, v)
        gnnprop[k], _ = _gnnprop(k, v)
        fst[k], fst_data[k] = _fst(k, treefiles[k], population_key)

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
    fig_gnnprop_all, gnnprop_all = _gnnprop(
        "All chromosomes", gnn_all, visible=True, metadata=individuals
    )
    fig_fst_all = _fst_heatmap("All chromosomes", fst_all, population_key, visible=True)
    fig_worldmap = gnnprop_all.world_map()

    doc.add_root(Div(text="""<h3>Average over all chromosomes</h3><p><br></p>"""))
    doc.add_root(row(fig_hm_all, fig_fst_all))
    if fig_worldmap is not None:
        doc.add_root(row(fig_worldmap))
    doc.add_root(row(fig_gnnprop_all))
    doc.add_root(Div(text="""<p><br></p>"""))

    # Single chromosomes
    doc.add_root(
        Div(
            text=(
                f"""<br></br><hr style="width:{docwidth}px;">"""  # noqa: E231,E702
                """</hr><br></br><h3>Single chromosomes</h3>"""
            )
        )
    )

    for k in list(hm.keys()):
        doc.add_root(Div(text=f"<h4>{k}</h4>"))
        doc.add_root(row(hm[k], fst[k]))
        doc.add_root(row(gnnprop[k]))

    html = file_html(doc, CDN, "Tsinfer EDA")
    if output_file is not None:
        with open(output_file, "w") as fh:  # noqa: F821
            fh.write(html)
    else:
        print(html)

    from bokeh.plotting import figure

    p = figure(width=400, height=400)

    # add a scatter circle renderer with a size, color, and alpha
    p.scatter([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], size=20, color="navy", alpha=0.5)

    ##############################
    # Save high-res versions of fst, gnnclust, map, and gnnprop
    ##############################

    try:
        for fig, fn in zip(
            (fig_hm_all, fig_fst_all, fig_worldmap, fig_gnnprop_all),
            ("gnnprop.png", "gnnfst.png", "worldmap.png", "gnnpropall.png"),
        ):
            export_png(fig, filename=fn)
    except RuntimeError as e:
        logging.error(e)
        logging.error("Failed to generate png %s", fn)
