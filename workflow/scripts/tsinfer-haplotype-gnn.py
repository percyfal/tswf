#!/usr/bin/env python3
import os
import sys
import collections
import tskit
import json
import numpy as np
import pandas as pd
import argparse

group_id = "population"
sample_id = "SM"

try:
    tree_sequence = snakemake.input.tree_sequence
    sample = snakemake.wildcards.sample
    window_size = snakemake.wildcards.window_size
    group_id = snakemake.params.group_id
    sample_id = snakemake.params.sample_id
    populations = snakemake.params.populations
    output_file = snakemake.output.txt
    plot = snakemake.params.plot
    start_pos = snakemake.wildcards.start_pos
    end_pos = snakemake.wildcards.end_pos
except NameError as e:
    tree_sequence = None
    sample = None
    window_size = None
    output_file = None
    populations = None
    plot = False
    start_pos = 0
    end_pos = None

parser = argparse.ArgumentParser(
    description="""
    Retrieve haplotype gnn data for a sample.
    """
)
parser.add_argument(
    "tree_sequence", help="The tree sequence to load in .trees format"
)
parser.add_argument(
    "sample", help="The sample for which to calculate data"
)
parser.add_argument(
    "--window-size", help="The window size in bp", default=None, type=int
)
parser.add_argument(
    "--group-id", help="The metadata group identifier to group samples by", default=group_id, type=str
)
parser.add_argument(
    "--sample-id", help="The metadata sample identifier", default=sample_id, type=str
)
parser.add_argument(
    "--populations", help="Populations to look at", default=None, nargs="*"
)
parser.add_argument(
    "--output-file", help="output file name", default=sys.stdout
)
parser.add_argument(
    "--plot", help="output plot file name", default=None, action='store'
)
parser.add_argument(
    "--start-pos", help="windows start position", default=start_pos, type=int
)
parser.add_argument(
    "--end-pos", help="windows end position", default=None, type=int
)
parser.add_argument(
    "--haplotype", help="windows end position", default=0, type=int
)
parser.add_argument(
    "--title", help="plot title", default=None, type=str
)

args = parser.parse_args()

# Copy from https://github.com/tskit-dev/tskit/pull/683/files#diff-e5e589330499b325320b2e3c205eaf350660b50691d3e1655f8789683e49dca6R399
def parse_time_windows(ts, time_windows):
    if time_windows is None:
        time_windows = [0.0, ts.max_root_time]
    return np.array(time_windows)

def windowed_genealogical_nearest_neighbours(
        ts,
        focal,
        reference_sets,
        windows=None,
        time_windows=None,
        span_normalise=True,
        time_normalise=True
):
    reference_set_map = np.full(ts.num_nodes, tskit.NULL, dtype=int)
    for k, reference_set in enumerate(reference_sets):
        for u in reference_set:
            if reference_set_map[u] != tskit.NULL:
                raise ValueError("Duplicate value in reference sets")
            reference_set_map[u] = k

    windows_used = windows is not None
    time_windows_used = time_windows is not None
    windows = ts.parse_windows(windows)
    num_windows = windows.shape[0] - 1
    time_windows = parse_time_windows(ts, time_windows)
    num_time_windows = time_windows.shape[0] - 1
    A = np.zeros((num_windows, num_time_windows, len(focal), len(reference_sets)))
    K = len(reference_sets)
    parent = np.full(ts.num_nodes, tskit.NULL, dtype=int)
    sample_count = np.zeros((ts.num_nodes, K), dtype=int)
    time = ts.tables.nodes.time
    norm = np.zeros((num_windows, num_time_windows, len(focal)))

    # Set the initial conditions.
    for j in range(K):
        sample_count[reference_sets[j], j] = 1

    window_index = 0
    ## Loop the tree sequence
    for (t_left, t_right), edges_out, edges_in in ts.edge_diffs():
        for edge in edges_out:
            parent[edge.child] = tskit.NULL
            v = edge.parent
            while v != tskit.NULL:
                sample_count[v] -= sample_count[edge.child]
                v = parent[v]
        for edge in edges_in:
            parent[edge.child] = edge.parent
            v = edge.parent
            while v != tskit.NULL:
                sample_count[v] += sample_count[edge.child]
                v = parent[v]

        # Update the windows
        assert window_index < num_windows
        while windows[window_index] < t_right and window_index + 1 <= num_windows:
            w_left = windows[window_index]
            w_right = windows[window_index + 1]
            left = max(t_left, w_left)
            right = min(t_right, w_right)
            span = right - left
            # Process this tree.
            for j, u in enumerate(focal):
                focal_reference_set = reference_set_map[u]
                delta = int(focal_reference_set != tskit.NULL)
                p = u
                while p != tskit.NULL:
                    total = np.sum(sample_count[p])
                    if total > delta:
                        break
                    p = parent[p]
                if p != tskit.NULL:
                    scale = span / (total - delta)
                    time_index = np.searchsorted(time_windows, time[p]) - 1
                    if 0 <= time_index < num_time_windows:
                        for k in range(len(reference_sets)):
                            n = sample_count[p, k] - int(focal_reference_set == k)
                            A[window_index, time_index, j, k] += n * scale
                        norm[window_index, time_index, j] += span
            assert span > 0
            if w_right <= t_right:
                window_index += 1
            else:
                # This interval crosses a tree boundary, so we update it again
                # in the next tree
                break

    # Reshape norm depending on normalization selected
    # Return NaN when normalisation value is 0
    if span_normalise and time_normalise:
        reshaped_norm = norm.reshape((num_windows, num_time_windows, len(focal), 1))
    elif span_normalise and not time_normalise:
        norm = np.sum(norm, axis=1)
        reshaped_norm = norm.reshape((num_windows, 1, len(focal), 1))
    elif time_normalise and not span_normalise:
        norm = np.sum(norm, axis=0)
        reshaped_norm = norm.reshape((1, num_time_windows, len(focal), 1))

    with np.errstate(invalid="ignore", divide="ignore"):
        A /= reshaped_norm
    A[np.all(A == 0, axis=3)] = np.nan

    # Remove dimension for windows and/or time_windows if parameter is None
    if not windows_used and time_windows_used:
        A = A.reshape((num_time_windows, len(focal), len(reference_sets)))
    elif not time_windows_used and windows_used:
        A = A.reshape((num_windows, len(focal), len(reference_sets)))
    elif not windows_used and not time_windows_used:
        A = A.reshape((len(focal), len(reference_sets)))
    return A

## Follow treeseq-inference/src/analyse_human_data.py::process_hg01933_local_gnn
ts = tskit.load(args.tree_sequence)

## FIXME: region_sample_set_map - partition samples into groups. In
## reality we don't know on what key to group; this should be an
## option
if args.populations is None:
    populations = [json.loads(pop.metadata.decode())[group_id] for pop in ts.populations()]
else:
    populations = args.populations
group_sample_set_map = collections.defaultdict(list)
for pop in ts.populations():
    md = json.loads(pop.metadata.decode())
    group = md[group_id]
    if group not in populations:
        continue
    group_sample_set_map[group].extend(list(ts.samples(
        population=pop.id)))
groups = list(group_sample_set_map.keys())
group_sample_sets = [group_sample_set_map[k] for k in groups]

## Define windows
windows = None
if args.window_size is not None:
    if args.end_pos is None:
        n = int((ts.sequence_length - args.start_pos) / args.window_size) + 1
        windows = np.linspace(start=args.start_pos, stop=ts.sequence_length, num=n)
    else:
        n = int((args.end_pos - args.start_pos) / args.window_size) + 1
        windows = np.linspace(start=args.start_pos, stop=args.end_pos, num=n)
    # Currently not possible to start within haplotype block
    n = int((ts.sequence_length - args.start_pos) / args.window_size) + 1
    windows = np.linspace(start=0.0, stop=ts.sequence_length, num=n)

for ind in ts.individuals():
    md = json.loads(ind.metadata.decode())
    if md[sample_id] == args.sample:
        # Focal node is haplotype node, *not* a particular variant
        # site; the sites are traveresed with ts.edge_diffs
        A = windowed_genealogical_nearest_neighbours(ts, ind.nodes, group_sample_sets, windows=windows)

dflist = []
if windows is None:
    for i in range(A.shape[0]):
        x = pd.DataFrame(A[i, :])
        x = x.T
        x.columns = populations
        x["haplotype"] = i
        x["start"] = 0
        x["end"] = ts.sequence_length
        dflist.append(x)
else:
    for i in range(A.shape[1]):
        x = pd.DataFrame(A[:, i, :])
        x.columns = populations
        x["haplotype"] = i
        x["start"] = windows[0:-1]
        x["end"] = windows[1:]
        dflist.append(x)

df = pd.concat(dflist)
df = df[df.haplotype == args.haplotype]
df.set_index(["haplotype", "start", "end"], inplace=True)
df.to_csv(args.output_file)

# Plot haplotype blocks
if args.plot is not None:
    import bokehutils
    from bokeh.plotting import show, figure, output_file
    from bokeh.models import ColumnDataSource, PrintfTickFormatter
    source = ColumnDataSource(df.reset_index())
    p = figure(title=args.title, plot_width=1800, plot_height=400,
               min_border=0, y_range=(0, 1), x_range=(0, ts.sequence_length))
    if args.window_size is None:
        width = ts.sequence_length
    else:
        width = args.window_size
    p.vbar_stack(populations, x='start', width=width,
                 color=bokehutils._get_palette(n=len(populations)), source=source,
                 legend_label=populations)
    p.add_layout(p.legend[0], "right")
    p.legend[0].label_text_font_size = "14pt"
    p.axis.major_tick_line_color = "black"
    p.axis.minor_tick_line_color = "black"
    p.xaxis.axis_label = "Base pairs"
    p.xaxis.axis_label_text_font_size = "18pt"
    p.xaxis.major_label_orientation = 1.0
    p.xaxis.major_label_text_font_size = "16pt"
    p.xaxis[0].formatter = PrintfTickFormatter(format="%4.1e")
    p.yaxis.major_label_text_font_size = "16pt"
    p.yaxis.axis_label = "GNN proportion"
    p.yaxis.axis_label_text_font_size = "18pt"
    p.axis.axis_line_color = "black"
    p.grid.grid_line_color = None
    p.outline_line_color = "black"
    p.title.text_font_size = "20pt"

    _, ext = os.path.splitext(args.plot)
    if ext == ".png":
        from bokeh.io import export_png
        p.toolbar_location = None
        export_png(p, filename=args.plot)
    elif ext == ".html":
        output_file(args.plot)
        show(p)
