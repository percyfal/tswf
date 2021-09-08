#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import scipy
from bokeh.transform import linear_cmap, transform
from bokeh.palettes import Set3, Viridis256
from bokeh.models import LinearColorMapper, ColumnDataSource, HoverTool, ColorBar, BasicTicker, CategoricalScale, CategoricalAxis
from bokeh.models.ranges import FactorRange
from bokeh.layouts import row
from bokeh import plotting
from stats import cluster_gnn_map
    

def _get_palette(cmap=Set3[12], n=12, start=0, end=1):
    import matplotlib, numpy as np
    linspace = np.linspace(start, end, n)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("customcmap", cmap)
    palette = cmap(linspace)
    hex_palette = [matplotlib.colors.rgb2hex(c) for c in palette]
    return hex_palette


class Figure:
    def __init__(self, *arg, **kw):
        self._fig = None
        self._kw = kw

    def _add_dendrogram(self, linkage, xmax=5):
        '''Add dendrogram'''
        results = scipy.cluster.hierarchy.dendrogram(linkage,
                                                     no_plot=True)
        ymax = float(linkage.shape[0]) + 0.5

        ycoord = pd.DataFrame(results['icoord'])
        ycoord = ycoord * (ymax / ycoord.max().max())
        ycoord = ycoord.values

        xcoord = pd.DataFrame(results['dcoord'])
        xcoord = xcoord * (xmax / xcoord.max().max()) - xmax
        xcoord = xcoord.values
        
        for x, y in zip(xcoord, ycoord):
            x = list(map(lambda z: -z, x))
            self._fig.line(x=x, y=y, line_color='black')
            


    def heatmap(self, data, low=None, high=None, palette="Viridis256",
                color_bar=True, dendrogram=True, dcols=5, *args, **kwargs):
        '''Make heatmap of data, possibly with dendrogram

        Args:
            data (:class:`~pd.DataFrame`) :
                A DataFrame n x n matrix with 'z' values to plot
 
            palette (str or seq[color], optional) :
                A palette to use to colormap z

        Any additional keyword arguments are passed to
        :func:`bokeh.plotting.rect`

        '''
        if data.index.name is None:
            data.index.name = "Row"
        rowname = data.index.name
        if data.columns.name is None:
            data.columns.name = "Column"
        colname = data.columns.name

        # Cluster -> recalculate data!
        if dendrogram:
            data, linkage = cluster_gnn_map(data, by=rowname)
        
        # Create source
        df = pd.DataFrame(data.stack(), columns=['z']).reset_index()
        source = ColumnDataSource(df)

        # Setup color mapper
        if low is None:
            low = df.z.min()
        if high is None:
            high = df.z.max()
        color_mapper = LinearColorMapper(palette=palette, low=low, high=high)

        # Setup ranges; need to redefine figure since default range is
        # linear.
        factors = [rowname] + list(data.index)
        if dendrogram:
            factors = list(map(lambda x: f"__{x}", range(dcols))) + factors

        # Setup ranges
        if 'x_range' not in self._kw.keys():
            self._kw['x_range'] = FactorRange(factors=factors, start=-dcols)
        if 'y_range' not in self._kw.keys():
           self._kw['y_range'] = list(reversed(data.index))
        self._fig = plotting.figure(**self._kw)


        # Make heatmap 
        self._fig.rect(x=colname, y=rowname,
                       width=1, height=1,
                       line_color='black', line_alpha=0.2, alpha=0.8,
                       source=source, fill_color=transform('z', color_mapper),
                       *args, **kwargs)

        # Add group colors
        ngroups = len(data.index)
        csource = pd.DataFrame({'x': dcols + 0.7, 'group': data.index, 'color': _get_palette(n=ngroups)})
        self._fig.rect(
            x='x', y='group', fill_color='color',
            line_color=None, line_alpha=0.2, alpha=0.8,
            source=csource, width=.6, height=1)

        # Add dendrogram
        if dendrogram:
            self._add_dendrogram(linkage, xmax=dcols)
            
        # Setup hover tooltips
        hover = HoverTool()
        hover.tooltips = [
            ("Query (row)", f"@{rowname}"),
            ("Subject (column)", f"@{colname}"),
            ("GNN Proportion (Z-scaled)", "@z")
        ]
        self._fig.add_tools(hover)

        # Add colorbar
        if color_bar:
            color_bar = ColorBar(color_mapper=color_mapper, location=(0,0),
                                 ticker=BasicTicker(), border_line_color=None,
                                 title="GNN Proportion")

            color_bar_fig = plotting.figure(height=self._fig.plot_height,
                                            width=int(self._fig.plot_width * 0.2),
                                            toolbar_location=None, min_border=0,
                                            outline_line_color=None)
            self._fig.add_layout(color_bar, "right")

        # Configure axes
        self._fig.axis.major_tick_line_color = None
        self._fig.axis.minor_tick_line_color = None
        self._fig.xaxis.major_label_overrides = dict(zip(factors, map(lambda x: '' if x.startswith("__") else x, factors)))
        self._fig.xaxis.major_label_orientation = 1.0
        self._fig.axis.axis_line_color = None
        self._fig.grid.grid_line_color = None
        self._fig.outline_line_color = None
        return self._fig


def figure(**kwargs):
    return Figure(**kwargs)


