#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import scipy
import collections
import itertools
import json
from bokeh.transform import linear_cmap, transform
from bokeh.palettes import Set3, Viridis256
from bokeh.models import LinearColorMapper, ColumnDataSource, HoverTool, ColorBar, BasicTicker, CategoricalScale, CategoricalAxis
from bokeh.models.ranges import FactorRange
from bokeh.layouts import row, column
from bokeh import plotting
from stats import cluster_gnn_map
    

def _get_palette(cmap=Set3[12], n=12, start=0, end=1):
    import matplotlib, numpy as np
    linspace = np.linspace(start, end, n)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("customcmap", cmap)
    palette = cmap(linspace)
    hex_palette = [matplotlib.colors.rgb2hex(c) for c in palette]
    return hex_palette

class GNNFigure:
    def __init__(self, data, *arg, **kw):
        self._data = data
        self._fig = None
        self._kw = kw
        if self._data.index.name is None:
            self._data.index.name = "Row"
        if self._data.columns.name is None:
            self._data.columns.name = "Column"
        # GNN data structure from tsinfer-gnn; we know we have columns
        # "Sample node", "Individual", "Species"
        self._data.set_index(["Individual", "Species", "Sample node"],
                             inplace=True)

    @property
    def data(self):
        return self._data

    
    def _add_dendrogram(self, linkage, xmax=5):
        '''Add dendrogram'''
        results = scipy.cluster.hierarchy.dendrogram(linkage,
                                                     no_plot=True)
        ymax = float(linkage.shape[0]) + 0.5

        ycoord = pd.DataFrame(results['icoord'])
        ycoord = ycoord * (ymax / ycoord.max().max())
        ycoord = - ycoord.values + ymax + 0.5

        xcoord = pd.DataFrame(results['dcoord'])
        xcoord = xcoord * (xmax / xcoord.max().max()) - xmax
        xcoord = xcoord.values
        
        for x, y in zip(xcoord, ycoord):
            x = list(map(lambda z: -z, x))
            self._fig.line(x=x, y=y, line_color='black')
            


    def heatmap(self, low=None, high=None, palette="Viridis256",
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
        rowname = self.data.index.name
        colname = self.data.columns.name
        data = self.data.copy()
        
        # Cluster -> recalculate data!
        if dendrogram:
            data, linkage = cluster_gnn_map(data, by=rowname)
        else:
            dcols = 0

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
            self._kw['x_range'] = FactorRange(factors=factors)
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
        csource = pd.DataFrame({'x': dcols + 0.6, 'group': sorted(data.index), 'color': _get_palette(n=ngroups)})
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
            color_bar = ColorBar(color_mapper=color_mapper,  location=(1,0),
                                 ticker=BasicTicker(), border_line_color=None,
                                 title="GNN Proportion")
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


    def bar_chart(self, groups):
        data = self.data.copy()
        data["Haplotype"] = list(map(lambda x: f"{x[0]}/{x[1]}",
                                           zip(data["Individual"],
                                               data[data.columns[0]] % 2 + 1)))
        data["factors"] = list(zip(data["Species"],
                                         data["Haplotype"]))
        source = ColumnDataSource(data)

        if 'x_range' not in self._kw.keys():
            self._kw['x_range'] = FactorRange(*data["factors"])
        self._fig = plotting.figure(**self._kw)
        self._fig.vbar_stack(groups, source=source, x='factors', color=_get_palette(n=len(groups)),
                             legend_label=sorted(groups), width=1,
                             line_color='black', line_alpha=0.2, line_width=0.5)
        self._fig.add_layout(self._fig.legend[0], 'right')

        hover = HoverTool()
        hover.tooltips = [
            ("Haplotype", "@Haplotype"),
            ("Species", "@Species")
        ]
        hover.tooltips.extend(list(map(lambda x: (x[0], f"@{x[1]}"), zip(groups, groups))))
        self._fig.add_tools(hover)

        self._fig.axis.major_tick_line_color = None
        self._fig.axis.minor_tick_line_color = None
        self._fig.xaxis.group_label_orientation = 1.0
        self._fig.xaxis.subgroup_label_orientation = 1.0
        self._fig.xaxis.major_label_orientation = 1.0
        self._fig.xaxis.major_label_text_font_size = '0pt'
        self._fig.axis.axis_line_color = None
        self._fig.grid.grid_line_color = None
        self._fig.outline_line_color = 'black'

        return self._fig


class TSFigure:
    def __init__(self, ts, *arg, **kw):
        self._ts = ts
        self._fig = None
        self._data = None
        self._kw = kw

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        self._data = data
        if self._data.index.name is None:
            self._data.index.name = "Row"
        if self._data.columns.name is None:
            self._data.columns.name = "Column"
    
    @property
    def ts(self):
        return self._ts

    @property
    def populations(self):
        return [json.loads(pop.metadata.decode())["population"] for pop in self.ts.populations()]

    def group_samples(self, by="population"):
        sample_group_set_map = collections.defaultdict(list)
        for population in ts.populations():
            md = json.loads(population.metadata.decode())
            key = md[by]
            sample_group_set_map[key].extend(list(self.ts.samples(
                population=population.id)))
        groups = list(sample_group_set_map.keys())
        sample_group_sets = [sample_group_set_map[k] for k in groups]
        return groups, sample_group_sets

    
    def fst(self, by="population", **kw):
        groups, sample_group_sets = self.group_samples(by=by)
        k = len(list(ts.populations()))
        i = list(itertools.product(list(range(k)), list(range(k))))
        self.ts.Fst(sample_group_sets, indexes=i, **kw)
        df = pd.DataFrame(np.reshape(fst, newshape=(k, k)),
                          columns=groups, index=groups)
        df.index.name = by
        df.columns.name = by
        self.set_data(df)
        

    def heatmap(self, **kw):
        data = self.data.copy()
        
        # Cluster -> recalculate data!
        if dendrogram:
            data, linkage = cluster_gnn_map(data, by=rowname)
        else:
            dcols = 0

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
            self._kw['x_range'] = FactorRange(factors=factors)
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
        csource = pd.DataFrame({'x': dcols + 0.6, 'group': sorted(data.index), 'color': _get_palette(n=ngroups)})
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
            color_bar = ColorBar(color_mapper=color_mapper,  location=(1,0),
                                 ticker=BasicTicker(), border_line_color=None,
                                 title="GNN Proportion")
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
        



class Figure:
    def __init__(self, data, *args, **kw):
        self._data = data
        self._fig = None
        self._kw = kw
        

class MatrixFigure(Figure):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)

    def _process_colors(self, colors, axis=0):
        labels = None
        if colors is not None:
            if colors is True:
                if axis == 0:
                    index = self._data.data.index
                else:
                    index = self._data.data.T.index
                d = {}
                for x in index.names:
                    levels = sorted(index.get_level_values(x).unique())
                    d[x] = list(map(str, _get_palette(n=len(levels))))
                colors = pd.DataFrame(d, index=levels)
            if isinstance(colors, (pd.Series, pd.DataFrame)):
                if axis == 0:
                    colors = colors.reindex(self._data.data.index)
                else:
                    colors = colors.reindex(self._data.data.columns)
                if isinstance(colors, pd.DataFrame):
                    labels = list(colors.columns)
                else:
                    if colors.name is None:
                        labels = [""]
                    else:
                        labels = [colors.name]
        return colors, labels

    def _add_dendrogram(self, padding, axis=0):
        '''Add dendrogram'''
        results = scipy.cluster.hierarchy.dendrogram(self._data.linkage(axis=axis),
                                                     no_plot=True)
        if axis == 0:
            ymax = float(len(self._data.linkage(axis=axis))) + 0.5
            xmax = padding
        else:
            ymax = padding
            xmax = float(len(self._data.linkage(axis=axis))) + 0.5
            
        ycoord = pd.DataFrame(results['icoord'])
        ycoord = ycoord * (ymax / ycoord.max().max())
        ycoord = - ycoord.values + ymax + 0.5

        xcoord = pd.DataFrame(results['dcoord'])
        xcoord = xcoord * (xmax / xcoord.max().max()) - xmax
        xcoord = xcoord.values
        
        for x, y in zip(xcoord, ycoord):
            x = list(map(lambda z: -z, x))
            self._fig.line(x=x, y=y, line_color='black')
            

        
    def heatmap(self, low=None, high=None, palette="Viridis256",
                row_colors=None, col_colors=None, color_bar=True,
                col_cluster=False, row_cluster=False, cbar_title=None,
                dendrogram=True, dendrogram_ratio=.2, *args, **kwargs):
        '''Make heatmap of data, possibly with dendrogram

        Args:
            data (:class:`~pd.DataFrame`) :
                A DataFrame matrix with 'z' values to plot
 
            palette (str or seq[color], optional) :
                A palette to use to colormap z

        Any additional keyword arguments are passed to
        :func:`bokeh.plotting.rect`

        '''
        data = self._data
        padding_cols, padding_rows = 0, 0
        row_colors, row_labels = self._process_colors(row_colors)
        col_colors, col_labels = self._process_colors(col_colors, axis=1)

        if row_cluster:
            data.cluster()
        if col_cluster:
            data.cluster(axis=1)
        def _factors(labels, cluster, axis=0):
            padding = 0
            factors = list(data.linkage_index(axis))
            if labels:
                factors = labels + factors
            if cluster and dendrogram:
                padding = int(dendrogram_ratio * len(factors))
                factors = list(map(lambda x: f"__{x}", range(padding))) + list(factors)
            return factors, padding

        row_factors, row_padding = _factors(row_labels, row_cluster)
        col_factors, col_padding = _factors(col_labels, col_cluster, axis=1)

        # Setup ranges
        if 'x_range' not in self._kw.keys():
            self._kw['x_range'] = FactorRange(factors=row_factors)
        if 'y_range' not in self._kw.keys():
           self._kw['y_range'] = FactorRange(factors=list(reversed(col_factors)))
        self._fig = plotting.figure(**self._kw)

        # Apparently reordering is necessary? I thought it would
        # suffice to set the factors on the x/y ranges
        df = pd.DataFrame(data.reorder().data.stack(), columns=['z']).reset_index()
        source = ColumnDataSource(df)

        # Setup color mapper
        if low is None:
            low = df.z.min()
        if high is None:
            high = df.z.max()
        color_mapper = LinearColorMapper(palette=palette, low=low, high=high)

        # Make heatmap 
        self._fig.rect(x=data.colname, y=data.rowname,
                       width=1, height=1,
                       line_color='black', line_alpha=0.2, alpha=0.8,
                       source=source, fill_color=transform('z', color_mapper),
                       *args, **kwargs)

                             
        # Add dendrogram
        if dendrogram:
            if row_cluster:
                self._add_dendrogram(row_padding, axis=0)
            if col_cluster:
                # FIXME: currently not correct
                # self._add_dendrogram(col_padding, axis=1)
                pass
            
            
        # Setup hover tooltips
        hover = HoverTool()
        hover.tooltips = [
            ("Row", f"@{data.rowname}"),
            ("Column", f"@{data.colname}"),
            ("Value", "@z")
        ]
        self._fig.add_tools(hover)

        # Add colorbar
        if color_bar:
            color_bar = ColorBar(color_mapper=color_mapper, location=(1,0),
                                 ticker=BasicTicker(), border_line_color=None,
                                 title=cbar_title)
            self._fig.add_layout(color_bar, "right")

        # Add group colors
        if row_colors is not None:
            row_colors = row_colors.stack()
            row_colors = row_colors.reset_index()
            row_colors.columns = ["y", "x", "color"]
            self._fig.rect(
                x='x', y='y', fill_color='color',
                line_color=None, line_alpha=0.2, alpha=0.8,
                source=row_colors, width=.6, height=1)


        # Configure axes
        self._fig.axis.major_tick_line_color = None
        self._fig.axis.minor_tick_line_color = None
        self._fig.xaxis.major_label_overrides = dict(zip(row_factors, map(lambda x: '' if x.startswith("__") else x, row_factors)))
        self._fig.xaxis.major_label_orientation = 1.0
        self._fig.axis.axis_line_color = None
        self._fig.grid.grid_line_color = None
        self._fig.outline_line_color = None
        return self._fig
        
    
class Matrix:
    def __init__(self, data, *args, **kw):
        self.data = data
        self._row_linkage = None
        self._col_linkage = None
        self._row_colors = None
        self._col_colors = None

    @property
    def is_square(self):
        return self.data.shape[0] == self.data.shape[1]
        
    def order(self, axis=0):
        """Return order of axis indices. Will change if clustering"""
        if axis == 1:
            return self.col_order
        return self.row_order

    def linkage_index(self, axis=0):
        """Return linkage index"""
        # By default return other axis if no linkage and square matrix
        if self.linkage(axis) is None and self.is_square:
            if self.linkage(1-axis) is not None:
                axis = 1 - axis
        order = self.order(axis=axis)
        if axis == 0:
            return self.data.index.values[order]
        elif axis == 1:
            return self.data.T.index.values[order]

    def linkage(self, axis=0):
        if axis == 0:
            return self.row_linkage
        elif axis == 1:
            return self.col_linkage
        
    @property
    def row_order(self):
        if self.row_linkage is None:
            return list(range(self.data.shape[0]))
        return scipy.cluster.hierarchy.leaves_list(self.row_linkage)

    @property
    def col_order(self):
        if self.col_linkage is None:
            return list(range(self.data.shape[1]))
        return scipy.cluster.hierarchy.leaves_list(self.col_linkage)

    @property
    def row_colors(self):
        return self._row_colors

    @row_colors.setter
    def row_colors(self, value):
        self._row_colors = value

    @property
    def col_colors(self):
        return self._col_colors

    @col_colors.setter
    def col_colors(self, value):
        self._col_colors = value

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        if not isinstance(data, pd.DataFrame):
            print("only data frames allowed")
            raise Exception
        self._data = data
        if self.rowname is None:
            self.rowname = "Row"
        if self.colname is None:
            self.colname = "Column"

    # Will fail on multiindices
    @property
    def rowname(self):
        return self._data.index.name

    @rowname.setter
    def rowname(self, value):
        self._data.index.name = value
    
    @property
    def colname(self):
        return self._data.columns.name

    @colname.setter
    def colname(self, value):
        self._data.columns.name = value
    
    @property
    def row_linkage(self):
        return self._row_linkage

    @row_linkage.setter
    def row_linkage(self, value):
        self._row_linkage = value

    @property
    def col_linkage(self):
        return self._col_linkage

    @col_linkage.setter
    def col_linkage(self, value):
        self._col_linkage = value

    def copy(self):
        obj = Matrix(self.data)
        # Check and set all other attributes?
        obj.row_linkage = self.row_linkage
        obj.col_linkage = self.col_linkage
        return obj

    # Unnecessary? In any case make use of linkage_order
    def reorder(self, axis=0, inplace=False, both=True):
        if inplace:
            obj = self
        else:
            obj = self.copy()
        data = obj.data
        if axis == 0:
            order = obj.row_order
        elif axis == 1:
            data = data.T
            order = obj.col_order
        data = data.reindex(obj.linkage_index(axis=axis))
        if both:
            data = data.T
            indexnames = data.index.names
            data.reindex(index=data.T.index)
            data.index.names = indexnames
            data = data.T
        if axis == 0:
            obj.data = data
        elif axis == 1:
            obj.data = data.T
        return obj


    def _zscore(self, axis=1):
        '''Standardize mean and variance. 0=rows, 1=columns'''
        if axis == 1:
            zscore = self.data.copy()
        else:
            zscore = self.data.copy().T
        for col in list(zscore):
            zscore[col] = scipy.stats.zscore(zscore[col])
        if axis == 1:
            self.data = zscore
        else:
            self.data = zscore.T
        return self

    def rescale(self, zscore=True, standardize=False, axis=1):
        if zscore:
            return self._zscore(axis=axis)
        return self


    def cluster(self, axis=0, method="average", metric="euclidean",
                **kw):
        optimal_ordering = kw.pop("optimal_ordering", True)
        return self._calculate_linkage(axis=axis, method=method,
                                       optimal_ordering=optimal_ordering, **kw)
        
    def _calculate_linkage(self, axis=0, method="average",
                          optimal_ordering=True, **kw):
        data = self.data.copy()
        if axis == 1:
            data = data.T
        linkage = scipy.cluster.hierarchy.linkage(data, method=method,
                                                  optimal_ordering=optimal_ordering, **kw)
        if axis == 0:
            self.row_linkage = linkage
        else:
            self.col_linkage = linkage
        return self
            

            

class GNNMatrix:
    pass

class TSMatrix(Matrix):
    def __init__(self, ts, *args, **kw):
        super().__init__(*args, **kw)
        self._ts = ts


    @property
    def ts(self):
        return self._ts


    def _group_samples(self, by="population"):
        sample_group_set_map = collections.defaultdict(list)
        for population in self.ts.populations():
            md = json.loads(population.metadata.decode())
            key = md[by]
            sample_group_set_map[key].extend(list(self.ts.samples(
                population=population.id)))
        groups = list(sample_group_set_map.keys())
        sample_group_sets = [sample_group_set_map[k] for k in groups]
        return groups, sample_group_sets

    
    def fst(self, by="population", **kw):
        groups, sample_group_sets = self._group_samples(by=by)
        k = len(list(self.ts.populations()))
        i = list(itertools.product(list(range(k)), list(range(k))))
        fst = self.ts.Fst(sample_group_sets, indexes=i, **kw)
        df = pd.DataFrame(np.reshape(fst, newshape=(k, k)),
                          columns=groups, index=groups)
        df.index.name = by
        df.columns.name = by
        self.data = df
        return self

    

def figure(data, **kwargs):
    if isinstance(data, Matrix):
        return MatrixFigure(data, **kwargs)
    # elif figtype == "ts":
    #     return TSFigure(data, **kwargs)


