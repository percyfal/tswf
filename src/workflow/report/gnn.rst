Genealogical Nearest Neighbours (GNN) plot for chromosome {{ snakemake.wildcards["chrom"] }}

Given a reference set (e.g population) of K classes, we calculate for
a focal node the reference class of the nearest neighbours, defined as
the nodes sharing the same parent as the focal node. The result is a
K-vector describing the proportion of each reference class in the
nearest neighbours set. See the `method description of GNN`_ for more
information.

The plot shows for each individual the GNN average over all nodes,
i.e. all sites, for chromosome {{ snakemake.wildcards["chrom"] }}. The
plot has been produced in R.

.. _method description of GNN: https://www.nature.com/articles/s41588-019-0483-y#Sec7
