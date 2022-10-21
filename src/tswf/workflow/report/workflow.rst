wg-genealogy-smk
================

The wg-genealogy-smk_ workflow runs applications for generation and
analysis of whole-genome genealogies. This analysis is based on commit
version {{ snakemake.config["__workflow_commit__"] }}.

The analysis can be rerun with the following command:

.. code-block:: shell

   cd {{ snakemake.config["__workflow_workdir__"] }}
{% if snakemake.config["__workflow_basedir__"] == snakemake.config["__workflow_workdir__"] %}
   snakemake -j 1 --use-conda
{% else %}
   snakemake -j 1 --use-conda -s {{ snakemake.config["__workflow_basedir__"] }}/Snakefile
{% endif %}

and the report

.. code-block:: shell

   cd {{ snakemake.config["__workflow_workdir__"] }}
{% if snakemake.config["__workflow_basedir__"] == snakemake.config["__workflow_workdir__"] %}
   snakemake -j 1 --report report.html
{% else %}
   snakemake -j 1 -s {{ snakemake.config["__workflow_basedir__"] }}/Snakefile --report report.html
{% endif %}



Supported applications
----------------------

tsinfer
~~~~~~~

- tree generation
- gnn plots

relate
~~~~~~~

WIP


Data organization
=================

.. code-block:: text

   {{ snakemake.config["__workflow_workdir__"] }}/                                <- top-level project folder
   |
   ├── config                   <- configuration directory for Snakemake and other things
   │
   ├── data
   │   ├── external             <- data from third party sources
   │   ├── interim              <- Intermediate data that can be safely deleted
   │   ├── metadata             <- metadata describing raw data files
   │   ├── processed            <- Final processed data used for analyses
   │   └── raw                  <- The original immutable data dump to be treated as read-only.
   │
   ├── logs                     <- Collection of log outputs, e.g. from cluster managers
   │
   ├── reports                  <- Generated analyses and articles as html, pdf and more, including multiqc.html
   │   └── figures              <- Graphics for use in reports.
   │
   └── results                  <- Final results for sharing with collaborators, typically derived


.. _wg-genealogy-smk: https://github.com/percyfal/wg-genealogy-smk
