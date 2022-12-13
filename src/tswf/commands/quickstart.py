"""# Quickstart (WIP).

See the repo README for additional information
(https://github.com/percyfal/tswf).


## Project initialization

Start by creating a project configuration file as follows. Create a
directory <project_name> for your analyses, cd to the directory and
initialize the tswf configuration file:
    \b
    tswf config init

By default, this will create a project configuration file named
<project_name>.yaml. Currently the workflow configuration defines
project name and snakemake profiles (see below).

## Snakemake workflow

The snakemake workflow can be run with the command
    \b
    tswf smk run

To get an idea what the workflow does, run the test
    \b
    tswf smk run --test

Running the snakemake workflow requires creating a number of
configuration and sample files:

1. a snakemake configuration file (default: config/config.yaml)
2. a samplesheet file (default: resources/samples.tsv)
3. a population definition file (default: resources/populations.tsv)
4. (OPTIONAL): a snakemake profile (e.g., config/<profile>/config.yaml)

The following sections describe in more detail how to generate example
files, what their purpose is, and where to install them.

### Configuration

The workflow configuration is typically located at config/config.yaml
and contains run parameter settings that you can customize for your
analyses. You can generate an example file with
    \b
    mkdir -p config
    tswf config example workflow > config/config.yaml

In addition to the configuration file, you need to define your input
samples and populations:
    \b
    mkdir -p resources
    tswf config example samples > resources/samples.tsv
    tswf config example populations > resources/populations.tsv

### Input variant files

Input variant files are placed one by one in separate directories
named data/raw/variants/{dataset}/{chrom}. Here, {dataset} is a
wildcard that is expanded from the configuration file setting
'dataset', and {chrom} is a wildcard that must match the entries in
the configuration section 'chromosomes'. The variant file must end
with '.vcf.gz' or '.bcf.'

### Snakemake profiles

Snakemake profiles
(https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
allow tailoring of snakemake options for specific compute
environments. In their most basic form, they consist of a directory
<profile> with a configuration file config.yaml that maps snakemake
options to values. You can also specify resource requirements for
rules, such as runtime and memory usage, which comes in handy when
submitting jobs to a cluster. To facilitate generation of profiles
tailored for HPC job managers, see the snakemake cookiecutter profile
page <https://github.com/Snakemake-Profiles>. For instance, to install
a SLURM configuration profile in config/slurm, install cookiecutter
and run
    \b
    mkdir -p config/slurm
    cookiecutter --output-dir config/slurm gh:Snakemake-Profiles/slurm

Although you can pass the path to the profile to snakemake, to save
typing, you can also define the location of your snakemake profiles in
the tswf configuration file itself. The snakemake profile section
allows defining multiple profiles that you can switch between using
the --profile option. Using the slurm profile above as an example, you
would set the 'snakemake-profile' configuration section in the tswf main
configuration file (<project_name>.yaml) to
    \b
    snakemake-profiles:
      slurm: config/slurm

Now you can activate the profile with
    \b
    tswf smk run --profile slurm

tswf provides a configuration example for a *local* profile (i.e., not
for cluster submission). To generate such a profile do
    \b
    mkdir -p config/local
    tswf config example profile > config/local/config.yaml

and as before, set the 'snakemake-profile' configuration section in
the tswf main configuration file:
    \b
    snakemake-profiles:
      local: config/local

Now you can activate the local profile with
    \b
    tswf smk run --profile local

### Environment modules

If you use the snakemake flag --use-envmodules (in profile
configuration set "use-envmodules: true"), snakemake will use
environment modules for deployment (see
https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-environment-modules).
All tswf rules can be configured to use environment modules.
Configurations are stored in an environment configuration file
(default: config/envmodules.yaml). An example can be generated with
    \b
    tswf config example envmodules

"""
import logging

import click
from tswf.cli import pass_environment

__shortname__ = __name__.split(".")[-1]


logger = logging.getLogger(__name__)


@click.command(name=__shortname__)
@pass_environment
def main(env):
    logger.debug(f"Running {__shortname__}")
    print()
    print(__doc__)
