"""# Quickstart (WIP).

## Project initialization

Create a directory for your analyses, cd to the directory and
initialize the tswf configuration file:
    \b
    tswf config init

Currently the workflow configuration defines project name and
snakemake profiles (see below).

## Snakemake workflow

The snakemake workflow is run with the command
    \b
    tswf smk run

To get an idea what the workflow does, run the test
    \b
    tswf smk run --test

### Configuration

Create a configuration file config/config.yaml. You can generate an
example file with
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

You can define snakemake profile locations in the tswf configuration
file. The snakemake profile section allows defining multiple profiles
that you can switch between using the --profile option.

To get started, you can generate a local profile with
    \b
    mkdir -p config/local
    tswf config example profile > config/local/config.yaml

and set the 'snakemake-profile' configuration section to
    \b
    snakemake-profiles:
      local: config/local

Now you can activate the profile with
    \b
    tswf smk run --profile local

To facilitate generation of profiles tailored for HPC job managers,
see the snakemake cookiecutter profile page
<https://github.com/Snakemake-Profiles>.

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
