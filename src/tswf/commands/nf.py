"""Run nextflow workflows.

Main entry point to run nextflow workflows.
"""

import logging

import click

import tswf.wrappers as wrappers
from tswf.cli import pass_environment


__shortname__ = __name__.split(".")[-1]


logger = logging.getLogger(__name__)


@click.group(help=__doc__, name=__shortname__)
@click.pass_context
def main(ctx):
    logger.debug(f"Running {__shortname__}")


@main.command(context_settings=dict(ignore_unknown_options=True), help="run help")
@click.argument("nextflow_args", nargs=-1, type=click.UNPROCESSED)
@pass_environment
def run(env, nextflow_args):
    """Run nextflow workflow"""
    options = list(nextflow_args)
    workflow = "hello"
    wrappers.nextflow(options=" ".join(options), workflow=workflow)
