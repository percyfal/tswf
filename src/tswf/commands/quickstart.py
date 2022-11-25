"""Quickstart

Documentation command.
"""
import logging

import click
from tswf.cli import pass_environment

__shortname__ = __name__.split(".")[-1]


logger = logging.getLogger(__name__)


@click.command(help=__doc__, name=__shortname__)
@pass_environment
def main(env):
    logger.debug(f"Running {__shortname__}")
    print(__doc__)
