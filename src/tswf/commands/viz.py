"""Visualization tools.

Visualization tools. The subcommands can be called as standalone
scripts with prefix tswf-.

"""
import logging

import click
from tswf.cli import tswf_CLI
from tswf.config import PKG_DIR


__shortname__ = __name__.split(".")[-1]


logger = logging.getLogger(__name__)


class tswf_viz_CLI(tswf_CLI):
    module = "tswf.tools.viz"
    cmd_folder = PKG_DIR / "tools" / "viz"


@click.command(cls=tswf_viz_CLI, help=__doc__, name="viz")
@click.pass_context
def main(ctx):
    logger.debug(f"Running {__shortname__}")
