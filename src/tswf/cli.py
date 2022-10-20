"""Console script for Tree Sequence Workflow (tswf)

Run tree sequence workflows.
"""
import logging
import os
import pathlib

import click
from tswf import decorators
from tswf.config import load_config
from tswf.env import Environment

from . import __version__

__author__ = "Per Unneberg"

logger = logging.getLogger(__name__)


ROOT_DIR = pathlib.Path(__file__).absolute().parent.parent.parent
CONTEXT_SETTINGS = dict(auto_envvar_prefix="TSWF", show_default=True)

pass_environment = click.make_pass_decorator(Environment, ensure=True)
cmd_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), "commands"))
tools_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), "tools"))


class tswf_CLI(click.MultiCommand):
    def list_commands(self, ctx):
        rv = []
        for filename in os.listdir(cmd_folder):
            if filename.endswith(".py") and not filename.startswith("__"):
                rv.append(filename[:-3])
        rv.sort()
        return rv

    def get_command(self, ctx, name):
        try:
            mod = __import__(f"tswf.commands.{name}", None, None, ["main"])
        except ImportError:
            raise
            return
        return mod.main


@click.command(
    cls=tswf_CLI,
    context_settings=CONTEXT_SETTINGS,
    help=__doc__,
    name="tswf",
)
@click.version_option(version=__version__)
@click.option("--config-file", help="configuration file", type=click.Path(exists=True))
@decorators.debug_option()
@pass_environment
def cli(env, config_file):
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s"
    )
    if env.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    env.home = ROOT_DIR
    if config_file is None:
        config_file = env.home / "tswf.yaml"
    if config_file.exists():
        config = load_config(file=config_file)
    else:
        config = load_config(data=dict(project_name="tswf"))
    env.config = config
