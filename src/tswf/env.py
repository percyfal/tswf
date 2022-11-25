import os
import pathlib
import sys

import click
from tswf.config import Config

try:
    from dotenv import dotenv_values
except Exception:
    pass

from . import __version__


class Environment:
    def __init__(self):
        self.verbose = False
        self.home = pathlib.Path(os.getcwd())
        self.debug = False
        self.dry_run = False
        self.version = __version__
        self.config = Config(data={"project_name": "tswf"})
        try:
            self.dotenv = dotenv_values(self.home / ".env")
        except Exception:
            self.dotenv = dict()

    def log(self, msg, *args):
        """Logs a message"""
        if args:
            msg %= args
        click.echo(msg, file=sys.stderr)

    def vlog(self, msg, *args):
        """Logs a message to stderr"""
        if self.verbose:
            self.log(msg, *args)
