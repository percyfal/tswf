"""Run snakemake workflows.

Main entry point to run snakemake workflows.
"""
import logging

import click
import tswf.snakemake.config as config
import tswf.wrappers as wrappers
from tswf.cli import pass_environment
from tswf.snakemake.utils import get_profile
from tswf.snakemake.utils import jobs_opt
from tswf.snakemake.utils import profile_opt
from tswf.snakemake.utils import test_opt

__shortname__ = __name__.split(".")[-1]


logger = logging.getLogger(__name__)


@click.group(help=__doc__, name=__shortname__)
@click.pass_context
def main(ctx):
    logger.debug(f"Running {__shortname__}")


@main.command(context_settings=dict(ignore_unknown_options=True), help="run help")
@profile_opt()
@jobs_opt()
@test_opt()
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
@pass_environment
def run(env, profile, jobs, test, snakemake_args):
    options = list(snakemake_args) + ["-j", str(jobs)]
    snakefile = config.SNAKEMAKE_ROOT / "Snakefile"
    if profile is not None:
        profile = get_profile(profile, env.config)
        options.extend(["--profile", profile])
    if test:
        options.extend(["--directory", str(config.PKG_DIR / "tests")])
    options = " ".join(options)
    target = ""
    wrappers.snakemake(options=options, snakefile=snakefile, targets=target)
