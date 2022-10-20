"""Snakemake utilities"""
import logging

import click

logger = logging.getLogger(__name__)


def get_profile(uri, config):
    """Retrieve snakemake profile from config"""
    try:
        uri = config["snakemake_profiles"][uri]
    except TypeError as e:
        logger.debug(f"TypeError: '{e}'")
    except KeyError as e:
        logger.debug(f"KeyError: no such snakemake-profiles key {e}")
    finally:
        logger.debug(f"trying snakemake profile at uri '{uri}'")
    return uri


def profile_opt(default_profile=None):
    """Add snakemake profile option"""
    func = click.option(
        "--profile",
        help=(
            "snakemake profile, either defined as key:value pair in config"
            " or a URI pointing to profile directory"
        ),
        default=default_profile,
    )
    return func


def jobs_opt():
    """Add jobs option"""
    return click.option("--jobs", "-j", type=int, default=1, help="snakemake jobs")


def test_opt():
    """Add test option"""
    return click.option(
        "--test", is_flag=True, help="run workflow on small test data set"
    )


def directory_opt():
    """Add snakemake directory option"""
    return click.option(
        "--directory/-d",
        help=(
            "Specify working directory (relative paths in the snakefile will "
            "use this as their origin). (default: project directory)"
        ),
    )
