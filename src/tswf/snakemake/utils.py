"""Snakemake utilities"""

import logging
from typing import Callable

import click
from click.decorators import FC


logger = logging.getLogger(__name__)


def get_profile(uri: str, config: dict[str, dict[str, str]]) -> None | str:
    """Retrieve snakemake profile from config"""
    uristr = None
    try:
        uristr = config["snakemake_profiles"][uri]
    except TypeError as err:
        logger.debug("TypeError: '%s'", err)
    except KeyError as err:
        logger.debug("KeyError: no such snakemake-profiles key %s", err)
    finally:
        logger.debug("trying snakemake profile at uri '%s'", uristr)
    return uristr


def profile_opt(default_profile: None | str = None) -> Callable[[FC], FC]:
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


def jobs_opt() -> Callable[[FC], FC]:
    """Add jobs option"""

    return click.option("--jobs", "-j", type=int, default=1, help="snakemake jobs")


def test_opt() -> Callable[[FC], FC]:
    """Add test option"""
    return click.option(
        "--test", is_flag=True, help="run workflow on small test data set"
    )


def testdir_opt() -> Callable[[FC], FC]:
    """Add test directory option option"""
    return click.option(
        "--testdir",
        type=click.Path(),
        help="setup workflow test in test directory",
    )


def directory_opt() -> Callable[[FC], FC]:
    """Add snakemake directory option"""
    return click.option(
        "--directory/-d",
        help=(
            "Specify working directory (relative paths in the snakefile will "
            "use this as their origin). (default: project directory)"
        ),
    )
