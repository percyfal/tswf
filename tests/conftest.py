#!/usr/bin/env python3
import os

import pytest
from click.testing import CliRunner


def pytest_configure(config):
    pytest.dname = os.path.dirname(__file__)
    pytest.project = os.path.dirname(pytest.dname)


def pytest_addoption(parser):
    parser.addoption(
        "--num-cores",
        action="store",
        default=4,
        type=int,
        help="Number of cores to use",
    )


@pytest.fixture(autouse=False)
def cd_tmp_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)


@pytest.fixture(scope="function")
def runner(request):
    return CliRunner(mix_stderr=False)


@pytest.fixture
def num_cores(request):
    ncores = request.config.getoption("--num-cores")
    if ncores < 1:
        raise ValueError("Number of cores must be greater than 0")
    return ncores
