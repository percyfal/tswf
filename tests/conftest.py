#!/usr/bin/env python3
import os

import pytest
from click.testing import CliRunner


def pytest_configure(config):
    pytest.dname = os.path.dirname(__file__)
    pytest.project = os.path.dirname(pytest.dname)


@pytest.fixture(autouse=False)
def cd_tmp_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)


@pytest.fixture(scope="function")
def runner(request):
    return CliRunner(mix_stderr=False)
