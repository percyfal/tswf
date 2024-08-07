"""Test snakemake workflow."""

import os
import shutil
import subprocess as sp
from pathlib import Path

import pytest


TESTDIR = Path(os.path.dirname(os.path.abspath(__file__)))
ROOTDIR = TESTDIR.parent
SNAKEFILE = ROOTDIR / "src/tswf/workflow/Snakefile"
DATADIR = ROOTDIR / "src/tswf/tests"


testfiles = [
    "src/tswf/tests/config/config.bio2zarr.yaml",
    "src/tswf/tests/config/config.derive_aa.yaml",
    "src/tswf/tests/config/config.yaml",
    "src/tswf/tests/config/envmodules.yaml",
    "src/tswf/tests/data/raw/variants/ooa/1/ooa_1_PASS.vcf.gz",
    "src/tswf/tests/data/raw/variants/ooa/2/ooa_2_PASS.vcf.gz",
    "src/tswf/tests/resources/populations.tsv",
    "src/tswf/tests/resources/ref_tsk_CHB_0.fa",
    "src/tswf/tests/resources/ref_tsk_CHB_0.fa.fai",
    "src/tswf/tests/resources/samples.tsv",
]


def parse_sem(version):
    """Parse semantic version."""
    semver = tuple(map(int, version.split(".")))
    if len(semver) == 1:
        semver = (semver[0], None, None)
    elif len(semver) == 2:
        semver = (semver[0], semver[1], None)
    return semver


def get_conda_frontend():
    """Determine the conda frontend."""
    major, minor, patch = map(
        int,
        sp.run(["conda", "--version"], capture_output=True, text=True)
        .stdout.split()[1]
        .split("."),
    )
    if major >= 24:
        return "conda"
    if (major == 23) and (minor >= 10) and (patch >= 0):
        return "conda"
    return "mamba"


@pytest.fixture
def conda_frontend():
    """Return conda frontend."""
    return get_conda_frontend()


@pytest.fixture
def snakemake_version():
    """Return Snakemake version."""
    version = sp.run(["snakemake", "--version"], capture_output=True, text=True)
    return version.stdout.strip()


@pytest.fixture
def python_version():
    """Return Python version."""
    version = sp.run(["python", "--version"], capture_output=True, text=True)
    return version.stdout.strip().split()[1]


@pytest.fixture
def conda_env_dir(snakemake_version, python_version):
    """Return conda environment directory."""
    pver = "-".join([str(x) for x in parse_sem(python_version)[0:2]])
    sver = "-".join([str(x) for x in parse_sem(snakemake_version)])
    envdir = TESTDIR / "envs" / f".nox/snakemake-python-{pver}-snakemake-{sver}"
    return envdir


@pytest.fixture(name="deps")
def install_deps(conda_env_dir, conda_frontend, num_cores):
    """Install conda environment dependencies."""
    configfile = DATADIR / "config/config.derive_aa.yaml"
    options = [
        "-s",
        str(SNAKEFILE),
        "--directory",
        str(DATADIR),
        "--configfile",
        str(configfile),
        "--conda-create-envs-only",
        "--conda-frontend",
        conda_frontend,
        "--use-conda",
        "--show-failed-logs",
        "--conda-prefix",
        str(conda_env_dir),
        "--cores",
        str(num_cores),
    ]
    cmd = ["snakemake"] + options

    sp.run(cmd, check=True, text=True)


@pytest.fixture()
def testdir(tmpdir_factory):
    """Copy test files."""
    dname = tmpdir_factory.mktemp("test")
    for testfile in testfiles:
        outfile = dname.join(testfile)
        outfile.dirpath().ensure_dir()
        shutil.copy(testfile, outfile)
    return dname / "src/tswf/tests"


def test_workflow(testdir, conda_env_dir, conda_frontend, num_cores, deps):
    """Test workflow."""
    configfile_aa = testdir / "config/config.derive_aa.yaml"
    options = [
        "--directory",
        str(testdir),
        "--conda-frontend",
        conda_frontend,
        "--conda-prefix",
        str(conda_env_dir),
        "--configfile",
        str(configfile_aa),
        "--use-conda",
        "--show-failed-logs",
        "--cores",
        str(num_cores),
        "-s",
        str(SNAKEFILE),
    ]
    cmd = ["snakemake"] + options
    sp.run(cmd, check=True, text=True)
