"""Nox sessions."""

import os
import shlex
import sys
import tempfile
from pathlib import Path
from textwrap import dedent

import nox


try:
    from nox_poetry import Session
    from nox_poetry import session
except ImportError:
    message = f"""\
    Nox failed to import the 'nox-poetry' package.

    Please install it using the following command:

    {sys.executable} -m pip install nox-poetry"""
    raise SystemExit(dedent(message)) from None


package = "tswf"
python_versions = ["3.12", "3.11"]

nox.needs_version = ">= 2021.6.6"
nox.options.sessions = (
    "mypy",
    "pre-commit",
    "tests",
)


def parse_sem(version):
    """Parse semantic version."""
    semver = tuple(map(int, version.split(".")))
    if len(semver) == 1:
        semver = (semver[0], None, None)
    elif len(semver) == 2:
        semver = (semver[0], semver[1], None)
    return semver


def pip_install_pulp(session, snakemake_version):
    """Pip install correct pulp version.

    cf https://github.com/snakemake/snakemake/issues/2607.
    """
    major, minor, patch = map(int, snakemake_version.split("."))
    if (major >= 8) and (minor >= 2):
        session.install("pulp>=2.8")
    elif (major == 8) and (minor == 1) and (patch >= 2):
        session.install("pulp>=2.8")
    else:
        session.install("pulp<2.8")


def install_poetry_groups(session, *groups: str) -> None:
    """Install dependencies from poetry groups."""
    with tempfile.NamedTemporaryFile() as requirements:
        session.run(
            "poetry",
            "export",
            *[f"--only={group}" for group in groups],
            "--format=requirements.txt",
            "--without-hashes",
            f"--output={requirements.name}",
            external=True,
        )
        session.install("-r", requirements.name)


def install_poetry_plugins_ci(session: Session) -> None:
    """Install poetry plugins for CI."""
    if "CI" in os.environ:
        session.run("pip", "install", "poetry-plugin-export")
        session.run("pip", "install", "poetry-dynamic-versioning")


def install_snakemake(session):
    """Install a version of snakemake compatible with Python version."""
    major, minor, _ = parse_sem(session.python)
    snakemake_version = "8.16.0"
    if major == 3 and minor <= 10:
        snakemake_version = "7.32.4"
    pip_install_pulp(session, snakemake_version)
    session.install(f"snakemake=={snakemake_version}")


def activate_virtualenv_in_precommit_hooks(session: Session) -> None:
    """Activate virtualenv in hooks installed by pre-commit.

    This function patches git hooks installed by pre-commit to activate the
    session's virtual environment. This allows pre-commit to locate hooks in
    that environment when invoked from git.

    Args:
        session: The Session object.
    """
    assert session.bin is not None  # noqa: S101

    # Only patch hooks containing a reference to this session's bindir. Support
    # quoting rules for Python and bash, but strip the outermost quotes so we
    # can detect paths within the bindir, like <bindir>/python.
    bindirs = [
        bindir[1:-1] if bindir[0] in "'\"" else bindir
        for bindir in (repr(session.bin), shlex.quote(session.bin))
    ]

    virtualenv = session.env.get("VIRTUAL_ENV")
    if virtualenv is None:
        return

    headers = {
        # pre-commit < 2.16.0
        "python": f"""\
            import os
            os.environ["VIRTUAL_ENV"] = {virtualenv!r}
            os.environ["PATH"] = os.pathsep.join((
                {session.bin!r},
                os.environ.get("PATH", ""),
            ))
            """,
        # pre-commit >= 2.16.0
        "bash": f"""\
            VIRTUAL_ENV={shlex.quote(virtualenv)}
            PATH={shlex.quote(session.bin)}"{os.pathsep}$PATH"
            """,
        # pre-commit >= 2.17.0 on Windows forces sh shebang
        "/bin/sh": f"""\
            VIRTUAL_ENV={shlex.quote(virtualenv)}
            PATH={shlex.quote(session.bin)}"{os.pathsep}$PATH"
            """,
    }

    hookdir = Path(".git") / "hooks"
    if not hookdir.is_dir():
        return

    for hook in hookdir.iterdir():
        if hook.name.endswith(".sample") or not hook.is_file():
            continue

        if not hook.read_bytes().startswith(b"#!"):
            continue

        text = hook.read_text()

        if not any(
            Path("A") == Path("a") and bindir.lower() in text.lower() or bindir in text
            for bindir in bindirs
        ):
            continue

        lines = text.splitlines()

        for executable, header in headers.items():
            if executable in lines[0].lower():
                lines.insert(1, dedent(header))
                hook.write_text("\n".join(lines))
                break


@session(name="pre-commit", python=python_versions[0])
def precommit(session: Session) -> None:
    """Lint using pre-commit."""
    args = session.posargs or [
        "run",
        "--all-files",
        "--hook-stage=manual",
        "--show-diff-on-failure",
    ]
    session.install(
        "black",
        "darglint",
        "flake8",
        "flake8-bandit",
        "flake8-bugbear",
        "flake8-docstrings",
        "flake8-rst-docstrings",
        "isort",
        "pep8-naming",
        "pre-commit",
        "pre-commit-hooks",
        "pyupgrade",
    )
    session.run("pre-commit", *args)
    if args and args[0] == "install":
        activate_virtualenv_in_precommit_hooks(session)


@session(python=python_versions)
def tests(session: Session) -> None:
    """Run the test suite."""
    install_poetry_plugins_ci(session)
    session.install(".")
    install_snakemake(session)
    session.install("coverage[toml]", "pytest", "pygments")
    session.run("ls", "src/tswf", external=True)
    try:
        session.run("coverage", "run", "--parallel", "-m", "pytest", *session.posargs)
    finally:
        if session.interactive:
            session.notify("coverage", posargs=[])


@session(python=python_versions[0])
def coverage(session: Session) -> None:
    """Produce the coverage report."""
    args = session.posargs or ["report"]

    session.install("coverage[toml]")

    if not session.posargs and any(Path().glob(".coverage.*")):
        session.run("coverage", "combine")

    session.run("coverage", *args)


@session(python=python_versions)
def mypy(session: Session) -> None:
    """Type-check using mypy."""
    args = session.posargs or ["src", "tests"]
    install_poetry_plugins_ci(session)
    session.install(".")
    session.install("mypy", "pytest")
    install_poetry_groups(session, "dev")
    session.run("mypy", *args)
    if not session.posargs:
        session.run("mypy", f"--python-executable={sys.executable}", "noxfile.py")
