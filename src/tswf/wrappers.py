"""Command wrappers"""
import logging
import shutil
import subprocess

logger = logging.getLogger(__name__)


class Wrapper:
    """Documentation wrapper base class."""

    def __init__(self):
        pass

    def run(self, *args, **kwargs):
        """Run wrapper using subprocess.run"""
        raise NotImplementedError


def snakemake(targets=None, options=None, snakefile=None):
    cmdlist = [
        "snakemake",
        f"{'-s ' + str(snakefile) if snakefile else ''}",
        f"{str(options) or ''}",
        f"{str(targets) or ''}",
    ]
    cmd = " ".join(cmdlist)
    if shutil.which("snakemake") is None:
        logger.error("snakemake not installed; cannot run command:")
        logger.error(f"  {cmd}")
        return

    try:
        logger.debug(f"running {cmd}")
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"{cmd} failed")
        print(e)
        raise


def nextflow(workflow, options=None):
    cmdlist = ["nextflow", "run", f"{str(options) or ''}", workflow]
    cmd = " ".join(cmdlist)

    if shutil.which("nextflow") is None:
        logger.error("nextflow not installed; cannot run command:")
        logger.error(f"    {cmd}")
        return

    try:
        logger.debug(f"running {cmd}")
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"{cmd} failed")
        print(e)
        raise
