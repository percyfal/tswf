"""Command wrappers"""

import logging
import shutil
import subprocess


logger = logging.getLogger(__name__)


def snakemake(
    targets: None | list[str] = None,
    options: None | list[str] = None,
    snakefile: None | str = None,
) -> None:
    """Wrap snakemake call."""
    cmdlist = [
        "snakemake",
        f"{'-s ' + str(snakefile) if snakefile else ''}",
        f"{str(options) or ''}",
        f"{str(targets) or ''}",
    ]
    cmd = " ".join(cmdlist)
    if shutil.which("snakemake") is None:
        logger.error("snakemake not installed; cannot run command:")
        logger.error("  %s", cmd)
        return

    try:
        logger.info("running %s", cmd)
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error("%s failed", cmd)
        print(e)
        raise


def nextflow(workflow: str, options: None | str = None) -> None:
    """Wrap Nextflow call."""
    cmdlist = ["nextflow", "run", f"{str(options) or ''}", workflow]
    cmd = " ".join(cmdlist)

    if shutil.which("nextflow") is None:
        logger.error("nextflow not installed; cannot run command:")
        logger.error("    %s", cmd)
        return

    try:
        logger.debug("running %s", cmd)
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error("%s failed", cmd)
        print(e)
        raise
