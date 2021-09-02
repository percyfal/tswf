from pathlib import Path
from snakemake.utils import validate
import pandas as pd

# Determine wrapper prefix since we mix local wrappers with wrappers
# from snakemake-wrappers
SMK_WRAPPER_VERSION = "0.67.0"
SMK_WRAPPER_PREFIX_RAW = "https://github.com/snakemake/snakemake-wrappers/raw"
SMK_WRAPPER_PREFIX = f"{SMK_WRAPPER_PREFIX_RAW}/{SMK_WRAPPER_VERSION}"
WRAPPER_PREFIX = workflow.wrapper_prefix.rstrip("/")
if WRAPPER_PREFIX == SMK_WRAPPER_PREFIX_RAW:
    # Change main to version number once we start sem-versioning
    WRAPPER_PREFIX = "https://raw.githubusercontent.com/NBISweden/manticore-smk/main/workflow/wrappers"


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


##############################
# Core configuration
##############################
include: "core/config.py"


##### load config and sample sheets #####
configfile: Path("config/config.yaml")


validate(config, schema="../schemas/config.schema.yaml")

populations = PopulationData(config["populations"])
samples = SampleData(config["samples"])
samples.merge(populations, left_on="population", right_index=True)

# Wrap config dictionary
cfg = Config(config, samples)

##############################
# Paths and globals
##############################
__INTERIM__ = Path("data/interim")
__RESULTS__ = Path("results")
__RAW__ = Path("data/raw")

try:
    with open(f"{cfg.ref}.fai") as fh:
        refdict = dict(
            [
                (x.split()[0], int(x.split()[1]))
                for x in fh.readlines()
                if x.startswith("chr")
            ]
        )
except Exception as e:
    logger.error("please run samtools faidx on the reference file!")
    raise


##############################
## Wildcard constraints
##############################
wildcard_constraints:
    analysis=wildcards_or(cfg.analysisnames),
    chrom=wildcards_or(cfg.chromosomes),
    dot="(.|)",
    interim=str(__INTERIM__),
    population=wildcards_or(cfg.samples.populations),
    results=str(__RESULTS__),
    sample=wildcards_or(cfg.samples.samples, empty=True),


##################################################
# Input collection functions
##################################################
def all(wildcards):
    d = {"tsinfer": all_tsinfer(wildcards)}
    return d


def all_tsinfer(wildcards):
    """Collect all tsinfer targets"""
    out = []
    for analysis in cfg.analyses("tsinfer"):
        logger.info(f"Collecting files for tsinfer: {analysis.name}")
        # FIXME: move to Analysis object
        fmt = (
            __RESULTS__
            / f"{analysis.name}/{analysis.dataset}/{analysis.fmt}{{plot}}.png"
        )
        d = {"chrom": analysis.chromosomes, "plot": [".gnn", ".R.gnn", ".mean"]}
        out.extend(expand(fmt, **d))
    return out
