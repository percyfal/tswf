import os
import pandas as pd
from pathlib import Path

from snakemake.utils import validate
from snakemake.io import load_configfile
from tswf.config import get_schema
from tswf.snakemake.config import add_gitinfo
from tswf.snakemake.config import PopulationData
from tswf.snakemake.config import SampleData
from tswf.snakemake.config import Config
from tswf.snakemake.config import wildcards_or


container: "docker://continuumio/miniconda3"


envmodules_file = Path(os.environ.get("TSWF_ENVMODULES", "config/envmodules.yaml"))
try:
    envmodules = load_configfile(envmodules_file)["envmodules"]
except FileNotFoundError:
    envmodules = {}
schema = get_schema("ENVMODULES_SCHEMA")
schema.validate(envmodules)

##############################
# Core configuration
##############################
# FIXME: remove
# include: "core/config.py"


##### load config and sample sheets #####
schema = get_schema("WORKFLOW_CONFIGURATION_SCHEMA")


configfile: Path("config/config.yaml")


schema.validate(config)

# validate(config, schema="../schemas/config.schema.yaml")
schema = get_schema("SAMPLES_SCHEMA")
populations = PopulationData(config["populations"])
samples = SampleData(config["samples"])
samples.merge(populations, left_on="population", right_index=True)

# Add git information
add_gitinfo(config, workflow)
# Wrap config dictionary
config["samples"] = samples
cfg = Config(config, samples)

##############################
# Paths and globals
##############################
__INTERIM__ = Path("data/interim")
__RESULTS__ = Path("results")
__RAW__ = Path("data/raw")

try:
    with open(f"{cfg.ref}.fai") as fh:
        refdict = dict([(x.split()[0], int(x.split()[1])) for x in fh.readlines()])
except Exception as e:
    logger.error("please run samtools faidx on the reference file!")
    raise


##############################
## Wildcard constraints
##############################
wildcard_constraints:
    bcf="(.vcf.gz|.bcf)",
    chrom=wildcards_or(cfg.chromosomes),
    dot="(.|)",
    interim=str(__INTERIM__),
    population=wildcards_or(cfg.samples.populations),
    prefix="(.+|)",
    results=str(__RESULTS__),
    sample=wildcards_or(cfg.samples.samples, empty=True),
    suffix="([_\-\.].+|)",
    vcf="(.vcf.gz|.vcf|.bcf)",


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
        eda = __RESULTS__ / f"{analysis.name}/{analysis.dataset}/population.eda.html"
        targz = __RESULTS__ / f"{analysis.name}/{analysis.dataset}/gnn.population.tar.gz"
        out.append(eda)
        out.append(targz)
    return out


##############################
# tsinfer
##############################
def tsinfer_eda_input(wildcards):
    """Return input files for tsinfer eda"""
    csv = expand(
        expand(
            "{{{{results}}}}/{{{{analysis}}}}/{{{{dataset}}}}/{fmt}.gnn.{{{{mode}}}}.csv",
            fmt=fmt(wildcards),
        ),
        chrom=cfg.get_analysis(wildcards.analysis).chromosomes,
    )
    trees = expand(
        expand(
            __INTERIM__ / "{{{{analysis}}}}/{{{{dataset}}}}/{fmt}.trees",
            fmt=fmt(wildcards),
        ),
        chrom=cfg.get_analysis(wildcards.analysis).chromosomes,
    )
    return {"csv": csv, "trees": trees}


##################################################
# Formatting functions
##################################################
def fmt(wildcards):
    value = cfg.get_analysis(wildcards.analysis).fmt
    if len(cfg.derive_aa.outgroups) != 0:
        value = value + f"_AA_{cfg.derive_aa.method}"
    return value
