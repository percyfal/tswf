import os
import re
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
# load config and sample sheets
schema = get_schema("WORKFLOW_CONFIGURATION_SCHEMA")


configfile: Path("config/config.yaml")


schema.validate(config)

schema = get_schema("SAMPLES_SCHEMA")
populations = PopulationData(config["populations"])
schema = get_schema("SAMPLES_SCHEMA")
samples = SampleData(config["samples"])
samples.merge(populations, left_on="population", right_index=True)

# Add git information
add_gitinfo(config, workflow)
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
        refdict = dict([(x.split()[0], int(x.split()[1])) for x in fh.readlines()])
except Exception as e:
    logger.error("please run samtools faidx on the reference file!")
    raise

# Map all files found in data/raw to chromosomes
dataset = __RAW__ / "variants" / cfg.dataset
variants = dict()
VCF = "(.vcf.gz|.bcf|.vcf)"
regex = re.compile(rf"(.*){VCF}$")
for p in dataset.iterdir():
    filt = filter(lambda x: regex.match(str(x)), list(p.iterdir()))
    filt = list(map(lambda x: regex.sub("\\1", str(x.name)), filt))
    variants[str(p.name)] = filt[0]

##############################
## Wildcard constraints
##############################
BCF = "(.vcf.gz|.bcf)"


wildcard_constraints:
    bcf=BCF,
    chrom=wildcards_or(cfg.chromosomes),
    dot="(.|)",
    interim=str(__INTERIM__),
    population=wildcards_or(cfg.samples.populations),
    prefix="(.+|)",
    results=str(__RESULTS__),
    sample=wildcards_or(cfg.samples.samples, empty=True),
    vcf=VCF,


##################################################
# Input collection functions
##################################################
def all(wildcards):
    d = {"tsinfer": all_tsinfer(wildcards)}
    return d


def all_tsinfer(wildcards):
    """Collect all tsinfer targets"""
    out = []
    out.append(
        __RESULTS__
        / "tsinfer"
        / f"{cfg.analysis}"
        / f"{cfg.dataset}"
        / "population.eda.html"
    )
    return out


##############################
# tsinfer
##############################
def tsinfer_eda_input(wildcards):
    """Return input files for tsinfer eda"""
    fmt = (
        Path(f"{wildcards.results}")
        / "tsinfer"
        / f"{wildcards.analysis}"
        / f"{wildcards.dataset}"
        / "{chrom}"
        / f"{wildcards.mode}"
        / "{prefix}"
    )
    gnn = expand(
        f"{str(fmt)}.gnn.csv",
        zip,
        chrom=cfg.chromosomes,
        prefix=[str(variants[c]) for c in cfg.chromosomes],
    )
    fmt = (
        __INTERIM__
        / "tsinfer"
        / f"{wildcards.analysis}"
        / f"{wildcards.dataset}"
        / "infer"
        / "{chrom}"
        / "{prefix}"
    )
    trees = expand(
        f"{str(fmt)}.trees",
        zip,
        chrom=cfg.chromosomes,
        prefix=[variants[c] for c in cfg.chromosomes],
    )
    return {"gnn": gnn, "trees": trees}


##################################################
# Formatting functions
##################################################
def fmt(wildcards):
    value = cfg.get_analysis(wildcards.analysis).fmt
    if len(cfg.derive_aa.outgroups) != 0:
        value = value + f"_AA_{cfg.derive_aa.method}"
    return value
