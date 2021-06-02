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
include: "core/config.smk"

##### load config and sample sheets #####
configfile: Path("config/config.yaml")
validate(config, schema="../schemas/config.schema.yaml")

populations = PopulationData(config["populations"])
samples = SampleData(config["samples"])
samples.merge(populations, left_on="population", right_index=True)

# Wrap config dictionary
cfg = Config(config, samples)

##############################
## Wildcard constraints
##############################
wildcard_constraints:
    population = wildcards_or(cfg.samples.populations),
    sample = wildcards_or(cfg.samples.samples)

##################################################
# Input collection functions
##################################################
def all(wildcards):
    d = {}
    return d


print(cfg.get_analysis("tsinfer/blackwhite").samples.data)
