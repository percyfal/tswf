from ruamel.yaml import YAML
from tswf.cli import ROOT_DIR
from tswf.config import load_config
from tswf.config import Schema

# Main config resides in tswf.yaml
config = load_config(file=ROOT_DIR / "tswf.yaml")

SNAKEMAKE_ROOT = ROOT_DIR / "src" / "workflow"


class SchemaFiles:
    CONFIGURATION_SCHEMA = ROOT_DIR / "schemas" / "config.schema.yaml"
    SAMPLES_SCHEMA = ROOT_DIR / "schemas" / "samples.schema.yaml"


def get_schema(schema="CONFIGURATION_SCHEMA"):
    schemafile = getattr(SchemaFiles, schema)
    with open(schemafile) as fh:
        schema = YAML().load(fh)
    return Schema(schema)
