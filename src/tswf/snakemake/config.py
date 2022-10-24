from ruamel.yaml import YAML
from tswf.cli import PKG_DIR
from tswf.config import Schema

SNAKEMAKE_ROOT = PKG_DIR / "workflow"


class SchemaFiles:
    CONFIGURATION_SCHEMA = PKG_DIR / "schemas" / "config.schema.yaml"
    SAMPLES_SCHEMA = PKG_DIR / "schemas" / "samples.schema.yaml"


def get_schema(schema="CONFIGURATION_SCHEMA"):
    schemafile = getattr(SchemaFiles, schema)
    with open(schemafile) as fh:
        schema = YAML().load(fh)
    return Schema(schema)
