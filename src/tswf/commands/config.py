"""Configuration administration utilities.

"""
import logging
import sys

import click
import pkg_resources
from tswf.cli import pass_environment
from tswf.config import Config
from tswf.config import get_schema
from tswf.config import SchemaFiles

logger = logging.getLogger(__name__)

__shortname__ = __name__.split(".")[-1]


@click.group(help=__doc__, name=__shortname__)
@pass_environment
def main(env):
    logger.debug(f"Running {__shortname__} subcommand.")


@main.command()
@click.option("--config-file", help="configuration file name")
@click.option("--show", help="render configuration", is_flag=True)
@pass_environment
def init(env, config_file, show):
    """Initialize a main tswf configuration file.

    By default will save to file PROJECT_NAME.yaml in the project home
    directory, where PROJECT_NAME is derived from the basename
    directory of the configuration file name.

    """
    logger.info("Initializing configuration file")
    project_name = env.home.name
    if config_file is None:
        config_file = env.home / f"{project_name}.yaml"
    schema = get_schema()
    if not config_file.exists():
        if show:
            Config.from_schema(schema, file=sys.stdout, project_name=str(project_name))
        else:
            with open(config_file, "w") as fh:
                Config.from_schema(schema, file=fh, project_name=project_name)
    else:
        logger.info(f"{config_file} exists; not overwriting")


@main.command()
@click.argument(
    "configuration_schema",
    type=click.Choice(("main", "samples", "populations", "envmodules", "workflow")),
    default="main",
)
@click.option(
    "--as-yaml",
    is_flag=True,
    help=(
        "print configuration as yaml. Only takes effect for "
        "samples and populations schemas"
    ),
)
@pass_environment
def example(env, configuration_schema, as_yaml):
    """Show example configuration files."""
    conf_map = {
        "main": "CONFIGURATION_SCHEMA",
        "samples": "SAMPLES_SCHEMA",
        "populations": "POPULATIONS_SCHEMA",
        "envmodules": "ENVMODULES_SCHEMA",
        "workflow": "WORKFLOW_CONFIGURATION_SCHEMA",
    }
    kwargs = {}
    tsv = False
    schema = get_schema(conf_map[configuration_schema])
    schemafile = pkg_resources.resource_filename(
        "tswf", str(getattr(SchemaFiles, conf_map[configuration_schema]))
    )

    required = schema._schema.get("required", None)
    if configuration_schema == "main":
        kwargs = {"project_name": env.home.name}
    if configuration_schema in ["samples", "populations"] and not as_yaml:
        tsv = True

    print()
    print(
        f"#\n# Showing example configuration for schema {conf_map[configuration_schema]}"
    )
    print(f"# See schema file {schemafile} for more details.\n#")
    if tsv:
        print("# Use option --as-yaml to see column descriptions.\n#")
    if required is not None:
        print(f"# Required fields: {','.join(required)}\n#")
    print()

    Config.from_schema(schema, file=sys.stdout, tsv=tsv, example=True, **kwargs)
    print()
