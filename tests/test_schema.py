import pytest
from tswf.config import get_schema


@pytest.fixture
def envmodules():
    data = {"envmodules": {"bcftools_index": ["module1", "module2"]}}
    return data


def test_schemafiles():
    schema = get_schema()
    assert sorted(list(schema.asdict()['properties'].keys())) == sorted(
        ["project_name", "snakemake_profiles"]
    )


def test_envmodules(envmodules):
    schema = get_schema("ENVMODULES_SCHEMA")
    schema.validate(envmodules)
