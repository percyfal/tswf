import jsonschema
import pandas as pd
import pytest
from tswf.config import get_schema
from tswf.config import Schema
from tswf.snakemake.config import PopulationData
from tswf.snakemake.config import SampleData


@pytest.fixture
def envmodules():
    data = {"envmodules": {"bcftools_index": ["module1", "module2"]}}
    return data


def _sample_read_tsv(args):
    return pd.DataFrame(args)


@pytest.fixture
def samples():
    data = pd.DataFrame([["s1", "p1"], ["s2", "p2"]])
    data.columns = ["SM", "population"]
    return data


@pytest.fixture
def populations():
    data = pd.DataFrame([["pop1", "spec1"], ["pop2", "spec2"]])
    data.columns = ["population", "species"]
    return data


def test_schemafiles():
    schema = get_schema()
    assert sorted(list(schema.asdict()['properties'].keys())) == sorted(
        ["project_name", "snakemake_profiles"]
    )


def test_envmodules(envmodules):
    schema = get_schema("ENVMODULES_SCHEMA")
    schema.validate(envmodules)


def test_sample_data(samples):
    SampleData(samples)


@pytest.mark.parametrize("header", [["SM2", "population"], ["SM", "populations"]])
def test_sample_data_missing(samples, header):
    samples.columns = header
    with pytest.raises(jsonschema.exceptions.ValidationError):
        SampleData(samples)


def test_population_data(populations):
    PopulationData(populations)


@pytest.mark.parametrize(
    "header", [["population2", "species"], ["population", "specise"]]
)
def test_population_data_missing(populations, header):
    populations.columns = header
    with pytest.raises(jsonschema.exceptions.ValidationError):
        PopulationData(populations)


def test_empty_schema():
    Schema(None)
