import sys

import pytest

from tswf.config import Config
from tswf.config import Schema
from tswf.config import SchemaFiles
from tswf.config import get_schema


params = ["object"] + [x for x in dir(SchemaFiles) if x.endswith("SCHEMA")]


@pytest.fixture
def objectschema():
    d = {
        'properties': {
            'envmodules': {
                'type': 'object',
                'default': {},
                'example': {'bwa': ['bwa/0.7.17']},
            }
        }
    }
    return Schema(d)


@pytest.fixture(params=params)
def schema(request, objectschema):
    if request.param == "object":
        return objectschema
    return get_schema(request.param)


def test_config_from_schema(schema):
    Config.from_schema(schema, file=sys.stdout, example=True)
