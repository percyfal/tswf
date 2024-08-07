#!/usr/bin/env python3
import json
import os
from typing import Any


def natural_earth(resolution: str = "50m") -> Any | dict[Any, Any]:
    """Return Natural Earth data"""
    _files = {
        "50m": os.path.join(os.path.dirname(__file__), "ne_50m_admin_0_countries.json"),
    }
    with open(_files[resolution]) as fh:
        data = "".join(fh.readlines())
    return json.loads(data)
