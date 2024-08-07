"""Lmod module for the environment module system."""

import copy
from typing import Any


def get_envmodules(config: dict[str, Any], envmodules: str | list[str]) -> list[str]:
    """Get the environment modules for a given environment."""
    retmodules = []
    if not isinstance(envmodules, list):
        envmodules = [envmodules]
    try:
        retmodules = copy.deepcopy(config["envmodules"]["__site__"])
    except KeyError:
        retmodules = []
    for mod in envmodules:
        try:
            retmodules.extend(config["envmodules"][mod])
        except KeyError:
            pass
    return retmodules
