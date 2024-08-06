"""Lmod module for the environment module system."""

import copy


def get_envmodules(config, envmodules):
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
