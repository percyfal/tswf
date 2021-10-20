#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import json


def natural_earth(resolution="50m"):
    _files = {
        "50m": os.path.join(os.path.dirname(__file__), "ne_50m_admin_0_countries.json"),
    }
    with open(_files[resolution], "r") as fh:
        data = "".join(fh.readlines())
    return json.loads(data)
