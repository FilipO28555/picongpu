"""
This file is part of PIConGPU.
Copyright 2025 PIConGPU contributors
Authors: Julian Lenz
License: GPLv3+
"""

from .layout import Layout

from .... import util


class Random(Layout):
    _name = "random"

    ppc = util.build_typesafe_property(int)
    """particles per cell (random layout), >0"""

    def _get_serialized(self) -> dict | None:
        return {"ppc": self.ppc}
