"""
This file is part of PIConGPU.
Copyright 2025 PIConGPU contributors
Authors: Julian Lenz
License: GPLv3+
"""

from ....rendering import SelfRegisteringRenderedObject

import typeguard


@typeguard.typechecked
class Layout(SelfRegisteringRenderedObject):
    """
    (abstract) parent class of all layouts

    A layout describes the particle layout within a cell.
    """

    pass
