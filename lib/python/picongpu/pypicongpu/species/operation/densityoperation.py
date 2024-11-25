"""
This file is part of PIConGPU.
Copyright 2023-2024 PIConGPU contributors
Authors: Brian Edward Marre
License: GPLv3+
"""

from .operation import Operation

import typeguard


@typeguard.typechecked
class DensityOperation(Operation):
    """
    common interface for all operations that create density
      and the not-placed operation
    """
