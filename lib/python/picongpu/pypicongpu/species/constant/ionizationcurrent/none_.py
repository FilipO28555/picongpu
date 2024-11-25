"""
This file is part of PIConGPU.
Copyright 2024-2024 PIConGPU contributors
Authors: Brian Edward Marre
License: GPLv3+
"""

from .ionizationcurrent import IonizationCurrent

import typeguard


@typeguard.typechecked
class None_(IonizationCurrent):
    PICONGPU_NAME: str = "None"
