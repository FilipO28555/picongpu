"""
This file is part of PIConGPU.
Copyright 2021-2024 PIConGPU contributors
Authors: Masoud Afshari
License: GPLv3+
"""

from .. import util
from ..species import Species

from .plugin import Plugin
from .timestepspec import TimeStepSpec


import typeguard
import typing
from enum import Enum


class EMFieldScaleEnum(Enum):
    AUTO = -1
    PLASMA_WAVE = 3
    CUSTOM = 6
    INCIDENT = 7

    @classmethod
    def _missing_(cls, value):
        """Ensure strings map correctly to Enum values."""
        for member in cls:
            if member.value == value:
                return member
        raise ValueError(f"{value} is not a valid {cls.__name__}")


class ColorScaleEnum(Enum):
    RED = "red"
    GREEN = "green"
    BLUE = "blue"
    GRAY = "gray"
    GRAY_INV = "grayInv"
    NONE = "none"

    @classmethod
    def _missing_(cls, value):
        """Ensure strings map correctly to Enum values."""
        for member in cls:
            if member.value == value:
                return member
        raise ValueError(f"{value} is not a valid {cls.__name__}")


@typeguard.typechecked
class Png(Plugin):
    species = util.build_typesafe_property(Species)
    period = util.build_typesafe_property(TimeStepSpec)
    axis = util.build_typesafe_property(str)
    slicePoint = util.build_typesafe_property(float)
    folder = util.build_typesafe_property(str)
    scale_image = util.build_typesafe_property(float)
    scale_to_cellsize = util.build_typesafe_property(bool)
    white_box_per_GPU = util.build_typesafe_property(bool)
    EM_FIELD_SCALE_CHANNEL1 = util.build_typesafe_property(EMFieldScaleEnum)
    EM_FIELD_SCALE_CHANNEL2 = util.build_typesafe_property(EMFieldScaleEnum)
    EM_FIELD_SCALE_CHANNEL3 = util.build_typesafe_property(EMFieldScaleEnum)
    preParticleDensCol = util.build_typesafe_property(ColorScaleEnum)
    preChannel1Col = util.build_typesafe_property(ColorScaleEnum)
    preChannel2Col = util.build_typesafe_property(ColorScaleEnum)
    preChannel3Col = util.build_typesafe_property(ColorScaleEnum)
    customNormalizationSI = util.build_typesafe_property(typing.List[float])
    preParticleDens_opacity = util.build_typesafe_property(float)
    preChannel1_opacity = util.build_typesafe_property(float)
    preChannel2_opacity = util.build_typesafe_property(float)
    preChannel3_opacity = util.build_typesafe_property(float)
    preChannel1 = util.build_typesafe_property(str)
    preChannel2 = util.build_typesafe_property(str)
    preChannel3 = util.build_typesafe_property(str)

    _name = "png"

    def __init__(self):
        "do nothing"

    def _get_serialized(self) -> typing.Dict:
        """Return the serialized representation of the object."""

        # Transform customNormalizationSI into a list of dictionaries
        custom_normalization_si_serialized = [{"value": val} for val in self.customNormalizationSI]

        return {
            "species": self.species.get_rendering_context(),
            "period": self.period.get_rendering_context(),
            "axis": self.axis,
            "slicePoint": self.slicePoint,
            "folder": self.folder,
            "scale_image": self.scale_image,
            "scale_to_cellsize": self.scale_to_cellsize,
            "white_box_per_GPU": self.white_box_per_GPU,
            "EM_FIELD_SCALE_CHANNEL1": self.EM_FIELD_SCALE_CHANNEL1.value,
            "EM_FIELD_SCALE_CHANNEL2": self.EM_FIELD_SCALE_CHANNEL2.value,
            "EM_FIELD_SCALE_CHANNEL3": self.EM_FIELD_SCALE_CHANNEL3.value,
            "preParticleDensCol": self.preParticleDensCol.value,
            "preChannel1Col": self.preChannel1Col.value,
            "preChannel2Col": self.preChannel2Col.value,
            "preChannel3Col": self.preChannel3Col.value,
            "customNormalizationSI": custom_normalization_si_serialized,
            "preParticleDens_opacity": self.preParticleDens_opacity,
            "preChannel1_opacity": self.preChannel1_opacity,
            "preChannel2_opacity": self.preChannel2_opacity,
            "preChannel3_opacity": self.preChannel3_opacity,
            "preChannel1": self.preChannel1,
            "preChannel2": self.preChannel2,
            "preChannel3": self.preChannel3,
        }
