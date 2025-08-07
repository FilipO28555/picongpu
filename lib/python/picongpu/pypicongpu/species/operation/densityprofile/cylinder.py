"""
This file is part of PIConGPU.
Copyright 2023-2024 PIConGPU contributors
Authors: Kristin Tippey, Brian Edward Marre
License: GPLv3+
"""

from .densityprofile import DensityProfile
from .plasmaramp import PlasmaRamp, None_
from .... import util

import typeguard
import math


@typeguard.typechecked
class Cylinder(DensityProfile):
    """
     Describes a cylindrical density distribution of particles with gaussian up-ramp
    with a constant density region in between. It can have an arbitrary orientation
    and position in space.

    Will create the following profile:
      n = density if r < reduced_radius
      n is 0 or follows the exponential ramp if r > reduced_radius
      n is 0 if r > reduced_radius + prePlasmaCutoff
      the reduced_radius is equal = @f[\\sqrt{R^2 -L^2} -L @f]
      with R - cylinder_radius and L - prePlasmaLength (scale length of the ramp)
      the reduced radius ensures mass conservation
    """

    _name = "cylinder"

    density_si = util.build_typesafe_property(float)
    """particle number density at at the foil plateau (m^-3)"""

    center_position_si = util.build_typesafe_property(tuple[float, float, float])
    """center of the cylinder [x, y, z], [m]"""

    radius_si = util.build_typesafe_property(float)
    """cylinder radius, [m]"""

    cylinder_axis = util.build_typesafe_property(tuple[float, float, float])
    """cylinder axis [x, y, z], [unitless]"""

    pre_plasma_ramp = util.build_typesafe_property(PlasmaRamp)
    """pre plasma ramp"""

    def __init__(self):
        # (nothing to do, overwrite from abstract parent)
        pass

    def check(self) -> None:
        if self.density_si <= 0.0:
            raise ValueError("density must be > 0")
        min_radius = (
            math.sqrt(2.0) * self.pre_plasma_ramp.PlasmaLength if type(self.pre_plasma_ramp) is not None_ else 0.0
        )
        if self.radius_si < min_radius:
            raise ValueError(
                f"radius must be > sqrt(2)*pre_plasma_length = {min_radius}, so that the reduced radius stays non negative. In case of no preplasma radius must be >= 0.0."
            )
        self.pre_plasma_ramp.check()

    def _get_serialized(self) -> dict:
        self.check()

        return {
            "density_si": self.density_si,
            "center_position_si": [{"component": x} for x in self.center_position_si],
            "radius_si": self.radius_si,
            "cylinder_axis": [{"component": x} for x in self.cylinder_axis],
            "pre_plasma_ramp": self.pre_plasma_ramp.get_rendering_context(),
        }
