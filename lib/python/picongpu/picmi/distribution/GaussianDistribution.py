"""
This file is part of PIConGPU.
Copyright 2024-2024 PIConGPU contributors
Authors: Brian Edward Marre
License: GPLv3+
"""

from ...pypicongpu import species
from ...pypicongpu import util

from .Distribution import Distribution

import typeguard
import typing
import numpy as np


@typeguard.typechecked
class GaussianDistribution(Distribution):
    """
    Describes a density distribution of particles with gaussian up- and down-ramps with a constant density region in
    between

    Will create the following profile:
    - for y < center_front:                density * exp(factor * abs(y - center_front/sigma_front)**power)
    - for center_front <= y <=center_rear: density
    - for y > center_rear:                 density * exp(factor * abs(y - center_rear/sigma_rear)**power)

    with y being the position in the simulation box
    """

    density: float
    """particle number density, [m^-3]"""

    center_front: float
    """center of gaussian ramp at the front, [m]"""
    center_rear: float
    """center of the gaussian ramp at the rear, [m]"""

    sigma_front: float
    """sigma of the gaussian ramp at the front, [m]"""
    sigma_rear: float
    """sigma of the gaussian ramp at the rear, [m]"""

    power: float
    """power used in exponential function, 2 will yield a gaussian, 4+ a super-gaussian, unitless"""
    factor: float
    """sign and scaling factor, must be < 0, unitless"""

    vacuum_front: float
    """size of the vacuum in front of density, gets rounded down to full cells, [m]"""

    lower_bound: typing.Tuple[float, float, float] | typing.Tuple[None, None, None] = (
        None,
        None,
        None,
    )
    upper_bound: typing.Tuple[float, float, float] | typing.Tuple[None, None, None] = (
        None,
        None,
        None,
    )

    cell_size: tuple[float, float, float] | None = None

    # @details pydantic provides an automatically generated __init__/constructor method which allows initialization off
    #   all attributes as keyword arguments

    # @note user may add additional attributes by hand, these will be available but not type verified

    def get_as_pypicongpu(self, grid) -> species.operation.densityprofile.DensityProfile:
        self.cell_size = grid.get_cell_size()
        util.unsupported("fill in not active", self.fill_in, True)

        # @todo support bounds, Brian Marre, 2024
        util.unsupported("lower bound", self.lower_bound, (None, None, None))
        util.unsupported("upper bound", self.upper_bound, (None, None, None))

        gaussian_profile = species.operation.densityprofile.Gaussian()

        if self.center_rear < self.center_front:
            raise ValueError("center_front must be <= center_rear")
        if self.density <= 0.0:
            raise ValueError("density must be > 0")

        # @todo change to constructor call once we switched PyPIConGPU to use pydantic, Brian Marre, 2024
        gaussian_profile.gas_center_front = self.center_front
        gaussian_profile.gas_center_rear = self.center_rear
        gaussian_profile.gas_sigma_front = self.sigma_front
        gaussian_profile.gas_sigma_rear = self.sigma_rear
        gaussian_profile.gas_factor = self.factor
        gaussian_profile.gas_power = self.power
        gaussian_profile.vacuum_cells_front = int(self.vacuum_front / grid.get_cell_size()[1])
        gaussian_profile.density = self.density

        return gaussian_profile

    def __call__(self, x, y, z):
        if self.cell_size is None:
            message = (
                "Due to inconsistencies in the backend, evaluation of this function requires information about the cell_size."
                " You can either set it manually "
                " or you can perform anything that includes writing the input files on your simulation object."
                " This is a temporary workaround and will be fixed in the future."
            )
            raise NotImplementedError(message)

        # The definition of this density uses the origin of the cell
        # while the call operator uses the center.
        x += -0.5 * self.cell_size[0]
        y += -0.5 * self.cell_size[1]
        z += -0.5 * self.cell_size[2]

        # The last term undoes the shift to the cell origin.
        vacuum_y = int(self.vacuum_front / self.cell_size[1]) * self.cell_size[1] - 0.5 * self.cell_size[1]

        # We do this to get the correct shape after broadcasting:
        exponent = 0 * (x + y + z)
        exponent[y < self.center_front] = np.abs((y - self.center_front) / self.sigma_front)[y < self.center_front]
        exponent[y >= self.center_rear] = np.abs((y - self.center_rear) / self.sigma_rear)[y >= self.center_rear]

        result = np.exp(self.factor * exponent**self.power)
        result[y < vacuum_y] = 0.0
        return self.density * result
