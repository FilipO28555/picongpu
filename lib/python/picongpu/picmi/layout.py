"""
This file is part of PIConGPU.
Copyright 2021-2024 PIConGPU contributors
Authors: Hannes Troepgen, Brian Edward Marre
License: GPLv3+
"""

import picmistandard
import typeguard

from ..pypicongpu.species.operation.layout import Random, OnePosition


@typeguard.typechecked
class PseudoRandomLayout(picmistandard.PICMI_PseudoRandomLayout):
    # note: is translated from outside, does not do any checks itself
    def check(self):
        """
        check validity of self

        if ok pass silently, raise on error
        """
        assert self.n_macroparticles_per_cell is not None, "macroparticles per cell must be given"
        assert self.n_macroparticles is None, "total number of macrosparticles not supported"

        assert self.n_macroparticles_per_cell > 0, "at least one particle per cell required"

        # Note: Call PICMI check interface once available upstream

    def get_as_pypicongpu(self):
        return Random()


@typeguard.typechecked
class GriddedLayout(picmistandard.PICMI_GriddedLayout):
    # note: is translated from outside, does not do any checks itself

    def __init__(self, *args, **kwargs):
        # The standard seems to have a typo here:
        # PseudoRandomLayout allows for n_macroparticles_per_cell (with an s)
        # while this one only has n_macroparticle_per_cell.
        # We fix this here:
        if "n_macroparticles_per_cell" in kwargs:
            if len(args) == 0:
                args = (kwargs.pop("n_macroparticles_per_cell"),)
            else:
                raise ValueError("You provided n_macroparticles_per_cell and an unnamed first arg.")
        super().__init__(*args, **kwargs)
        # The standard apparently has a typo in its interface.
        self.n_macroparticles_per_cell = self.n_macroparticle_per_cell
        if self.grid is not None:
            raise NotImplementedError("Non-default grid is not implemented.")

    def check(self):
        """
        check validity of self

        if ok pass silently, raise on error
        """
        assert self.n_macroparticles_per_cell is not None, "macroparticles per cell must be given"
        assert self.n_macroparticles_per_cell > 0, "at least one particle per cell required"

        # Note: Call PICMI check interface once available upstream

    def get_as_pypicongpu(self):
        return OnePosition()
