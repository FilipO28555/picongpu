"""
This file is part of PIConGPU.
Copyright 2025 PIConGPU contributors
Authors: Julian Lenz
License: GPLv3+
"""

from .densityprofile import DensityProfile
import typeguard
from ....rendering.pmaccprinter import PMAccPrinter
from sympy import Expr, lambdify


@typeguard.typechecked
class FreeFormula(DensityProfile):
    _name = "freeformula"

    def __init__(self, density_expression: Expr) -> None:
        self.density_expression = density_expression

    def check(self):
        pass

    def _get_serialized(self) -> dict | None:
        return {"function_body": PMAccPrinter().doprint(self.density_expression)}

    def __call__(self, x, y, z):
        return lambdify("x,y,z", self.density_expression, "numpy")(x, y, z)
