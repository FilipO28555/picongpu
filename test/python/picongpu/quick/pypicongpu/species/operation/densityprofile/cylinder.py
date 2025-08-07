"""
This file is part of PIConGPU.
Copyright 2024-2024 PIConGPU contributors
Authors: Hannes Troepgen, Brian Edward Marre
License: GPLv3+
"""

import unittest
import typeguard

from picongpu.pypicongpu.species.operation.densityprofile import DensityProfile
from picongpu.pypicongpu.species.operation.densityprofile.cylinder import Cylinder
from picongpu.pypicongpu.species.operation.densityprofile.plasmaramp import None_


class TestCylinder(unittest.TestCase):
    def test_inheritance(self):
        """cylinder is a density profile"""
        self.assertTrue(isinstance(Cylinder(), DensityProfile))

    def test_basic(self):
        """simple scenario works, no ramp"""
        c = Cylinder()
        c.density_si = 1.0e20
        c.center_position_si = (0.0, 0.0, 0.0)
        c.radius_si = 1.0e-3
        c.cylinder_axis = (0.0, 1.0, 0.0)
        c.pre_plasma_ramp = None_()  # valid ramp

        # passes
        c.check()

    def test_value_pass_through(self):
        """values are passed through correctly"""
        c = Cylinder()
        c.density_si = 2.5e23
        c.center_position_si = (1.0, 2.0, 3.0)
        c.radius_si = 5.0e-4
        c.cylinder_axis = (1.0, 0.0, 0.0)
        ramp = None_()
        c.pre_plasma_ramp = ramp

        self.assertAlmostEqual(2.5e23, c.density_si)
        self.assertEqual((1.0, 2.0, 3.0), c.center_position_si)
        self.assertAlmostEqual(5.0e-4, c.radius_si)
        self.assertEqual((1.0, 0.0, 0.0), c.cylinder_axis)
        self.assertEqual(ramp, c.pre_plasma_ramp)

    def test_typesafety(self):
        """typesafety is ensured"""
        c = Cylinder()

        for invalid in [None, "str", [], {}]:
            with self.assertRaises(typeguard.TypeCheckError):
                c.density_si = invalid

        for invalid in [None, "str", [], (1.0), (1.0, 2.0)]:
            with self.assertRaises(typeguard.TypeCheckError):
                c.center_position_si = invalid

        for invalid in [None, "str", [], {}]:
            with self.assertRaises(typeguard.TypeCheckError):
                c.radius_si = invalid

        for invalid in [None, "str", [], (1.0, 2.0)]:
            with self.assertRaises(typeguard.TypeCheckError):
                c.cylinder_axis = invalid

        for invalid in [None, "str", [], 0]:
            with self.assertRaises(typeguard.TypeCheckError):
                c.pre_plasma_ramp = invalid

    def test_check_unsetParameters(self):
        """validity check on self for no parameters set"""
        c = Cylinder()
        with self.assertRaises(Exception):
            c.check()

    def test_check_density(self):
        """invalid density"""
        c = Cylinder()
        c.pre_plasma_ramp = None_()
        c.center_position_si = (0.0, 0.0, 0.0)
        c.radius_si = 1.0e-3
        c.cylinder_axis = (0.0, 1.0, 0.0)

        for invalid in [0, -1, -1.0e20]:
            c.density_si = invalid
            with self.assertRaisesRegex(ValueError, ".*density must be > 0.*"):
                c.check()

    def test_check_radius(self):
        """invalid radius"""
        c = Cylinder()
        c.pre_plasma_ramp = None_()
        c.center_position_si = (0.0, 0.0, 0.0)
        c.cylinder_axis = (0.0, 1.0, 0.0)
        c.density_si = 1.0e20

        for invalid in [-1e-5, -123.0]:
            c.radius_si = invalid
            with self.assertRaisesRegex(ValueError, ".*radius must be > sqrt(2)*"):
                c.check()

    def test_rendering(self):
        """check rendering context"""
        c = Cylinder()
        c.density_si = 42.17
        c.center_position_si = (1.0, 2.0, 3.0)
        c.radius_si = 1.0e-3
        c.cylinder_axis = (0.5, 0.5, 0.707)
        c.pre_plasma_ramp = None_()
        c.check()

        context = c.get_rendering_context()
        # optional: verify some structure
        self.assertTrue(context.get("typeID", {}).get("cylinder", False))
        data = context["data"]
        self.assertAlmostEqual(c.density_si, data["density_si"])
        self.assertEqual([{"component": pos} for pos in c.center_position_si], data["center_position_si"])
        self.assertAlmostEqual(c.radius_si, data["radius_si"])
        self.assertEqual([{"component": ax} for ax in c.cylinder_axis], data["cylinder_axis"])

        # expected "no ramp" structure
        expectedContextNoRamp = {
            "typeID": {"exponential": False, "none": True},
            "data": None,
        }
        self.assertEqual(expectedContextNoRamp, data["pre_plasma_ramp"])
