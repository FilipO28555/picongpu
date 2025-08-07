"""
This file is part of PIConGPU.
Copyright 2025 PIConGPU contributors
Authors: Julian Lenz
License: GPLv3+
"""

import unittest
import sympy
from sympy.abc import x, y
from picongpu.pypicongpu.rendering import PMAccPrinter

CODE_SNIPPETS = {
    "abs": (sympy.Abs(x), "pmacc::math::abs(x)"),
    "min": (sympy.Min(x, y), "pmacc::math::min(x, y)"),
    "max": (sympy.Max(x, y), "pmacc::math::max(x, y)"),
    "erf": (sympy.erf(x), "pmacc::math::erf(x)"),
    "exp": (sympy.exp(x), "pmacc::math::exp(x)"),
    "log": (sympy.log(x), "pmacc::math::log(x)"),
    "pow": (sympy.Pow(x, y), "pmacc::math::pow(x, y)"),
    "sqrt": (sympy.sqrt(x), "pmacc::math::sqrt(x)"),
    "cbrt": (sympy.cbrt(x), "pmacc::math::cbrt(x)"),
    "floor": (sympy.floor(x), "pmacc::math::floor(x)"),
    "ceil": (sympy.ceiling(x), "pmacc::math::ceil(x)"),
    "sin": (sympy.sin(x), "pmacc::math::sin(x)"),
    "cos": (sympy.cos(x), "pmacc::math::cos(x)"),
    "tan": (sympy.tan(x), "pmacc::math::tan(x)"),
    "asin": (sympy.asin(x), "pmacc::math::asin(x)"),
    "acos": (sympy.acos(x), "pmacc::math::acos(x)"),
    "atan": (sympy.atan(x), "pmacc::math::atan(x)"),
    "atan2": (sympy.atan2(x, y), "pmacc::math::atan2(x, y)"),
    "sinh": (sympy.sinh(x), "pmacc::math::sinh(x)"),
    "cosh": (sympy.cosh(x), "pmacc::math::cosh(x)"),
    "tanh": (sympy.tanh(x), "pmacc::math::tanh(x)"),
    "asinh": (sympy.asinh(x), "pmacc::math::asinh(x)"),
    "acosh": (sympy.acosh(x), "pmacc::math::acosh(x)"),
    "atanh": (sympy.atanh(x), "pmacc::math::atanh(x)"),
    "fmod2_native": (x % y, "pmacc::math::fmod(x, y)"),
    "fmod_sympy": (sympy.Mod(x, y), "pmacc::math::fmod(x, y)"),
    "rsqrt": (1 / sympy.sqrt(x), "pmacc::math::pow(x, -1.0/2.0)"),
    "int_division": (x // y, "pmacc::math::floor(x/y)"),
    "pi": (
        sympy.S.Pi,
        "pmacc::math::Pi<float_X>::Value",
    ),
    "2_pi": (
        sympy.S.Pi / 2,
        "pmacc::math::Pi<float_X>::halfValue",
    ),
    "pi_half": (
        sympy.S.Pi / 4,
        "pmacc::math::Pi<float_X>::quarterValue",
    ),
    "pi_quarter": (
        2 / sympy.S.Pi,
        "pmacc::math::Pi<float_X>::doubleReciprocalValue",
    ),
    "log8": (sympy.log(x, 8), "pmacc::math::log(x)/pmacc::math::log(8)"),
    # for these two there's technically specialised functions in PMAcc
    # but it's kind of a pain to catch these expressions
    # before they got through the printer
    # because sympy treats them as log(x)/log(2), ...
    # so it's very hard to actually identify such expressions.
    "log2": (sympy.log(x, 2), "pmacc::math::log(x)/pmacc::math::log(2)"),
    "log10": (sympy.log(x, 10), "pmacc::math::log(x)/pmacc::math::log(10)"),
    # PMAcc functionality that's not supported/expressible with sympy:
    #        "round": (round(x), "pmacc::math::round(x)"),
    #    subtly different from fmod (in PMAcc at least):
    #        "remainder": (???, "pmacc::math::remainder(x, y)"),
    #        "trunc": (???, "pmacc::math::trunc(x)"),
    #        "lround": (???, "pmacc::math::lround(x)"),
    #        "llround": (???, "pmacc::math::llround(x)"),
    # Sinc can be used but isn't translated into the pmacc equivalent
    # but the defining formula.
    # Bessel functions are available upon request...
}


class TestPMAccPrinterMeta(type):
    def __new__(cls, name, bases, dictionary):
        # Generate one test for each example in the examples folder
        for code_name, (expression, cpp_code) in CODE_SNIPPETS.items():
            code_name = "test_" + code_name
            dictionary[code_name] = (
                # This is slightly convoluted:
                # Python's semantics around variables implement
                # "sharing" semantics (not even quite reference semantics).
                # Also, lambdas capture the variable and not the value.
                # So after the execution of a loop all lambdas refer to
                # the last value of the loop variable
                # if they tried to capture it.
                # So, we need to eagerly evaluate the `example` variable
                # which we achieve via an immediately evaluated lambda expression.
                # Please excuse my C++ dialect.
                lambda expression, cpp_code: lambda self: self.generic_test(expression, cpp_code)
            )(expression, cpp_code)
        return type.__new__(cls, name, bases, dictionary)


class TestPMAccPrinter(unittest.TestCase, metaclass=TestPMAccPrinterMeta):
    def generic_test(self, expression, cpp_code):
        self.assertEqual(PMAccPrinter().doprint(expression), cpp_code)


if __name__ == "__main__":
    unittest.main()
