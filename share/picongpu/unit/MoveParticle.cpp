/* Copyright 2023 Julian Lenz
 *
 * This file is part of PIConGPU.
 *
 * PIConGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIConGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PIConGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#include <pmacc/boost_workaround.hpp>

#include <picongpu/simulation_defines.hpp>

#include <functional>
#include <numeric>

#include <catch2/catch_approx.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <picongpu/particles/Particles.hpp>

using ::picongpu::floatD_X;
using ::pmacc::localCellIdx;
using ::pmacc::multiMask;
using position = ::picongpu::position<>;
using ::picongpu::simDim;
using ::picongpu::particles::moveParticle;

const auto superCellSize = ::picongpu::SuperCellSize::toRT();
constexpr const auto superCellVolume = ::pmacc::math::CT::volume<::picongpu::SuperCellSize>::type::value;

template<typename T>
static bool isApproxEqual(T const& a, T const& b)
{
    return std::transform_reduce(
        &a[0], // should be std::cbegin(a),
        (&a[simDim - 1]) + 1, // should be std::cend(a),
        &b[0], // should be std::cbegin(b),
        true,
        std::logical_and{},
        [](auto const& lhs, auto const& rhs)
        { return lhs == Catch::Approx(rhs).margin(std::numeric_limits<typename T::type>::epsilon()); });
}


/** A tiny stub of a particle implementing the interface expected by moveParticle()
 */
struct ParticleStub
{
    floatD_X pos = floatD_X::create(0.);
    int localCellIdxValue = 0;
    int multiMaskValue = 1;

    int& operator[](localCellIdx const& index)
    {
        return localCellIdxValue;
    }

    floatD_X& operator[](position const& index)
    {
        return pos;
    }

    int& operator[](multiMask const& index)
    {
        return multiMaskValue;
    }

    bool operator==(ParticleStub const& other) const
    {
        return isApproxEqual(pos, other.pos) and localCellIdxValue == other.localCellIdxValue
            and multiMaskValue == other.multiMaskValue;
    }

    std::string toString() const
    {
        std::stringstream ss;
        ss << "pos: " << pos << ", localCellIdxValue: " << localCellIdxValue << ", multiMaskValue: " << multiMaskValue;
        return ss.str();
    }
};

TEST_CASE("unit::moveParticle", "[moveParticle test]")
{
    ParticleStub particle;
    auto expectedParticle = particle;
    floatD_X newPos = particle.pos;

    SECTION("does nothing for unchanged position")
    {
        REQUIRE(newPos == particle.pos);

        moveParticle(particle, newPos);

        CHECK(particle == expectedParticle);
    }

    SECTION("moves trivially inside cell")
    {
        auto i = GENERATE(range(0u, simDim));

        newPos[i] = .42;
        expectedParticle.pos = newPos;

        REQUIRE(newPos != particle.pos);

        moveParticle(particle, newPos);

        CHECK(particle == expectedParticle);
    }


    SECTION("moves outside of cell in positive direction")
    {
        auto i = GENERATE(range(0u, simDim));

        const std::array<int, 3> neighbouringLocalCellIdxPositive{1, 8, 64};
        newPos[i] = 1.1;
        expectedParticle.pos[i] = .1;
        expectedParticle.localCellIdxValue = neighbouringLocalCellIdxPositive[i];

        moveParticle(particle, newPos);

        CHECK(particle == expectedParticle);
    }

    SECTION("moves outside of cell in negative direction")
    {
        auto i = GENERATE(range(0u, simDim));

        const auto lastCell = superCellVolume - 1;
        // last number can be hard-coded because if this is used, we know we're 3D:
        const std::array<int, 3> neighbouringLocalCellIdxNegative{
            lastCell - 1,
            lastCell - superCellSize[simDim - 2],
            191};

        // make sure that we don't cross supercell borders just here:
        particle.localCellIdxValue = lastCell;

        newPos[i] = -.3;
        expectedParticle.pos[i] = .7;
        expectedParticle.localCellIdxValue = neighbouringLocalCellIdxNegative[i];

        moveParticle(particle, newPos);

        CHECK(particle == expectedParticle);
    }

    SECTION("moves diagonally out of cell")
    {
        newPos[0] = 1.4;
        newPos[1] = 1.6;
        expectedParticle.pos[0] = .4;
        expectedParticle.pos[1] = .6;
        expectedParticle.localCellIdxValue = 9;

        moveParticle(particle, newPos);
        std::cout << particle.toString() << "\n";

        CHECK(particle == expectedParticle);
    }
}
