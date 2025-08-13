/* Copyright 2025-2025 Filip Optolowicz
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PIConGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "picongpu/defines.hpp"
#include "picongpu/particles/fusion/param.hpp"


#if (BOOST_LANG_CUDA || BOOST_COMP_HIP)
#    include <mallocMC/mallocMC.hpp>
#endif

namespace picongpu::particles::fusion
{
    namespace detail
    {
        struct CreationFusion{

            template<
                typename T_Worker,
                typename T_ParAccessor0,
                typename T_ParAccessor1,
                typename T_ParAccessor2,
                typename T_ParAccessor3
                >
            DINLINE decltype(auto) createParticles(
                T_Worker const& worker,
                    IdGenerator& idGen,
                T_ParAccessor0 const& r1,
                T_ParAccessor1 const& r2,
                float_X const& productWeighting,
                float3_X const& mom1,
                float3_X const& mom2,
                T_ParAccessor2& p1r1, // product 1 at pos 1
                T_ParAccessor2& p1r2, // product 1 at pos 2
                T_ParAccessor3& p2r1, // product 2 at pos 1
                T_ParAccessor3& p2r2 // product 2 at pos 2
                ) const
            {
                /** for not mixing operations::assign up with the nvidia functor assign */
                namespace partOp = pmacc::particles::operations;
                /** Set all product particles multimask to 0 - we will set it to 1 at the end of fusion stage */
                p1r1[multiMask_] = 1u;
                p1r2[multiMask_] = 1u;
                p2r1[multiMask_] = 1u;
                p2r2[multiMask_] = 1u;

                /** each thread initializes a clone of the parent particle but leaving out
                 * some attributes:
                 * - multiMask: reading from global memory takes longer than just setting it again explicitly
                 * - momentum: we have the momentum
                 */
                auto targetClone2 = partOp::deselect<pmacc::mp_list<multiMask, momentum, weighting>>(p1r1);
                auto targetClone3 = partOp::deselect<pmacc::mp_list<multiMask, momentum, weighting>>(p1r2);
                auto targetClone4 = partOp::deselect<pmacc::mp_list<multiMask, momentum, weighting>>(p2r1);
                auto targetClone5 = partOp::deselect<pmacc::mp_list<multiMask, momentum, weighting>>(p2r2);

                targetClone2.derive(worker, idGen, r1);
                targetClone3.derive(worker, idGen, r2);

                targetClone4.derive(worker, idGen, r1);
                targetClone5.derive(worker, idGen, r2);

                p1r1[momentum_] = mom1;
                p1r1[weighting_] = productWeighting/2.;
                p1r2[momentum_] = mom1;
                p1r2[weighting_] = productWeighting/2.;
                p2r1[momentum_] = mom2;
                p2r1[weighting_] = productWeighting/2.;
                p2r2[momentum_] = mom2;
                p2r2[weighting_] = productWeighting/2.;
            }
        };

    } // namespace detail
} // namespace picongpu::particles::fusion

