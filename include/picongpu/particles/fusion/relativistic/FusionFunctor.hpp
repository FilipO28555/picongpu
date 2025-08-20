/* Copyright 2022-2024 Rene Widera, Pawel Ordyna
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

#include "picongpu/particles/fusion/kernels.def"
#include "picongpu/particles/fusion/relativistic/FusionAlgorithm.hpp"


#include <string>

namespace picongpu
{
    namespace particles
    {
        namespace fusion
        {
            namespace relativistic
            {
                namespace acc
                {
                    //! Coulomb logarithm functor for a fixed logarithm defined at compile time
                    template<typename T_Param>
                    struct CalcCrossSection
                    {
                        DINLINE float_COLL operator()(float_COLL const& Energy) const
                        {
                            float_COLL const Snum = ((((Energy*T_Param::A5)+T_Param::A4)*Energy
                                + T_Param::A3)*Energy + T_Param::A2)*Energy + T_Param::A1;
                            float_COLL const Sden = (((((Energy*T_Param::B4)+T_Param::B3)*Energy
                                + T_Param::B2)*Energy + T_Param::B1)*Energy + 1._COLL);
                            float_COLL const S = Snum / Sden;
                            float_COLL const Eexp = Energy * math::exp(T_Param::BG/math::sqrt(Energy));
                            return S / Eexp;
                        }
                    };

                } // namespace acc

                template<typename T_Param>
                struct FusionFunctorImpl
                {
                    template<typename T_Species0, typename T_Species1, typename T_Species2, typename T_Species3>
                    struct apply
                    {
                        using type = FusionFunctorImpl<T_Param>;
                    };

                    HINLINE FusionFunctorImpl(uint32_t currentStep) {};

                    using AccFunctorImpl = acc::FusionAlg<acc::CalcCrossSection<T_Param>>;
                    using AccFunctor = fusion::acc::IBinary<AccFunctorImpl>;
                    // define kernel that should be used to call this functor
                    using CallingInterKernel = InterCollision;
                    using CallingIntraKernel = IntraCollision;

                    /** create device manipulator functor
                     *
                     * @param worker lockstep worker
                     * @param offset (in supercells, without any guards) to the origin of the local domain
                     * @param density0 cell density of the 1st species
                     * @param density1 cell density of the 2nd species
                     * @param potentialPartners number of potential collision partners for a macro particle in
                     *   the cell.
                     * @param coulombLog Coulomb logarithm
                     */
                    
                    HDINLINE auto operator()()
                    {
                        using namespace picongpu::particles::collision::precision;
                        return AccFunctor{AccFunctorImpl{}};
                    }
                    
                    //! get the name of the functor
                    HINLINE static std::string getName()
                    {
                        return "FusionFunctor";
                    }
                };
            } // namespace relativistic
        } // namespace collision
    } // namespace particles
} // namespace picongpu
