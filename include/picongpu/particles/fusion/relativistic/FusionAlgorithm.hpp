/* Copyright 2015-2024 Rene Widera, Pawel Ordyna
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

#include <pmacc/random/distributions/Uniform.hpp>

#include <cmath>
#include <cstdio>
#include <type_traits>
#include <utility>

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
                    using namespace pmacc;
                    using namespace picongpu::particles::fusion::precision;
                    constexpr float_COLL c = static_cast<float_COLL>(sim.pic.getSpeedOfLight());

                    /* Calculate @f[ \gamma^* m @f]
                     *
                     * Returns particle mass times its Lorentz factor in the COM frame.
                     *
                     * @param labMomentum particle momentum in the labFrame
                     * @param mass particle mass
                     * @param gamma
                     * @param gammaComs Lorentz factor of the COM frame in the lab frame
                     * @param comsVelocity COM system velocity in the lab frame
                     */
                    DINLINE float_COLL coeff(
                        float3_COLL labMomentum,
                        float_COLL mass,
                        float_COLL gamma,
                        float_COLL gammaComs,
                        float3_COLL comsVelocity)
                    {
                        float3_COLL labVelocity = labMomentum / gamma / mass;
                        float_COLL dot = pmacc::math::dot(comsVelocity, labVelocity);

                        float_COLL val = gammaComs * gamma - dot * (gammaComs * gamma / (c * c));
                        val *= mass;
                        return val;
                    }

                    /* Convert momentum from the lab frame into the COM frame.
                     *
                     * @param labMomentum momentum in the lab frame
                     * @param mass particle mass
                     * @param gamma particle Lorentz factor in the lab frame
                     * @param gammaComs @f[ \gamma_C @f] Lorentz factor of the COM frame in the lab frame
                     * @param factorA @f[ \frac{\gamma_C - 1}{v_C^2} @f]
                     * @param comsVelocity @f[ v_C @f] COM system velocity in the lab frame
                     */
                    DINLINE float3_COLL labToComs(
                        float3_COLL labMomentum,
                        float_COLL mass,
                        float_COLL gamma,
                        float_COLL gammaComs,
                        float_COLL factorA,
                        float3_COLL comsVelocity)
                    {
                        float3_COLL labVelocity = labMomentum / gamma / mass;
                        float_COLL dot = pmacc::math::dot(comsVelocity, labVelocity);
                        float_COLL factor = (factorA * dot - gammaComs);
                        factor *= mass * gamma;
                        float3_COLL diff = factor * comsVelocity;
                        return labMomentum + diff;
                    }

                    /* Calculate relative velocity in the COM system
                     *
                     * @param comsMementumMag0 1st particle momentum (in the COM system) magnitude
                     * @param mass0 1st particle mass
                     * @param mass1 2nd particle mass
                     * @param gamma0 1st particle Lorentz factor
                     * @param gamma1 2nd particle Lorentz factor
                     * @param gammaComs Lorentz factor of the COM frame in the lab frame
                     */
                    DINLINE float_COLL calcRelativeComsVelocity(
                        float_COLL comsMomentumMag0,
                        float_COLL mass0,
                        float_COLL mass1,
                        float_COLL gamma0,
                        float_COLL gamma1,
                        float_COLL coeff0,
                        float_COLL coeff1,
                        float_COLL gammaComs)
                    {
                        float_COLL val = (mass0 * gamma0 + mass1 * gamma1) * comsMomentumMag0;
                        val = val / (coeff0 * coeff1 * gammaComs); // TODO: null division?
                        // TODO:
                        // this is actually not the relative velocity! since we are missing (1 + v1v2/c^2).
                        return val;
                    }

                    /* Convert momentum from the COM frame into the lab frame.
                     *
                     * @param labMomentum momentum in the COM frame
                     * @param mass particle mass
                     * @param gamma particle Lorentz factor in the lab frame
                     * @param gammaComs @f[ \gamma_C @f] Lorentz factor of the COM frame in the lab frame
                     * @param factorA @f[ \frac{\gamma_C - 1}{v_C^2} @f]
                     * @param comsVelocity @f[ v_C @f] COM system velocity in the lab frame
                     */
                    DINLINE float3_COLL comsToLab(
                        float3_COLL comsMomentum,
                        float_COLL mass,
                        float_COLL coeff,
                        float_COLL gammaComs,
                        float_COLL factorA,
                        float3_COLL comsVelocity)
                    {
                        // (13) in [Perez 2012]
                        float_COLL dot = pmacc::math::dot(comsVelocity, comsMomentum);
                        float_COLL factor = (factorA * dot + coeff * gammaComs);
                        float3_COLL diff = factor * comsVelocity;

                        return comsMomentum + diff;
                    }


                    /* Calculate the momentum after the collision in the COM frame
                     *
                     * @param p momentum in the COM frame
                     * @param cosXi cosine of the scattering angle
                     * @param phi azimuthal scattering angle from [0, 2pi]
                     */
                    DINLINE float3_COLL
                    calcFinalComsMomentum(float3_COLL const p, float_COLL const cosXi, float_COLL const phi)
                    {
                        float_COLL sinPhi, cosPhi;
                        pmacc::math::sincos(phi, sinPhi, cosPhi);
                        float_COLL sinXi = math::sqrt(1.0_COLL - cosXi * cosXi);

                        // (12) in [Perez 2012]
                        float3_COLL finalVec;
                        float_COLL const pNorm2 = math::sqrt(pmacc::math::l2norm2(p));
                        float_COLL const pPerp = math::sqrt(p.x() * p.x() + p.y() * p.y());
                        // TODO chose a better limit?
                        // limit px->0 py=0. this also covers the pPerp = pAbs = 0 case. An alternative would
                        // be to let the momentum unchanged in that case.
                        if(pPerp <= math::max(std::numeric_limits<float_COLL>::epsilon(), 1.0e-10_COLL) * pNorm2)
                        {
                            finalVec[0] = pNorm2 * sinXi * cosPhi;
                            finalVec[1] = pNorm2 * sinXi * sinPhi;
                            finalVec[2] = pNorm2 * cosXi;
                        }
                        else // normal case
                        {
                            finalVec[0] = (p.x() * p.z() * sinXi * cosPhi - p.y() * pNorm2 * sinXi * sinPhi) / pPerp
                                          + p.x() * cosXi;
                            finalVec[1] = (p.y() * p.z() * sinXi * cosPhi + p.x() * pNorm2 * sinXi * sinPhi) / pPerp
                                          + p.y() * cosXi;
                            finalVec[2] = -1.0_COLL * pPerp * sinXi * cosPhi + p.z() * cosXi;
                        }
                        return finalVec;
                    }

                    //! Stores some precalculated values used in the collision algorithm
                    struct Variables
                    {
                        PMACC_ALIGN(normalizedWeight0, float_COLL);
                        PMACC_ALIGN(normalizedWeight1, float_COLL);
                        PMACC_ALIGN(labMomentum0, float3_COLL);
                        PMACC_ALIGN(labMomentum1, float3_COLL);
                        PMACC_ALIGN(mass0, float_COLL);
                        PMACC_ALIGN(mass1, float_COLL);
                        PMACC_ALIGN(charge0, float_COLL);
                        PMACC_ALIGN(charge1, float_COLL);
                        PMACC_ALIGN(gamma0, float_COLL);
                        PMACC_ALIGN(gamma1, float_COLL);
                        PMACC_ALIGN(comsVelocity, float3_COLL);

                        PMACC_ALIGN(comsMomentum0, float3_COLL);
                        PMACC_ALIGN(comsMomentum0Norm2, float_COLL);
                        PMACC_ALIGN(gammaComs, float_COLL);
                        PMACC_ALIGN(factorA, float_COLL);
                        PMACC_ALIGN(coeff0, float_COLL);
                        PMACC_ALIGN(coeff1, float_COLL);

                        template<typename T_Par0, typename T_Par1>
                        DINLINE Variables(T_Par0 const& par0, T_Par1 const& par1)
                            : normalizedWeight0(precisionCast<float_COLL>(par0[weighting_]) / WEIGHT_NORM_COLL)
                            , normalizedWeight1(precisionCast<float_COLL>(par1[weighting_]) / WEIGHT_NORM_COLL)
                            , labMomentum0(precisionCast<float_COLL>(par0[momentum_]) / normalizedWeight0)
                            , labMomentum1(precisionCast<float_COLL>(par1[momentum_]) / normalizedWeight1)
                            , mass0(
                                  precisionCast<float_COLL>(
                                      picongpu::traits::attribute::getMass(WEIGHT_NORM_COLL, par0)))
                            , mass1(
                                  precisionCast<float_COLL>(
                                      picongpu::traits::attribute::getMass(WEIGHT_NORM_COLL, par1)))
                            , charge0(
                                  precisionCast<float_COLL>(
                                      picongpu::traits::attribute::getCharge(WEIGHT_NORM_COLL, par0)))
                            , charge1(
                                  precisionCast<float_COLL>(
                                      picongpu::traits::attribute::getCharge(WEIGHT_NORM_COLL, par1)))
                            , gamma0(picongpu::gamma<float_COLL>(labMomentum0, mass0))
                            , gamma1(picongpu::gamma<float_COLL>(labMomentum1, mass1))
                            , comsVelocity((labMomentum0 + labMomentum1) / (mass0 * gamma0 + mass1 * gamma1))
                        {
                            float_COLL const comsVelocityNorm2 = pmacc::math::l2norm2(comsVelocity);

                            if(comsVelocityNorm2 != 0.0_COLL)
                            {
                                float_COLL const comsVelocityAbs = math::sqrt(comsVelocityNorm2);
                                // written as (1-v)(1+v) rather than (1-v^2) for better performance when v close to
                                // c
                                gammaComs = 1.0_COLL
                                            / math::sqrt(
                                                (1.0_COLL - comsVelocityAbs / c) * (1.0_COLL + comsVelocityAbs / c));
                                // used later for comsToLab:
                                factorA = (gammaComs - 1.0_COLL) / comsVelocityNorm2;

                                // Stared gamma times mass, from [Perez 2012].
                                coeff0 = coeff(labMomentum0, mass0, gamma0, gammaComs, comsVelocity);
                                // gamma^* . mass
                                coeff1 = coeff(labMomentum1, mass1, gamma1, gammaComs, comsVelocity);
                                // (2) in [Perez 2012]
                                comsMomentum0
                                    = labToComs(labMomentum0, mass0, gamma0, gammaComs, factorA, comsVelocity);
                            }
                            else
                            {
                                // Lab frame is the same as the COMS frame
                                gammaComs = 1.0_COLL;
                                // used later for comsToLab:
                                // lim v_coms-->0 for (gamma_coms -1 / v_coms^2) is 1/(2c^2)
                                factorA = 1.0_COLL / (2.0_COLL * c * c);
                                // Stared gamma times mass, from [Perez 2012].
                                coeff0 = mass0 * gamma0;
                                // gamma^* . mass
                                coeff1 = mass1 * gamma1;
                                comsMomentum0 = labMomentum0;
                            }
                            comsMomentum0Norm2 = pmacc::math::l2norm2(comsMomentum0);
                        }
                    };


                    template<typename T_CrossSection, bool ifDebug>
                    struct FusionAlg
                    {
                        HDINLINE FusionAlg(){};
                        PMACC_ALIGN(crossSection, T_CrossSection);


                    public:
                        template<typename T_Worker, typename T_Par0, typename T_Par1, typename T_RngHandle>
                        DINLINE void fuse(T_Worker const& worker, T_Par0 par0, T_Par1 par1, uint32_t duplicationCorrection, float_X probabilityFactor, float3_X &mom1, float3_X &mom2, T_RngHandle& rngHandle){
                            // if((par0[momentum_] == float3_X{0.0_X, 0.0_X, 0.0_X})
                            //    && (par1[momentum_] == float3_X{0.0_X, 0.0_X, 0.0_X}))
                                // return;
                            // Variables const v{par0, par1};
                            // if(v.comsMomentum0Norm2 == 0.0_COLL)
                                // return;
                                
                            // Get a random float value from 0,1
                            using UniformFloat = pmacc::random::distributions::Uniform<
                                pmacc::random::distributions::uniform::ExcludeOne<float_COLL>::Reduced>;
                            auto rng = rngHandle.template applyDistribution<UniformFloat>();
                            float_COLL rngValue1 = rng(worker);
                            float_COLL rngValue2 = rng(worker);
                            float_COLL rngValue3 = rng(worker);

                            float_X someEnergy = math::dot(par0[momentum_], par0[momentum_]);
                            float_X test_sigma = crossSection(someEnergy);
                            float_X P = 0.01_X * probabilityFactor;
                            test_sigma *= (rngValue1 < P);

                            // rngValues 2 and 3 are used to generate the scattering angle
                            float_COLL x1 = 2.0_COLL * rngValue2 - 1.0_COLL; // [-1,1]
                            float_COLL x2 = 2.0_COLL * rngValue3 - 1.0_COLL; // [-1,1]
                            while(x1 * x1 + x2 * x2 > 1.0_COLL)
                            {
                                // rejection sampling
                                rngValue2 = rng(worker);
                                rngValue3 = rng(worker);
                                x1 = 2.0_COLL * rngValue2 - 1.0_COLL; // [-1,1]
                                x2 = 2.0_COLL * rngValue3 - 1.0_COLL; // [-1,1]
                            }
                            float_COLL s = math::sqrt(1.0_COLL - x1 * x1 - x2 * x2);
                            float_COLL x = 2.0_COLL * x1 * s;
                            float_COLL y = 2.0_COLL * x2 * s;
                            float_COLL z = 1.0_COLL - 2*(x * x + y * y);

                            float3_X dir = float3_X(x, y, z);
                            mom1 = dir*test_sigma;
                            mom2 = -dir*test_sigma;
                        }
                    };
                } // namespace acc
            } // namespace relativistic
        } // namespace collision
    } // namespace particles
} // namespace picongpu
