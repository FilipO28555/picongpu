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
                    constexpr float_COLL c2 = c * c;
                    constexpr float_COLL c3 = c2 * c;
                    constexpr float_COLL c4 = c3 * c;

                    //! Stores some precalculated values used in the collision algorithm
                    struct Variables
                    {
                        PMACC_ALIGN(labMomentum0, float3_COLL);
                        PMACC_ALIGN(labMomentum1, float3_COLL);
                        PMACC_ALIGN(mass0, float_COLL);
                        PMACC_ALIGN(mass1, float_COLL);
                        PMACC_ALIGN(massProduct1, float_COLL);
                        PMACC_ALIGN(massProduct2, float_COLL);
                        PMACC_ALIGN(gamma0, float_COLL);
                        PMACC_ALIGN(gamma1, float_COLL);
                        PMACC_ALIGN(V_cm, float3_COLL);
                        
                        PMACC_ALIGN(gamma_cm, float_COLL);
                        PMACC_ALIGN(gamma0_cm, float_COLL);
                        PMACC_ALIGN(gamma1_cm, float_COLL);
                        PMACC_ALIGN(factorA, float_COLL);
                        PMACC_ALIGN(E_r, float_COLL);
                        PMACC_ALIGN(V_rel, float_COLL);

                        template<typename T_Par0, typename T_Par1>
                        DINLINE Variables(T_Par0 const& par0, T_Par1 const& par1, float_X const weight0, float_X const weight1)
                            : 
                             labMomentum0(precisionCast<float_COLL>(par0[momentum_]) / weight0 )
                            , labMomentum1(precisionCast<float_COLL>(par1[momentum_]) / weight1 )
                            , mass0(
                                // previously WEIGHT_NORM_COLL was used here 
                                  precisionCast<float_COLL>(
                                      picongpu::traits::attribute::getMass(1, par0))) // check if the weighting is right
                            , mass1(
                                  precisionCast<float_COLL>(
                                      picongpu::traits::attribute::getMass(1, par1)))
                            , gamma0(picongpu::gamma<float_COLL>(labMomentum0, mass0))
                            , gamma1(picongpu::gamma<float_COLL>(labMomentum1, mass1))
                            , V_cm((labMomentum0 + labMomentum1) / (mass0 * gamma0 + mass1 * gamma1))
                        {
                            float3_COLL const u0 = labMomentum0/mass0;
                            float3_COLL const u1 = labMomentum1/mass1;
                            // calculate CM velocity and gamma factor
                            float_COLL const V_cm_mag = pmacc::math::l2norm(V_cm);
                            gamma_cm = 1.0_COLL / math::sqrt(
                                        (1.0_COLL - V_cm_mag / c) * (1.0_COLL + V_cm_mag / c));

                            gamma0_cm = gamma_cm*(gamma0-pmacc::math::dot(V_cm, u0)/c2);
                            gamma1_cm = gamma_cm*(gamma1-pmacc::math::dot(V_cm, u1)/c2);

                            factorA = (gamma_cm - 1.0_COLL) / (V_cm_mag*V_cm_mag);
                            using pmacc::math::dot;
                            // Boost the momenta into the CM frame
                            float3_COLL const u0_cm = u0 + (dot(V_cm, u0) * factorA - gamma_cm*gamma0)*V_cm;
                            float3_COLL const u1_cm = u1 + (dot(V_cm, u1) * factorA - gamma_cm*gamma1)*V_cm;
                            float3_COLL const V0_cm = u0_cm/gamma0_cm;
                            float3_COLL const V1_cm = u1_cm/gamma1_cm;
                            // calculate relative velocity in the CM frame
                            float_COLL const m_r = mass0 * mass1 / (mass0 + mass1);
                            V_rel = pmacc::math::l2norm(
                                (V0_cm - V1_cm)/(1.0_COLL - pmacc::math::dot(V0_cm, V1_cm)/c/c));

                            float_COLL const gamma_r = 1.0_COLL / math::sqrt(
                                        (1.0_COLL - V_rel / c) * (1.0_COLL + V_rel / c));
                            // calculate the relative energy in the CM frame
                            E_r = m_r * c * c * (gamma_r-1.0_COLL);
                        }

                        // void P(float3_COLL const& dir)
                        // {
                        //     float_COLL const Q = 0.0_COLL;
                        //     float_COLL const mP0 = 1.0_COLL;
                        //     float_COLL const mP1 = 2.0_COLL;
                        //     // energy of product 0
                        //     float_COLL const Ep0 = (E_r + Q + mP1 * c * c) / (E_r + Q + (mP0+mP1)*c*c) / 2.0_COLL;
                        //     float_COLL const p0mag = math::sqrt((Ep0 + mP0*c*c) * (Ep0 + mP0*c*c) - mP0 * mP0 * c * c)/c;
 
                        //     // Momentum vectors in the CM frame
                        //     float3_COLL const p0_cm = p0mag * dir;
                        //     float3_COLL const p1_cm = -p0_cm;


                        //     // --- Inverse Lorentz Boost back to Lab Frame ---
                        //     // We apply the reverse transformation using the pre-calculated V_cm and gamma_cm.
                        //     // The structure is similar to the forward boost, but the sign of the velocity-dependent term is flipped.

                        //     // For Product 0:
                        //     float_COLL  const gamma_p0_cm = E_p0_tot_cm / (mP0 * c2);
                        //     float3_COLL const u0_cm = p0_cm / mP0;
                        //     float3_COLL const u0_lab = u0_cm + (math::dot(V_cm, u0_cm) * factorA + gamma_cm * gamma_p0_cm) * V_cm;
                        //     labMomentum0 = u0_lab * mP0;

                        //     // For Product 1:
                        //     float_COLL  const gamma_p1_cm = E_p1_tot_cm / (mP1 * c2);
                        //     float3_COLL const u1_cm = p1_cm / mP1;
                        //     float3_COLL const u1_lab = u1_cm + (math::dot(V_cm, u1_cm) * factorA + gamma_cm * gamma_p1_cm) * V_cm;
                        //     labMomentum1 = u1_lab * mP1;
                        // }
                        template<typename T_Product0Box, typename T_Product1Box>
                        DINLINE void P_gemini(float3_COLL const& dir)
                        {
                            // --- Define reaction properties ---
                            // Q-value of the reaction (energy released)
                            // float_COLL const Q = 0.0_COLL;
                            // Rest masses of the two product particles - masses of one real particle
                            float_COLL const mP0 = picongpu::traits::frame::getMass<typename T_Product0Box::FrameType>();
                            float_COLL const mP1 = picongpu::traits::frame::getMass<typename T_Product1Box::FrameType>();

                            // --- Corrected Energy Calculation (CM Frame) ---
                            // The total energy in the CM frame is the sum of the initial particles' kinetic and rest mass energies.
                            // E_r is the kinetic energy of the system, so E_cm_tot = E_r + (mass0 + mass1) * c^2
                            // which also equals E_r + Q + (mP0 + mP1) * c^2
                            float_COLL const E_cm_tot = E_r + (mass0 + mass1) * c2;

                            // Correct relativistic formula for the TOTAL energy of product 0 in the CM frame
                            float_COLL const E_p0_tot_cm = (E_cm_tot * E_cm_tot + (mP0 * mP0 - mP1 * mP1) * c4) / (2.0_COLL * E_cm_tot);
                            // TOTAL energy of product 1
                            float_COLL const E_p1_tot_cm = E_cm_tot - E_p0_tot_cm;


                            // --- Calculate Product Momenta (CM Frame) ---
                            // Magnitude of momentum for product 0, from p = sqrt(E_tot^2 - (mc^2)^2) / c
                            float_COLL const p0_mag_cm = math::sqrt(E_p0_tot_cm * E_p0_tot_cm - mP0 * mP0 * c4) / c;

                            // Momentum vectors in the CM frame
                            float3_COLL const p0_cm = p0_mag_cm * dir;
                            float3_COLL const p1_cm = -p0_cm;


                            // --- Inverse Lorentz Boost back to Lab Frame ---
                            // We apply the reverse transformation using the pre-calculated V_cm and gamma_cm.
                            // The structure is similar to the forward boost, but the sign of the velocity-dependent term is flipped.

                            // For Product 0:
                            float_COLL  const gamma_p0_cm = E_p0_tot_cm / (mP0 * c2);
                            float3_COLL const u0_cm = p0_cm / mP0;
                            float3_COLL const u0_lab = u0_cm + (math::dot(V_cm, u0_cm) * factorA + gamma_cm * gamma_p0_cm) * V_cm;
                            labMomentum0 = u0_lab * mP0;

                            // For Product 1:
                            float_COLL  const gamma_p1_cm = E_p1_tot_cm / (mP1 * c2);
                            float3_COLL const u1_cm = p1_cm / mP1;
                            float3_COLL const u1_lab = u1_cm + (math::dot(V_cm, u1_cm) * factorA + gamma_cm * gamma_p1_cm) * V_cm;
                            labMomentum1 = u1_lab * mP1;
                        }
                        DINLINE float3_X P0() const
                        {
                            return precisionCast<float_X>(labMomentum0);
                        }

                        DINLINE float3_X P1() const
                        {
                            return precisionCast<float_X>(labMomentum1);
                        }

                    };


                    template<typename T_CrossSection>
                    struct FusionAlg
                    {
                        HDINLINE FusionAlg(){};
                        PMACC_ALIGN(crossSection, T_CrossSection);


                    public:
                        template<typename T_Product0Box, typename T_Product1Box, typename T_Worker, typename T_Par0, typename T_Par1, typename T_RngHandle>
                        DINLINE void fuse(T_Worker const& worker, T_Par0 par0, T_Par1 par1, float_X weightingR1, float_X weightingR2, float_X probabilityFactor, float3_X &mom0, float3_X &mom1, T_RngHandle& rngHandle){
                            if((par0[momentum_] == float3_X{0.0_X, 0.0_X, 0.0_X})
                               && (par1[momentum_] == float3_X{0.0_X, 0.0_X, 0.0_X}))
                                return;
                            // calculate boost and relative energy
                            Variables v{par0, par1, weightingR1, weightingR2};

                            // Convert energy from PIC units to keV
                            constexpr float_COLL picEnergy_to_Joule = sim.unit.energy();
                            constexpr float_COLL joule_to_eV = 1.0 / sim.si.get_eV();
                            constexpr float_COLL eV_to_keV = 1e-3;
                            constexpr float_COLL convToKeV = picEnergy_to_Joule * joule_to_eV * eV_to_keV;

                            // Convert cross section from millibarns to PIC area units
                            constexpr float_COLL millibarn_to_m2 = 1e-31;  // 1 millibarn = 1e-31 m²
                            constexpr float_COLL picLength_to_m = sim.unit.length();
                            constexpr float_COLL m2_to_picArea = 1.0 / (picLength_to_m * picLength_to_m);
                            constexpr float_COLL millibarn_to_picArea = millibarn_to_m2 * m2_to_picArea;

                            // Apply conversions
                            float_X sigma_picArea = crossSection(v.E_r * convToKeV) * millibarn_to_picArea;
                            float_X P = probabilityFactor * sigma_picArea * v.V_rel * v.gamma_cm;

                            // Get a random float value from 0,1
                            using UniformFloat = pmacc::random::distributions::Uniform<
                                pmacc::random::distributions::uniform::ExcludeOne<float_COLL>::Reduced>;
                            auto rng = rngHandle.template applyDistribution<UniformFloat>();
                            float_COLL rngValue1 = rng(worker);

                            // print with probability 1e-8
                            if (worker.workerIdx() == 0 && rng(worker) < 1e-9)
                            printf("Worker %d,millibarn_to_picArea: %e, sigma_picArea: %e, probabilityFactor: %e, v.V_rel: %e, v.gamma_cm: %e, P: %e\n",
                                   worker.workerIdx(), millibarn_to_picArea, sigma_picArea, probabilityFactor, v.V_rel, v.gamma_cm, P);

                            if(rngValue1 < P){

                                float_COLL rngValue2 = rng(worker);
                                float_COLL rngValue3 = rng(worker);

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
                                float_COLL z = 1.0_COLL - 2*(x1 * x1 + x2 * x2);


                                // returns momentum of one particle - not multiplied by weighting.
                                // Multiplication by weighting is later in creation of particles
                                float3_COLL const dir = float3_COLL(x, y, z);
                                v.P_gemini<T_Product0Box, T_Product1Box>(dir);
                                mom0 = v.P0();
                                mom1 = v.P1();

                            }
                            else
                            {
                                mom0 = float3_X{0.0_X, 0.0_X, 0.0_X};
                                mom1 = float3_X{0.0_X, 0.0_X, 0.0_X};
                            }
                        }
                    };
                } // namespace acc
            } // namespace relativistic
        } // namespace collision
    } // namespace particles
} // namespace picongpu
