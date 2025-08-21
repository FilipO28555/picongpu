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
                    /**
                    * @brief Calculates the total relativistic energy of a particle.
                    *
                    * @tparam T_Float The floating point type for the calculation.
                    * @tparam T_Vec The type of the 3-momentum vector.
                    * @tparam T_Mass The type of the rest mass.
                    * @param momentum The relativistic 3-momentum vector of the particle.
                    * @param mass The rest mass of the particle.
                    * @return The total relativistic energy E.
                    */
                    template<typename T_Float, typename T_Vec, typename T_Mass>
                    DINLINE T_Float energy(T_Vec const& momentum, T_Mass const& mass)
                    {
                        // Using the formula E = sqrt((pc)^2 + (mc^2)^2)
                        // which is E = sqrt(p^2 * c^2 + m^2 * c^4)
                        // where p is the magnitude of the 3-momentum vector.

                        T_Float const p_sq = pmacc::math::l2norm2(momentum);
                        return math::sqrt(p_sq * c2 + mass * mass * c4);
                    }

                    //! Stores some precalculated values used in the collision algorithm
                    struct Variables
                    {
                        PMACC_ALIGN(labMomentum0, float3_COLL);
                        PMACC_ALIGN(labMomentum1, float3_COLL);
                        PMACC_ALIGN(mass0, float_COLL);
                        PMACC_ALIGN(mass1, float_COLL);
                        // PMACC_ALIGN(massProduct1, float_COLL); // These are not used
                        // PMACC_ALIGN(massProduct2, float_COLL); // These are not used
                        PMACC_ALIGN(V_cm, float3_COLL);
                        PMACC_ALIGN(gamma_cm, float_COLL);
                        PMACC_ALIGN(factorA, float_COLL);
                        PMACC_ALIGN(E_cm_tot, float_COLL); // Total energy in CM frame
                        PMACC_ALIGN(V_rel_mag, float_COLL); // Magnitude of relative velocity for cross-section

                        template<typename T_Par0, typename T_Par1>
                        DINLINE Variables(T_Par0 const& par0, T_Par1 const& par1, float_X const weight0, float_X const weight1)
                            : labMomentum0(precisionCast<float_COLL>(par0[momentum_]) / weight0)
                            , labMomentum1(precisionCast<float_COLL>(par1[momentum_]) / weight1)
                            , mass0(precisionCast<float_COLL>(picongpu::traits::attribute::getMass(1, par0)))
                            , mass1(precisionCast<float_COLL>(picongpu::traits::attribute::getMass(1, par1)))
                        {
                            // --- Calculate total 4-momentum in Lab Frame ---
                            float_COLL const E0_lab = energy<float_COLL>(labMomentum0, mass0);
                            float_COLL const E1_lab = energy<float_COLL>(labMomentum1, mass1);
                            float_COLL const E_tot_lab = E0_lab + E1_lab;
                            float3_COLL const p_tot_lab = labMomentum0 + labMomentum1;

                            // --- Calculate Invariant CM Energy and CM Velocity ---
                            // The square of the total CM energy is the invariant mass squared of the system.
                            float_COLL const E_cm_tot_sq = E_tot_lab * E_tot_lab - pmacc::math::l2norm2(p_tot_lab) * c2;
                            E_cm_tot = math::sqrt(E_cm_tot_sq);

                            // CM velocity is needed for the inverse boost.
                            V_cm = p_tot_lab * c2 / E_tot_lab;

                            // --- Calculate parameters for the inverse boost ---
                            float_COLL const V_cm_mag_sq = pmacc::math::l2norm2(V_cm);
                            if (V_cm_mag_sq > 1.e-32_COLL) // Use a safe epsilon
                            {
                                gamma_cm = 1.0_COLL / math::sqrt(1.0_COLL - V_cm_mag_sq / c2);
                                factorA = (gamma_cm - 1.0_COLL) / V_cm_mag_sq;
                            }
                            else
                            {
                                gamma_cm = 1.0_COLL;
                                factorA = 0.5_COLL / c2; // Non-relativistic limit: (gamma-1)/v^2 -> 1/(2c^2)
                            }
                            // instead of if:
                            // // Note: A small V_cm_mag_sq could still lead to gamma_cm being exactly 1.0
                            // // due to precision limits, which is safe in the formula below.
                            // gamma_cm = 1.0_COLL / math::sqrt(1.0_COLL - V_cm_mag_sq / c2);

                            // // This numerically stable, branchless formula avoids the 0/0 problem
                            // // for small V_cm. It is equivalent to (gamma_cm - 1.0) / V_cm_mag_sq.
                            // factorA = gamma_cm * gamma_cm / (c2 * (gamma_cm + 1.0_COLL));



                            // --- Calculate relative velocity for cross-section ---
                            // This is needed for the probability calculation.
                            // s = (p0_4 + p1_4)^2 = E_cm_tot^2
                            float_COLL const s = E_cm_tot_sq;
                            float_COLL const p_cm_mag_sq =
                                (s - (mass0 + mass1) * (mass0 + mass1) * c4) *
                                (s - (mass0 - mass1) * (mass0 - mass1) * c4) / (4.0_COLL * s);
                            V_rel_mag = math::sqrt(p_cm_mag_sq) * s / (E_cm_tot * mass0 * mass1 * c3);
                        }

                        template<typename T_Product0Box, typename T_Product1Box>
                        DINLINE void P_gemini(float3_COLL const& dir)
                        {
                            // --- Define reaction properties ---
                            float_COLL const mP0 = picongpu::traits::frame::getMass<typename T_Product0Box::FrameType>();
                            float_COLL const mP1 = picongpu::traits::frame::getMass<typename T_Product1Box::FrameType>();

                                                        // --- FINAL DEBUG: Manually define product masses in PIC units ---
                            // This test overrides the suspected faulty compile-time mass retrieval.
                            // constexpr float_COLL amu = 1.0_COLL / 5.48579909065e-4_COLL;
                            // constexpr float_COLL mP0 = 1.00866491595_COLL * amu; // Neutron mass
                            // constexpr float_COLL mP1 = 4.00260325413_COLL * amu; // He4 mass

                            // --- DEBUG: Manually inject the Q-value in correct PIC units ---
                            // This is a test to verify the rest of the kinematic calculations.

                            // // 1. Define masses in SI units (kg)
                            // constexpr float_64 U_SI = 1.66053906660e-27; // atomic mass unit in kg
                            // constexpr float_64 C_SI = 299792458.0;       // speed of light in m/s
                            // constexpr float_64 M_D_SI  = 2.01410177812 * U_SI;
                            // constexpr float_64 M_T_SI  = 3.0160492779  * U_SI;
                            // constexpr float_64 M_N_SI  = 1.00866491595 * U_SI;
                            // constexpr float_64 M_HE4_SI= 4.00260325413 * U_SI;

                            // // 2. Calculate Q-value in SI units (Joules)
                            // constexpr float_64 MASS_DEFECT_SI = (M_D_SI + M_T_SI) - (M_N_SI + M_HE4_SI);
                            // constexpr float_64 Q_VALUE_SI = MASS_DEFECT_SI * C_SI * C_SI;

                            // // 3. Convert Q-value from Joules to PIC energy units
                            // // sim.unit.energy() gives the value of 1 PIC energy unit in Joules.
                            // constexpr float_COLL Q_value_pic = Q_VALUE_SI / sim.unit.energy();

                            // // 4. Add the correctly-scaled Q-value to the total CM energy
                            // float_COLL const E_cm_tot_with_Q = E_cm_tot + Q_value_pic;
                            
                            // // Relativistic formula for the TOTAL energy of product 0 in the CM frame
                            // // We use the energy WITH the manually added Q-value.
                            // float_COLL const E_p0_tot_cm = (E_cm_tot_with_Q * E_cm_tot_with_Q + (mP0 * mP0 - mP1 * mP1) * c4) / (2.0_COLL * E_cm_tot_with_Q);
                            // // TOTAL energy of product 1
                            // float_COLL const E_p1_tot_cm = E_cm_tot_with_Q - E_p0_tot_cm;



                            // Relativistic formula for the TOTAL energy of product 0 in the CM frame
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
                            // if((par0[momentum_] == float3_X{0.0_X, 0.0_X, 0.0_X})
                            //    && (par1[momentum_] == float3_X{0.0_X, 0.0_X, 0.0_X}))
                            //     return;
                            // calculate boost and relative energy
                            Variables v{par0, par1, weightingR1, weightingR2};

                            // Convert energy from PIC units to keV
                            // The E_r calculation has been removed. We need to calculate the kinetic energy in the CM frame.
                            float_COLL const E_kin_cm = v.E_cm_tot - (v.mass0 + v.mass1) * c2;
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
                            float_X sigma_picArea = crossSection(E_kin_cm * convToKeV) * millibarn_to_picArea;
                            float_X P = probabilityFactor * sigma_picArea * v.V_rel_mag * v.gamma_cm;

                            // Get a random float value from 0,1
                            using UniformFloat = pmacc::random::distributions::Uniform<
                                pmacc::random::distributions::uniform::ExcludeOne<float_COLL>::Reduced>;
                            auto rng = rngHandle.template applyDistribution<UniformFloat>();
                            float_COLL rngValue1 = rng(worker);
                            if constexpr(alwaysFuseQ) P=1.0_COLL; // always fuse if this is set to true

                            // print with probability 1e-2
                            if (debugFusion || (worker.workerIdx() == 0 && rng(worker) < 1e-2))
                            printf("Worker %d,millibarn_to_picArea: %e, sigma_picArea: %e, probabilityFactor: %e, v.V_rel: %e, v.gamma_cm: %e, P: %e\n",
                                   worker.workerIdx(), millibarn_to_picArea, sigma_picArea, probabilityFactor, v.V_rel_mag, v.gamma_cm, P);

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
