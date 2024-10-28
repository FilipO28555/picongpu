/* Copyright 2024 Brian Marre
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

/** @file FillIPDSumFields ionization potential depression(IPD) sub-stage or an ion species
 *
 * implements filling of IPD sum fields from reduction of all macro particles of the specified **ion** species
 */

#pragma once

#include "picongpu/defines.hpp"
#include "picongpu/particles/atomicPhysics/ionizationPotentialDepression/kernel/FillIPDSumFields_Ion.kernel"
#include "picongpu/particles/atomicPhysics/ionizationPotentialDepression/stage/FillIPDSumFields_Ion.def"
#include "picongpu/particles/atomicPhysics/localHelperFields/TimeRemainingField.hpp"
#include "picongpu/particles/param.hpp"

#include <pmacc/particles/meta/FindByNameOrType.hpp>

#include <string>

namespace picongpu::particles::atomicPhysics::ionizationPotentialDepression::stage
{
    //! short hand for IPD namespace
    namespace s_IPD = picongpu::particles::atomicPhysics::ionizationPotentialDepression;

    //! call of kernel for every superCell
    template<typename T_IonSpecies, typename T_TemperatureFunctional>
    HINLINE void FillIPDSumFields_Ion<T_IonSpecies, T_TemperatureFunctional>::operator()(
        picongpu::MappingDesc const mappingDesc) const
    {
        // might be alias, from here on out no more
        //! resolved type of alias T_ParticleSpecies
        using IonSpecies = pmacc::particles::meta::FindByNameOrType_t<VectorAllSpecies, T_IonSpecies>;

        static_assert(
            pmacc::traits::HasIdentifiers<typename IonSpecies::FrameType, MakeSeq_t<weighting, momentum>>::type::value,
            "atomic physics: species is missing one of the following attributes: weighting, momentum");

        // full local domain, no guards
        pmacc::AreaMapping<CORE + BORDER, MappingDesc> mapper(mappingDesc);
        pmacc::DataConnector& dc = pmacc::Environment<>::get().DataConnector();

        auto& timeRemainingField = *dc.get<
            picongpu::particles::atomicPhysics::localHelperFields::TimeRemainingField<picongpu::MappingDesc>>(
            "TimeRemainingField");

        // pointer to memory, we will only work on device, no sync required
        // init pointer to particles and localSumFields
        auto& ions = *dc.get<IonSpecies>(IonSpecies::FrameType::getName());

        auto& localSumWeightAllField
            = *dc.get<s_IPD::localHelperFields::SumWeightAllField<picongpu::MappingDesc>>("SumWeightAllField");
        auto& localSumTemperatureFunctionalField
            = *dc.get<s_IPD::localHelperFields::SumTemperatureFunctionalField<picongpu::MappingDesc>>(
                "SumTemperatureFunctionalField");

        auto& localSumChargeNumberIonsField
            = *dc.get<s_IPD::localHelperFields::SumChargeNumberIonsField<picongpu::MappingDesc>>(
                "SumChargeNumberIonsField");
        auto& localSumChargeNumberSquaredIonsField
            = *dc.get<s_IPD::localHelperFields::SumChargeNumberSquaredIonsField<picongpu::MappingDesc>>(
                "SumChargeNumberSquaredIonsField");

        // macro for call of kernel on every superCell, see pull request #4321
        PMACC_LOCKSTEP_KERNEL(s_IPD::kernel::FillIPDSumFieldsKernel_Ion<T_TemperatureFunctional>())
            .config(mapper.getGridDim(), ions)(
                mapper,
                timeRemainingField.getDeviceDataBox(),
                ions.getDeviceParticlesBox(),
                localSumWeightAllField.getDeviceDataBox(),
                localSumTemperatureFunctionalField.getDeviceDataBox(),
                localSumChargeNumberIonsField.getDeviceDataBox(),
                localSumChargeNumberSquaredIonsField.getDeviceDataBox());
    }
} // namespace picongpu::particles::atomicPhysics::ionizationPotentialDepression::stage
