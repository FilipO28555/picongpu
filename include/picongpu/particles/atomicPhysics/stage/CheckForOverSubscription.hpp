/* Copyright 2023-2024 Brian Marre
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

/** @file check for overSubscription of histogram bins and cells and calculate rejectionProbability for each*/

#pragma once

// need picongpu::atomicPhysics::ElectronHistogram from atomicPhysics.param
#include "picongpu/defines.hpp"
#include "picongpu/fields/FieldE.hpp"
#include "picongpu/particles/atomicPhysics/electronDistribution/LocalHistogramField.hpp"
#include "picongpu/particles/atomicPhysics/kernel/CheckForOverSubscription.kernel"
#include "picongpu/particles/atomicPhysics/localHelperFields/RejectionProbabilityCacheField.hpp"
#include "picongpu/particles/atomicPhysics/localHelperFields/SharedResourcesOverSubscribedField.hpp"
#include "picongpu/particles/atomicPhysics/localHelperFields/TimeRemainingField.hpp"

#include <pmacc/Environment.hpp>
#include <pmacc/mappings/kernel/AreaMapping.hpp>
#include <pmacc/type/Area.hpp>

#include <cstdint>

namespace picongpu::particles::atomicPhysics::stage
{
    /** CheckForAndRejectOversubscription atomic physics sub-stage
     *
     * check each histogram bin for deltaWeight > weight0, if yes mark bin as over subscribed
     *
     * @tparam T_numberAtomicPhysicsIonSpecies specialization template parameter used to prevent compilation of all
     *  atomicPhysics kernels if no atomic physics species is present.
     */
    template<uint32_t T_numberAtomicPhysicsIonSpecies>
    struct CheckForOverSubscription
    {
        //! call of kernel for every superCell
        HINLINE void operator()(picongpu::MappingDesc const mappingDesc) const
        {
            // full local domain, no guards
            pmacc::AreaMapping<CORE + BORDER, MappingDesc> mapper(mappingDesc);
            pmacc::DataConnector& dc = pmacc::Environment<>::get().DataConnector();

            auto& timeRemainingField = *dc.get<
                picongpu::particles::atomicPhysics::localHelperFields::TimeRemainingField<picongpu::MappingDesc>>(
                "TimeRemainingField");

            auto& electronHistogramField
                = *dc.get<picongpu::particles::atomicPhysics::electronDistribution::
                              LocalHistogramField<picongpu::atomicPhysics::ElectronHistogram, picongpu::MappingDesc>>(
                    "Electron_HistogramField");

            auto& sharedResourcesOverSubscribedField
                = *dc.get<picongpu::particles::atomicPhysics::localHelperFields::SharedResourcesOverSubscribedField<
                    picongpu::MappingDesc>>("SharedResourcesOverSubscribedField");

            auto& rejectionProbabilityCacheField
                = *dc.get<picongpu::particles::atomicPhysics::localHelperFields::RejectionProbabilityCacheField<
                    picongpu::MappingDesc>>("RejectionProbabilityCacheField");

            auto& eField = *dc.get<FieldE>(FieldE::getName());

            using FieldEnergyUseCacheField
                = picongpu::particles::atomicPhysics::localHelperFields::FieldEnergyUseCacheField<
                    picongpu::MappingDesc>;
            auto& fieldEnergyUseCacheField = *dc.get<FieldEnergyUseCacheField>("FieldEnergyUseCacheField");

            // macro for call of kernel for every superCell, see pull request #4321
            PMACC_LOCKSTEP_KERNEL(picongpu::particles::atomicPhysics::kernel::CheckForOverSubscriptionKernel<
                                      picongpu::atomicPhysics::ElectronHistogram,
                                      FieldEnergyUseCacheField::ElementType,
                                      T_numberAtomicPhysicsIonSpecies>())
                .template config<picongpu::atomicPhysics::ElectronHistogram::numberBins>(mapper.getGridDim())(
                    mapper,
                    timeRemainingField.getDeviceDataBox(),
                    electronHistogramField.getDeviceDataBox(),
                    eField.getDeviceDataBox(),
                    fieldEnergyUseCacheField.getDeviceDataBox(),
                    sharedResourcesOverSubscribedField.getDeviceDataBox(),
                    rejectionProbabilityCacheField.getDeviceDataBox());

            /// @todo implement photon histogram, Brian Marre, 2023
        }
    };
} // namespace picongpu::particles::atomicPhysics::stage
