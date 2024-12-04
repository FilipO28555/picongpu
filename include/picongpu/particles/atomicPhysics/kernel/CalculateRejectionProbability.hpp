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

#pragma once

#include "picongpu/defines.hpp"

#include <cstdint>

namespace picongpu::particles::atomicPhysics::kernel
{
    //! methods for calculating the rejection probability
    struct CalculateRejectionProbability
    {
        using VectorIdx = DataSpace<picongpu::simDim>;

        /** store histogram bin rejection probability for specified bin in the passed rejectionProbabilityCacheBin
         *
         * @param binIndex
         * @param histogram
         * @param rejectionProbabilityCacheBin entry for bin will be set to -1 if no rejection necessary or rejection
         *  probability, >= 0, otherwise
         * @param sharedResourcesOverSubscribed previous state of sharedResourcesOverSubscribed
         *
         * @return bin is over subscribed
         */
        template<typename T_Histogram, typename T_RejectionProbabilityCache_Bin>
        HDINLINE static bool ofHistogramBin(
            uint32_t const binIndex,
            T_Histogram const& histogram,
            T_RejectionProbabilityCache_Bin& rejectionProbabilityCacheBin)
        {
            float_X const weight0 = histogram.getBinWeight0(binIndex);
            float_X const deltaWeight = histogram.getBinDeltaWeight(binIndex);

            float_X rejectionProbability = -1._X;
            bool sharedResourcesOverSubscribed = false;
            if(weight0 < deltaWeight)
            {
                // bin is oversubscribed by suggested changes

                // calculate rejection probability
                rejectionProbability = math::max(
                    // proportion of weight we want to reject
                    (deltaWeight - weight0) / deltaWeight,
                    // but at least one average one macro ion should be rejected
                    picongpu::sim.unit.typicalNumParticlesPerMacroParticle() / deltaWeight);

                // set flag that we found at least one over subscribed resource
                sharedResourcesOverSubscribed = true;
            }

            rejectionProbabilityCacheBin.setBin(binIndex, rejectionProbability);
            return sharedResourcesOverSubscribed;
        }

        /** store cell rejection probability for specified linearCellIndex in the passed rejectionProbabilityCacheCell
         *
         * @param linearCellIndex 1D index of the cell
         * @param superCellCellOffset offset of the superCell in cells
         * @param eFieldBox dataBox giving access to the eField Values of all local cells
         * @param eFieldEnergyUseCacheCell cache of the EField energy use for each cell
         * @param rejectionProbabilityCacheBin entry for bin will be set to -1 if no rejection necessary or rejection
         *  probability, >= 0, otherwise
         * @param sharedResourcesOverSubscribed previous state of sharedResourcesOverSubscribed
         *
         * @return cell is oversubscribed
         */
        template<typename T_EFieldBox, typename T_EFieldEnergyUseCacheCell, typename T_RejectionProbabilityCache_Cell>
        HDINLINE static bool ofCell(
            uint32_t const linearCellIndex,
            VectorIdx const& superCellCellOffset,
            T_EFieldBox const eFieldBox,
            T_EFieldEnergyUseCacheCell const& eFieldEnergyUseCacheCell,
            T_RejectionProbabilityCache_Cell& rejectionProbabilityCacheCell)
        {
            VectorIdx const localCellIndex
                = pmacc::math::mapToND(picongpu::SuperCellSize::toRT(), static_cast<int>(linearCellIndex));
            VectorIdx const cellIndex = localCellIndex + superCellCellOffset;

            /* unit: unit_charge^2 * unit_time^2 / (unit_mass * unit_length^3)
             *  * unit_length^3
             * = unit_charge^2 * unit_time^2 / unit_mass */
            constexpr float_X eps0HalfTimesCellVolume
                = (picongpu::sim.pic.getEps0() / 2._X) * picongpu::sim.pic.getCellSize().productOfComponents();

            /* unit: unit_charge^2 * unit_time^2 / unit_mass * ((unit_mass * unit_length)/(unit_time^2 *
             * unit_charge))^2 = unit_charge^2 * unit_time^2 * unit_mass^2 * unit_length^2 / (unit_mass * unit_time^4 *
             * unit_charge^2) = unit_mass * unit_length^2/ (unit_time^2 * unit_length) = unit_energy */
            float_X const eFieldEnergy = eps0HalfTimesCellVolume * pmacc::math::l2norm2(eFieldBox(cellIndex));

            // unit: eV * 1 = eV * unit_energy/unit_energy = (ev / unit_energy) * unit_energy = unit_energy
            float_X const eFieldEnergyUse
                = picongpu::sim.pic.get_eV() * eFieldEnergyUseCacheCell.energyUsed(linearCellIndex);

            float_X rejectionProbability = -1._X;
            bool sharedResourcesOverSubscribed = false;
            if(eFieldEnergyUse > eFieldEnergy)
            {
                // cell is oversubscribed by suggested changes

                // calculate rejection probability
                rejectionProbability = pmacc::math::max(
                    // proportion of weight we want to reject
                    (eFieldEnergyUse - eFieldEnergy) / eFieldEnergyUse,
                    // but approximately at least one average one macro ion per cell should be rejected
                    1._X / static_cast<float_X>(sim.getTypicalNumParticlesPerCell()));

                // set flag that we found at least one over subscribed resource
                sharedResourcesOverSubscribed = true;
            }

            rejectionProbabilityCacheCell.setCell(linearCellIndex, rejectionProbability);
            return sharedResourcesOverSubscribed;
        }
    };
} // namespace picongpu::particles::atomicPhysics::kernel
