/* Copyright 2013-2024 Axel Huebl, Felix Schmitt, Heiko Burau, Rene Widera,
 *                     Richard Pausch, Alexander Debus, Marco Garten,
 *                     Benjamin Worpitz, Alexander Grund, Sergei Bastrakov,
 *                     Filip Optolowicz
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

#pragma once

#include "picongpu/param/synchrotron.param"

#include <pmacc/meta/ForEach.hpp>
#include <pmacc/particles/traits/FilterByFlag.hpp>

#include <cstdint>

namespace picongpu::simulation::stage
{
    /** Functor for the stage of the PIC loop performing synchrotron radiation
     *
     * Only affects particle species with the Synchrotron attribute.
     */
    class SynchrotronRadiation
    {
    private:
        /** exponential integration
         *
         * @brief
         * 1) approximate function with exponential function by using the xLeft and xMiddle point to fit the exponent
         * 2) integrate the exponential function on log scale,
         *
         * @note overestimates function below xMiddle but underestimates function above xMiddle leading to better fit
         * overal (for bessel function).
         *
         * @param xLeft left point
         * @param xMiddle middle point
         * @param xRight right point
         * @param yLeft function value at xLeft
         * @param yMiddle function value at xMiddle
         */
        template<typename T_Number>
        static T_Number integrateAsExponential(
            T_Number const xLeft,
            T_Number const xMiddle,
            T_Number const xRight,
            T_Number const yLeft,
            T_Number const yMiddle)
        {
            if constexpr(picongpu::particles::synchrotron::T_Debug)
            {
                if(xLeft == xMiddle)
                {
                    throw "xLeft and xMiddle cannot be the same.";
                }
                if(yLeft <= 0 || yMiddle <= 0)
                {
                    throw "yLeft and yMiddle must be positive.";
                }
                if(yLeft == yMiddle)
                {
                    throw "yLeft and yMiddle must not be equal.";
                }
            }

            // fitting function: y = a * e^(b * x)
            T_Number const b = (math::log(yMiddle) - math::log(yLeft)) / (xMiddle - xLeft);
            T_Number const a = yLeft / math::exp(b * xLeft);

            T_Number const integral = (a / b) * (math::exp(b * xRight) - math::exp(b * xLeft));
            return integral;
        }

        /// @todo insure intervals consistent
        /** compute first synchrotorn function. See paper: "Extended particle-in-cell schemes for physics in
         * ultrastrong laser fields: Review and developments" by A. Gonoskov et.Al.
         *
         * @param zq
         *
         * @returns zq * (integral of 2nd kind bessel function from zq to ~infinity)
         */
        static float_64 firstSynchrotronFunction(float_64 const zq)
        {
            using namespace picongpu::particles::synchrotron::params; // for "FirstSynchrotronFunctionParams" struct

            float_64 const log_start = std::log2(zq); // from zq to FirstSynchrotronFunctionParams::logEnd
            float_64 const log_step = (FirstSynchrotronFunctionParams::logEnd - log_start)
                / (FirstSynchrotronFunctionParams::numberSamplePoints - 1);

            float_64 integral = 0;

            float_64 xLeft;
            float_64 xRight = zq;

            for(unsigned int i = 0; i < FirstSynchrotronFunctionParams::numberSamplePoints - 1; ++i)
            {
                xLeft = xRight;
                xRight = math::pow(2., log_start + log_step * (i + 1));
                float_64 xMiddle = (xLeft + xRight) / 2.0;

                // try and catch errors in the bessel function
                try
                {
                    float_64 yLeft = std::cyl_bessel_k(5.0 / 3.0, xLeft);
                    float_64 yMiddle = std::cyl_bessel_k(5.0 / 3.0, xMiddle);
                    integral += integrateAsExponential(
                        xLeft,
                        xMiddle,
                        xRight,
                        yLeft,
                        yMiddle); // computes the integral over one interval: [xLeft, xRight] using the exponential
                                  // approximation between: [xLeft, xMiddle]
                }
                catch(std::exception& e)
                {
                    std::cout << "Caught exception when precomputing firstSynchrotronFunction: " << e.what()
                              << std::endl;
                }
            }
            return zq * integral;
        }

        //! @note see paper: "Extended particle-in-cell schemes for physics in ultrastrong laser fields: Review and
        //! developments" by A. Gonoskov et.Al.
        static float_64 secondSynchrotronFunction(float_64 const x)
        {
            return x * std::cyl_bessel_k(2. / 3, x);
        }

    public:
        /** Create a synchrotron radiation instance: precompute first and second synchrotron functions
         *
         * @param cellDescription mapping for kernels
         */
        SynchrotronRadiation(MappingDesc const cellDescription) : cellDescription(cellDescription)
        {
            using namespace picongpu::particles::synchrotron::params; // for "InterpolationParams" struct
            auto data_space = DataSpace<2>{InterpolationParams::numberTableEntries, 2};
            auto grid_layout = GridLayout<2>{data_space};

            // capture the space for table and variables
            tableValuesF1F2 = std::make_shared<GridBuffer<float_X, 2>>(grid_layout);
            failedRequirementQ = std::make_shared<GridBuffer<int32_t, 1>>(DataSpace<1>{1});
            // set the values of variables to false
            failedRequirementQ->getHostBuffer().getDataBox()(DataSpace<1>{0}) = false;
            failedRequirementPrinted = false;

            constexpr float_64 minZqExp = InterpolationParams::minZqExponent;
            constexpr float_64 maxZqExp = InterpolationParams::maxZqExponent;
            constexpr float_64 tableEntries = InterpolationParams::numberTableEntries;

            // precompute F1 and F2 on log scale
            for(uint32_t iZq = 0; iZq < tableEntries; iZq++)
            {
                float_64 zq = std::pow(2, minZqExp + (maxZqExp - minZqExp) * iZq / (tableEntries - 1));
                // inverse function for index retrieval:
                // index = (log2(zq) - minZqExp) / (maxZqExp - minZqExp) * (tableEntries-1);

                float_64 const F1 = firstSynchrotronFunction(zq);
                float_64 const F2 = secondSynchrotronFunction(zq);

                tableValuesF1F2->getHostBuffer().getDataBox()(DataSpace<2>{iZq, u32(Accessor::f1)})
                    = static_cast<float_X>(F1);
                tableValuesF1F2->getHostBuffer().getDataBox()(DataSpace<2>{iZq, u32(Accessor::f2)})
                    = static_cast<float_X>(F2);
            }
            // move the data to the device
            tableValuesF1F2->hostToDevice();
            failedRequirementQ->hostToDevice();

            //! debug only, print F1 and F2, @todo remove
            if constexpr(picongpu::particles::synchrotron::T_Debug)
            {
                int timesMore = 1;
                for(uint32_t zq = 0; zq < InterpolationParams::numberTableEntries * timesMore; zq++)
                {
                    float_64 zq_ = std::pow(
                        2,
                        minZqExp
                            + (maxZqExp - minZqExp) * zq
                                / float_64(InterpolationParams::numberTableEntries * timesMore));
                    uint32_t index
                        = (log2(zq_) - minZqExp) * (InterpolationParams::numberTableEntries) / (maxZqExp - minZqExp);


                    std::cout << "fractional index = " << zq / float(timesMore) << " index = " << index << std::endl;
                    std::cout << "zq = " << zq_
                              << " F1 = " << tableValuesF1F2->getHostBuffer().getDataBox()(DataSpace<2>{index, 0})
                              << " F2 = " << tableValuesF1F2->getHostBuffer().getDataBox()(DataSpace<2>{index, 1})
                              << std::endl;
                }
            }
        }

        /** Radiation happens here
         *
         * @param step index of time iteration
         */
        void operator()(uint32_t const step)
        {
            using pmacc::particles::traits::FilterByFlag;
            using SpeciesWithSynchrotron = typename FilterByFlag<VectorAllSpecies, picongpu::synchrotron<>>::type;

            // call the synchrotron radiation for each particle species with the synchrotron attribute
            pmacc::meta::ForEach<SpeciesWithSynchrotron, particles::CallSynchrotron<boost::mpl::_1>>
                synchrotronRadiation;

            synchrotronRadiation(
                cellDescription,
                step,
                tableValuesF1F2->getDeviceBuffer().getDataBox(),
                failedRequirementQ);

            // check if the requirements are met
            if constexpr(particles::synchrotron::params::supressRequirementWarning == false)
            {
                // retrieve the failedRequirementQ from the device
                failedRequirementQ->deviceToHost();
                if(failedRequirementQ->getHostBuffer().getDataBox()(DataSpace<1>{0}) == true)
                {
                    if((failedRequirementPrinted) == false)
                    {
                        printf("Synchrotron Extension requirement1 or requirement2 failed; should be less than 0.1 -> "
                               "reduce the timestep. \n\tCheck the requrement by specifying the predicted maxHeff and "
                               "maxGamma in"
                               "\n\tpicongpu/docs/source/usage/workflows/synchrotronRequirements.py\n");
                        printf("This warning is printed only once per simulation. Next warnings are dots.\n");
                        failedRequirementPrinted = true;
                    }
                    printf(".");
                    failedRequirementQ->getHostBuffer().getDataBox()(DataSpace<1>{0})
                        = static_cast<int32_t>(false); // reset the requirement flag
                    failedRequirementQ->hostToDevice();
                }
            }
        }

    private:
        //! Mapping for kernels
        MappingDesc cellDescription;

        /// <CACHE F1F2>
        // precomputed first and second synchrotron functions:
        // -> 2d grid of floats_64 -> tableValuesF1F2[zq][0/1] ; 0/1 = F1/F2
        std::shared_ptr<GridBuffer<float_X, 2>> tableValuesF1F2;
        // used for table access
        enum struct Accessor : uint32_t
        {
            f1 = 0u,
            f2 = 1u
        };
        static constexpr uint32_t u32(Accessor const t)
        {
            return static_cast<uint32_t>(t);
        }
        /// <END CACHE F1F2>

        // flag to check if the requirements 1 and 2 are met -> in
        // picongpu/include/picongpu/particles/synchrotron/AlgorithmSynchrotron.hpp <- in class SynchrotronIdea(); we
        // check the requirements
        std::shared_ptr<GridBuffer<int32_t, 1>> failedRequirementQ;
        bool failedRequirementPrinted; // this means that the operator() can't be const anymore
                                       //@idea: maby using shared pointer for this will be faster?
    };
} // namespace picongpu::simulation::stage
