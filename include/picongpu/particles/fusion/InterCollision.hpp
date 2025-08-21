/*
 * Copyright 2019-2024 Rene Widera, Pawel Ordyna, Filip Optolowicz
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
 * along with PIConGPU. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// PIConGPU Includes
#include "picongpu/defines.hpp"
#include "picongpu/fields/FieldTmp.hpp"
#include "picongpu/particles/fusion/detail/Creation.hpp"
#include "picongpu/particles/fusion/detail/FusionContext.hpp"
#include "picongpu/particles/fusion/detail/ListEntry.hpp"
#include "picongpu/particles/fusion/detail/cellDensity.hpp"
#include "picongpu/particles/fusion/fieldSlots.hpp"
#include "picongpu/particles/filter/IUnary.def"

// PMacc Includes
#include <pmacc/lockstep.hpp>
#include <pmacc/mappings/kernel/AreaMapping.hpp>
#include <pmacc/math/Vector.hpp>
#include <pmacc/math/operation.hpp>
#include <pmacc/memory/shared/Allocate.hpp>
#include <pmacc/mpi/MPIReduce.hpp>
#include <pmacc/mpi/reduceMethods/Reduce.hpp>
#include <pmacc/particles/algorithm/ForEach.hpp>
#include <pmacc/random/RNGProvider.hpp>
#include <pmacc/random/distributions/Uniform.hpp>

// Standard Library Includes
#include <array>
#include <cstddef>
#include <cstdio>
#include <utility>

namespace picongpu::particles::fusion
{
    /**
     * @brief Handles inter-species particle collisions within a supercell.
     *
     * This functor orchestrates the binary collision process between two
     * reactant particle species, resulting in the creation of two product species.
     * The process involves preparing particle lists, calculating densities,
     * shuffling for randomness, executing the collision logic in chunks,
     * and managing memory for new particles.
     */
    struct InterCollision
    {
    public:
        HINLINE InterCollision() = default;

        
        template<typename T_worker, typename T_arr>
        DINLINE void zeroArray(T_worker const& worker, T_arr* arr, uint32_t const& size) const
        {
            for (int i = worker.workerIdx();i < size; i += worker.numWorkers())
            {
                arr[i] = 0;
            }
            worker.sync();
        }

        template<bool debug = false, typename T_worker, typename T_arr>
        DINLINE void maxArrayDestroy(T_worker const& worker, T_arr& arr, int const& size) const
        {
            uint32_t pow = 1;
            while(pow < size){
                for(uint32_t i = worker.workerIdx(); pow*(2*i+1) < size; i += 2*pow*worker.numWorkers())
                {
                    arr[2*i*pow] = std::max(arr[2*i*pow],arr[pow*(2*i+1)]);
                }
                pow <<= 1; //*2
                worker.sync();
                if constexpr (debug){
                    if(worker.workerIdx() == 0)
                        printArray(arr);
                    worker.sync();
                }
            }
            // max is now at arr[0];
        }

        template<std::size_t... Is, std::size_t N>
        DINLINE void printArrayImpl(memory::Array<uint32_t, N>& arr, std::index_sequence<Is...>) const
        {
            printf("array: ");
            ((printf("%u, ", arr[Is])), ...);
            printf("\n");
        }

        template<std::size_t N>
        DINLINE void printArray(memory::Array<uint32_t, N>& arr) const
        {
            printArrayImpl(arr, std::make_index_sequence<N>{});
        }


        /**
         * @brief Main operator to execute the inter-species collision kernel.
         *
         * @tparam T_Reactant1ParBox Particle box for the first reactant.
         * @tparam T_Reactant2ParBox Particle box for the second reactant.
         * @tparam T_Product1ParBox Particle box for the first product.
         * @tparam T_Product2ParBox Particle box for the second product.
         * @tparam T_Mapping Maps grid indices to data.
         * @tparam T_Worker The parallel worker (e.g., a CUDA thread).
         * @tparam T_DeviceHeapHandle Handle for dynamic memory allocation on the device.
         * @tparam T_RngHandle Handle for the random number generator.
         * @tparam T_SrcCollisionFunctor Functor containing the physics of the fusion.
         * @tparam T_Filter0 Filter for the first reactant species.
         * @tparam T_Filter1 Filter for the second reactant species.
         */
        template<
            typename T_Reactant1ParBox,
            typename T_Reactant2ParBox,
            typename T_Product1ParBox,
            typename T_Product2ParBox,
            typename T_Mapping,
            typename T_Worker,
            typename T_DeviceHeapHandle,
            typename T_RngHandle,
            typename T_SrcCollisionFunctor,
            typename T_Filter0,
            typename T_Filter1>
        DINLINE void operator()(
            T_Worker const& worker,
            T_Reactant1ParBox reactant1Box,
            T_Reactant2ParBox reactant2Box,
            T_Product1ParBox product1Box,
            T_Product2ParBox product2Box,
            IdGenerator& idGen,
            T_Mapping const mapper,
            T_DeviceHeapHandle deviceHeapHandle,
            T_RngHandle rngHandle,
            T_SrcCollisionFunctor const collisionFunctor,
            T_Filter0 filter0,
            T_Filter1 filter1) const
        {
            // Type aliases for clarity
            using namespace pmacc::particles::operations;
            constexpr auto numCellsPerSuperCell = pmacc::math::CT::volume<SuperCellSize>::type::value;

            // --- 1. Initialization ---

            DataSpace<simDim> const superCellIdx = mapper.getSuperCellIndex(worker.blockDomIdxND());
            DataSpace<simDim> const localSuperCellOffset = superCellIdx - mapper.getGuardingSuperCells();

            auto& reactant1SuperCell = reactant1Box.getSuperCell(superCellIdx);
            auto& reactant2SuperCell = reactant2Box.getSuperCell(superCellIdx);

            // Early exit if there's nothing to collide.
            if (reactant1SuperCell.getNumParticles() == 0 || reactant2SuperCell.getNumParticles() == 0)
            {
                return;
            }
            
            auto onlyMaster = lockstep::makeMaster(worker);

            // --- 2. Shared Memory Allocation ---

            PMACC_SMEM(worker, nppc, memory::Array<uint32_t, numCellsPerSuperCell>);

            PMACC_SMEM(worker, reactant1CellList, detail::ListEntry<T_Reactant1ParBox, numCellsPerSuperCell>);
            PMACC_SMEM(worker, reactant2CellList, detail::ListEntry<T_Reactant2ParBox, numCellsPerSuperCell>);
            PMACC_SMEM(worker, reactant1Density, memory::Array<float_X, numCellsPerSuperCell>);
            PMACC_SMEM(worker, reactant2Density, memory::Array<float_X, numCellsPerSuperCell>);
            

            // --- 3. Prepare Particle Data ---

            // Initialize RNG for this supercell.
            initializeRNG(worker, mapper, superCellIdx, rngHandle, localSuperCellOffset);

            // Prepare filtered lists of particles in each cell of the supercell.
            auto accFilter0 = filter0(worker, localSuperCellOffset);
            auto accFilter1 = filter1(worker, localSuperCellOffset);
            auto forEachCell = lockstep::makeForEach<numCellsPerSuperCell>(worker);
            prepareList(
                worker,
                forEachCell,
                reactant1Box,
                superCellIdx,
                deviceHeapHandle,
                reactant1CellList,
                nppc,
                accFilter0);

            prepareList(
                worker,
                forEachCell,
                reactant2Box,
                superCellIdx,
                deviceHeapHandle,
                reactant2CellList,
                nppc,
                accFilter1);


            
            // Calculate particle densities.
            detail::cellDensity<typename T_Reactant1ParBox::FramePtr>(
                worker,
                forEachCell,
                reactant1CellList,
                reactant1Density,
                accFilter0);
            detail::cellDensity<typename T_Reactant2ParBox::FramePtr>(
                worker,
                forEachCell,
                reactant2CellList,
                reactant2Density,
                accFilter1);

            worker.sync();

            // Find the maximum number of particles (either species) per cell.
            PMACC_SMEM(worker, maxNumParticlesInCell, uint32_t);
            for(uint32_t i = worker.workerIdx(); i<numCellsPerSuperCell; i+=worker.numWorkers()){
                nppc[i] = std::max(reactant1CellList.numParticles[i], reactant2CellList.numParticles[i]);
            }
            // now in nppc[i] we have the maximum number of particles in each cell
            worker.sync();
            maxArrayDestroy<false>(worker, nppc, numCellsPerSuperCell);
            // now in nppc[0] we have the maximum number of particles in the supercell
            onlyMaster([&]() {
                maxNumParticlesInCell = nppc[0];
                
                if constexpr (debugFusion){
                    printf("worker %d: maxNumParticlesInCell = %d\n", worker.workerIdx(), maxNumParticlesInCell);
                }
            });
            // don't need sync

            // --- 4. Shuffle Particle Lists ---
            // To ensure random pairing, shuffle the longer list in each cell.
            forEachCell([&](uint32_t const linearIdx) {
                    uint32_t maxListLength = math::max(reactant1CellList.size(linearIdx), reactant2CellList.size(linearIdx));

                    uint32_t* parIdListLong = reactant1CellList.size(linearIdx) == maxListLength
                                                  ? reactant1CellList.particleIds(linearIdx)
                                                  : reactant2CellList.particleIds(linearIdx);
                    detail::shuffle(worker, parIdListLong, maxListLength, rngHandle);
            });



            // allocate memory for the list where we store how many times did we use the weighting
            // After processing each cell we update the reactant particles using this info
            // We need to subtract the number of times it underwent fusion*minWeighting*something else
            PMACC_SMEM(worker, weightArray, float_X*);
            onlyMaster([&]() {
                constexpr uint32_t chunkSizePerCell = cellListChunkSize * sizeof(float_X);
                weightArray = (float_X*)
                (reactant1CellList.template allocMem<chunkSizePerCell>(worker, sizeof(float_X) * maxNumParticlesInCell, deviceHeapHandle));
            });

            worker.sync();

            // --- 5. Collision Loop ---
            processCollisionsInChunks(
                worker,
                idGen,
                collisionFunctor,
                superCellIdx,
                reactant1CellList,
                reactant2CellList,
                product1Box,
                product2Box,
                reactant1Density,
                reactant2Density,
                weightArray,
                maxNumParticlesInCell,
                rngHandle);

            //! @todo check if this is required
            worker.sync();

            // --- 6. Finalization ---
            reactant1CellList.finalize(worker, deviceHeapHandle);
            reactant2CellList.finalize(worker, deviceHeapHandle);
            // Free the memory allocated for the weighting array
            finalizeWeightArray(worker, deviceHeapHandle, weightArray);
        }


    private:
            /**
             * @brief Frees the memory allocated for the temporary weighting array.
             */
            template<typename T_Worker, typename T_DeviceHeapHandle>
            DINLINE void finalizeWeightArray(
                T_Worker const& worker,
                T_DeviceHeapHandle& deviceHeapHandle,
                float_X*& weightArray) const
            {
                // The master thread that allocated the memory is responsible for freeing it.
                auto onlyMaster = lockstep::makeMaster(worker);
                onlyMaster(
                    [&]()
                    {
                        if(weightArray != nullptr)
                        {
        #if (BOOST_LANG_CUDA || BOOST_COMP_HIP)
                            // Free memory on the GPU device
                            deviceHeapHandle.free(worker.getAcc(), static_cast<void*>(weightArray));
        #else
                            // Free memory on the CPU
                            delete[] weightArray;
        #endif
                            weightArray = nullptr;
                        }
                    });
            }
        /**
         * @brief Corrects for uneven particle list sizes by duplicating particles from the shorter list.
         *
         * This ensures that every particle in the longer list has a collision partner.
         * The formula is from Higginson et al. 2020, DOI: 10.1016/j.jcp.2020.109450.
         *
         * @param idx Index in the longer list.
         * @param sizeShort Size of the shorter particle list.
         * @param sizeLong Size of the longer particle list.
         * @return The duplication factor for the particle from the shorter list.
         */
        DINLINE static uint32_t duplicationCorrection(
            uint32_t const idx,
            uint32_t const sizeShort,
            uint32_t const sizeLong)
        {
            if (sizeLong == sizeShort)
                return 1u;

            uint32_t correction = sizeLong / sizeShort;
            uint32_t modulo = sizeLong % sizeShort;
            if ((idx % sizeShort) < modulo)
            {
                correction += 1u;
            }
            return correction;
        }

        /**
         * @brief Initializes the Random Number Generator for the current supercell.
         */
        template<typename T_Worker, typename T_Mapping, typename T_RngHandle>
        DINLINE void initializeRNG(
            T_Worker const& worker,
            T_Mapping const& mapper,
            DataSpace<simDim> const& superCellIdx,
            T_RngHandle& rngHandle,
            DataSpace<simDim> const localSuperCellOffset) const
        {
            auto rngOffset = DataSpace<simDim>::create(0);
            rngOffset.x() = worker.workerIdx();

            auto numRNGsPerSuperCell = DataSpace<simDim>::create(1);
            numRNGsPerSuperCell.x() = numFrameSlots;

            rngHandle.init(localSuperCellOffset * numRNGsPerSuperCell + rngOffset);
        }

        /**
         * @brief Processes particle collisions in manageable chunks to handle memory allocation.
         */
        template<
            typename T_Worker,
            typename T_SrcCollisionFunctor,
            typename T_Reactant1List,
            typename T_Reactant2List,
            typename T_Product1ParBox,
            typename T_Product2ParBox,
            typename T_RngHandle,
            typename T_DensityArray>
        DINLINE void processCollisionsInChunks(
            T_Worker const& worker,
            IdGenerator& idGen,
            T_SrcCollisionFunctor const& collisionFunctor,
            DataSpace<simDim> const& superCellIdx,
            T_Reactant1List& reactant1CellList,
            T_Reactant2List& reactant2CellList,
            T_Product1ParBox& product1Box,
            T_Product2ParBox& product2Box,
            T_DensityArray& reactant1Density,
            T_DensityArray& reactant2Density,
            float_X* weightingArray,
            uint32_t weightingArraySize,
            T_RngHandle& rngHandle) const
        {
            // Create a small buffer of target frames for new particles
            // Two empty frames at the end because for each fusion reaction we will create two product particles
            // We use 3 frames: [current_partially_filled, next_empty, next_empty]
            constexpr uint32_t NUM_PRODUCT_FRAMES = 3;
            using ProductFramePtr1 = typename T_Product1ParBox::FramePtr;
            using ProductFramePtr2 = typename T_Product2ParBox::FramePtr;
            using FrameArray1 = memory::Array<ProductFramePtr1, NUM_PRODUCT_FRAMES>;
            using FrameArray2 = memory::Array<ProductFramePtr2, NUM_PRODUCT_FRAMES>;
            constexpr uint32_t numCellsPerSuperCell = pmacc::math::CT::volume<SuperCellSize>::type::value;

            PMACC_SMEM(worker, product1Frames, FrameArray1);
            PMACC_SMEM(worker, product2Frames, FrameArray2);
            PMACC_SMEM(worker, particlesCreatedInChunk, uint32_t);
            PMACC_SMEM(worker, product1FillLevel, uint32_t);
            PMACC_SMEM(worker, product2FillLevel, uint32_t);
            // correction factor from Wu et al. 2022, DOI: 10.1063/5.0051178
            PMACC_SMEM(worker, correctionFactor, memory::Array<float_X, numCellsPerSuperCell>); //scalar used for correction factor. n_a/n_ba=n_a/min(w1,w2)
            worker.sync(); // do we need this sync after declaring shared memory?
            // for every cell sum the minimum weighting of the reactants
            for(uint32_t i = worker.workerIdx(); i < numCellsPerSuperCell; i += worker.numWorkers())
            {
                correctionFactor[i] = 0._X; // initialize to zero
                uint32_t const size1 = reactant1CellList.numParticles[i];
                uint32_t const size2 = reactant2CellList.numParticles[i];
                if (size1 == 0 || size2 == 0)
                    continue;

                bool const isList1Longer = (size1 >= size2);
                uint32_t const maxNumParticles = isList1Longer ? size1 : size2;
                uint32_t const minNumParticles = isList1Longer ? size2 : size1;

                auto accessor1 = reactant1CellList.getParticlesAccessor(i);
                auto accessor2 = reactant2CellList.getParticlesAccessor(i);
                for(uint32_t j = 0; j < maxNumParticles; ++j)
                {
                    auto reactant1 = accessor1[j % size1];
                    auto reactant2 = accessor2[j % size2];
                    auto duplicationFactor = duplicationCorrection(j, minNumParticles, maxNumParticles);
                    
                    float_X weightingR1 = reactant1[weighting_] / (isList1Longer ? 1 : duplicationFactor);
                    float_X weightingR2 = reactant2[weighting_] / (isList1Longer ? duplicationFactor : 1);
                    bool const isWeightingR1Greater = (weightingR1 >= weightingR2);
                    correctionFactor[i] += (isWeightingR1Greater ? weightingR2 : weightingR1);
                }
                float_X const densityLonger = isList1Longer ? reactant1Density[i] : reactant2Density[i];
                float_X constexpr cellVolume = sim.pic.getCellSize().productOfComponents();
                correctionFactor[i] = densityLonger * cellVolume / correctionFactor[i]; // n_a/n_ba = n_a/min(w1,w2); for Na > Nb
            }

            
            // Master thread pre-allocates the next two empty frames for each product.
            auto onlyMaster = lockstep::makeMaster(worker);
            onlyMaster([&]() {
                // Get current fill levels and frames.
                product1FillLevel = product1Box.getSuperCell(superCellIdx).getSizeLastFrame();
                product2FillLevel = product2Box.getSuperCell(superCellIdx).getSizeLastFrame();

                product1Frames[0] = product1Box.getLastFrame(superCellIdx);
                product2Frames[0] = product2Box.getLastFrame(superCellIdx);
                // if lastFrame is null allocate a new empty frame
                if (product1Frames[0] == nullptr)
                {
                    product1Frames[0] = product1Box.getEmptyFrame(worker);
                    product1Box.setAsLastFrame(worker, product1Frames[0], superCellIdx);
                }
                if (product2Frames[0] == nullptr)
                {
                    product2Frames[0] = product2Box.getEmptyFrame(worker);
                    product2Box.setAsLastFrame(worker, product2Frames[0], superCellIdx);
                }

                product1Frames[1] = product1Box.getEmptyFrame(worker);
                product1Box.setAsLastFrame(worker, product1Frames[1], superCellIdx);
                product1Frames[2] = product1Box.getEmptyFrame(worker);
                product1Box.setAsLastFrame(worker, product1Frames[2], superCellIdx);

                product2Frames[1] = product2Box.getEmptyFrame(worker);
                product2Box.setAsLastFrame(worker, product2Frames[1], superCellIdx);
                product2Frames[2] = product2Box.getEmptyFrame(worker);
                product2Box.setAsLastFrame(worker, product2Frames[2], superCellIdx);

                particlesCreatedInChunk = 0u;
            });
            worker.sync();

            constexpr auto particlesPerFrame1 = T_Product1ParBox::frameSize;
            constexpr auto particlesPerFrame2 = T_Product2ParBox::frameSize;
            constexpr uint32_t numPairsAtOnce = (particlesPerFrame1 <= particlesPerFrame2) ? particlesPerFrame1 : particlesPerFrame2;

            static_assert(numPairsAtOnce > 0, "Frame size for product species must be greater than zero.");

            // Iterate over all cells in the supercell
            for (int cellIdx = 0; cellIdx < numCellsPerSuperCell ; ++cellIdx)
            {
                // sync() inside
                zeroArray(worker, weightingArray, weightingArraySize); 

                uint32_t const size1 = reactant1CellList.numParticles[cellIdx];
                uint32_t const size2 = reactant2CellList.numParticles[cellIdx];
                if (size1 == 0 || size2 == 0) continue;

                auto accessor1 = reactant1CellList.getParticlesAccessor(cellIdx);
                auto accessor2 = reactant2CellList.getParticlesAccessor(cellIdx);

                // Determine which particle list is longer to iterate over it.
                bool const isList1Longer = (size1 >= size2);
                uint32_t const maxNumParticles = isList1Longer ? size1 : size2;
                uint32_t const minNumParticles = isList1Longer ? size2 : size1;

                bool const isDensity1Greater = (reactant1Density[cellIdx] >= reactant2Density[cellIdx]);
                float_X const minReactantDensity = isDensity1Greater ? reactant2Density[cellIdx] : reactant1Density[cellIdx];


                //! @todo remove
                if(maxNumParticles>weightingArraySize) printf("weightingArraySize is too SMALL!!!");

                // Thread collective loop
                // Process particles in chunks to manage memory frame allocations
                for (uint32_t chunkStart = 0; chunkStart < maxNumParticles; chunkStart += numPairsAtOnce)
                {
                    // Parallel grid-stride loop over the current chunk
                    constexpr uint32_t step = std::min(worker.numWorkers(), numPairsAtOnce);
                    for (int i = chunkStart + worker.workerIdx(); i < chunkStart + numPairsAtOnce && i < maxNumParticles; i += step)
                    {
                        auto reactant1 = accessor1[i % size1];
                        auto reactant2 = accessor2[i % size2];
                        auto duplicationFactor = duplicationCorrection(i, minNumParticles, maxNumParticles);
                        
                        
                        float_X weightingR1 = reactant1[weighting_] / (isList1Longer ? 1 : duplicationFactor);
                        float_X weightingR2 = reactant2[weighting_] / (isList1Longer ? duplicationFactor : 1);

                        bool const isWeightingR1Greater = (weightingR1 >= weightingR2);
                        float_X const minWeighting = isWeightingR1Greater ? weightingR2 : weightingR1;
                        float_X Fmult = maxFmult; 
                        float_X productWeighting = minWeighting/Fmult;
                        if(productWeighting<productMinWeighting){
                            Fmult = std::max(1._X,minWeighting/productMinWeighting);
                            productWeighting = minWeighting/Fmult;
                        }
                        if constexpr (debugFusion){
                            printf("Worker %d, cell %d, i: %d, Fmult: %f, productWeighting: %f, minWeighting: %f\n",
                                worker.workerIdx(), cellIdx, i, Fmult, productWeighting, minWeighting);
                        }

                        float3_X product1Momentum{0._X};
                        float3_X product2Momentum{0._X};

                        // P = n_min * n_a / n_ba * Fmult * minWeighting * dt * (sigma*v_rel*gamma_cm) <- this inside fuse()
                        float_X const probabilityCorrectionFactor = minReactantDensity * correctionFactor[cellIdx] * Fmult * sim.pic.getDt();
                        // print probabilityCorrectionFactor;
                        if constexpr (debugFusion){
                            printf("Worker %d, cell %d, duplicationFactor: %u, probabilityCorrectionFactor: %f\n",
                                worker.workerIdx(), cellIdx, duplicationFactor, probabilityCorrectionFactor);
                            }
                        // The actual fusion physics calculation
                        T_SrcCollisionFunctor fuser = collisionFunctor;
                                                fuser().template fuse<T_Product1ParBox, T_Product2ParBox>(worker, reactant1, reactant2, weightingR1, weightingR2, probabilityCorrectionFactor, product1Momentum, product2Momentum, rngHandle);
                                                
                        // If a reaction occurred, create the product particles
                        if (product1Momentum != float3_X{0._X} || product2Momentum != float3_X{0._X})
                        {
                            weightingArray[i] += productWeighting; // no atomic needed

                            uint32_t freeIndex = alpaka::atomicAdd(
                                worker.getAcc(),
                                &particlesCreatedInChunk,
                                2u, // two particles are created per reaction per product
                                ::alpaka::hierarchy::Threads{});

                            // Calculate indices into the target frames
                            auto idx1 = (product1FillLevel + freeIndex);
                            auto idx2 = (product2FillLevel + freeIndex);
                            
                            auto product1AtR1Pos = product1Frames[idx1 / particlesPerFrame1][idx1 % particlesPerFrame1];
                            auto product2AtR1Pos = product2Frames[idx2 / particlesPerFrame2][idx2 % particlesPerFrame2];

                            idx1++;
                            idx2++;

                            auto product1AtR2Pos = product1Frames[idx1 / particlesPerFrame1][idx1 % particlesPerFrame1];
                            auto product2AtR2Pos = product2Frames[idx2 / particlesPerFrame2][idx2 % particlesPerFrame2];
                            detail::CreationFusion creator;
                            creator.createParticles(
                                worker, idGen,
                                reactant1, reactant2,
                                productWeighting,
                                product1Momentum, product2Momentum,
                                product1AtR1Pos, product1AtR2Pos,
                                product2AtR1Pos, product2AtR2Pos
                            );
                        }

                    } // end grid-stride loop for chunk
                    worker.sync();

                    // Master thread checks if new frames are needed and allocates them.
                    if (worker.workerIdx() == 0)
                    {
                        product1FillLevel = manageFrameAllocation<T_Product1ParBox>(
                            worker, superCellIdx, product1Frames,product1Box, product1FillLevel, particlesCreatedInChunk);

                        product2FillLevel = manageFrameAllocation<T_Product2ParBox>(
                            worker, superCellIdx, product2Frames,product2Box, product2FillLevel, particlesCreatedInChunk);
                        
                        particlesCreatedInChunk = 0u;
                    }
                    worker.sync();
                } // end chunk loop


                // do the magic with weighting and changing momenta
                uint32_t step = std::min(worker.numWorkers(), minNumParticles);
                for (uint32_t chunkStart = 0; chunkStart < maxNumParticles; chunkStart += minNumParticles){
                    for(int i = chunkStart + worker.workerIdx(); i<chunkStart+minNumParticles && i < maxNumParticles; i += step){
                        float_X const oldWeighting1 = accessor1[i % size1][weighting_];
                        float_X const oldWeighting2 = accessor2[i % size2][weighting_];
                        // change the reactant particles according to the weighting array
                        // accessor1[i % size1][weighting_] =std::max(0, accessor1[i % size1][weighting_] - weightingArray[i]);
                        // accessor2[i % size2][weighting_] =std::max(0, accessor2[i % size2][weighting_] - weightingArray[i]);
                        accessor1[i % size1][weighting_] -= weightingArray[i];
                        accessor2[i % size2][weighting_] -= weightingArray[i];

                        accessor1[i % size1][multiMask_] = (accessor1[i % size1][weighting_] > 1e-6);
                        accessor2[i % size2][multiMask_] = (accessor2[i % size2][weighting_] > 1e-6);

                        // change the momenta as well
                        accessor1[i % size1][momentum_] *= accessor1[i % size1][weighting_] / oldWeighting1;
                        accessor2[i % size2][momentum_] *= accessor2[i % size2][weighting_] / oldWeighting2;
                        
                        if( accessor1[i % size1][weighting_]<1e-6){
                            // printf("Deleting Particle: Worker %d, cell %d, i: %d, weighting1: %f, momentum1: %f, %f, %f\n",
                            //     worker.workerIdx(), cellIdx, i,
                            //     accessor1[i % size1][weighting_],
                            //     accessor1[i % size1][momentum_][0],
                            //     accessor1[i % size1][momentum_][1],
                            //     accessor1[i % size1][momentum_][2]);
                            accessor1[i % size1][multiMask_] = 0; // mark for deletion
                        }
                        if( accessor2[i % size2][weighting_]<1e-6){
                            // printf("Deleting Particle: Worker %d, cell %d, i: %d, weighting2: %f, momentum2: %f, %f, %f\n",
                            //     worker.workerIdx(), cellIdx, i,
                            //     accessor2[i % size2][weighting_],
                            //     accessor2[i % size2][momentum_][0],
                            //     accessor2[i % size2][momentum_][1],
                            //     accessor2[i % size2][momentum_][2]);
                            accessor2[i % size2][multiMask_] = 0; // mark for deletion
                        }
                    }
                }
                worker.sync();

            } // end cell loop
        }
        
        /**
         * @brief Manages the allocation of new particle frames when the current ones are full.
         *
         * @return The new fill level for the current frame.
         */
        template<typename T_ProductParBox, typename T_Worker, size_t N>
        DINLINE uint32_t manageFrameAllocation(
            T_Worker const& worker,
            DataSpace<simDim> const& superCellIdx,
            memory::Array<typename T_ProductParBox::FramePtr, N>& productFrames,
            T_ProductParBox productBox,
            uint32_t currentFillLevel,
            uint32_t particlesCreated) const
        {
            constexpr auto particlesPerFrame = T_ProductParBox::frameSize;
            uint32_t newFillLevel = currentFillLevel + particlesCreated;
            
            if (newFillLevel > particlesPerFrame)
            {
                // First new frame is needed
                productFrames[0] = productFrames[1];
                productFrames[1] = productFrames[2];
                productFrames[2] = productBox.getEmptyFrame(worker);
                productBox.setAsLastFrame(worker, productFrames[2], superCellIdx);
                newFillLevel -= particlesPerFrame;

                if (newFillLevel > particlesPerFrame)
                {
                    // Second new frame is also needed
                    productFrames[0] = productFrames[1];
                    productFrames[1] = productFrames[2];
                    productFrames[2] = productBox.getEmptyFrame(worker);
                    productBox.setAsLastFrame(worker, productFrames[2], superCellIdx);
                    newFillLevel -= particlesPerFrame;

                    if (newFillLevel > particlesPerFrame)
                    {
                        printf("Error: Your logic is flawed - too many particles created for frame management to handle.\n");
                    }
                }
            }
            
            if constexpr (debugFusion){
                printf("WorkerIDx: %u, Current fill level: %u, particles created: %u, particles per frame: %u, New fill level: %u\n",
                    worker.workerIdx(), currentFillLevel, particlesCreated, particlesPerFrame, newFillLevel);
            }

            // print fill level after allocation
            return newFillLevel;
        }
    };


    /**
     * @brief Kernel launcher for inter-species collisions.
     *
     * This struct sets up the environment and launches the main `InterCollision`
     * kernel for a specific pair of reactant and product species.
     */
    template<
        typename T_CollisionFunctor,
        typename T_FilterPair,
        typename T_ReactantSpecies0,
        typename T_ReactantSpecies1,
        typename T_ProductSpecies1,
        typename T_ProductSpecies2,
        uint32_t colliderId,
        uint32_t pairId>
    struct DoInterCollision
    {
        /**
         * @brief Runs the collision kernel.
         *
         * @param deviceHeap A pointer to device heap for dynamic memory.
         * @param currentStep The current simulation step.
         * @param idGen The unique ID generator for new particles.
         */
        HINLINE void operator()(std::shared_ptr<DeviceHeap> const& deviceHeap, uint32_t currentStep, IdGenerator idGen)
        {
            // --- Type Aliases for Readability ---
            using Species0 = T_ReactantSpecies0;
            using Filter0 = typename T_FilterPair::first::template apply<Species0>::type;

            using Species1 = T_ReactantSpecies1;
            using Filter1 = typename T_FilterPair::second::template apply<Species1>::type;

            // --- Data Access ---
            auto& dc = Environment<>::get().DataConnector();
            auto species0 = dc.get<Species0>(Species0::FrameType::getName());
            auto species1 = dc.get<Species1>(Species1::FrameType::getName());
            auto productSpecies1 = dc.get<T_ProductSpecies1>(T_ProductSpecies1::FrameType::getName());
            auto productSpecies2 = dc.get<T_ProductSpecies2>(T_ProductSpecies2::FrameType::getName());

            // --- Kernel Configuration and Launch ---
            auto const mapper = makeAreaMapper<CORE + BORDER>(species0->getCellDescription());
            using RNGFactory = pmacc::random::RNGProvider<simDim, random::Generator>;
            using Kernel = InterCollision; // The refactored kernel functor

            PMACC_LOCKSTEP_KERNEL(Kernel{})
                .config(mapper.getGridDim(), *species0)(
                    species0->getDeviceParticlesBox(),
                    species1->getDeviceParticlesBox(),
                    productSpecies1->getDeviceParticlesBox(),
                    productSpecies2->getDeviceParticlesBox(),
                    idGen,
                    mapper,
                    deviceHeap->getAllocatorHandle(),
                    RNGFactory::createHandle(),
                    T_CollisionFunctor(currentStep),
                    particles::filter::IUnary<Filter0>{currentStep, idGen},
                    particles::filter::IUnary<Filter1>{currentStep, idGen});
                
            species0->fillAllGaps();
            species1->fillAllGaps();
            productSpecies1->fillAllGaps();
            productSpecies2->fillAllGaps();
        }
    };
} // namespace picongpu::particles::fusion
