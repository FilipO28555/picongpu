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
#include "picongpu/particles/fusion/detail/arrayHelpers.hpp"
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
    struct IntraCollision
    {
    public:
        HINLINE IntraCollision() = default;

        

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
         * @tparam T_Filter Filter for the first reactant species.
         */
        template<
            typename T_ReactantParBox,
            typename T_Product1ParBox,
            typename T_Product2ParBox,
            typename T_Mapping,
            typename T_Worker,
            typename T_DeviceHeapHandle,
            typename T_RngHandle,
            typename T_SrcCollisionFunctor,
            typename T_Filter>
        DINLINE void operator()(
            T_Worker const& worker,
            T_ReactantParBox reactantBox,
            T_Product1ParBox product1Box,
            T_Product2ParBox product2Box,
            IdGenerator& idGen,
            T_Mapping const mapper,
            T_DeviceHeapHandle deviceHeapHandle,
            T_RngHandle rngHandle,
            T_SrcCollisionFunctor const collisionFunctor,
            T_Filter filter) const
        {
            // Type aliases for clarity
            using namespace pmacc::particles::operations;
            constexpr auto numCellsPerSuperCell = pmacc::math::CT::volume<SuperCellSize>::type::value;

            // --- 1. Initialization ---

            DataSpace<simDim> const superCellIdx = mapper.getSuperCellIndex(worker.blockDomIdxND());
            DataSpace<simDim> const localSuperCellOffset = superCellIdx - mapper.getGuardingSuperCells();

            auto& reactantSuperCell = reactantBox.getSuperCell(superCellIdx);

            // Early exit if there's nothing to collide.
            if (reactantSuperCell.getNumParticles() == 0)
            {
                return;
            }
            
            auto onlyMaster = lockstep::makeMaster(worker);

            // --- 2. Shared Memory Allocation ---

            PMACC_SMEM(worker, nppc, memory::Array<uint32_t, numCellsPerSuperCell>);

            PMACC_SMEM(worker, reactantCellList, detail::ListEntry<T_ReactantParBox, numCellsPerSuperCell>);
            PMACC_SMEM(worker, reactantDensity, memory::Array<float_X, numCellsPerSuperCell>);

            // --- 3. Prepare Particle Data ---

            // Initialize RNG for this supercell.
            initializeRNG(worker, mapper, superCellIdx, rngHandle, localSuperCellOffset);

            // Prepare filtered lists of particles in each cell of the supercell.
            auto accFilter = filter(worker, localSuperCellOffset);
            auto forEachCell = lockstep::makeForEach<numCellsPerSuperCell>(worker);
            prepareList(
                worker,
                forEachCell,
                reactantBox,
                superCellIdx,
                deviceHeapHandle,
                reactantCellList,
                nppc,
                accFilter);

            // Calculate particle densities.
            detail::cellDensity<typename T_ReactantParBox::FramePtr>(
                worker,
                forEachCell,
                reactantCellList,
                reactantDensity,
                accFilter);

            worker.sync();

            PMACC_SMEM(worker, maxNumParirsInCell, uint32_t);
            // now in nppc[i] we have the number of particles in each cell
            worker.sync();
            detail::maxArrayDestroy<false>(worker, nppc, numCellsPerSuperCell);
            // now in nppc[0] we have the maximum number of particles in the supercell
            onlyMaster([&]() {
                // we only need half as much pairs as particles
                maxNumParirsInCell = nppc[0]/2+1;
                
                if constexpr (debugFusion){
                    printf("worker %d: maxNumParirsInCell = %d\n", worker.workerIdx(), maxNumParirsInCell);
                }
            });
            // don't need sync

            // --- 4. Shuffle Particle Lists ---
            // To ensure random pairing, shuffle the longer list in each cell.
            forEachCell(
                [&](uint32_t const linearIdx)
                {
                    detail::shuffle(
                        worker,
                        reactantCellList.particleIds(linearIdx),
                        reactantCellList.size(linearIdx),
                        rngHandle);
                });


            // allocate memory for the list where we store how many times did we use the weighting
            // After processing each cell we update the reactant particles using this info
            // We need to subtract the number of times it underwent fusion*minWeighting*something else
            PMACC_SMEM(worker, weightArray, float_X*);
            onlyMaster([&]() {
                constexpr uint32_t chunkSizePerCell = cellListChunkSize * sizeof(float_X);
                weightArray = (float_X*)
                (reactantCellList.template allocMem<chunkSizePerCell>(worker, sizeof(float_X) * maxNumParirsInCell, deviceHeapHandle));
            });

            worker.sync();

            // --- 5. Collision Loop ---
            processCollisionsInChunks(
                worker,
                idGen,
                collisionFunctor,
                superCellIdx,
                reactantCellList,
                product1Box,
                product2Box,
                reactantDensity,
                weightArray,
                maxNumParirsInCell,
                rngHandle);

            //! @todo check if this is required
            worker.sync();

            // --- 6. Finalization ---
            reactantCellList.finalize(worker, deviceHeapHandle);
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
            
            
        /* Get the duplication correction for a collision
         *
         * A particle duplication is how many times a particle collides in the current time step.
         * The duplication correction is equal to max(D_0, D_1) where D_0, D_1 are the duplications
         * of the colliding particles.
         * In a case of internal collisions all particles collide once expect for the 1st one if the total
         * number of particles is odd.
         *
         * @param idx the index of the particle in the particle list of the cell
         * @param sizeAll the length of the particle list of the cell
         */
        DINLINE static uint32_t duplicationCorrection(uint32_t const idx, uint32_t const sizeAll)
        {
            // All particle collide once.
            if(sizeAll % 2u == 0u)
                return 1u;
            // The first particle collides twice and the last one collides with the first one.
            // for the rest the correction is 1.
            else
                return (idx == 0u || idx == sizeAll - 1u) ? 2u : 1u;
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
            typename T_ReactantList,
            typename T_Product1ParBox,
            typename T_Product2ParBox,
            typename T_RngHandle,
            typename T_DensityArray>
        DINLINE void processCollisionsInChunks(
            T_Worker const& worker,
            IdGenerator& idGen,
            T_SrcCollisionFunctor const& collisionFunctor,
            DataSpace<simDim> const& superCellIdx,
            T_ReactantList& reactantCellList,
            T_Product1ParBox& product1Box,
            T_Product2ParBox& product2Box,
            T_DensityArray& reactantDensity,
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
            // for every cell sum the minimum weighting of the reactant
            for(uint32_t i = worker.workerIdx(); i < numCellsPerSuperCell; i += worker.numWorkers())
            {
                uint32_t const size = reactantCellList.numParticles[i];
                if (size < 2)
                    continue;

                bool const isOdd = (size % 2 == 1);
                auto accessor = reactantCellList.getParticlesAccessor(i);
                
                correctionFactor[i] = 0._X; // initialize to zero
                for(uint32_t j = 0; j < size; j+=2)
                {
                    //parOdd becomes 0 if j+1 == size and size is odd
                    auto parEven = accessor[j];
                    auto parOdd = accessor[(j + 1) % size]; 
                    // we divide by 2 if the number of particles is odd and we are at the first or last pair
                    uint32_t duplicationCorrection1 = 1+(isOdd && j == 0u); 
                    uint32_t duplicationCorrection2 = 1+(isOdd && j == size - 1u); 

                    float_X weightingR1 = parEven[weighting_] / duplicationCorrection1;
                    float_X weightingR2 = parOdd[weighting_] / duplicationCorrection2;
                    bool const isWeightingR1Greater = (weightingR1 >= weightingR2);
                    correctionFactor[i] += (isWeightingR1Greater ? weightingR2 : weightingR1);
                }
                float_X const densityLonger = reactantDensity[i];
                float_X constexpr cellVolume = sim.pic.getCellSize().productOfComponents();
                correctionFactor[i] = densityLonger * cellVolume / correctionFactor[i]; // n_a/n_aa = n_a/min(w1,w2)
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
                detail::zeroArray(worker, weightingArray, weightingArraySize); 

                uint32_t const size = reactantCellList.numParticles[cellIdx];
                if (size < 2) continue;
                bool const isOdd = (size % 2 == 1);

                auto accessor = reactantCellList.getParticlesAccessor(cellIdx);
                float_X const reactantDensity_local = reactantDensity[cellIdx];

                // Thread collective loop
                // Process particles in chunks to manage memory frame allocations
                for (uint32_t chunkStart = 0; chunkStart < size; chunkStart += numPairsAtOnce)
                {
                    // Parallel grid-stride loop over the current chunk
                    constexpr uint32_t step = std::min(worker.numWorkers(), numPairsAtOnce);
                    for (int i = chunkStart + worker.workerIdx(); i < chunkStart + numPairsAtOnce && i < size/2; i += step)
                    {
                        uint32_t j = i * 2;
                        //parOdd becomes 0 if j+1 == size and size is odd
                        auto parEven = accessor[j];
                        auto parOdd = accessor[(j + 1) % size]; 
                        // we divide by 2 if the number of particles is odd and we are at the first or last pair
                        uint32_t duplicationCorrection1 = 1+(isOdd && j == 0u); 
                        uint32_t duplicationCorrection2 = 1+(isOdd && j == size - 1u); 

                        // for debug
                        uint32_t duplicationFactor = duplicationCorrection1+duplicationCorrection2-1;

                        float_X weightingR1 = parEven[weighting_] / duplicationCorrection1;
                        float_X weightingR2 = parOdd[weighting_] / duplicationCorrection2;

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

                        // WU:
                        // P = n_min * n_a / n_ba * Fmult * minWeighting * dt * (sigma*v_rel*gamma_cm) <- this inside fuse()
                        float_X const probabilityCorrectionFactor = reactantDensity_local * correctionFactor[cellIdx] * Fmult * sim.pic.getDt();

                        // Higginson:
                        // P = minWeighting/V * maxNumParticles * Fmult * sim.pic.getDt() * (sigma*v_rel*gamma_cm) <- this inside fuse()
                        // float_X const maxWeighting = isWeightingR1Greater ? weightingR1 : weightingR2;
                        // float_X const maxReactantDensity = isDensity1Greater ? reactant1Density[cellIdx] : reactant2Density[cellIdx];
                        // float_X constexpr cellVolume = sim.pic.getCellSize().productOfComponents();
                        // float_X const probabilityCorrectionFactor = maxWeighting * maxNumParticles * Fmult * sim.pic.getDt() / cellVolume;

                        // print probabilityCorrectionFactor;
                        if constexpr (debugFusion){
                            printf("Worker %d, cell %d, duplicationFactor: %u, probabilityCorrectionFactor: %f\n",
                                worker.workerIdx(), cellIdx, duplicationFactor, probabilityCorrectionFactor);
                            }
                        // The actual fusion physics calculation
                        T_SrcCollisionFunctor fuser = collisionFunctor;
                                                fuser().template fuse<T_Product1ParBox, T_Product2ParBox>(worker, parEven, parOdd, weightingR1, weightingR2, probabilityCorrectionFactor, product1Momentum, product2Momentum, rngHandle);
                                                
                        // If a reaction occurred, create the product particles
                        if (product1Momentum != float3_X{0._X} || product2Momentum != float3_X{0._X})
                        {
                            weightingArray[i] = productWeighting;

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
                                parEven, parOdd,
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

                // --- 5b. Update Reactant Particles ---
                // loop over all particles, and the first one may have to add weighting from weightingArray[size/2]
                for(int i = worker.workerIdx(); i<size; i += worker.numWorkers()){
                    // necessary for weighting array
                    uint32_t pairIdx = i/2;

                    float_X const oldWeighting = accessor[i][weighting_];
                    // change the reactant particles according to the weighting array
                    // weightingArrays is indexed by the pair index
                    accessor[i][weighting_] -= weightingArray[pairIdx];
                    // if the number of particles is odd the first and last particle could interact twice
                    if(isOdd && (i==0 || i == size-1)) accessor[i][weighting_] -= weightingArray[size/2];
                    
                    // delete the particle
                    accessor[i][multiMask_] = (accessor[i][weighting_] > 1e-6);
                    if (debugFusion){
                        // print weighting before and after
                        printf("Worker %d, cell %d, i: %d, oldWeighting: %f, newWeighting: %f\n",
                            worker.workerIdx(), cellIdx, i, oldWeighting, accessor[i][weighting_]);
                        if(accessor[i][weighting_] < 1e-6){
                            printf("Deleting particle with index %d in cell %d\n", i, cellIdx);
                        }
                    }

                    // change the momenta as well
                    accessor[i][momentum_] *= accessor[i][weighting_] / oldWeighting;
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
     * This struct sets up the environment and launches the main `IntraCollision`
     * kernel for a specific pair of reactant and product species.
     */
    template<
        typename T_CollisionFunctor,
        typename T_FilterPair,
        typename T_ReactantSpecies,
        typename T_ProductSpecies1,
        typename T_ProductSpecies2,
        uint32_t colliderId,
        uint32_t pairId>
    struct DoIntraCollision
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
            if(debugFusion){
                printf("Starting IntraCollision for colliderId: %u, pairId: %u at step %u\n", colliderId, pairId, currentStep);
            }
            // --- Type Aliases for Readability ---
            using Species = T_ReactantSpecies;
            using Filter = typename T_FilterPair::first::template apply<Species>::type;

            // --- Data Access ---
            auto& dc = Environment<>::get().DataConnector();
            auto species = dc.get<Species>(Species::FrameType::getName());
            auto productSpecies1 = dc.get<T_ProductSpecies1>(T_ProductSpecies1::FrameType::getName());
            auto productSpecies2 = dc.get<T_ProductSpecies2>(T_ProductSpecies2::FrameType::getName());

            // --- Kernel Configuration and Launch ---
            auto const mapper = makeAreaMapper<CORE + BORDER>(species->getCellDescription());
            using RNGFactory = pmacc::random::RNGProvider<simDim, random::Generator>;
            using Kernel = IntraCollision; // The refactored kernel functor

            PMACC_LOCKSTEP_KERNEL(Kernel{})
                .config(mapper.getGridDim(), *species)(
                    species->getDeviceParticlesBox(),
                    productSpecies1->getDeviceParticlesBox(),
                    productSpecies2->getDeviceParticlesBox(),
                    idGen,
                    mapper,
                    deviceHeap->getAllocatorHandle(),
                    RNGFactory::createHandle(),
                    T_CollisionFunctor(currentStep),
                    particles::filter::IUnary<Filter>{currentStep, idGen});
                
            species->fillAllGaps();
            productSpecies1->fillAllGaps();
            productSpecies2->fillAllGaps();
        }
    };
} // namespace picongpu::particles::fusion
