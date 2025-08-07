// /* Copyright 2019-2024 Rene Widera, Pawel Ordyna
//  *
//  * This file is part of PIConGPU.
//  *
//  * PIConGPU is free software: you can redistribute it and/or modify
//  * it under the terms of the GNU General Public License as published by
//  * the Free Software Foundation, either version 3 of the License, or
//  * (at your option) any later version.
//  *
//  * PIConGPU is distributed in the hope that it will be useful,
//  * but WITHOUT ANY WARRANTY; without even the implied warranty of
//  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  * GNU General Public License for more details.
//  *
//  * You should have received a copy of the GNU General Public License
//  * along with PIConGPU.
//  * If not, see <http://www.gnu.org/licenses/>.
//  */

// #pragma once

// #include "picongpu/defines.hpp"
// #include "picongpu/fields/FieldTmp.hpp"
// #include "picongpu/particles/fusion/detail/Creation.hpp"
// #include "picongpu/particles/fusion/detail/FusionContext.hpp"
// #include "picongpu/particles/fusion/detail/ListEntry.hpp"
// #include "picongpu/particles/fusion/detail/cellDensity.hpp"
// #include "picongpu/particles/fusion/fieldSlots.hpp"
// #include "picongpu/particles/filter/IUnary.def"

// #include <pmacc/lockstep.hpp>
// #include <pmacc/mappings/kernel/AreaMapping.hpp>
// #include <pmacc/math/Vector.hpp>
// #include <pmacc/math/operation.hpp>
// #include <pmacc/memory/shared/Allocate.hpp>
// #include <pmacc/mpi/MPIReduce.hpp>
// #include <pmacc/mpi/reduceMethods/Reduce.hpp>
// #include <pmacc/particles/algorithm/ForEach.hpp>
// #include <pmacc/random/RNGProvider.hpp>
// #include <pmacc/random/distributions/Uniform.hpp>

// #include <array>
// #include <cstddef>
// #include <cstdio>

// namespace picongpu::particles::fusion
// {
//     struct InterCollision
//     {
//         HINLINE InterCollision()
//         {
//         }

//     private:
//         DINLINE static uint32_t duplicationCorrection(
//             uint32_t const idx,
//             uint32_t const sizeShort,
//             uint32_t const sizeLong)
//         {
//             uint32_t duplication_correction(1u);
//             if(sizeLong > sizeShort) // no need for duplications when sizeLong = sizeShort
//             {
//                 // Taken from Higginson 2020 DOI: 10.1016/j.jcp.2020.109450
//                 duplication_correction = sizeLong / sizeShort;
//                 uint32_t modulo = sizeLong % sizeShort;
//                 if((idx % sizeShort) < modulo)
//                     duplication_correction += 1u;
//             }
//             return duplication_correction;
//         }

//     public:
//         template<
//             typename T_Reactant1ParBox,
//             typename T_Reactant2ParBox,
//             typename T_Product1ParBox,
//             typename T_Product2ParBox,
//             typename T_Mapping,
//             typename T_Worker,
//             typename T_DeviceHeapHandle,
//             typename T_RngHandle,
//             typename T_SrcCollisionFunctor,
//             typename T_Filter0,
//             typename T_Filter1,
//             typename T_SumCoulombLogBox,
//             typename T_SumSParamBox,
//             typename T_TimesCollidedBox>
//         DINLINE void operator()(
//             T_Worker const& worker,
//             T_Reactant1ParBox pbReactant1,
//             T_Reactant2ParBox pbReactant2,
//             T_Product1ParBox pbProduct1,
//             T_Product2ParBox pbProduct2,
//             IdGenerator& idGen,
//             T_Mapping const mapper,
//             T_DeviceHeapHandle deviceHeapHandle,
//             T_RngHandle rngHandle,
//             T_SrcCollisionFunctor const srcCollisionFunctor,
//             T_Filter0 filter0,
//             T_Filter1 filter1,
//             T_SumCoulombLogBox sumCoulombLogBox,
//             T_SumSParamBox sumSParamBox,
//             T_TimesCollidedBox timesCollidedBox) const
//         {
//             using namespace pmacc::particles::operations;

//             constexpr uint32_t numCellsPerSuperCell = pmacc::math::CT::volume<SuperCellSize>::type::value;

//             PMACC_SMEM(worker, nppc, memory::Array<uint32_t, numCellsPerSuperCell>);

//             PMACC_SMEM(worker, reactant1CellList, detail::ListEntry<T_Reactant1ParBox, numCellsPerSuperCell>);
//             PMACC_SMEM(worker, reactant2CellList, detail::ListEntry<T_Reactant2ParBox, numCellsPerSuperCell>);
//             PMACC_SMEM(worker, reactant1DensityArray, memory::Array<float_X, numCellsPerSuperCell>);
//             PMACC_SMEM(worker, reactant2DensityArray, memory::Array<float_X, numCellsPerSuperCell>);

            
//             DataSpace<simDim> const superCellIdx = mapper.getSuperCellIndex(worker.blockDomIdxND());

//             auto& reactant1SuperCell = pbReactant1.getSuperCell(superCellIdx);
//             uint32_t numReactant1Particles = reactant1SuperCell.getNumParticles();

//             auto& reactant2SuperCell = pbReactant2.getSuperCell(superCellIdx);
//             uint32_t numReactant2Particles = reactant2SuperCell.getNumParticles();

//             // if we have no particles in one species there is no need to perform any calculations
//             if(numReactant1Particles == 0 || numReactant2Particles == 0)
//                 return;

//             // offset of the superCell (in cells, without any guards) to the
//             // origin of the local domain
//             DataSpace<simDim> const localSuperCellOffset = superCellIdx - mapper.getGuardingSuperCells();
//             auto rngOffset = DataSpace<simDim>::create(0);
//             rngOffset.x() = worker.workerIdx();
//             auto numRNGsPerSuperCell = DataSpace<simDim>::create(1);
//             numRNGsPerSuperCell.x() = numFrameSlots;

//             rngHandle.init(localSuperCellOffset * numRNGsPerSuperCell + rngOffset);

//             auto accFilter0 = filter0(worker, localSuperCellOffset);
//             auto accFilter1 = filter1(worker, localSuperCellOffset);

//             auto forEachCell = lockstep::makeForEach<numCellsPerSuperCell>(worker);

//             prepareList(worker, forEachCell, pbReactant1, superCellIdx, deviceHeapHandle, reactant1CellList, nppc, accFilter0);

//             prepareList(worker, forEachCell, pbReactant2, superCellIdx, deviceHeapHandle, reactant2CellList, nppc, accFilter1);

//             using FramePtr0 = typename T_Reactant1ParBox::FramePtr;
//             using FramePtr1 = typename T_Reactant2ParBox::FramePtr;
//             detail::cellDensity<FramePtr0>(worker, forEachCell, reactant1CellList, reactant1DensityArray, accFilter0);
//             detail::cellDensity<FramePtr1>(worker, forEachCell, reactant2CellList, reactant2DensityArray, accFilter1);
//             worker.sync();

//             // shuffle indices list of the longest particle list
//             forEachCell(
//                 [&](uint32_t const linearIdx)
//                 {
//                     uint32_t maxListLength = math::max(reactant1CellList.size(linearIdx), reactant2CellList.size(linearIdx));

//                     uint32_t* parIdListLong = reactant1CellList.size(linearIdx) == maxListLength
//                                                   ? reactant1CellList.particleIds(linearIdx)
//                                                   : reactant2CellList.particleIds(linearIdx);
//                     detail::shuffle(worker, parIdListLong, maxListLength, rngHandle);
//                 });

//             worker.sync();


//             // Run the collision functor for each pair and create particles if needed
//             // create 2 target frames per target species
//             using ProductFramePtr1 = typename T_Product1ParBox::FramePtr;
//             using ProductFramePtr2 = typename T_Product2ParBox::FramePtr;
            
//             using FrameArray1 = memory::Array<ProductFramePtr1, 3>;
//             PMACC_SMEM(worker, product1TargetFrames, FrameArray1);
//             using FrameArray2 = memory::Array<ProductFramePtr2, 3>;
//             PMACC_SMEM(worker, product2TargetFrames, FrameArray2);

//             PMACC_SMEM(worker, particles_created, uint32_t);
//             // get last frames of target species
//             product1TargetFrames[0] = pbProduct1.getLastFrame(superCellIdx);
//             product2TargetFrames[0] = pbProduct2.getLastFrame(superCellIdx);


//             // fill level for Product 1
//             auto& superCellProduct1 = pbProduct1.getSuperCell(superCellIdx);
//             auto& superCellProduct2 = pbProduct2.getSuperCell(superCellIdx);

//             // numParticles % frameSize gives us the number of particles that are in the current partially filled frame
//             uint32_t product1FrameFillLvl = superCellProduct1.getSizeLastFrame();
//             uint32_t product2FrameFillLvl = superCellProduct2.getSizeLastFrame();
//             // only master creates two new frames for each target species - in total we have six frames (two partially filled and two empty)
//             onlyMaster(
//                 [&]()
//                 {
//                     product1TargetFrames[1] = T_Product1ParBox::getEmptyFrame(worker); // can I do it on type and not instance of pb?
//                     T_Product1ParBox::setAsLastFrame(worker, product1TargetFrames[1], superCellIdx);
//                     product1TargetFrames[2] = T_Product1ParBox::getEmptyFrame(worker); 
//                     T_Product1ParBox::setAsLastFrame(worker, product1TargetFrames[2], superCellIdx);
//                     product2TargetFrames[1] = T_Product2ParBox::getEmptyFrame(worker); 
//                     T_Product2ParBox::setAsLastFrame(worker, product2TargetFrames[1], superCellIdx);
//                     product2TargetFrames[2] = T_Product2ParBox::getEmptyFrame(worker); 
//                     T_Product2ParBox::setAsLastFrame(worker, product2TargetFrames[2], superCellIdx);
//                     // reset particles_created counter
//                     particles_created = 0u;
//                 });
//             worker.sync();

            
//             // we check how many particles we can create untill we need a new frame - that's why ProductParBox::frameSize
//             auto ppf1 = T_Product1ParBox::frameSize; // particles per frame
//             auto ppf2 = T_Product2ParBox::frameSize;
//             constexpr uint32_t numPairsAtOnce = math::min(ppf1, ppf2); // we will be createing pairs of particles in chunks of numPairsAtOnce
//             if(numPairsAtOnce==1) printf("Not good");
//             if(numPairsAtOnce<1) printf("You're fucked");
            
            
//             for(int cellIdx = 0; cellIdx<numCellsPerSuperCell; cellIdx++){
//                 // Get the particle counts for the current cell once.
//                 uint32_t const size1 = reactant1CellList.numParticles[cellIdx];
//                 uint32_t const size2 = reactant2CellList.numParticles[cellIdx];

//                 // Get the accessors for each list once.
//                 auto accessor1 = reactant1CellList.getParticlesAccessor(cellIdx);
//                 auto accessor2 = reactant2CellList.getParticlesAccessor(cellIdx);

//                 // Declare the final variables we will use.
//                 uint32_t maxNumParticles;
//                 uint32_t minNumParticles;
//                 decltype(accessor1) pAcc1; // Use decltype to match the type automatically
//                 decltype(accessor2) pAcc2;

//                 // Determine which list is longer and assign all variables in one go.
//                 if (size1 >= size2)
//                 {
//                     maxNumParticles = size1;
//                     minNumParticles = size2;
//                     pAcc1 = accessor1; // pAcc1 is the longer list's accessor
//                     pAcc2 = accessor2; // pAcc2 is the shorter list's accessor
//                 }
//                 else
//                 {
//                     maxNumParticles = size2;
//                     minNumParticles = size1;
//                     pAcc1 = accessor2; // pAcc1 is the longer list's accessor (swapped)
//                     pAcc2 = accessor1; // pAcc2 is the shorter list's accessor (swapped)
//                 }
//                 uint32_t const step = worker.num();
//                 // ## FOR LOOP 1: Iterate over the data in chunks of size ppf. ##
//                 for (uint32_t chunk_start = 0; chunk_start < maxNumParticles; chunk_start += numPairsAtOnce)
//                 {

//                     // 1. FOR LOOP: The parallel grid-stride loop over the CURRENT CHUNK. ##
//                     // Each worker starts at its ID and processes elements within the chunk.
//                     for (int i = chunk_start + worker.id(); i < chunk_start + numPairsAtOnce; i += step)
//                     {
//                         // Boundary check for the very last chunk, which might be smaller than ppf.
//                         if (i >= maxNumParticles) {
//                             break;
//                         }

//                         auto reactant1 = pAcc1(i);
//                         auto reactant2 = pAcc2(i%minNumParticles);
//                         auto duplicationFactor = duplicationCorrection(i,minNumParticles,maxNumParticles);

//                         float3_X product1Momentum = float3_X(0._X);
//                         float3_X product2Momentum = float3_X(0._X);
//                         T_SrcCollisionFunctor fuser = srcCollisionFunctor();
//                         fuser.fuse(reactant1, reactant2, duplicationFactor, product1Momentum, product2Momentum, rngHandle);
                        
//                         if(product1Momentum != float3_X(0._X) || product2Momentum != float3_X(0._X)){

//                             uint32_t freeIndex = alpaka::atomicAdd(
//                                         worker.getAcc(),
//                                         particles_created, // pointer?
//                                         2u,
//                                         ::alpaka::hierarchy::Threads{}); // what other hierarchys are there?

//                             auto indexFromFirstFrame1 = (product1FrameFillLvl + freeIndex);
//                             auto indexFromFirstFrame2 = (product2FrameFillLvl + freeIndex);
//                             auto product1AtReactant1Pos = product1TargetFrames[indexFromFirstFrame1/ppf1][indexFromFirstFrame1%ppf1];
//                             auto product2AtReactant1Pos = product2TargetFrames[indexFromFirstFrame2/ppf2][indexFromFirstFrame2%ppf2];

//                             indexFromFirstFrame1 += 1;
//                             indexFromFirstFrame2 += 1;
//                             auto product1AtReactant2Pos = product1TargetFrames[indexFromFirstFrame1/ppf1][indexFromFirstFrame1%ppf1];
//                             auto product2AtReactant2Pos = product2TargetFrames[indexFromFirstFrame2/ppf2][indexFromFirstFrame2%ppf2];

//                             detail::CreationFusion::createParticles(
//                                 worker,
//                                 idGen,
//                                 reactant1, // reactant 1
//                                 reactant2, // reactant 2
//                                 product1Momentum,// momentum of product 1
//                                 product2Momentum,// momentum of product 2
//                                 product1AtReactant1Pos, // product 1 at location of reactant 1
//                                 product1AtReactant2Pos, // product 1 at location of reactant 2
//                                 product2AtReactant1Pos, // product 2 at location of reactant 1
//                                 product2AtReactant2Pos  // product 2 at location of reactant 2
//                             );
//                         }
//                     }

//                     // Wait for all workers to finish the current chunk before starting ALLOCATION.
//                     worker.sync();
                    
//                     // 2. ALLOCATION: Performed ONCE by the master thread for the new chunk.
//                     if (worker.id() == 0) {
//                         // create frames for product 1
//                         if(particles_created > (ppf1 - product1FrameFillLvl)){
//                             // we need to create a new frame for the first product species
//                             product1TargetFrames[0] = product1TargetFrames[1];
//                             product1TargetFrames[1] = product1TargetFrames[2];
//                             product1TargetFrames[2] = T_Product1ParBox::getEmptyFrame(worker);
//                             T_Product1ParBox::setAsLastFrame(worker, product1TargetFrames[2], superCellIdx);
//                             product1FrameFillLvl = (particles_created*2) - ppf1; // particles_created*2 - number of product1 particles created

//                             if(product1FrameFillLvl > ppf1){
//                                 // we need to create a second new frame
//                                 product1TargetFrames[0] = product1TargetFrames[1];
//                                 product1TargetFrames[1] = product1TargetFrames[2];
//                                 product1TargetFrames[2] = T_Product1ParBox::getEmptyFrame(worker);
//                                 T_Product1ParBox::setAsLastFrame(worker, product1TargetFrames[2], superCellIdx);
//                                 product1FrameFillLvl = product1FrameFillLvl - ppf1; 
//                                 if(product1FrameFillLvl>ppf1){
//                                     printf("Error: too many particles created for product1 frame management");
//                                 }
//                             }
//                         }

//                         // create frames for product 2
//                         if(particles_created > (ppf2 - product2FrameFillLvl)){
//                             // we need to create a new frame for the second product species
//                             product2TargetFrames[0] = product2TargetFrames[1];
//                             product2TargetFrames[1] = product2TargetFrames[2];
//                             product2TargetFrames[2] = T_Product2ParBox::getEmptyFrame(worker);
//                             T_Product2ParBox::setAsLastFrame(worker, product2TargetFrames[2], superCellIdx);
//                             product2FrameFillLvl = (particles_created*2) - ppf2; // particles_created*2 - number of product2 particles created

//                             if(product2FrameFillLvl > ppf2){
//                                 // we need to create a second new frame
//                                 product2TargetFrames[0] = product2TargetFrames[1];
//                                 product2TargetFrames[1] = product2TargetFrames[2];
//                                 product2TargetFrames[2] = T_Product2ParBox::getEmptyFrame(worker);
//                                 T_Product2ParBox::setAsLastFrame(worker, product2TargetFrames[2], superCellIdx);
//                                 product2FrameFillLvl = product2FrameFillLvl - ppf2; 
//                                 if(product2FrameFillLvl>ppf2){
//                                     printf("Error: too many particles created for product2 frame management");
//                                 }
//                             }
//                         }
                    
//                         particles_created = 0u;
//                     }
//                     // All other workers wait here until the allocation is complete.
//                     worker.sync();
//                 }
//             }

//             reactant1CellList.finalize(worker, deviceHeapHandle);
//             reactant2CellList.finalize(worker, deviceHeapHandle);

//         }

//         template<
//             typename T_ParAccessor0,
//             typename T_ParAccessor1,
//             typename T_PairAccessor>
//         DINLINE decltype(auto) populatePairList(
//             T_ParAccessor0 const& parAccessor0,
//             T_ParAccessor1 const& parAccessor1,
//             T_PairAccessor pairList,
//             uint32_t startIndex) const
//         {
//             uint32_t const size0 = parAccessor0.size();
//             uint32_t const size1 = parAccessor1.size();
//             uint32_t const minListLength = math::min(size0, size1);
//             uint32_t const maxListLength = math::max(size0, size1);

//             if(minListLength != 0u)
//             {
//                 for(uint32_t i = 0; i < maxListLength; ++i)
//                 {
//                     auto par0 = parAccessor0[i % size0];
//                     auto par1 = parAccessor1[i % size1];
//                     pairList[startIndex + i] = T_PairAccessor(par0, par1, duplicationCorrection(i, minListLength, maxListLength));
                    
//                 }
//             }
//         }
//     };

//     /* Run kernel for collisions between two species.
//      *
//      * @tparam T_CollisionFunctor A binary particle functor defining a single macro particle collision in
//      * the binary-collision algorithm.
//      * @tparam T_FilterPair A pair of particle filters, each for each species
//      *     in the colliding pair.
//      * @tparam T_ReactantSpecies0 1st colliding species.
//      * @tparam T_ReactantSpecies1 2nd colliding species.
//      */
//     template<
//         typename T_CollisionFunctor,
//         typename T_FilterPair,
//         typename T_ReactantSpecies0,
//         typename T_ReactantSpecies1,
//         typename T_ProductSpecies1,
//         typename T_ProductSpecies2,
//         uint32_t colliderId,
//         uint32_t pairId>
//     struct DoInterCollision
//     {
//         /* Run kernel
//          *
//          * @param deviceHeap A pointer to device heap for allocating particle lists.
//          * @param currentStep The current simulation step.
//          */
//         HINLINE void operator()(std::shared_ptr<DeviceHeap> const& deviceHeap, uint32_t currentStep, IdGenerator idGen)
//         {
//             using Species0 = T_ReactantSpecies0;
//             using FrameType0 = typename Species0::FrameType;
//             using Filter0 = typename T_FilterPair::first ::template apply<Species0>::type;

//             using Species1 = T_ReactantSpecies1;
//             using FrameType1 = typename Species1::FrameType;
//             using Filter1 = typename T_FilterPair::second ::template apply<Species1>::type;

//             using ProductSpecies1 = T_ProductSpecies1;
//             using ProductFrameType1 = typename ProductSpecies1::FrameType;

//             using ProductSpecies2 = T_ProductSpecies2;
//             using ProductFrameType2 = typename ProductSpecies2::FrameType;

//             using CollisionFunctor = T_CollisionFunctor;

//             // Access particle data:
//             DataConnector& dc = Environment<>::get().DataConnector();
//             auto species0 = dc.get<Species0>(FrameType0::getName());
//             auto species1 = dc.get<Species1>(FrameType1::getName());

//             auto productSpecies1 = dc.get<ProductSpecies1>(ProductFrameType1::getName());
//             auto productSpecies2 = dc.get<ProductSpecies2>(ProductFrameType2::getName());

//             // Use mapping information from the first species:
//             auto const mapper = makeAreaMapper<CORE + BORDER>(species0->getCellDescription());

//             //! random number generator
//             using RNGFactory = pmacc::random::RNGProvider<simDim, random::Generator>;
//             using Kernel = typename CollisionFunctor::CallingInterKernel;

//             PMACC_LOCKSTEP_KERNEL(Kernel{}).config(mapper.getGridDim(), *species0)(
//                 species0->getDeviceParticlesBox(),
//                 species1->getDeviceParticlesBox(),
//                 productSpecies1->getDeviceParticlesBox(),
//                 productSpecies2->getDeviceParticlesBox(),
//                 idGen,
//                 mapper,
//                 deviceHeap->getAllocatorHandle(),
//                 RNGFactory::createHandle(),
//                 CollisionFunctor(currentStep),
//                 particles::filter::IUnary<Filter0>{currentStep, idGen},
//                 particles::filter::IUnary<Filter1>{currentStep, idGen},
//                 nullptr,
//                 nullptr,
//                 nullptr);
//         }
//     };
// } // namespace picongpu::particles::collision

