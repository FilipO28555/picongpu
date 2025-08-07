# Handle tests with dimensionality
file(GLOB_RECURSE DIMENSIONAL_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
foreach(dim 2 3)
    foreach(testCaseFilepath ${DIMENSIONAL_TESTS})
        get_filename_component(testCaseFilename ${testCaseFilepath} NAME)
        string(REPLACE "UT.cpp" "" testCase ${testCaseFilename})
        set(testExe "${PROJECT_NAME}-${testCase}-${dim}D")
        alpaka_add_executable(${testExe} ${testCaseFilepath})
        target_compile_definitions(${testExe} PRIVATE TEST_DIM=${dim})
        target_link_libraries(${testExe} PUBLIC Catch2 Catch2WithMain)
        target_link_libraries(${testExe} PRIVATE pmacc::pmacc)
        if(${testExe} MATCHES ".*CLI.*")
          target_link_libraries(${testExe} PUBLIC argparser)
        endif()
        target_include_directories(
            ${testExe}
            SYSTEM PRIVATE
            $<TARGET_PROPERTY:nlohmann_json::nlohmann_json,INTERFACE_INCLUDE_DIRECTORIES>)
        target_link_libraries(${testExe} PRIVATE nlohmann_json::nlohmann_json)
        add_test(NAME "${testCase}-${dim}D" COMMAND mpiexec ${MPI_RUNTIME_FLAGS} -n 1 ./${testExe})
    endforeach()
    string(REPLACE "-DTEST_DIM=${dim}" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endforeach()
