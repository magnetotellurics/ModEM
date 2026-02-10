cmake_minimum_required(VERSION 3.17.0)

project(ModEM VERSION 0
        DESCRIPTION "ModEM"
        LANGUAGES Fortran C)

enable_language(Fortran C)

option(BUILD_MPI "Build MPI version of MPI" ON)
option(FG "Build the Fine-Grained ModEM Version - (Only compatiable with SP2)" OFF)
option(USE_C_TIMERS "Use the C timers in ModEM_timers" OFF)

set(FORWARD_FORMULATION "SP2" CACHE STRING "Set the forward formulation: < MF | SP | SP2 >")
set(MODEM_BUILD_DIMS "3D" CACHE STRING "Specify whether to build either the 2D or 3D versino of ModEM < 2D | 3D >")
set(BUILD_GPU "OFF" CACHE STRING "Specify whether to build GPU options: < OFF | CUDA | HIP >")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(SUPPORTED_FORWARD_FORMULATIONS "MF;SP;SP2")
set(SUPPORTED_MODEM_DIMS "2D;3D")
set(SUPPORTED_GPU_OPTS "off;OFF;CUDA;HIP")

if (USE_C_TIMERS)
        enable_language(C)
endif()

if (NOT BUILD_GPU IN_LIST SUPPORTED_GPU_OPTS)
        message(FATAL_ERROR "-DBUILD_GPU option ('${BUILD_GPU}') was not valid - choose from: ${SUPPORTED_GPU_OPTS}")
endif()

if (BUILD_GPU STREQUAL CUDA)
        enable_language(CUDA)

        # Turning CMAKE_CUDA_ARCHITECTURES to OFF will make
        # nvcc target the lowest possible CUDA Compute Capability.
        # But this will allow users to pass in their own desired Compute
        # Capability via: -DCMAKE_CUDA_ARCHITECTURES=<option>.
        # Options here can be: all, all-major, native, or a specific version(s): "35;50;72"
        # Note: The removal of decimal places.
        #
        # One option here would be to use 'native', but on HPC systems, there might not be a 
        # GPU available where one compiles the code as they are only available on compute nodes,
        # thus, we shouldn't use it here, but users can pass it in themselves.
        if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES OFF)
        else()
            set(CUDA_ARCHITECTURES $CMAKE_CUDA_ARCHITECTURES)
        endif()

elseif (BUILD_GPU STREQUAL HIP)
        enable_language(HIP)
endif()

if (NOT FORWARD_FORMULATION IN_LIST SUPPORTED_FORWARD_FORMULATIONS )
        message(FATAL_ERROR "-DFORWARD_FORMULATION option ('${FORWARD_FORMULATION}') was not valid choose from: ${SUPPORTED_FORWARD_FORMULATIONS}")
endif()

if (NOT MODEM_BUILD_DIMS IN_LIST SUPPORTED_MODEM_DIMS)
        message(FATAL_ERROR "-DMODEM_BUILD_DIMS option ('${MODEM_BUILD_DIMS}') was not valid choose from ${SUPPORTED_MODEM_DIMS}")
endif()

if (FG AND NOT FORWARD_FORMULATION STREQUAL SP2)
        message(FATAL_ERROR "You cannot build 'FG' with '${FORWARD_FORMULATION}' - Only SP2: -DFG=on -DFORWARD_FORMULATION=SP2")
endif()

if (BUILD_MPI)
        add_compile_definitions(MPI)
endif()

if (BUILD_GPU STREQUAL "CUDA")
        add_compile_definitions(CUDA)
elseif (BUILD_GPU STREQUAL "HIP")
        add_compile_definitions(HIP)
endif()

if (FG)
        add_compile_definitions(FG)
endif()

# Source Code
add_subdirectory(f90)
