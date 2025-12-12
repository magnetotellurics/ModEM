cmake_minimum_required(VERSION 3.17.0)

project(ModEM VERSION 0
        DESCRIPTION "ModEM"
        LANGUAGES Fortran)

enable_language(Fortran)

option(MPI "Build MPI version of MPI" ON)
option(FG "Build the Fine-Grained ModEM Version - (Only compatiable with SP2)" OFF)
option(USE_C_TIMERS "Use the C timers in ModEM_timers" OFF)

set(FORWARD_FLAVOR "SP2" CACHE STRING "Set the forward flavor: < MF | SP | SP2 >")
set(MODEM_BUILD_DIMS "3D" CACHE STRING "Specify whether to build either the 2D or 3D versino of ModEM < 2D | 3D >")
set(BUILD_GPU "OFF" CACHE STRING "Specify whether to build GPU options: < OFF | CUDA | HIP >")

set(SUPPORTED_FORWARD_FLAVORS "MF;SP;SP2")
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
elseif (BUILD_GPU STREQUAL HIP)
        enable_language(HIP)
endif()

if (NOT FORWARD_FLAVOR IN_LIST SUPPORTED_FORWARD_FLAVORS )
        message(FATAL_ERROR "-DFORWARD_FLAVOR option ('${FORWARD_FLAVOR}') was not valid choose from: ${SUPPORTED_FORWARD_FLAVORS}")
endif()

if (NOT MODEM_BUILD_DIMS IN_LIST SUPPORTED_MODEM_DIMS)
        message(FATAL_ERROR "-DMODEM_BUILD_DIMS option ('${MODEM_BUILD_DIMS}') was not valid choose from ${SUPPORTED_MODEM_DIMS}")
endif()

if (FG AND NOT FORWARD_FLAVOR STREQUAL SP2)
        message(FATAL_ERROR "You cannot build 'FG' with '${FORWARD_FLAVOR}' - Only SP2: -DFG=on -DFORWARD_FLAVOR=SP2")
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-cpp>
        $<$<COMPILE_LANGUAGE:Fortran>:-O3>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>
        $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>
    )
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-O3>
        $<$<COMPILE_LANGUAGE:Fortran>:-fpp>
        $<$<COMPILE_LANGUAGE:Fortran>:-w>
    )
    add_link_options(
        "-L${MKLROOT}/lib/intel64"
        "-lmkl_intel_lp64"
        "-lmkl_intel_thread"
        "-lmkl_core"
        "-lpthread"
        "-qopenmp"
        "-parallel"
    )
else()
    messaing(FATAL_ERROR "Unknown/Unsupported Fortran Compiler: '${CMake_Fortran_COMPILER_ID}' - Supported compilers are < GNU | Ifort")
endif()


# Source Code
add_subdirectory(f90)
