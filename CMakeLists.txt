cmake_minimum_required(VERSION 4.1.0)

project(ModEM VERSION 0
        DESCRIPTION "ModEM"
        LANGUAGES Fortran)

enable_language(Fortran)

option(MPI "Build MPI version of MPI" ON)
option(FG "Build the Fine-Grained ModEM Version - (Only compatiable with SP2)" OFF)
set(FORWARD_FLAVOR "SP2" CACHE STRING "Set the forward flavor: < MF | SP | SP2 >")
set(MODEM_BUILD_DIMS "3D" CACHE STRING "Specify whether to build either the 2D or 3D versino of ModEM < 2D | 3D >")

set(SUPPORTED_FORWARD_FLAVORS "MF;SP;SP2")
set(SUPPORTED_MODEM_DIMS "2D;3D")

if (NOT FORWARD_FLAVOR IN_LIST SUPPORTED_FORWARD_FLAVORS )
        message(FATAL_ERROR "-DFORWARD_FLAVOR option ('${FORWARD_FLAVOR}') was not valid choose from: ${SUPPORTED_FORWARD_FLAVORS}")
endif()

if (NOT MODEM_BUILD_DIMS IN_LIST SUPPORTED_MODEM_DIMS)
        message(FATAL_ERROR "-DMODEM_BUILD_DIMS option ('${MODEM_BUILD_DIMS}') was not valid choose from ${SUPPORTED_MODEM_DIMS}")
endif()

if (FG AND NOT FORWARD_FLAVOR STREQUAL SP2)
        message(FATAL_ERROR "You cannot build 'FG' with '${FORWARD_FLAVOR}' - Only SP2: -DFG=on -DFORWARD_FLAVOR=SP2")
endif()


# Source Code
add_subdirectory(f90)