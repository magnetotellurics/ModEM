cmake_minimum_required(VERSION 4.1.0)

project(ModEM VERSION 0
        DESCRIPTION "ModEM"
        LANGUAGES Fortran)

enable_language(Fortran)

option(BUILD_2D "Build 2D version of ModEM" OFF)
option(BUILD_3D "Build 3D version of ModEM" ON)
# Code to ensure both options are not specified at once

option(BUILD_MF "Build Matrix-Free version of ModEM" ON)
option(BUILD_SP "Build SP version of ModEM" OFF)
option(BUILD_SP2 "Build SP2 version of ModEM" OFF)
# Code to ensure both options are not specified at once

option(BUILD_MPI "Build MPI version of MPI" ON)

# Source Code
add_subdirectory(f90)