### f90 CMakeLists.txt

set (MY_MODS Mod3DMT.f90
            ModEM_Tester.f90
            UserCtrl.f90)

if (MODEM_BUILD_DIMS STREQUAL 2D)
    set(MODEM_EXE Mod2DMT_${FORWARD_FLAVOR})
    set(MODEM_MAIN_FILE Mod2DMT.f90)
    set(DIM_DIR 2D_MT)
elseif (MODEM_BUILD_DIMS STREQUAL 3D)
    set(MODEM_EXE Mod3DMT_${FORWARD_FLAVOR})
    set(MODEM_MAIN_FILE Mod3DMT.f90)
    set(DIM_DIR 3D_MT)
endif()

add_executable(${MODEM_EXE} ${MODEM_MAIN_FILE})
set_target_properties(${MODEM_EXE} PROPERTIES
    Fortran_PREPROCESS ON
)

if (BUILD_MPI)
    add_subdirectory("MPI")
endif()

add_subdirectory(${DIM_DIR})
add_subdirectory("FIELDS/FiniteDiff3D")
add_subdirectory("INV")
add_subdirectory("SENS")
add_subdirectory("UTILS")

target_sources(${MODEM_EXE} PRIVATE UserCtrl.f90)

include(FindBLAS)
if (BLAS_FOUND)
    message(STATUS " WE FOUND BLAS")
endif()

include(FindLAPACK)
if (LAPACK_FOUND)
    message(STATUS " WE FOUND LAPACK")
endif()

target_link_options(${MODEM_EXE} PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})