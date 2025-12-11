### f90 CMakeLists.txt

set (MY_MODS Mod3DMT.f90
            ModEM_Tester.f90
            UserCtrl.f90)

if (BUILD_2D)
    set(MODEM_EXE Mod2DMT)
    add_executable(${MODEM_EXE} Mod2DMT.f90)
    add_subdirectory(2D_MT)
endif()

if (BUILD_3D)
    set(MODEM_EXE Mod3DMT)
    add_executable(${MODEM_EXE} Mod3DMT.f90)
    add_subdirectory(3D_MT)
endif()

if (BUILD_MPI)
    add_subdirectory("MPI")
endif()


set_target_properties(${MODEM_EXE} PROPERTIES
    Fortran_PREPROCESS ON
)

message(STATUS "Building ModEM ExE: ${MODEM_EXE}")

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