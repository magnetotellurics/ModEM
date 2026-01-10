# f90 CMakeLists.txt

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
    Fortran_MODULE_DIRECTORY obj
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

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    include(FindBLAS)
    if (NOT BLAS_FOUND)
        message(FATAL_ERROR "We could not find the BLAS Library")
    endif()

    include(FindLAPACK)
    if (NOT LAPACK_FOUND)
        message(FATAL_ERROR "We could not find the LAPACK library")
    endif()

    target_link_libraries(${MODEM_EXE} PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()

if (BUILD_GPU STREQUAL "CUDA")
    find_package(CUDAToolkit REQUIRED COMPONENTS cublas cudart cusparse nvml)
    target_link_libraries(${MODEM_EXE} PRIVATE CUDA::cublas CUDA::cudart CUDA::cusparse CUDA::nvml)
elseif (BUILD_GPU STREQUAL "HIP")
    target_link_libraries(${MODEM_EXE} PRIVATE "-lamdhip64 -lhipblas -lhipsparse -lstdc++")
endif()
