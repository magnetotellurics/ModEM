# f90 CMakeLists.txt

if (MODEM_BUILD_DIMS STREQUAL 2D)
    set(MODEM_EXE Mod2DMT_${FORWARD_FORMULATION})
    set(MODEM_MAIN_FILE Mod2DMT.f90)
    set(DIM_DIR 2D_MT)
elseif (MODEM_BUILD_DIMS STREQUAL 3D)
    set(MODEM_EXE Mod3DMT_${FORWARD_FORMULATION})
    set(MODEM_MAIN_FILE Mod3DMT.f90)
    set(DIM_DIR 3D_MT)
endif()

add_executable(${MODEM_EXE} ${MODEM_MAIN_FILE})
set_target_properties(${MODEM_EXE} PROPERTIES
    Fortran_PREPROCESS ON
    Fortran_MODULE_DIRECTORY obj
)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
        ${MODEM_EXE}
        PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-cpp>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>
        $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>
        $<$<CONFIG:RELEASE>:-O3>
        $<$<CONFIG:DEBUG>:-O0>
        $<$<CONFIG:DEBUG>:-g>
        $<$<CONFIG:DEBUG>:-fbacktrace>
        $<$<CONFIG:DEBUG>:-fbounds-check>
    )
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    enable_language(C)
    find_library(MKL_INTEL_LP64 mkl_intel_lp64 HINTS $ENV{MKLROOT}/lib/intel64 ...)
    find_library(MKL_INTEL_THREAD mkl_intel_thread HINTS $ENV{MKLROOT}/lib/intel64 ...)
    find_library(MKL_CORE_LIB mkl_core HINTS $ENV{MKLROOT}/lib/intel64 ...)
    find_package(Threads REQUIRED)

    target_compile_options(
        ${MODEM_EXE}
        PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>-fpp>
        $<$<COMPILE_LANGUAGE:Fortran>-parallel>
        $<$<CONFIG:RELEASE>:-O3>
        $<$<CONFIG:RELEASE>:-w>
        $<$<CONFIG:DEBUG>:-O0>
        $<$<COMPILE_LANGUAGE:Fortran>:-fp-stack-check>
        $<$<COMPILE_LANGUAGE:Fortran>:-g>
        $<$<COMPILE_LANGUAGE:Fortran>:-traceback>
        $<$<COMPILE_LANGUAGE:Fortran>:-debug extend>
        $<$<COMPILE_LANGUAGE:Fortran>:-check bounds>
        $<$<COMPILE_LANGUAGE:Fortran>:-heap-array>
    )
    target_link_libraries(
        ${MODEM_EXE}
        PRIVATE
        ${MKL_INTEL_LP64}
        ${MKL_INTEL_THREAD}
        ${MKL_CORE_LIB}
        Threads::Threads
    )
    target_link_options(
        ${MODEM_EXE}
        PRIVATE
        -qopenmp
    )
else()
    message(FATAL_ERROR "Unknown/Unsupported Fortran Compiler: '${CMake_Fortran_COMPILER_ID}' - Supported compilers are < GNU | Ifort")
endif()


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
