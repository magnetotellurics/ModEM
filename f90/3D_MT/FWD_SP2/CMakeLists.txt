# f90/3D_MT/FWD_SP2 CMakeLists.txt

target_sources(${MODEM_EXE}
    PRIVATE EMsolve3D.f90
    PRIVATE modelOperator3D.f90
    PRIVATE solver.f90
)

if (BUILD_GPU STREQUAL CUDA)
    target_sources(${MODEM_EXE}
        PRIVATE kernel_c.cu
        PRIVATE cudaFortMap.f90
    )
elseif (BUILD_GPU STREQUAL HIP)
    target_sources(${MODEM_EXE}
        PRIVATE kernel_c.hip
        PRIVATE hipFortMap.f90
    )
endif()
