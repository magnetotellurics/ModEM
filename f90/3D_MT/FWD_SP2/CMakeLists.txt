### f90/3D_MT/FWD_SP2 CMakeLists.txt

target_sources(${MODEM_EXE}
    #PRIVATE kernal_c.cu
    #PRIVATE kernal_c.hip
    #PRIVATE cudaFortMap.f90
    PRIVATE EMsolve3D.f90
    PRIVATE modelOperator3D.f90
    PRIVATE solver.f90
    #PRIVATE spModOpTest.f90
)