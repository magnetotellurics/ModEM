# f90/3D_MT/FWD_SP CMakeLists.txt

target_sources(${MODEM_EXE}
    PRIVATE EMsolve3D.f90
    PRIVATE modelOperator3D.f90
    PRIVATE solver.f90
    # PRIVATE spModOpTest.f90
)