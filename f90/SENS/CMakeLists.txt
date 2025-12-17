# f90/SENS CMakeLists.txt

target_sources(${MODEM_EXE}
    PRIVATE DataSens.f90
    PRIVATE DataSpace.f90
    PRIVATE SensComp.f90
    PRIVATE SensMatrix.f90
    PRIVATE SymmetryTest.f90
)