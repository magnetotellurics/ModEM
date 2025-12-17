# f90/INV CMakeLists.txt

target_sources(${MODEM_EXE}
    PRIVATE DCG.f90
    PRIVATE INVcore.f90
    PRIVATE LBFGS.f90
    PRIVATE NLCG.f90
)