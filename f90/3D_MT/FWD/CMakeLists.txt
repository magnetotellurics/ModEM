# 3D_MT/FWD CMakeLists.txt

add_subdirectory(Mod2d)

target_sources(${MODEM_EXE}
    PRIVATE boundary_ws.f90
    #PRIVATE boundary_wsS.f90
    PRIVATE nestedEM.f90
)

if (FORWARD_FORMULATION STREQUAL MF)
    target_sources(${MODEM_EXE}
        PRIVATE EMsolve3D.f90
        PRIVATE modelOperator3D.f90
        PRIVATE solver.f90
    )
endif()