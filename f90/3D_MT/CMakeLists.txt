# 3D_MT CMakeLists.txt

if (FORWARD_FORMULATION STREQUAL MF)
    message(STATUS "Building BUILD_MF")
    add_subdirectory(FWD)
elseif(FORWARD_FORMULATION STREQUAL SP)
    message(STATUS "Building BUILD_SP")
    add_subdirectory(FWD)
    add_subdirectory(FWD_SP)
    add_subdirectory(SP_Topology)
elseif(FORWARD_FORMULATION STREQUAL SP2)
    message(STATUS "Building BUILD_SP2")
    add_subdirectory(FWD)
    add_subdirectory(FWD_SP2)
    add_subdirectory(SP_Topology)
endif()

add_subdirectory(DICT)
add_subdirectory(ioMod)
add_subdirectory(modelParam)

target_sources(${MODEM_EXE} 
    PRIVATE DataFunc.f90
    #PRIVATE DataIO_HDF5.f90 
    PRIVATE DataIO.f90
    PRIVATE EMfieldInterp.f90
    PRIVATE ForwardSolver.f90
    PRIVATE GridCalc.f90
    #PRIVATE GridCalcS.f90
    PRIVATE GridDef.f90
    PRIVATE Main.f90
    #PRIVATE mtInvSetup.f90
    PRIVATE SolnSpace.f90
    PRIVATE SolverSens.f90
    PRIVATE Sub_MPI.f90
)

if(BUILD_MPI)
    target_sources(${MODEM_EXE}
        PRIVATE Sub_MPI.f90
    )
endif()
