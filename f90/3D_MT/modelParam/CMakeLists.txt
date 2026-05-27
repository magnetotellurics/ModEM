# f90/3D_MT/modelParam CMakeLists.txt

target_sources(${MODEM_EXE}
    PRIVATE ModelSpace.f90
    PRIVATE modelparamIO.f90
    PRIVATE modelParamIO/modelParam_IO_WS.f90
    PRIVATE modelParamIO/modelParam_IO_HDF5.f90
    PRIVATE modelParamIO/modelParam_IO_Binary.f90
    PRIVATE modelParamIO/modelParam_IO_Mackie.f90
)
