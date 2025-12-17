# f90/3D_MT/ioMod CMakeLists.txt

target_sources(${MODEM_EXE}
    PRIVATE ioAscii.f90
    #PRIVATE ioBinary.f90
    #PRIVATE ioHDF5.f90
)