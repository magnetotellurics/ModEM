# f90/3D_MT/ioMod CMakeLists.txt

if (BUILD_HDF5)
    target_sources(${MODEM_EXE}
        PRIVATE ioHDF5.f90
    )
else()
    target_sources(${MODEM_EXE}
        PRIVATE ioAscii.f90
    )
endif()
