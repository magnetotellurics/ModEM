### DICT CMakeLists.txt

target_sources(${MODEM_EXE}
    PRIVATE dataTypes.f90
    PRIVATE receivers.f90
    PRIVATE transmitters.f90
    #PRIVATE txTypes.f90
)