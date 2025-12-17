# 3D_MT/FWD/Mod2D CMakeLists.txt

target_sources(${MODEM_EXE}
    PRIVATE FwdTEmod.f90
    #PRIVATE FWDtest.f90
    PRIVATE FwdTMmod.f90
    PRIVATE WSfwd1Dmod.f90
    PRIVATE WSfwd2Dmod.f90
    PRIVATE WSfwd2Dpar.f90
    PRIVATE WSutils.f90
)