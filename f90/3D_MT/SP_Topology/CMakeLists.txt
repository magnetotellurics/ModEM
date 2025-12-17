# f90/3D_MT/SP_Topology CMakeLists.txt

target_sources(${MODEM_EXE}
    PRIVATE MetricElements_CSG.f90
    PRIVATE spOpTools.f90
    PRIVATE spOpTopology_SG.f90
    PRIVATE vecTranslate.f90
)