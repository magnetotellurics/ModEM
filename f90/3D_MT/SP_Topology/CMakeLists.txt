# f90/3D_MT/SP_Topology CMakeLists.txt

target_sources(${MODEM_EXE}
    PRIVATE spOpTools.f90
    PRIVATE spOpTopology_SG.f90
    PRIVATE vecTranslate.f90
)

if (NOT BUILD_SPHERICAL)
    target_sources(${MODEM_EXE}
        PRIVATE MetricElements_CSG.f90
    )
else()
    target_sources(${MODEM_EXE}
        PRIVATE MetricElements_SSG.f90
    )
endif()
