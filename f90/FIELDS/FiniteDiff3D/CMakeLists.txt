# f90/FIELDS/FiniteDiff3D

target_sources(${MODEM_EXE}
    PRIVATE sg_boundary.f90
    PRIVATE sg_diff_oper.f90
    PRIVATE sg_scalar.f90
    PRIVATE sg_sparse_vector.f90
    PRIVATE sg_spherical.f90
    PRIVATE sg_vector.f90
)