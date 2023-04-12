!
!> Module briefing
!
module HDF5_3D_SG_GridReader
    !
    use hdf5
    use Constants
    use Grid3D_SG
    !
    type :: HDF5_3D_SG_GridReader_t
        !
        !> No derived properties
        !
        contains
            !
            procedure, public :: Read
            !
    end type HDF5_3D_SG_GridReader_t
    !
contains
    !
    !> No subroutine briefing
	!
    function Read(self, fileName) result(newGrid)
        !
        class(HDF5_3D_SG_GridReader_t), intent( in ) :: self
        character(*), intent( in ) :: fileName
        !
        class(Grid3D_SG_t), pointer :: newGrid
        !> ** HDF5 specific auxiliary variables
        type(C_PTR) :: hdf_cptr
        integer(HID_T) :: hdf_file_id
        integer(HID_T) :: hdf_grp_id
        integer(HID_T) :: hdf_att_id
        integer(HID_T) :: hdf_space_id
        integer :: hdf_error
        integer, target :: nx_t(1), ny_t(1), nz_t(1)
        integer(HSIZE_T) :: hdf_dim(1), hdf_maxdim(1)
        real, target :: hdf_grid_origin(3)
        real, target :: hdf_grid_rotation(1)
        !> Other variables
        integer :: nx, ny, nzair, nzEarth
        real( kind=prec ), dimension(:), allocatable, target :: dx, dy, dz
        real( kind=prec ) :: grid_origin(3)
        real( kind=prec ) :: grid_rotation
        
        !> Initializes HDF5 library
        call h5open_f(hdf_error)

        !> Open hdf5 file
        call h5fopen_f(fileName, H5F_ACC_RDONLY_F, hdf_file_id, hdf_error)

        ! NX        
        call h5aopen_by_name_f(hdf_file_id, 'Model', 'nx', &
                 &                                 hdf_att_id, hdf_error)
        call h5aget_space_f(hdf_att_id, hdf_space_id, hdf_error)
        hdf_cptr = C_LOC(nx_t)
        call h5aread_f(hdf_att_id, H5T_NATIVE_INTEGER, hdf_cptr, hdf_error)
        call h5aclose_f(hdf_att_id, hdf_error);
        nx = nx_t(1)

        ! NY
        call h5aopen_by_name_f(hdf_file_id, 'Model', 'ny', &
                 &                                 hdf_att_id, hdf_error)
        call h5aget_space_f(hdf_att_id, hdf_space_id, hdf_error)
        hdf_cptr = C_LOC(ny_t)
        call h5aread_f(hdf_att_id, H5T_NATIVE_INTEGER, hdf_cptr, hdf_error)
        call h5aclose_f(hdf_att_id, hdf_error);
        ny = ny_t(1)

        ! NZ
        call h5aopen_by_name_f(hdf_file_id, 'Model', 'nz', &
                 &                                 hdf_att_id, hdf_error)
        call h5aget_space_f(hdf_att_id, hdf_space_id, hdf_error)
        hdf_cptr = C_LOC(nz_t)
        call h5aread_f(hdf_att_id, H5T_NATIVE_INTEGER, hdf_cptr, hdf_error)
        call h5aclose_f(hdf_att_id, hdf_error);
        nzEarth = nz_t(1)

        !> Initialize and read in the grid
        nzAir = 10
        allocate(newGrid, source = Grid3D_SG_t(nx, ny, nzAir, nzEarth))
        !
        ! Read grid cells dimensions    **
        !
        allocate(dx(nx))
        allocate(dy(ny))
        allocate(dz(nzAir + nzEarth))

        dx = 0.; dy = R_ZERO; dz = 0.
        dx = 0.; dy = 0.; dz = 0.
        
        !> DX
        hdf_dim = 0; hdf_maxdim = 0;
        call h5aopen_by_name_f(hdf_file_id, 'Model', 'dx', &
                 &                                 hdf_att_id, hdf_error)
        call h5aget_space_f(hdf_att_id, hdf_space_id, hdf_error)
        call h5sget_simple_extent_dims_f(hdf_space_id, &
                 &                                                     hdf_dim, hdf_maxdim, hdf_error)
        hdf_cptr = C_LOC(dx)
        call h5aread_f(hdf_att_id, H5T_IEEE_F64LE, hdf_cptr, hdf_error)
        call h5aclose_f(hdf_att_id, hdf_error)

        !> DY
        hdf_dim = 0; hdf_maxdim = 0;
        call h5aopen_by_name_f(hdf_file_id, 'Model', 'dy', &
                 &                                 hdf_att_id, hdf_error)
        call h5aget_space_f(hdf_att_id, hdf_space_id, hdf_error)
        call h5sget_simple_extent_dims_f(hdf_space_id, &
                 &                                                     hdf_dim, hdf_maxdim, hdf_error)
        hdf_cptr = C_LOC(dy)
        call h5aread_f(hdf_att_id, H5T_IEEE_F64LE, hdf_cptr, hdf_error)
        call h5aclose_f(hdf_att_id, hdf_error)

        !> DZ
        hdf_dim = 0; hdf_maxdim = 0;
        call h5aopen_by_name_f(hdf_file_id, 'Model', 'dz', &
                 &                                 hdf_att_id, hdf_error)
        call h5aget_space_f(hdf_att_id, hdf_space_id, hdf_error)
        call h5sget_simple_extent_dims_f(hdf_space_id, &
                 &                                                     hdf_dim, hdf_maxdim, hdf_error)
        hdf_cptr = C_LOC(dz(nzAir + 1 : nzAir + nzEarth))
        call h5aread_f(hdf_att_id, H5T_IEEE_F64LE, hdf_cptr, hdf_error)
        call h5aclose_f(hdf_att_id, hdf_error)
        !
        !dx = hdf_dx
        !dy = hdf_dy        
        !dz(nzAir + 1 : nzAir + nzEarth) = hdf_dz
        !
        call newGrid%setCellSizes(dx, dy, dz)
        !
        !> ** Read grid origin **
        !
        call h5aopen_by_name_f(hdf_file_id, 'Model', 'Origin', &
                 &                                 hdf_att_id, hdf_error)
        call h5aget_space_f(hdf_att_id, hdf_space_id, hdf_error)
        hdf_cptr = C_LOC(hdf_grid_origin)
        call h5aread_f(hdf_att_id, H5T_IEEE_F32LE, hdf_cptr, hdf_error)
        call h5aclose_f(hdf_att_id, hdf_error)

        grid_origin = hdf_grid_origin
        call newGrid%setOrigin(grid_origin(1), grid_origin(2), grid_origin(3))

        !
        !> ** Read grid rotation **
        !
        call h5aopen_by_name_f(hdf_file_id, 'Model', 'Rotation', &
                 &                                 hdf_att_id, hdf_error)
        call h5aget_space_f(hdf_att_id, hdf_space_id, hdf_error)
        hdf_cptr = C_LOC(hdf_grid_rotation)
        call h5aread_f(hdf_att_id, H5T_IEEE_F32LE, hdf_cptr, hdf_error)
        call h5aclose_f(hdf_att_id, hdf_error)

        grid_rotation = hdf_grid_rotation(1)
        call newGrid%setGridRotation(grid_rotation)
        !
        !> Clean up
        !
        call h5close_f(hdf_error)
        
    end function Read
    
end module HDF5_3D_SG_GridReader
