submodule (ModelSpace:ModelSpaceIO) modelParam_IO_HDF5

#ifdef HDF5

use griddef
use hdf5
use ModEM_HDF5

implicit none

character (len=*), parameter :: METERS = 'meters'

contains

module subroutine write_modelParam_hdf5(m,cfile,comment)
    ! opens cfile on unit ioPrm, writes the object of
    ! type modelParam in HDF5/NetCDF4+ format, closes file

    type(modelParam_t), intent(in)	   :: m
    character(*), intent(in)             :: cfile
    character(*), intent(in), optional   :: comment

    character (len=512) :: hdf5_version_str

    integer(kind=HID_T) :: file_id 

    if (gridCoords .eq. SPHERICAL) then
        write(0,*) 'ERROR: HDF5 in spherical coordinates is not currently implemtned'
        call ModEM_abort()
    else
        write(0,*) 'Will be writing the model output in cartesian HDF5 format...'
    end if

    ! Open file here
    call ModEM_HDF5_create_file(cfile, H5F_ACC_TRUNC_F, file_id)

    call ModEM_HDF5_get_version(hdf5_version_str)
    call ModEM_HDF5_add_attr(file_id, '_NCProperties', "version=2,hdf5="// trim(hdf5_version_str))
    call ModEM_HDF5_add_attr(file_id, '_HDF5_version', trim(hdf5_version_str))
    call ModEM_HDF5_add_attr(file_id, '_nc3_strict', 1)

    ! Metadata attributes
    call ModEM_HDF5_add_attr(file_id, 'author_name', 'none')
    call ModEM_HDF5_add_attr(file_id, 'author_institution', 'none')
    call ModEM_HDF5_add_attr(file_id, 'author_email', 'none')
    call ModEM_HDF5_add_attr(file_id, 'author_url', 'none')
    call ModEM_HDF5_add_attr(file_id, 'id', 'none')
    call ModEM_HDF5_add_attr(file_id, 'model', 'none')
    call ModEM_HDF5_add_attr(file_id, 'reference', 'none')
    call ModEM_HDF5_add_attr(file_id, 'reference_pid', 'none')
    call ModEM_HDF5_add_attr(file_id, 'summary', 'none')
    call ModEM_HDF5_add_attr(file_id, 'title', 'none')
    call ModEM_HDF5_add_attr(file_id, 'data_revision', 'none')
    call ModEM_HDF5_add_attr(file_id, 'version', 'none')

    ! Write Model attributes
    call ModEM_HDF5_add_attr(file_id, 'Conventions', 'CF-1.0')
    call ModEM_HDF5_add_attr(file_id, 'Metadata_Conventions', 'Unidata Dataset Discovery v1.0')

    !! Origin
    call ModEM_HDF5_add_attr(file_id, 'model_origin_x', m % grid % ox)
    call ModEM_HDF5_add_attr(file_id, 'model_origin_y', m % grid % oy)
    call ModEM_HDF5_add_attr(file_id, 'model_origin_z', m % grid % oz)

    call ModEM_HDF5_add_attr(file_id, 'model_origin_location', 'model_center_ground_level')
    call ModEM_HDF5_add_attr(file_id, 'model_origin_description', 'defined in meters from upper southwest corner')
    call ModEM_HDF5_add_attr(file_id, 'model_primary_coords', 'xy')

    !! Rotation Angle
    call ModEM_HDF5_add_attr(file_id, 'model_rotation_angle', m % grid % rotdeg)
    call ModEM_HDF5_add_attr(file_id, 'model_rotation_units', 'degrees')

    ! Grid Dimension and refence
    call ModEM_HDF5_add_attr(file_id, 'grid_nx', m % grid % nx)

    call write_geometry_hdf5(file_id, m % grid)

    ! call write_gridSpacing_hdf5(file_id, m)

    call write_sigma_hdf5(file_id, m)

    ! Close file here
    call ModEM_HDF5_close_file(file_id)

end subroutine write_modelParam_hdf5

module subroutine write_geometry_hdf5(loc_id, grid)

    integer(kind=HID_T), intent(in)      :: loc_id
    type(grid_t), intent(in)	     :: grid

    integer                              :: Nx, Ny, NzEarth

    integer (kind=HID_T) :: nx_dspace_id, ny_dspace_id, nz_dspace_id
    integer (kind=HID_T) :: nx_dset_id, ny_dset_id, nz_dset_id
    integer (kind=HID_T) :: x_dspace_id, y_dspace_id, z_dspace_id
    integer (kind=HID_T) :: x_scale_id, y_scale_id, z_scale_id
    integer (kind=HID_T) :: x_dset_id, y_dset_id, z_dset_id

    integer :: hdferr

    Nx=grid%nx !this defines the length of the data array
    Ny=grid%ny
    NzEarth=grid%nz - grid%nzAir

    ! Write grid geometry definitions
    ! Create and write the X variable dimension
    call ModEM_HDF5_create_dataspace(rank(grid % xCenter), (/size(grid % xCenter, kind=HSIZE_T)/), x_dspace_id)
    call ModEM_HDF5_create_dataset(loc_id, 'x', H5T_NATIVE_DOUBLE, x_dspace_id, x_dset_id)
    call ModEM_HDF5_write_dataset(x_dset_id, grid % xCenter)
    call ModEM_HDF5_make_dataset_dimension(x_dset_id, hdferr)

    call ModEM_HDF5_add_attr(x_dset_id, 'long_name', 'Latitude; positive north')
    call ModEM_HDF5_add_attr(x_dset_id, 'standard_name', 'x')
    call ModEM_HDF5_add_attr(x_dset_id, 'units', METERS)
    call ModEM_HDF5_close_dataset(x_dset_id)
    call ModEM_HDF5_close_dataspace(x_dspace_id)

    ! Create and write the Y variable dimensions  
    call ModEM_HDF5_create_dataspace(rank(grid % yCenter), (/size(grid % yCenter, kind=HSIZE_T)/), y_dspace_id)
    call ModEM_HDF5_create_dataset(loc_id, 'y', H5T_NATIVE_DOUBLE, y_dspace_id, y_dset_id)
    call ModEM_HDF5_write_dataset(y_dset_id, grid % yCenter)
    call ModEM_HDF5_make_dataset_dimension(y_dset_id, hdferr)

    call ModEM_HDF5_add_attr(y_dset_id, 'long_name', 'Longitude; positive east')
    call ModEM_HDF5_add_attr(y_dset_id, 'standard_name', 'y')
    call ModEM_HDF5_add_attr(y_dset_id, 'units', METERS)
    call ModEM_HDF5_close_dataset(y_dset_id)
    call ModEM_HDF5_close_dataspace(y_dspace_id)

    ! Create and write the Z variable dimension
    call ModEM_HDF5_create_dataspace(rank(grid % zCenter), (/int(NzEarth, kind=HSIZE_T)/), z_dspace_id)
    call ModEM_HDF5_create_dataset(loc_id, 'z', H5T_NATIVE_DOUBLE, z_dspace_id, z_dset_id)
    call ModEM_HDF5_write_dataset(z_dset_id, grid % zCenter(grid % NzAir+1:grid%nz))
    call ModEM_HDF5_make_dataset_dimension(z_dset_id, hdferr)

    call ModEM_HDF5_add_attr(z_dset_id, 'long_name', 'depth below earth surface')
    call ModEM_HDF5_add_attr(z_dset_id, 'positive', 'down')
    call ModEM_HDF5_add_attr(z_dset_id, 'units', METERS)
    call ModEM_HDF5_add_attr(z_dset_id, 'standard_name', 'elevation')
    call ModEM_HDF5_close_dataset(z_dset_id)
    call ModEM_HDF5_close_dataspace(z_dspace_id)

end subroutine write_geometry_hdf5

subroutine write_gridSpacing_hdf5(file_id, m)
    ! Note: This subroutine is not currently used - we re-calculate
    ! the grid spacing from the cell centers and file attributes. 
    ! But perhaps this will be helpful in the future. We can save some
    ! space but not writing them out.

    integer(kind=HID_T), intent(in)       :: file_id 
    type(modelParam_t), intent(in)	      :: m

    ! local variables
    type(grid_t)                          :: grid
    type(rscalar)                         :: ccond
    character(80)                         :: paramType

    integer                               :: Nx, Ny, NzEarth
    integer (kind=HID_T) :: grid_spacing_group_id
    integer (kind=HID_T) :: root_group_id
    integer (kind=HID_T) :: dx_dset_id, dy_dset_id, dz_dset_id
    integer (kind=HID_T) :: dx_dspace_id, dy_dspace_id, dz_dspace_id 

    paramType = userParamType

    call getValue_modelParam(m,paramType,ccond)

    grid = ccond%grid
    Nx=grid%nx !this defines the length of the data array
    Ny=grid%ny
    NzEarth=grid%nz - grid%nzAir

    ! Write grid geometry definitions
!    call ModEM_HDF5_open_group(file_id, "/", root_group_id)

    call ModEM_HDF5_create_group(file_id, "GridSpacing", grid_spacing_group_id)
    ! Dx
    call ModEM_HDF5_create_dataspace(rank(grid % dx), (/size(grid % dx, kind=HSIZE_T)/), dx_dspace_id)
    call ModEM_HDF5_create_dataset(grid_spacing_group_id, 'Dx', H5T_NATIVE_DOUBLE, dx_dspace_id, dx_dset_id)
    call ModEM_HDF5_write_dataset(dx_dset_id, grid % dx)

    call ModEM_HDF5_close_dataset(dx_dset_id)
    call ModEM_HDF5_close_dataspace(dx_dspace_id)

    ! Dy
    call ModEM_HDF5_create_dataspace(rank(grid % dy), (/size(grid % dy, kind=HSIZE_T)/), dy_dspace_id)
    call ModEM_HDF5_create_dataset(grid_spacing_group_id, 'Dy', H5T_NATIVE_DOUBLE, dy_dspace_id, dy_dset_id)
    call ModEM_HDF5_write_dataset(dy_dset_id, grid % dy)


    call ModEM_HDF5_close_dataset(dy_dset_id)
    call ModEM_HDF5_close_dataspace(dy_dspace_id)

    ! Dz
    call ModEM_HDF5_create_dataspace(1, (/int(NzEarth, kind=HSIZE_T)/), dz_dspace_id)
    call ModEM_HDF5_create_dataset(grid_spacing_group_id, 'Dz', H5T_NATIVE_DOUBLE, dz_dspace_id, dz_dset_id)
    call ModEM_HDF5_write_dataset(dz_dset_id, grid % dz(grid%NzAir+1:grid%nz))

    call ModEM_HDF5_close_dataset(dz_dset_id)
    call ModEM_HDF5_close_dataspace(dz_dspace_id)

    call ModEM_HDF5_close_group(grid_spacing_group_id)
!    call ModEM_HDF5_close_group(root_group_id)

end subroutine write_gridSpacing_hdf5

!******************************************************************
subroutine write_sigma_hdf5(file_id, m, dimension_location, var_name)

    integer (kind=HID_T), intent(in)                  :: file_id
    type (modelParam_t), intent(in)	     :: m
    integer (kind=HID_T), optional :: dimension_location
    character (len=*), optional :: var_name

    ! local variables
    type(grid_t)                          :: grid
    type(rscalar)                         :: ccond, rho
    character(80)                         :: paramType =''

    integer                               :: Nx, Ny, NzEarth
    character (len=*), parameter :: prop = "sigma"

    integer (kind=HID_T) :: sigma_dset_id, sigma_dspace_id
    integer (kind=HID_T) :: x_dset_id, y_dset_id, z_dset_id
    integer :: i, j, k

    real(kind=prec), dimension(:,:,:), pointer :: sigma_out

    integer :: err 

    integer (kind=HID_T) :: dims_location_lcl
    character (len=512) :: var_name_lcl

    if (present(dimension_location)) then
        dims_location_lcl = dimension_location
    else
        dims_location_lcl = file_id
    end if

    if (present(var_name)) then
        var_name_lcl = trim(var_name)
    else
        var_name_lcl = trim(prop)
    end if


    ! Convert modelParam to natural log or log10 for output
    !paramType = userParamType
    paramType = 'LOG10' 
    call getValue_modelParam(m,paramType,ccond)


    grid = ccond%grid
    Nx=grid%nx !this defines the length of the data array
    Ny=grid%ny
    NzEarth=grid%nz - grid%nzAir
    

    call ModEM_HDF5_create_dataspace(3, (/int(Ny, kind=HSIZE_T), int(Nx, kind=HSIZE_T), int(NzEarth, kind=HSIZE_T)/), sigma_dspace_id)

    call ModEM_HDF5_create_dataset(file_id, trim(var_name_lcl), H5T_NATIVE_DOUBLE, sigma_dspace_id, sigma_dset_id)


    ! Attache the x, y, and z datasets (dimensions) to log_10_sigma.
    ! 
    ! Even without attaching the dimensions below, ncdump is able to recognize
    ! what dimensions go where; however, NOT doing this will leave off some 
    ! HDF5 dimension scale attributes, which may or maynot be necessary. 
    ! So, we perform this attachment, just to be sure other NetCDF applications 
    ! have an easy time reading this file.
    !
    ! For those interested, this will create a DIMENSION_LIST attribute on the HDF5 
    ! log_10_sigma variable and a REFERENCE_LIST on attribute on the dimensions (x, y, z).
    call ModEM_HDF5_attach_dim(dims_location_lcl, 'z', sigma_dset_id, 1)
    call ModEM_HDF5_attach_dim(dims_location_lcl, 'x', sigma_dset_id, 2)
    call ModEM_HDF5_attach_dim(dims_location_lcl, 'y', sigma_dset_id, 3)

    ! Assign attribute values
    call ModEM_HDF5_add_attr(sigma_dset_id, 'paramType', 'LOG10')
    call ModEM_HDF5_add_attr(sigma_dset_id, 'display_name', 'log(10) electrical conductivity, in Seimens/meter')

    ! Transpose
    allocate(sigma_out(Ny, Nx, NzEarth))
    do k = 1, NzEarth
        do j = 1, Ny
            do i = 1, Nx
                sigma_out(j, i, k) = ccond%v(i, j, k)
            end do
        end do
    end do

    ! Convert from conductivity to resitivity for output
    if ((index(paramType,'LOGE')>0) .or. (index(paramType,'LOG10')>0)) then
        sigma_out = - sigma_out
    else if (index(paramType,'LINEAR')>0) then
        sigma_out = ONE/sigma_out
    end if

    call ModEM_HDF5_add_attr(sigma_dset_id, 'long_name', 'electrical conductivity')
    call ModEM_HDF5_add_attr(sigma_dset_id, 'units', 'S/m')
    call ModEM_HDF5_add_attr(sigma_dset_id, 'missing_value', 99999.0_prec)

    call ModEM_HDF5_write_dataset(sigma_dset_id, sigma_out)

    call ModEM_HDF5_close_dataset(sigma_dset_id)
    call ModEM_HDF5_close_dataspace(sigma_dspace_id)

end subroutine write_sigma_hdf5

subroutine calculate_grid_spacing_from_centers(coord_origin, coord_centers, dVar)
    ! This might be better in GridCalc...

    implicit none

    real (kind=prec), intent(in) :: coord_origin 
    real (kind=prec), dimension(:), intent(in) :: coord_centers
    real (kind=prec), dimension(:), intent(out) :: dVar

    real (kind=prec) :: cursor

    integer :: i

    cursor = coord_origin
    do i = 1, size(coord_centers)
        dVar(i) = 2.0_prec * (coord_centers(i) - cursor)
        cursor = cursor + dVar(i)
    end do

end subroutine calculate_grid_spacing_from_centers

module subroutine read_modelParam_hdf5(grid,airLayers,m,cfile)

    ! opens cfile on unit ioPrm, reads the object of
    ! type modelParam in HDF5/NetCDF4+ format, closes file
    ! we can update the grid here, but the grid as an input is critical
    ! for setting the pointer to the grid in the modelParam

    type(grid_t), target, intent(inout)  :: grid
    type(airLayers_t), intent(inout)	   :: airLayers
    type(modelParam_t), intent(out)	   :: m
    character(*), intent(in)             :: cfile
    ! local variables

    integer(kind=HID_T)                 :: file_id
    type(rscalar)                        :: ccond
    character(80)                        :: paramType=''

    if (gridCoords .eq. SPHERICAL) then
        write(0,*) 'Will be reading the model input in spherical HDF5 format...'
    else
        write(0,*) 'Will be reading the model input in cartesian HDF5 format...'
    end if

    paramType = ''

    call ModEM_HDF5_open(cfile, file_id, H5F_ACC_RDONLY_F)

    ! First read grid geometry from HDF5 file
    call read_geometry_hdf5(file_id, grid, airLayers)
    call create_rscalar(grid, ccond, CELL_EARTH)
    call read_sigma_hdf5(file_id, grid, ccond, paramType)

    call ModEM_HDF5_close_file(file_id)

    ! Finally create the model parameter
    call create_modelParam(grid,paramType,m,ccond)

    ! In ModelSpace, save the user paramType for output
    userParamType = paramType

    ! ALWAYS convert modelParam to natural log for computations
    paramType = 'LOGE'
    call setType_modelParam(m,paramType)
    call deall_rscalar(ccond)

end subroutine read_modelParam_hdf5

subroutine read_geometry_hdf5(file_id, grid, airlayers)

    ! Given an open HDF5 file_id, and a unallocated grid and initalized airLayers,
    ! read in the cell centers of the HDF5 file it will also
    ! attempt to read the model origin (x, y, z) coordinates and then
    ! it will calculate dx, dy, dz based on the cell centers and model origin.
    !
    ! At the end of this routine, the grid will be fully allocated and
    ! initailzed and ready to perfrom computations.

    implicit none

    integer(kind=HID_T), intent(in) :: file_id
    type (grid_t) , intent(out) :: grid
    type (airLayers_t), intent(inout) :: airLayers

    integer(kind=HSIZE_T), dimension(1) :: dim1d

    character(80) :: someChar=''
    character(80) :: paramType=''
    real(kind=prec), dimension(:), allocatable :: xctr
    real(kind=prec), dimension(:), allocatable :: yctr
    real(kind=prec), dimension(:), allocatable :: zctr

    integer(kind=HSIZE_T) :: nx, ny, nz
    integer  :: grid_x, grid_y, grid_z, NzAir, i, j, k
    real(kind=prec) :: model_origin_x, model_origin_y, model_origin_z
    real(kind=prec) :: model_rotation_angle

    character (len=512) :: model_rotation_units

    integer (kind=HID_T) :: x_dset_id, y_dset_id, z_dset_id
    integer (kind=HID_T) :: x_dspace_id, y_dspace_id, z_dspace_id

    ! Read x
    call ModEM_HDF5_open_dataset(file_id, "x", x_dset_id)
    call ModEM_HDF5_get_dataspace(x_dset_id, x_dspace_id)
    call ModEM_HDF5_get_dataspace_size(x_dspace_id, nx)

    allocate(xctr(nx))

    call ModEM_HDF5_read_dataset(x_dset_id, H5T_NATIVE_DOUBLE, xctr)
    call ModEM_HDF5_close_dataset(x_dset_id)
    call ModEM_HDF5_close_dataspace(x_dspace_id)

    ! Read y
    call ModEM_HDF5_open_dataset(file_id, "y", y_dset_id)
    call ModEM_HDF5_get_dataspace(y_dset_id, y_dspace_id)
    call ModEM_HDF5_get_dataspace_size(y_dspace_id, ny)

    allocate(yctr(ny))

    call ModEM_HDF5_read_dataset(y_dset_id, H5T_NATIVE_DOUBLE, yctr)
    call ModEM_HDF5_close_dataset(y_dset_id)
    call ModEM_HDF5_close_dataspace(y_dspace_id)

    ! Read z
    call ModEM_HDF5_open_dataset(file_id, "z", z_dset_id)
    call ModEM_HDF5_get_dataspace(z_dset_id, z_dspace_id)
    call ModEM_HDF5_get_dataspace_size(z_dspace_id, nz)

    allocate(zctr(nz))

    call ModEM_HDF5_read_dataset(z_dset_id, H5T_NATIVE_DOUBLE, zctr)
    call ModEM_HDF5_close_dataset(z_dset_id)
    call ModEM_HDF5_close_dataspace(z_dspace_id)

    ! Allocate the grid
    call create_grid(int(nx, kind=4), int(ny, kind=4), airLayers % nz, int(nz, kind=4), grid)

    ! Read model origin and rotation attributes from the HDF5 file
    if (ModEM_HDF5_attr_exists(file_id, 'model_origin_x')) then
        call ModEM_HDF5_read_attr(file_id, 'model_origin_x', model_origin_x)
    else
        call errStop('Attribute model_origin_x does not exist in the HDF5 file.')
    end if

    if (ModEM_HDF5_attr_exists(file_id, 'model_origin_y')) then
        call ModEM_HDF5_read_attr(file_id, 'model_origin_y', model_origin_y)
    else
        call errStop('Attribute model_origin_y does not exist in the HDF5 file.')
    end if

    if (ModEM_HDF5_attr_exists(file_id, 'model_origin_z')) then
        call ModEM_HDF5_read_attr(file_id, 'model_origin_z', model_origin_z)
    else
        call errStop('Attribute model_origin_z does not exist in the HDF5 file.')
    end if

    if (ModEM_HDF5_attr_exists(file_id, 'model_rotation_angle')) then
        call ModEM_HDF5_read_attr(file_id, 'model_rotation_angle', model_rotation_angle)
    else
        call errStop('Attribute model_rotation does not exist in the HDF5 file.')
    end if

    if (ModEM_HDF5_attr_exists(file_id, 'model_rotation_units')) then
        call ModEM_HDF5_read_attr(file_id, 'model_rotation_units', model_rotation_units)
    else
        call errStop('Attribute model_rotation_units does not exist in the HDF5 file.')
    end if

    grid%rotdeg = model_rotation_angle

    ! Calculate grid spacing from cell centers
    call calculate_grid_spacing_from_centers(model_origin_x, xctr, grid%dx)
    call calculate_grid_spacing_from_centers(model_origin_y, yctr, grid%dy)
    call calculate_grid_spacing_from_centers(model_origin_z, zctr, &
        grid%dz(grid%nzAir+1:grid%nzAir+grid%nzEarth))

    ! Now, finish set up of the air layers structure using the grid
    call setup_airlayers(airLayers,grid)

    ! Finally, insert correct air layers in the grid
    call update_airlayers(grid,airLayers%Nz,airLayers%Dz)

    ! Preserve origin exactly as stored in file attributes.
    grid%ox = model_origin_x
    grid%oy = model_origin_y
    grid%oz = model_origin_z

    ! Recompute centers/edges using the stored grid origin.
    call setup_grid(grid)

end subroutine read_geometry_hdf5

subroutine read_sigma_hdf5(file_id, grid, sigma, paramType)

    ! Read in electrical conductivity from the file and store sigma in a rscalar variable
    ! note that this rscalar variable after reading this variable will still be conductivity

    implicit none

    integer (kind=HID_T), intent(in) :: file_id
    type (grid_t), intent(in) :: grid
    type (rscalar), intent(inout) :: sigma
    character(80), intent(out) :: paramType

    real (kind=prec), dimension(:,:,:), pointer :: sigma_read

    integer (kind=HID_T) :: dset_id, dspace_id
    integer :: i, j, k

    ! Read sigma from the HDF5 file
    ! TODO we will need a way to handle different variable names for sigma, as I have
    ! seen log10sigma, sigma, log10_sigma, etc.
    call ModEM_HDF5_open_dataset(file_id, "sigma", dset_id)
    call ModEM_HDF5_get_dataspace(dset_id, dspace_id)

    allocate(sigma_read(grid % ny, grid % nx, grid%nzEarth))
    call ModEM_HDF5_read_dataset(dset_id, H5T_NATIVE_DOUBLE, sigma_read)

    ! Read sigma into a temporary array, then transpose it into sigma%v
    ! Inverse mapping of write_sigma_hdf5 loop: (Ny, Nx, NzEarth) -> (Nx, Ny, NzEarth)
    do k = 1, grid % NzEarth
        do j = 1, grid % ny
            do i = 1, grid % nx
                sigma % v(i,j,k) = sigma_read(j,i,k)
            end do
        end do
    end do

    ! Set ParamType to be used when reading conductivity values
    ! Read the paramater type from the HDF5 file attributes, if available. If not, default to LOG10.
    if (ModEM_HDF5_attr_exists(dset_id, 'paramType')) then
        call ModEM_HDF5_read_attr(dset_id, 'paramType', paramType)
        call clean_null_term_string(paramType)
    else
        paramType = 'LOG10'
    end if

    if ((index(paramType,'LOGE')>0) .or. (index(paramType,'LOG10')>0)) then
        sigma%v = - sigma%v
    else if (index(paramType,'LINEAR')>0) then
        sigma%v = ONE/sigma%v
    end if

    call ModEM_HDF5_close_dataspace(dspace_id)
    call ModEM_HDF5_close_dataset(dset_id)

    deallocate(sigma_read)

end subroutine read_sigma_hdf5

#endif

end submodule ModelParam_IO_HDF5 
