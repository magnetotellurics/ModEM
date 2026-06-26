submodule (ModelSpace:ModelSpaceIO) modelParam_IO_HDF5

#ifdef HDF5
use hdf5
use griddef

implicit none

integer(HID_T) :: attr_id, dset_id, dspace_id, atype_id, aspace_id ! file, data set, and dataspace handles

contains

	!******************************************************************
	subroutine write_modelParam_hdf5(m,cfile,comment)

	  ! opens cfile on unit ioPrm, writes the object of
	  ! type modelParam in HDF5/NetCDF4+ format, closes file

	  type(modelParam_t), intent(in)	   :: m
	  character(*), intent(in)             :: cfile
	  character(*), intent(in), optional   :: comment

      integer(kind=HID_T) :: file_id
      integer :: hdferr

	  if (gridCoords .eq. SPHERICAL) then
	  	  write(0,*) 'Will be writing the model output in spherical HDF5 format...'
	  else
	  	  write(0,*) 'Will be writing the model output in cartesian HDF5 format...'
	  end if

      ! Open file here
      call open_hdf5(cfile, file_id, H5F_ACC_TRUNC_F)
	  call write_geometry_hdf5(file_id, m)
	  call write_sigma_hdf5(file_id, m)

      ! Close file here
      call close_hdf5(file_id)

	end subroutine write_modelParam_hdf5


	!******************************************************************
	subroutine read_modelParam_hdf5(grid,airLayers,m,cfile)

      ! opens cfile on unit ioPrm, reads the object of
      ! type modelParam in HDF5/NetCDF4+ format, closes file
      ! we can update the grid here, but the grid as an input is critical
      ! for setting the pointer to the grid in the modelParam

      type(grid_t), target, intent(inout)  :: grid
	  type(airLayers_t), intent(inout)	   :: airLayers
      type(modelParam_t), intent(out)	   :: m
      character(*), intent(in)             :: cfile
      integer		                       :: istat
      ! local variables

      integer (kind=HID_T)                 :: file_id
      type(rscalar)                        :: rho, ccond
      character(80)                        :: someChar='',paramType=''
      integer                              :: Nx, Ny, NzEarth, NzAir, i, j, k

      INTEGER(HSIZE_T), DIMENSION(3) :: dim3d ! Datasets dimensions for 3D arrays
      REAL(KIND=8), DIMENSION(3), allocatable :: Sigma(:,:,:) ! Read data buffer for 3D Arrays
      integer           :: hdferr


	  if (gridCoords .eq. SPHERICAL) then
	  	  write(0,*) 'Will be reading the model input in spherical HDF5 format...'
	  else
	  	  write(0,*) 'Will be reading the model input in cartesian HDF5 format...'
	  end if

      call open_read_hdf5(cfile, file_id)

	  ! First read grid geometry from HDF5 file
	  call read_geometry_hdf5(file_id, grid, airLayers)

	  paramType = ''

	  ! Now, reopen HDF5 to read the conductivity
	    write(0,*) cfile,' is open and ready to read electrical conductivity'

		!!!!!!!!! READ SIGMA DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! Open Sigma dataset and get dataspace dimensions
		CALL h5dopen_f(file_id, "log10sigma", dset_id, hdferr)

		! The npoints will be a differents size for each dataset for sigma use ny,nx,and nz
		call h5dget_space_f(dset_id, dspace_id, hdferr)

		!allocate the space for the local variable
		Ny = grid%ny
		Nx = grid%nx
		NzEarth = grid%nzEarth
		allocate(Sigma(Ny,Nx,NzEarth), STAT = istat)

	   ! Read  grid geometries from Hdf5
		CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, Sigma, dim3d, hdferr)

        call close_hdf5(file_id)

	  ! Set ParamType to be used when reading conductivity values
	  paramType = 'LOG10'

	  ! Create and read in the resistivity values
      call create_rscalar(grid,ccond,CELL_EARTH)

 	  ! Save the conductivity
		 do i=1,Nx
			do j=1,Ny
				do k=1,NzEarth
			      ccond%v(i,j,k) = Sigma(j,i,k)
				end do
			end do
		end do
		write(0,*) ccond%v(1:10,1,1)

	  ! Finally create the model parameter
	  call create_modelParam(grid,paramType,m,ccond)

 	  ! In ModelSpace, save the user paramType for output
	  userParamType = paramType

	  ! ALWAYS convert modelParam to natural log for computations
	  paramType = 'LOGE'
	  call setType_modelParam(m,paramType)

	  ! now done with the rscalars, so deallocate
	  deallocate(Sigma, STAT = istat)
	  call deall_rscalar(ccond)

	end subroutine read_modelParam_hdf5

    !*******************************************************************
    !This subroutine will either create a new hdf5 file based on the name
    !given in the input or open an already exisiting file to be appended to
    subroutine open_hdf5(cfile, file_id, access_flags)
        character(*), intent(in)                :: cfile
        integer(kind=HID_T), intent(out)        :: file_id
        integer, optional, intent(in)           :: access_flags

        integer                                 :: hdferr
        logical                                 :: lexist

        CALL h5open_f(hdferr)

        call h5fcreate_f(cfile, H5F_ACC_TRUNC_F, file_id, hdferr)

        if (hdferr < 0) then
            call h5eprint_f(hdferr)
            write(0,*) 'ERROR: In open_hdf5()'
            call ModEM_Abort()
        end if

    end subroutine open_hdf5

    !**********************************************************************

    !Subroutine to open the hdf5 for reading
    subroutine open_read_hdf5(cfile, file_id)
        character(*), intent(in)                :: cfile
        integer(kind=HID_T), intent(out)        :: file_id
        integer                                 :: hdferr
        logical                                 :: lexist

        CALL h5open_f(hdferr)
        inquire(file = cfile, exist = lexist)
        if (lexist) then
            CALL h5fopen_f(cfile, H5F_ACC_RDONLY_F, file_id, hdferr)
        else
            write(0,*) 'No HDF5 file to read'
        end if

    end subroutine open_read_hdf5

    !**********************************************************************
    subroutine close_hdf5(file_id)

        integer(kind=HID_T), intent(in)         :: file_id
        integer                                 :: hdferr

        ! TODO: Do we need to close all groups 
        CALL h5fclose_f(file_id, hdferr)

    end subroutine close_hdf5

    !**********************************************************************
    subroutine write_hdf5_attr(a_type, attr_name, attr_obj, path_id)
       INTEGER(HID_T), intent(in)                   :: path_id
        character(len=*), intent(in)                :: attr_name
        type(custom_att),pointer,  intent(in)       :: attr_obj
        character(len=3), intent(in)                :: a_type
        integer                                     :: hdferr
        INTEGER(HSIZE_T), DIMENSION(1)              :: dimsc = 1
        INTEGER(HSIZE_T)                            :: attrlen  ! Length of the attribute string

        ! Determine the HDF5 datatype of the value

          select case (a_type)
             case ('str')
                attrlen = len(attr_obj%att_string)
                CALL h5screate_simple_f(1, dimsc, aspace_id, hdferr)
                CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
                CALL h5tset_size_f(atype_id, attrlen, hdferr)
                CALL h5acreate_f(path_id, attr_name , atype_id, aspace_id, attr_id, hdferr)
                CALL h5awrite_f(attr_id, atype_id, attr_obj%att_string, dimsc, hdferr)
             case ('int')
                ! atype_id = H5T_NATIVE_INTEGER
                CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, hdferr)
                call H5Screate_simple_f(1, dimsc, aspace_id, hdferr)
                CALL h5acreate_f(path_id, attr_name , atype_id, aspace_id, attr_id, hdferr)
                CALL h5awrite_f(attr_id, atype_id, attr_obj%att_int, dimsc, hdferr)
             case('dbl')
                CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, hdferr)
                call H5Screate_simple_f(1, dimsc, aspace_id, hdferr)
                CALL h5acreate_f(path_id, attr_name , atype_id, aspace_id, attr_id, hdferr)
                CALL h5awrite_f(attr_id, atype_id, attr_obj%att_real, dimsc, hdferr)
                call h5aclose_f(attr_id, hdferr)
                call h5sclose_f(aspace_id, hdferr)
           end select

      end subroutine write_hdf5_attr

    !**********************************************************************
    function read_hdf5_attr(att_name) result( att_data)

        character(*), intent(in)                :: att_name
        integer                                 :: hdferr
        INTEGER(HSIZE_T), DIMENSION(1)              :: maxdims !Read buffer dimension
        INTEGER(HSIZE_T)  :: attrlen  ! Length of the attribute string
        INTEGER(SIZE_T)   :: dimsc(1) ! Scalar or single value string
        REAL(dp)             :: att_data !attribute read buffer
        integer (kind=HID_T) :: group_id

        call h5aopen_f(group_id, att_name, attr_id, hdferr)

        ! Get dataspace and allocate memory for read buffer.
        call H5Aget_space_f(attr_id, dspace_id, hdferr)
        call H5Sget_simple_extent_dims_f(dspace_id, dimsc, maxdims, hdferr)

        !Read the attribute
        call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, att_data, dimsc, hdferr )
        call h5aclose_f( attr_id, hdferr)

    end function read_hdf5_attr


    !**********************************************************************

    subroutine write_geometry_hdf5(file_id, m)

        integer(kind=HID_T), intent(in)      :: file_id
        type(modelParam_t), intent(in)	     :: m
        ! type(airLayers_t), intent(in)	     :: airLayers
        integer           		             :: istat
        ! local variables
        type(custom_att), target             :: attr_obj
        type(grid_t)                         :: grid
        integer                              :: Nx, Ny, NzEarth,  i, j, k

        character(len=2)                            :: xy


        INTEGER(8), DIMENSION(1) :: dim1d ! Datasets dimensions for 1D arrays
        integer           :: hdferr, ii
        integer (kind=HID_T) :: group_id

        grid = m%grid

        Nx=grid%nx !this defines the length of the data array
        Ny=grid%ny
        NzEarth=grid%nz - grid%nzAir

        ! nzAir = airLayers%Nz

        ! ! Write grid geometry definitions
        CALL h5gopen_f(file_id, "/", group_id, hdferr)

        ! Assign attribute values
        attr_obj%att_real = -grid%ox
        call write_hdf5_attr('dbl','model_origin_x', attr_obj, group_id)
        attr_obj%att_real = -grid%oy
        call write_hdf5_attr('dbl','model_origin_y', attr_obj, group_id)
        attr_obj%att_real = -grid%oz
        call write_hdf5_attr('dbl','model_origin_z', attr_obj, group_id)
        attr_obj%att_real = 0.0
        call write_hdf5_attr('dbl','model_rotation_angle', attr_obj, group_id)

        attr_obj%att_string = 'xy'
        call write_hdf5_attr('str','model_primary_coords', attr_obj, group_id)
        attr_obj%att_string = 'degrees'
        call write_hdf5_attr('str','model_rotation_units', attr_obj, group_id)

        dim1d(1) = nx   !might need to adjust this so that each celsize has the correct dimension length
        ! write the linear data array for NX
        CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)
        CALL h5dcreate_f(group_id, 'x', H5T_NATIVE_DOUBLE , dspace_id, dset_id, hdferr)
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , grid%xCenter, dim1d, hdferr)

        !write attributes for x dataset
        attr_obj%att_string = 'DIMENSION_SCALE'
        call write_hdf5_attr('str','CLASS', attr_obj, dset_id)
        attr_obj%att_string = 'x'
        call write_hdf5_attr('str','NAME', attr_obj, dset_id)
        attr_obj%att_string = 'Latitude; positive north'
        call write_hdf5_attr('str','long_name', attr_obj, dset_id)
        attr_obj%att_string = 'x'
        call write_hdf5_attr('str','standard_name', attr_obj, dset_id)
        attr_obj%att_string = 'meters'
        call write_hdf5_attr('str','units', attr_obj, dset_id)

        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(dspace_id, hdferr)

        ! write the linear data array for NY
        dim1d(1) = ny
        CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)
        CALL h5dcreate_f(group_id, 'y', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, grid%yCenter, dim1d, hdferr)

        !write attributes for y dataset
        attr_obj%att_string = 'DIMENSION_SCALE'
        call write_hdf5_attr('str','CLASS', attr_obj, dset_id)
        attr_obj%att_string = 'y'
        call write_hdf5_attr('str','NAME', attr_obj, dset_id)
        attr_obj%att_string = 'Latitude; positive east'
        call write_hdf5_attr('str','long_name', attr_obj, dset_id)
        attr_obj%att_string = 'y'
        call write_hdf5_attr('str','standard_name', attr_obj, dset_id)
        attr_obj%att_string = 'meters'
        call write_hdf5_attr('str','units', attr_obj, dset_id)

        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(dspace_id, hdferr)

        ! write the linear data array for NZ
        dim1d(1) = nzEarth
        CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)
        CALL h5dcreate_f(group_id, 'z', H5T_NATIVE_DOUBLE , dspace_id, dset_id, hdferr)
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , grid%zCenter(grid%NzAir+1:NzEarth+grid%NzAir), dim1d, hdferr)

        !write attributes for z dataset
        attr_obj%att_string = 'DIMENSION_SCALE'
        call write_hdf5_attr('str','CLASS', attr_obj, dset_id)
        attr_obj%att_string = 'z'
        call write_hdf5_attr('str','NAME', attr_obj, dset_id)
        attr_obj%att_string = 'depth below earth surface'
        call write_hdf5_attr('str','long_name', attr_obj, dset_id)
        attr_obj%att_string = 'down'
        call write_hdf5_attr('str','positive', attr_obj, dset_id)
        attr_obj%att_string = 'meters'
        call write_hdf5_attr('str','units', attr_obj, dset_id)

        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(dspace_id, hdferr)

    end subroutine write_geometry_hdf5

    !******************************************************************
    !write code for the nodes
    subroutine write_gridSpacing_hdf5(m,file_id)
        type(modelParam_t), intent(in)	      :: m
        integer(kind=HID_T), intent(in)       :: file_id 

        ! local variables
        type(grid_t)                          :: grid
        type(rscalar)                         :: rho,ccond
        character(80)                         :: paramType

        INTEGER(HSIZE_T), DIMENSION(1)        :: dim1d ! Datasets dimensions for 1D arrays
        integer                               :: hdferr, istat
        integer                               :: Nx, Ny, NzEarth, i, j, k
        integer (kind=HID_T) :: group_id

        paramType = userParamType

        call getValue_modelParam(m,paramType,ccond)

        grid = ccond%grid
        Nx=grid%nx !this defines the length of the data array
        Ny=grid%ny
        NzEarth=grid%nz - grid%nzAir

        CALL h5gcreate_f(file_id, "GridSpacing", group_id, hdferr)
        dim1d(1) = Nx
        ! write the linear data array for NodesX
        CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)
        CALL h5dcreate_f(group_id, 'Dx', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, grid%Dx, dim1d, hdferr)

        ! write the linear data array for NodesY
        dim1d(1) = Ny
        CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)
        CALL h5dcreate_f(group_id, 'Dy', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, grid%Dy, dim1d, hdferr)

        ! write the linear data array for NodesZ
        dim1d(1) = NzEarth
        CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)
        CALL h5dcreate_f(group_id, 'Dz', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, grid%Dz, dim1d, hdferr)

    end subroutine write_gridSpacing_hdf5

    !******************************************************************
    subroutine write_sigma_hdf5(file_id, m)
        integer(kind=HID_T), intent(in)                  :: file_id
        type(modelParam_t), intent(in)	     :: m

        ! local variables
        type(custom_att), target        :: attr_obj
        type(grid_t)                          :: grid
        type(rscalar)                         :: rho,ccond
        character(80)                         :: paramType =''

        INTEGER(HSIZE_T), DIMENSION(3)        :: dim3d ! Datasets dimensions for 1D arrays
        integer                               :: hdferr, istat
        integer                               :: Nx, Ny, NzEarth, i, j, k
        CHARACTER(LEN=10), parameter :: prop = "log10sigma"
        real(8)                               :: origin_x

        ! Convert modelParam to natural log or log10 for output
        !paramType = userParamType
        paramType = 'LOG10'
        call getValue_modelParam(m,paramType,ccond)

        grid = ccond%grid
        Nx=grid%nx !this defines the length of the data array
        Ny=grid%ny
        NzEarth=grid%nz - grid%nzAir

        dim3d(1) = Ny
        dim3d(2) = Nx
        dim3d(3) = NzEarth
        CALL h5screate_simple_f(3, dim3d, dspace_id, hdferr)
        call h5dcreate_f(file_id, prop, H5T_NATIVE_DOUBLE , dspace_id, dset_id, hdferr)

        ! Assign attribute values
        attr_obj%att_string = 'log(10) electrical conductivity, in S/m'
        call write_hdf5_attr('str','display_name', attr_obj, dset_id)

        attr_obj%att_string = 'electrical conductivity'
        call write_hdf5_attr('str','long_name', attr_obj, dset_id)

        attr_obj%att_string = 'S/m'
        call write_hdf5_attr('str','units', attr_obj, dset_id)

        attr_obj%att_real = 99999.0
        call write_hdf5_attr('dbl','missing_value', attr_obj, dset_id)

        ! Write the resistivity
        CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ccond%v, dim3d, hdferr)

    end subroutine write_sigma_hdf5

    !******************************************************************
    subroutine read_geometry_hdf5(file_id, grid, airlayers)

        integer(kind=HID_T), intent(in)                 :: file_id
        type (grid_t) , intent(inout)		:: grid
        type(airLayers_t), intent(inout)    :: airLayers

        integer		                         :: istat

        INTEGER(HSIZE_T), DIMENSION(1) :: dim1d ! Datasets dimensions for 1D arrays
        integer           :: hdferr

        character(80)                           :: someChar='',paramType=''
        REAL(KIND=8), DIMENSION(1), allocatable :: xctr(:), yctr(:), zctr(:) !Read data buffers for 1D Arrays
        INTEGER(HSIZE_T)                        :: nx, ny, nz !create memory buffers for X, Y, Z, and sigma to be used for getting the # of elements in each array
        integer                                 :: grid_x, grid_y, grid_z, NzAir, i, j, k
        real(8)                                    :: x, y, z
        real(8)                                    :: origin(3)

        integer (kind=HID_T) :: group_id

        !!!!!!!!! READ X DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Open X dataset and get dataspace dimensions
        call h5dopen_f(file_id, "x", dset_id, hdferr)

        ! The npoints will be a differents size for each dataset
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_npoints_f(dspace_id, nx, hdferr)

        !allocate the space for the local variable
        allocate(xctr(nx), STAT = istat)

        ! Read  grid geometries from Hdf5
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xctr, dim1d, hdferr)

        call h5dopen_f(file_id, "y", dset_id, hdferr)

        ! The npoints will be a differents size for each dataset
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_npoints_f(dspace_id, ny, hdferr)

        !allocate the space for the local variable
        allocate(yctr(ny), STAT = istat)

       ! Read  grid geometries from Hdf5
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, yctr, dim1d, hdferr)

        !!!!!!!!! READ Z DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Open Z dataset and get dataspace dimensions
        call h5dopen_f (file_id, "z", dset_id, hdferr)

        ! The npoints will be a differents size for each dataset
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        call h5sget_simple_extent_npoints_f(dspace_id, nz, hdferr)

        !allocate the space for the local variable
        allocate(zctr(nz), STAT = istat)

        ! Read  grid geometries from Hdf5
        CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, zctr, dim1d, hdferr)

        ! Setup grid
        grid_x = nx
        grid_y = ny
        grid_z = nz
        nzAir = airLayers%Nz
        call create_grid(grid_x, grid_y, nzAir, grid_z, grid)

        write(0,*) nzAir
        !Read Origin Attributes from the file to be used in X, Y, and Z centers

        call h5gopen_f (file_id, "/", group_id, hdferr) !need to grab group id for attribute function to work porperly
        x = read_hdf5_attr('model_origin_x')
        y = read_hdf5_attr('model_origin_y')
        z = read_hdf5_attr('model_origin_z')

        origin(1) = -x
        origin(2) = -y
        origin(3) = -z

        x = -x
        y = -y
        z = -z
        ! x = origin(1)
        do  i = 1, size(xctr)
         grid%dx(i) = 2*(xctr(i)-x)
         x = x + grid%dx(i)
        end do
        ! y= origin(2)
        do  i = 1, size(yctr)
         grid%dy(i) = 2*(yctr(i)-y)
         y = y + grid%dy(i)
        end do
        ! z = origin(3)
        do  i = 1, size(zctr)
         grid%dz(nzAir + i) = 2*(zctr(i)-z)
         z = z + grid%dz(nzAir +i)
        end do

        ! Finally, insert correct air layers in the grid and run setup_grid
        call setup_airlayers(airLayers,grid)
        call update_airlayers(grid,nzAir,airLayers%Dz)
        call setup_grid(grid, origin)

    end subroutine read_geometry_hdf5

#endif

end submodule ModelParam_IO_HDF5 
