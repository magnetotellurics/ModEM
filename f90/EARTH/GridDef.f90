! *****************************************************************************
! * 3D global spherical grid definitions and basic subroutines

module GridDef

  ! all the modules being used are being listed explicitly (no
  ! inheritance)
  use math_constants
  use file_units
  use utilities
  implicit none

  ! Don't forget to overload the '=' sign: depending on the compiler, might
  ! run into trouble with the default assignment, since that sometimes doesn't
  ! copy allocatable or pointer arrays
  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE copy_grid
  END INTERFACE

  ! Possible grid types for EMfield, storing the intention of use for types
  ! such as cvector, cscalar, rvector, rscalar, sparsevecc.
  character(len=80), parameter		:: FACE = 'FACE'
  character(len=80), parameter		:: EDGE = 'EDGE'
  character(len=80), parameter		:: CENTER = 'CELL'
  character(len=80), parameter		:: CORNER = 'NODE'
  character(len=80), parameter		:: CELL_EARTH = 'CELL EARTH'

  ! ***************************************************************************
  ! * BOP
  ! type grid_param consists of parameters that define the basic grid geometry
  ! used for three dimensional numerical modeling
  type :: grid_t
  ! * EOP

     ! Grid coordinate system; important - used in EMfield
     character (len=80)			:: coords = Spherical

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nzEarth is the number of earth layers used in the model parametrization
     ! nzCrust is the number of layers to define the thinsheet
     ! nzAir is the number of air layers
     ! nz is grid dimension (number of cells) in the z-direction:
     ! nz = nzAir + nzCrust + nzEarth
	 ! If no crust information is provided, nzCrust=0
     integer               :: nx, ny, nz, nzEarth, nzCrust, nzAir

     ! Grid geometry:
	 ! No grid geometry is currently defined in the grid type definition
	 ! This may (and should be) added in the future, when Spherical and
	 ! Cartesian codes are united.

     ! Book-keeping on cumulative distances
	 ! No cumulative distances are currently defined in the grid type

 	 ! storing the (spherical) grid in Randie Mackie's format
	 ! when allocated, dimensions will be x(nx), y(ny+1), z(nz+1)
	 ! (this should only be used for input, set up and output)
	 real(8), pointer, dimension(:)			    :: x,y,z

	 ! storing the cell node coordinates, in radians (ph,th) and km (r)
	 ! dimensions are ph(nx+1), th(ny+1), r(nz+1)
	 ! the nx+1'st longitude ph(nx+1) = 360.0 * d2r (for interpolation)
	 real(8), pointer, dimension(:)  			:: ph,th,r
	 ! finally, storing the distances between cell nodes
     real(8), pointer, dimension(:)             :: dp,dt,dr

     ! allocated:  .true.  all the arrays have been allocated
     logical		                             :: allocated = .false.

  end type grid_t

  public         :: create_grid, deall_grid, copy_grid

Contains

     !************************************************************************
     subroutine create_grid(nx,ny,nz,grid)
       !  creates finite differences grid_t structure of
       !  size Nz x Ny, allocates arrays
       !
       implicit none
       integer, intent(in)		        :: nx,ny,nz
       type (grid_t), intent(inout)	    :: grid
       integer                          :: istat

       grid%nx = nx
       grid%ny = ny
       grid%nz = nz
       allocate(grid%x(nx), grid%y(ny+1), grid%z(nz+1), STAT=istat)
       allocate(grid%ph(nx+1), grid%th(ny+1), grid%r(nz+1), STAT=istat)
       allocate(grid%dp(nx), grid%dt(ny), grid%dr(nz), STAT=istat)
       grid%allocated = .true.

     end subroutine create_grid

     !************************************************************************
     subroutine deall_grid(grid)
       !  deallocates finite differences grid_t structure
       !
       implicit none
       type (grid_t), intent(inout)		:: grid
       integer                          :: istat

       if (associated(grid%x)) deallocate(grid%x, STAT=istat)
       if (associated(grid%y)) deallocate(grid%y, STAT=istat)
       if (associated(grid%z)) deallocate(grid%z, STAT=istat)
       if (associated(grid%ph)) deallocate(grid%ph, STAT=istat)
       if (associated(grid%th)) deallocate(grid%th, STAT=istat)
       if (associated(grid%r)) deallocate(grid%r, STAT=istat)
       if (associated(grid%dp)) deallocate(grid%dp, STAT=istat)
       if (associated(grid%dt)) deallocate(grid%dt, STAT=istat)
       if (associated(grid%dr)) deallocate(grid%dr, STAT=istat)
       grid%allocated = .false.

     end subroutine deall_grid

     !************************************************************************
     subroutine copy_grid(gridOut,gridIn)
       !  overloads the '=' sign
       !
       implicit none
       type (grid_t), intent(in)		:: gridIn
       type (grid_t), intent(out)		:: gridOut
       integer		        			:: nx,ny,nz

       nx = gridIn%nx
       ny = gridIn%ny
       nz = gridIn%nz

       call deall_grid(gridOut)
       call create_grid(nx,ny,nz,gridOut)

       gridOut%coords = gridIn%coords
       gridOut%nzEarth = gridIn%nzEarth
       gridOut%nzCrust = gridIn%nzCrust
       gridOut%nzAir = gridIn%nzAir

       gridOut%x = gridIn%x
       gridOut%y = gridIn%y
       gridOut%z = gridIn%z
       gridOut%ph = gridIn%ph
       gridOut%th = gridIn%th
       gridOut%r = gridIn%r
       gridOut%dp = gridIn%dp
       gridOut%dt = gridIn%dt
       gridOut%dr = gridIn%dr

       gridOut%allocated = gridIn%allocated

     end subroutine copy_grid

  ! ***************************************************************************
  ! * read_grid reads the modelfile cfile to store the grid only.
  ! * Traditionally, for global spherical grid, assume the following directions:
  ! * x -> phi (longitude, varies from 0 to 360)
  ! * y -> theta (co-latitude = 90 - latitude; varies from 0 to 180)
  ! * z -> -r (radius from Earth's centre, r=2871.0 at CMB,
  ! *                                        6371.0 at Earth/air interface)
  ! *
  ! * For consistency with the original Randy Mackie, 06-27-85, 3-D model, and
  ! * also for consistency with the current forward solver subroutines, we are
  ! * keeping the following grid structure in this forward solver:
  ! *      line 1: dimensions of model (x,y,z)
  ! *      line 2: x(*) in degrees (interval)
  ! *      line 3: y(*) in degrees (position from n-pole)
  ! *      line 4: z(*) in km (distance from center of the earth, decreasing)

  subroutine read_grid(grid,cfile)

    character(*), intent(in)                        :: cfile
    type (grid_t) , intent(out)                     :: grid
    integer                                         :: nx,ny,nz,nzAir,nzCrust,nzEarth
    real(8), dimension(:), allocatable              :: x,y,z,ph,th,r
    integer                                         :: ios,istat,i
    logical                                         :: exists

    inquire(FILE=trim(cfile),EXIST=exists)
    if(exists) then
      open(ioGrd,file=cfile,status='old',form='formatted',iostat=ios)
      write(6,*) node_info,'Reading from the grid file ',trim(cfile)
    else
      write(0,*) node_info,'Error: (read_grid) input file does not exist'
      stop
    end if

    read(ioGrd,*) nx,ny,nzAir,nzCrust,nzEarth

    nz = nzAir + nzCrust + nzEarth

    ! model grid and resistivity memory allocation
    allocate(x(nx+1),y(ny+1),z(nz+1), STAT=istat)

    ! read the x-intervals and y in degrees, z in km from the top of the air layer down
    read(ioGrd,*) x(1:nx)
    read(ioGrd,*) y(1:ny+1)
    read(ioGrd,*) z(1:nz+1)

    close(ioGrd)

    ! round vertical grid values to nearest meter to prevent precision errors
    do i=1,nz+1
      z(i)=nearest_meter(z(i))
    end do

    ! check that the Earth boundary in the grid file matches the internal Earth radius
    if ((clean(z(nzAir+1)) > EARTH_R + EPS_GRID) .or. (clean(z(nzAir+1)) < EARTH_R - EPS_GRID)) then
        write(6,*) node_info,'Warning: Earth radius in the grid file ',z(nzAir+1),' does not match EARTH_R'
    end if

    ! fill in the grid structure
    call create_grid(nx,ny,nz,grid)
    grid%nx = nx
    grid%ny = ny
    grid%nz = nz
    grid%nzAir = nzAir
    grid%nzCrust = nzCrust
    grid%nzEarth = nzEarth
    grid%x(1:nx)   = x(1:nx)*d2r
    grid%y(1:ny+1) = y(1:ny+1)*d2r
    grid%z(1:nz+1) = z(1:nz+1)*1000.0D0
    grid%allocated = .true.

    ! now, define ph,th,r in radians and km
    allocate(ph(nx+1),th(ny+1),r(nz+1), STAT=istat)
    ph(1) = 0.0d0
    do i=1,nx
      ph(i+1) = ph(i)+x(i)*d2r
    end do
    !ph(nx+1) = ph(1) ! don't do that - problems with interpolation!!!!
    th(1:ny+1) = y(1:ny+1)*d2r
    r(1:nz+1) = z(1:nz+1)

    ! save the cell nodes and distances in radians and km, respectively
    grid%ph(1) = ph(1)
    do i=2,nx+1
      grid%ph(i) = ph(i)
      grid%dp(i-1) = ph(i) - ph(i-1)
    end do

    grid%th(1) = th(1)
    do i=2,ny+1
      grid%th(i) = th(i)
      grid%dt(i-1) = th(i) - th(i-1)
    end do

    grid%r(1) = r(1) ! note: r is decreasing from top to bottom
    do i=2,nz+1
      grid%r(i) = r(i)
      grid%dr(i-1) = r(i-1) - r(i)
    end do

    deallocate(x,y,z)

    return

  end subroutine read_grid

end module GridDef
