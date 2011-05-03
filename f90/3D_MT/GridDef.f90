! *****************************************************************************
! initializes and does basic calculations for the grid. Computes basic derived
! grid parameters and are used repeatedly by other routines.
! Belongs to SG_Basics class: staggered cartesian grid, data
! types defined on this grid, and operations defined on these data types. Not
! specific to EM problem, no dependency on outside (from other classes) modules.
module GridDef

  use math_constants
  implicit none

  ! Very important: '=' sign has to be overloaded, since by default it is
  ! legal in fortran to say y = x for data types, but that doesn't copy
  ! allocatable or pointer arrays
  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE copy_grid
  END INTERFACE

  ! Initialization routines
  public                            :: create_grid, deall_grid
  public                            :: copy_grid, setup_grid

  ! Possible grid types for EMfield, storing the intention of use for types
  ! such as cvector, cscalar, rvector, rscalar, sparsevecc.
  character(len=80), parameter		:: FACE = 'FACE'
  character(len=80), parameter		:: EDGE = 'EDGE'
  character(len=80), parameter		:: CENTER = 'CELL'
  character(len=80), parameter		:: CORNER = 'NODE'
  character(len=80), parameter		:: CELL_EARTH = 'CELL EARTH'

  ! ***************************************************************************
  ! type grid_param consists of parameters that define the basic grid geometry
  ! used for three dimensional numerical modeling
  type :: grid_t

     ! Grid coordinate system; important - used in EMfield
     character (len=80)			:: coords = Cartesian

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nzEarth is number of earth layers used in the grid modeling
     ! nzAir is number of air layers
     ! nz is grid dimension (number of cells) in the z-direction:
     ! nz = nzAir + nzEarth
     integer               :: nx, ny, nz, nzEarth, nzAir

     ! the origin of the model, by default set to zero
     real (kind=prec)			:: ox = 0.0, oy = 0.0, oz = 0.0
     !  the rotation angle in degrees, by defualt set to zero
     real (kind=prec)			:: rotdeg = 0.0

     ! Grid geometry:
     ! dx,delX are arrays of grid spacings in x-direction
     ! dx denotes spacings betgween cell edges: dimension: dx(nx)
     ! dxinv  = 1/ dx
     ! delX denotes spacings between cell centers: dimension: delX(nx+1)
     ! why dimensions: delX(nx+1) (delX(2) is distance between centers of cells
     ! 2 and 1)
     ! delXinv = 1/ delX
     ! dy,delY, dz, delZ are analagous arrays for other directions
     ! similarly, are dyinv and delYinv
     ! Note that the arrays are allocated dynamically
     real (kind=prec), pointer, dimension(:)        :: dx, dy, dz
     real (kind=prec), pointer, dimension(:)        :: dxinv, dyinv, dzinv
     real (kind=prec), pointer, dimension(:)        :: delX, delY, delZ
     real (kind=prec), pointer, dimension(:)        :: delXinv, delYinv, &
     delZinv

     ! Book-keeping on cumulative distances
     ! xEdge is the array for cumulative distance of the edge faces from the
     ! coordinate axis with dimensions nx+1
     ! xCenter is the array for cumulative distance of the edge center from
     ! the coordinate axis with dimensions nx
     ! yEdge, yCenter, zEdge, zCenter are analagous arrays for other directions
     ! Note that the arrays are allocated dynamically
     real (kind=prec), pointer, dimension(:)        :: xEdge, yEdge, zEdge
     real (kind=prec), pointer, dimension(:)        :: xCenter, yCenter, &
     zCenter

     ! total thickness of the air above
     real (kind=prec)				     :: zAirThick

     ! allocated:  .true.  all the arrays have been allocated
     logical		                             :: allocated = .false.

  end type grid_t

Contains

  !************************************************************************
  subroutine create_grid(Nx,Ny,NzAir,NzEarth,grid)
    !  creates finite differences grid_t structure of
    !  size  Nx x Ny Nz, allocates arrays
    !
    implicit none
    integer, intent(in)				    :: Nx,Ny,NzAir,NzEarth
    type (grid_t) , intent(inout)		:: grid

    !local variables
    integer					:: Nz

    Nz = NzEarth+NzAir
    grid%NzAir = NzAir
    grid%Nx = Nx
    grid%Ny = Ny
    grid%NzEarth = NzEarth
    grid%Nz = Nz
    allocate(grid%Dx(Nx))
    allocate(grid%Dy(Ny))
    allocate(grid%Dz(Nz))

    ! dxinv  = 1/ dx and similarly for dyinv and dzinv
    allocate(grid%dxinv(Nx))
    allocate(grid%dyinv(Ny))
    allocate(grid%dzinv(Nz))

    ! delX, delY, and delZ are the distances between the electrical field
    ! defined on the center of the edges in x, y, and z axes, respectively.
    allocate(grid%delX(Nx+1))
    allocate(grid%delY(Ny+1))
    allocate(grid%delZ(Nz+1))

    ! delXinv = 1/ delX and similarly for delYinv and delZinv
    allocate(grid%delXinv(Nx+1))
    allocate(grid%delYinv(Ny+1))
    allocate(grid%delZinv(Nz+1))

    ! xEdge is the array for cumulative distance of the edge for each
    ! grid (starting from the coordinate axes) with dimensions nx+1
    ! xCenter is the array for cumulative distance of the center for each
    ! grid (starting from the coordinate axes) with dimension n
    ! yEdge, yCenter, zEdge, zCenter are analagous arrays for other directions
    allocate(grid%xEdge(Nx+1))
    allocate(grid%yEdge(Ny+1))
    allocate(grid%zEdge(Nz+1))
    allocate(grid%xCenter(Nx))
    allocate(grid%yCenter(Ny))
    allocate(grid%zCenter(Nz))

	grid%coords = Cartesian
    grid%allocated = .true.

  end subroutine create_grid

  ! **************************************************************************
  subroutine copy_grid(gridOut,gridIn)

  !  copies gridIn to gridOut; cannot overwrite, of course!

  type(grid_t),intent(in)		    :: gridIn
  type(grid_t),intent(inout)		:: gridOut

     if(gridOut%allocated) then
        !  just deallocate, and start over cleanly
        call deall_grid(gridOut)
     endif

     call create_grid(gridIn%Nx,gridIn%Ny,gridIn%NzAir, &
             gridIn%NzEarth,gridOut)

     gridOut%Dz = gridIn%Dz
     gridOut%Dy = gridIn%Dy
     gridOut%Dx = gridIn%Dx
     gridOut%ox = gridIn%ox
     gridOut%oy = gridIn%oy
     gridOut%oz = gridIn%oz

     gridOut%rotdeg = gridIn%rotdeg
     gridOut%coords = gridIn%coords

     call setup_grid(gridOut)

  end subroutine copy_grid

  ! **************************************************************************
  subroutine deall_grid(grid)

    type (grid_t) , intent(inout)	:: grid

    deallocate(grid%Dx)
    deallocate(grid%Dy)
    deallocate(grid%Dz)

    deallocate(grid%dxinv)
    deallocate(grid%dyinv)
    deallocate(grid%dzinv)

    deallocate(grid%delX)
    deallocate(grid%delY)
    deallocate(grid%delZ)

    deallocate(grid%delXinv)
    deallocate(grid%delYinv)
    deallocate(grid%delZinv)

    deallocate(grid%xEdge)
    deallocate(grid%yEdge)
    deallocate(grid%zEdge)
    deallocate(grid%xCenter)
    deallocate(grid%yCenter)
    deallocate(grid%zCenter)

    grid%allocated = .false.
    grid%NzAir = 0
    grid%Nx = 0
    grid%Ny = 0
    grid%NzEarth = 0
    grid%Nz = 0

 end subroutine deall_grid

  ! **************************************************************************
  ! setup_grid does calculations for grid geometry, which cannot be done
  ! until dx, dy, dz, and the origin are set.
  ! Normal usage is to first call create_grid to set grid dimensions
  ! and allocate arrays, read dx, dy, dz and set these elements of the grid,
  ! then call setup_grid to do all other computations.  By including the optional
  ! origin argument, grid centers and edges are given in absolute coordinates
  ! (i.e., the origin of the grid at the Earth surface is set to the origin,
  ! and variables like xCenter, yEdge, etc. are given in the same coordinate
  ! system).  If argument origin is not present, whatever is set already in the grid
  ! origin is used; by default this is initialized to zero.
  subroutine setup_grid(grid, origin)

    implicit none
    type(grid_t), target, intent(inout)    :: grid
    real(kind=prec), intent(in), optional	  :: origin(3)

    integer                               :: ix,iy,iz,i,j
    integer                               :: status
    real (kind=prec)                         :: xCum, yCum, zCum
    real(kind=prec)                     :: alpha = 3.

    !   Following is Kush's approach to setting air layers:
    ! mirror imaging the dz values in the air layer with respect to
    ! earth layer as far as we can using the following formulation
    ! air layer(bottom:top) = (alpha)^(j-1) * earth layer(top:bottom)
    if (minval(grid%dz) .le. R_ZERO) then
	    i = grid%nzAir+1
	    j = 0
	    do iz = grid%nzAir, 1, -1
	        j = j + 1
	        grid%dz(iz) = ((alpha)**(j-1))*grid%dz(i)
	        i = i + 1
	    end do
    end if

    ! the topmost air layer has to be atleast 30 km
    if (grid%dz(1).lt.30000) then
        grid%dz(1) = 30000
    end if

    grid%dxinv = 1/ grid%dx
    grid%dyinv = 1/ grid%dy
    grid%dzinv = 1/ grid%dz

    grid%rotdeg = grid%rotdeg
    if (present(origin)) then
    	grid%ox = origin(1)
    	grid%oy = origin(2)
    	grid%oz = origin(3)
    end if

    grid%xEdge(1) = grid%ox
    grid%yEdge(1) = grid%oy
    grid%zEdge(1) = 0.0
    xCum = 0.0
    yCum = 0.0
    zCum = 0.0
    do ix = 1, grid%nx
       xCum = xCum + grid%dx(ix)
       grid%xEdge(ix+1) = xCum + grid%ox
    enddo
    do iy = 1, grid%ny
       yCum = yCum + grid%dy(iy)
       grid%yEdge(iy+1) = yCum + grid%oy
    enddo
    !  NOTE: adjust for origin later to get airthickness, reference to origin
    !    at Earth's surface correct!
    do iz = 1, grid%nz
       zCum = zCum + grid%dz(iz)
       grid%zEdge(iz+1) = zCum
    enddo
    grid%zAirThick = grid%zEdge(grid%nzAir+1)

    ! distance between center of the grids
    grid%delX(1) = grid%dx(1)
    DO ix = 2,grid%nx
       grid%delX(ix) = grid%dx(ix-1) + grid%dx(ix)
    ENDDO
    grid%delX(grid%nx+1) = grid%dx(grid%nx)
    grid%delX = grid%delX/2.0

    grid%delY(1)    = grid%dy(1)
    DO iy = 2,grid%ny
       grid%delY(iy) = grid%dy(iy-1) + grid%dy(iy)
    ENDDO
    grid%delY(grid%ny+1) = grid%dy(grid%ny)
    grid%delY = grid%delY/2.0

    grid%delZ(1)    = grid%dz(1)
    DO iz = 2,grid%nz
       grid%delZ(iz) = grid%dz(iz-1) + grid%dz(iz)
    ENDDO
    grid%delZ(grid%nz+1) = grid%dz(grid%nz)
    grid%delZ = grid%delZ/ 2.0

    grid%delXinv = 1/ grid%delX
    grid%delYinv = 1/ grid%delY
    grid%delZinv = 1/ grid%delZ

    xCum = 0.0
    yCum = 0.0
    zCum = 0.0
    ! cumulative distance between the centers, adjusted to model origin
    do ix = 1, grid%nx
       xCum = xCum + grid%delX(ix)
       grid%xCenter(ix) = xCum + grid%ox
    enddo
    do iy = 1, grid%ny
       yCum = yCum + grid%delY(iy)
       grid%yCenter(iy) = yCum + grid%oy
    enddo
    do iz = 1, grid%nz
       zCum = zCum + grid%delZ(iz)
       grid%zCenter(iz) = zCum
    enddo

    !  need to be careful here ... grid origin is given at Earth's surface,
    !   not top of model domain!
    do iz = 1, grid%nz
       grid%zCenter(iz) = grid%zCenter(iz)-grid%zAirThick+grid%oz
       grid%zEdge(iz) = grid%zEdge(iz)-grid%zAirThick+grid%oz
    enddo
    grid%zEdge(grid%nz+1) = grid%zEdge(grid%nz+1)-grid%zAirThick+grid%oz


  end subroutine setup_grid

end module GridDef
