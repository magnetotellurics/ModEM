    ! Initializes and does basic calculations for the grid. Computes basic derived
    ! grid parameters and are used repeatedly by other routines.
    ! Belongs to SG_Basics class: staggered cartesian grid, data
    ! types defined on this grid, and operations defined on these data types. Not
    ! specific to EM problem, no dependency on outside (from other classes) modules.
    ! modified by Cherevatova (Jan,2012) for the multi-grid
    module griddef

    use math_constants
    implicit none

    ! Very important: '=' sign has to be overloaded, since by default it is
    ! legal in fortran to say y = x for data types, but that doesn't copy
    ! allocatable or pointer arrays
    interface assignment (=)
      module procedure copy_mgrid
    end interface


    ! Initialization routines
    public  :: create_grid, create_mgrid ! allocates grid and multigrid(subgrids)
    public  :: deall_grid, deall_mgrid ! deallocates grid and multigrid
    public  :: setup_grid, setup_mgrid  ! setup grid geometry and compute subgrids dx,dy,...
    public  :: copy_grid, copy_mgrid ! copy grid, multigrid(subgrids)

    ! Possible grid types for EMfield, storing the intention of use for types
    ! such as cvector, cscalar, rvector, rscalar, sparsevecc.
    character(len=80), parameter		:: FACE = 'FACE'
    character(len=80), parameter		:: EDGE = 'EDGE'
    character(len=80), parameter		:: CENTER = 'CELL'
    character(len=80), parameter		:: CORNER = 'NODE'
    character(len=80), parameter		:: CELL_EARTH = 'CELL EARTH'

    ! ***************************************************************************
    ! Type grid_param consists of parameters that define the basic grid geometry
    ! used for three dimensional numerical modeling
    ! (Original grid type)
    type :: grid_orig

     ! Grid coordinate system; important - used in EMfield
     character (len=80)  :: coords = Cartesian

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nzEarth is number of earth layers used in the grid modeling
     ! nzAir is number of air layers
     ! nz is grid dimension (number of cells) in the z-direction:
     ! nz = nzAir + nzEarth
     integer  :: nx, ny, nz, nzEarth, nzAir

     ! the origin of the model, by default set to zero
     real (kind=prec)  :: ox = 0.0, oy = 0.0, oz = 0.0
     !  the rotation angle in degrees, by defualt set to zero
     real (kind=prec)  :: rotdeg = 0.0

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
     real (kind=prec), pointer, dimension(:)  :: dx, dy, dz
     real (kind=prec), pointer, dimension(:)  :: dxinv, dyinv, dzinv
     real (kind=prec), pointer, dimension(:)  :: delX, delY, delZ
     real (kind=prec), pointer, dimension(:)  :: delXinv, delYinv, delZinv

     ! Book-keeping on cumulative distances
     ! xEdge is the array for cumulative distance of the edge faces from the
     ! coordinate axis with dimensions nx+1
     ! xCenter is the array for cumulative distance of the edge center from
     ! the coordinate axis with dimensions nx
     ! yEdge, yCenter, zEdge, zCenter are analagous arrays for other directions
     ! Note that the arrays are allocated dynamically
     real (kind=prec), pointer, dimension(:)  :: xEdge, yEdge, zEdge
     real (kind=prec), pointer, dimension(:)  :: xCenter, yCenter, zCenter

     ! total thickness of the air above
     real (kind=prec)  :: zAirThick
     ! allocated:  .true.  all the arrays have been allocated
     logical  :: allocated = .false.

    end type grid_orig

    ! Extension to grid_orig type. For multigrid

    type, extends(grid_orig) :: grid_t ! multigrid type inherits all grid_orig
    ! flag = Surface; Air; Earth ! just to know where we are!!!!
    integer  :: flag
    ! to define interface Type
    ! coarse to fine / fine to coarse/ fine to fine or, possibly, coarse to coarse
    character(len = 10)  :: interfaceType(20)
    ! number of 'layers' with different resolutions
    integer  :: mgridSize! from file
    ! coarseness
    integer,allocatable  :: coarseness(:) ! from file
    ! number of Z layers in each subgrid
    integer,allocatable :: nzGrid(:)
    ! this is array of the subgrids
    type(grid_t), pointer :: gridArray(:)
    end type grid_t

    Contains

    !************************************************************************

    subroutine create_grid(Nx,Ny,NzAir,NzEarth,grid)
    !  creates finite differences grid_t structure of
    !  size  Nx x NyxNz, allocates arrays

    implicit none
    integer, intent(in)  :: Nx,Ny,NzAir,NzEarth
    type (grid_t) , intent(inout)  :: grid
    !local variables
    integer  :: Nz

    Nz = NzEarth+NzAir
    grid%NzAir = NzAir
    grid%Nx = Nx
    grid%Ny = Ny
    grid%NzEarth = NzEarth
    grid%Nz = Nz
    allocate(grid%Dx(Nx))
    allocate(grid%Dy(Ny))
    allocate(grid%Dz(Nz+1))

    ! dxinv  = 1/ dx and similarly for dyinv and dzinv
    allocate(grid%dxinv(Nx))
    allocate(grid%dyinv(Ny))
    allocate(grid%dzinv(Nz+1))

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
    allocate(grid%zCenter(Nz+1))

    grid%coords = Cartesian
    grid%allocated = .true.

    end subroutine create_grid

    ! **************************************************************************

    subroutine create_mgrid(mgrid)
    !  creates finite differences mgrid_t structure
    !   allocates sub grid arrays
    implicit none

    type(grid_t), intent(inout)  :: mgrid
    !  local variables
    integer  :: imgrid, ig
    integer  :: nzCum, nzAirGrid, nzEarthGrid

      ! Initialize multigrid structure

      ! allocate grid arrays
       allocate(mgrid%gridArray(mgrid%mgridSize))

      ! allocate subgrids
      ! start from the top atmosphere

      nzCum = 0
      do imgrid = 1,mgrid%mgridSize ! main do loop on subgrids

       mgrid%gridArray(imgrid)%nx = mgrid%nx/2**(mgrid%coarseness(imgrid)) ! nx in each subgrid
       mgrid%gridArray(imgrid)%ny = mgrid%ny/2**(mgrid%coarseness(imgrid)) ! ny in each subgrid

       nzAirGrid = 0
       nzEarthGrid = 0

         do ig = 1, mgrid%nzGrid(imgrid) ! do loop on nz in ezch subgrid

           if (nzCum.lt.mgrid%nzAir) then ! Air
             nzAirGrid = nzAirGrid+1
           else
             nzEarthGrid = nzEarthGrid+1  ! Earth
           endif
          nzCum = nzCum + 1
         enddo

          mgrid%gridArray(imgrid)%nzAir = nzAirGrid
          mgrid%gridArray(imgrid)%nzEarth = nzEarthGrid

          call create_grid(mgrid%gridArray(imgrid)%nx,mgrid%gridArray(imgrid)%ny,mgrid%gridArray(imgrid)%nzAir,&
                       mgrid%gridArray(imgrid)%nzEarth,mgrid%gridArray(imgrid))
          ! define flag
          if(mgrid%gridArray(imgrid)%nzAir.ne.0.and.mgrid%gridArray(imgrid)%nzEarth.ne.0) then ! Surface
            mgrid%gridArray(imgrid)%flag = 0  ! Surface
          else if(mgrid%gridArray(imgrid)%nzAir.ne.0.and.mgrid%gridArray(imgrid)%nzEarth.eq.0) then ! Air
            mgrid%gridArray(imgrid)%flag = 1  ! Air
          else if(mgrid%gridArray(imgrid)%nzAir.eq.0.and.mgrid%gridArray(imgrid)%nzEarth.ne.0) then ! Earth
            mgrid%gridArray(imgrid)%flag = 2  ! Earth
          endif
       enddo

       do imgrid = 1, mGrid%mgridSize-1

         if (mGrid%coarseness(imgrid) .gt. mGrid%coarseness(imgrid+1)) then
           mGrid%interfaceType(imgrid) = 'c2f'
         else if (mGrid%coarseness(imgrid) .lt. mGrid%coarseness(imgrid+1)) then
           mGrid%interfaceType(imgrid) = 'f2c'
         else if (mGrid%coarseness(imgrid) .eq. mGrid%coarseness(imgrid+1)) then
           mGrid%interfaceType(imgrid) = 'f2f'
         end if

       enddo ! imgrid

    mgrid%coords = Cartesian
    mgrid%allocated = .true.

    end subroutine create_mgrid

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
      ! copies multigrid
      subroutine copy_mgrid(mgridOut,mgridIn)

      implicit none
        type(grid_t),intent(in)  :: mgridIn
        type(grid_t),intent(inout)  :: mgridOut

        integer  :: errAll

        if(mgridOut%allocated) then
        ! just deallocate, and start over cleanly
         call deall_mgrid(mgridOut)
        endif

        allocate(mgridOut%nzGrid(mgridIn%mgridSize),STAT=errAll)
        allocate(mgridOut%coarseness(mgridIn%mgridSize), STAT=errall)

        call create_grid(mgridIn%Nx,mgridIn%Ny,mgridIn%NzAir, &
                 mgridIn%NzEarth,mgridOut)
          mgridOut%mgridSize = mgridIn%mgridSize
          mgridOut%nzGrid = mgridIn%nzGrid
          mgridOut%coarseness = mgridIn%coarseness
          mgridOut%interfaceType = mgridIn%interfaceType

        call create_mgrid(mgridOut)

          mgridOut%dz = mgridIn%dz
          mgridOut%dy = mgridIn%dy
          mgridOut%dx = mgridIn%dx
          mgridOut%ox = mgridIn%ox
          mgridOut%oy = mgridIn%oy
          mgridOut%oz = mgridIn%oz
          mgridOut%rotdeg = mgridIn%rotdeg
          mgridOut%coords = mgridIn%coords

          call setup_mgrid(mgridOut)

    end subroutine copy_mgrid

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

    subroutine deall_mgrid(mgrid)

    type (grid_t) , intent(inout)  :: mgrid
    integer  :: imgrid

    do imgrid = 1,mgrid%mgridSize ! deallocate subgrids
      call deall_grid(mgrid%gridArray(imgrid))
    enddo

    call deall_grid(mgrid) ! deallocate mgrid structure

    deallocate(mgrid%coarseness)
    deallocate(mgrid%nzGrid)
    deallocate(mgrid%gridArray)

    mgrid%allocated = .false.
    mgrid%interfaceType = 'zeroes'
    mgrid%NzAir = 0
    mgrid%Nx = 0
    mgrid%Ny = 0
    mgrid%NzEarth = 0
    mgrid%Nz = 0
    mgrid%mgridSize = 0

    end subroutine deall_mgrid

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
    DO iz = 2,grid%nz+1
      grid%delZ(iz) = grid%dz(iz-1) + grid%dz(iz)
    ENDDO
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
    do iz = 1, grid%nz+1
      zCum = zCum + grid%delZ(iz)
      grid%zCenter(iz) = zCum
    enddo

    !  need to be careful here ... grid origin is given at Earth's surface,
    !   not top of model domain!
    do iz = 1, grid%nz+1
      grid%zCenter(iz) = grid%zCenter(iz)-grid%zAirThick+grid%oz
      grid%zEdge(iz) = grid%zEdge(iz)-grid%zAirThick+grid%oz
    enddo

    end subroutine setup_grid

      ! *****************************************************************************
      ! sets up grid parameters at each subgrid

      subroutine setup_mgrid(mgrid,origin)
      ! set up subgrids geometry. its dx, dy, dz, delx ....

      implicit none
        ! input -- original grid
        ! output -- contains subgrid arrays of grid_t
        type(grid_t), target, intent(inout)  :: mgrid
        ! don't know for origin!!!
        real(kind=prec), intent(in), optional  :: origin(3)
        ! local variables

        integer  :: ccoeff_current, nzCum, nx, ny, nz
        integer  :: imgrid,ic,ix,iy, iz, izv


        call setup_grid(mgrid) ! setup original (usually finest) grid
        mgrid%dz(mgrid%nz+1) = 1.3*mgrid%dz(mgrid%nz)

        ! dx, dy,dz to each sub grid

        nzCum = 0  ! cumulative nz

        do imgrid = 1, mgrid%mgridSize ! main do loop on grids
          nx = mgrid%gridArray(imgrid)%nx
          ny = mgrid%gridArray(imgrid)%ny
          nz = mgrid%gridArray(imgrid)%nz
          ccoeff_current = 2**mgrid%coarseness(imgrid) ! 2^coarseness
          mgrid%gridArray(imgrid)%dx = .0
          do ix = 1, nx
            do ic = 1, ccoeff_current
              mgrid%gridArray(imgrid)%dx(ix) = mgrid%dx(ccoeff_current*(ix-1)+ic) & ! recompute dx
                                               +mgrid%gridArray(imgrid)%dx(ix)
            enddo
          enddo
          mgrid%gridArray(imgrid)%dy= .0
          do iy = 1, ny
            do ic = 1, ccoeff_current
              mgrid%gridArray(imgrid)%dy(iy) = mgrid%dy(ccoeff_current*(iy-1)+ic) & ! dy
                                               +mgrid%gridArray(imgrid)%dy(iy)
            enddo
          enddo

          do iz = 1, nz+1
            izv = iz + nzCum
             mgrid%gridArray(imgrid)%dz(iz) = mgrid%dz(izv)
          enddo

             mgrid%gridArray(imgrid)%ox = mgrid%ox
             mgrid%gridArray(imgrid)%oy = mgrid%oy
             mgrid%gridArray(imgrid)%oz = mgrid%oz
             mgrid%gridArray(imgrid)%rotdeg = mgrid%rotdeg
             mgrid%gridArray(imgrid)%coords = mgrid%coords

          call setup_grid(mgrid%gridArray(imgrid)) ! setup subgrid parameters. It will change dz, etc.
                                                  ! for Air dz
                                                  ! we need copy it again from the finest grid

          do iz = 1, nz+1
            izv = iz + nzCum

             mgrid%gridArray(imgrid)%dz(iz) = mgrid%dz(izv)
             mgrid%gridArray(imgrid)%dzinv(iz) = mgrid%dzinv(izv)
             mgrid%gridArray(imgrid)%Zedge(iz) = mgrid%Zedge(izv)
             mgrid%gridArray(imgrid)%delZ(iz) = mgrid%delZ(izv)
             mgrid%gridArray(imgrid)%delZinv(iz) = mgrid%delZinv(izv)
             mgrid%gridArray(imgrid)%Zcenter(iz) = mgrid%Zcenter(izv)
          enddo

             mgrid%gridArray(imgrid)%zAirThick = mgrid%zAirThick

        nzCum = nzCum + nz
        enddo


      end subroutine  setup_mgrid

    end module griddef
