!
!> Implementation of standard 3D Cartesian grid.
!
!> dx_inv = 1/ dx and similarly for dy_inv and dz_inv
!
!> del_x, del_y, and del_z are the distances between
!> the electrical field defined on the center of the
!> edges in x, y, and z axes, respectively.
!
!> x_edge is the array for cumulative distance of the edge
!> for each grid (starting from the coordinate axes) with
!> dimensions nx + 1.
!
!> x_center is the array for cumulative distance of the center
!> for each grid (starting from the coordinate axes) with
!> dimension n.
!
!> y_edge, y_center, z_edge, z_center are analagous arrays for
!> other directions.
!
module Grid3D_SG
    !
    use Grid
    !
    type, extends( Grid_t ) :: Grid3D_SG_t
        !
        real( kind=prec ), allocatable, dimension(:) :: dx_inv, dy_inv, dz_inv
        !
        real( kind=prec ), allocatable, dimension(:) :: del_x, del_y, del_z
        !
        real( kind=prec ), allocatable, dimension(:) :: del_x_inv, del_y_inv, del_z_inv
        !
        real( kind=prec ), allocatable, dimension(:) :: x_center, y_center, z_center
        !
        real( kind=prec ), allocatable, dimension(:) :: x_edge, y_edge, z_edge
        !
        real( kind=prec ) :: zAirThick
        !
        contains
            !
            final :: Grid3D_SG_dtor
            !
            procedure, public :: numberOfEdges => numberOfEdges_Grid3D_SG
            procedure, public :: numberOfFaces => numberOfFaces_Grid3D_SG
            procedure, public :: numberOfNodes => numberOfNodes_Grid3D_SG
            procedure, public :: numberOfCells => numberOfCells_Grid3D_SG
            !
            procedure, public :: length => length_Grid3D_SG
            !
            procedure, public :: allocateDim => allocateDim_Grid3D_SG
            !
            procedure, public :: deallocateDim => deallocateDim_Grid3D_SG
            !
            procedure, public :: setup => setup_Grid3D_SG
            !
            procedure, public :: slice1D => slice1D_Grid3D_SG
            procedure, public :: slice2D => slice2D_Grid3D_SG
            !
            procedure, public :: write => write_Grid3D_SG
            !
    end type Grid3D_SG_t
    !
    interface Grid3D_SG_t
        module procedure Grid3D_SG_t_ctor_n
        module procedure Grid3D_SG_t_ctor_n_d
    end interface Grid3D_SG_t
    !
contains
    !
    !> Class constructor for simple tensor product grid
    !> Usage obj = Grid_t3D(Dx,Dy,Dz,Nza)
    !> Dx, Dy, Dz are cell dimensions for x, y, z direction
    !> Nza is number of air layers to allow(included in Dz)
    !
    function Grid3D_SG_t_ctor_n( nx, ny, nzAir, nzEarth ) result( self )
        implicit none
        !
        integer, intent( in ) :: nx, ny, nzAir, nzEarth
        !
        type( Grid3D_SG_t ) :: self
        !
        !write( *, * ) "Constructor_n Grid3D_SG_t"
        !
        self%n_grids = 1
        !
        call self%baseInit
        !
        call self%create( nx, ny, nzAir, nzEarth )
        !
    end function Grid3D_SG_t_ctor_n
    !
    !> Class constructor for simple tensor product grid
    !> Usage obj = Grid_t3D(Dx,Dy,Dz,Nza)
    !> Dx, Dy, Dz are cell dimensions for x, y, z direction
    !> Nza is number of air layers to allow(included in Dz)
    !
    function Grid3D_SG_t_ctor_n_d( nx, ny, nzAir, nzEarth, dx, dy, dz ) result( self )
        implicit none
        !
        integer, intent( in ) :: nx, ny, nzAir, nzEarth
        real( kind=prec ), dimension(:), intent( in ) :: dx, dy, dz
        !
        type( Grid3D_SG_t ) :: self
        !
        !write( *, * ) "Constructor_n_d Grid3D_SG_t"
        !
        self%n_grids = 1
        !
        call self%baseInit
        !
        call self%create( nx, ny, nzAir, nzEarth )
        !
        call self%setCellSizes( dx, dy, dz )
        !
        call self%setup
        !
    end function Grid3D_SG_t_ctor_n_d
    !
    !> Deconstructor routine:
    !>     Calls the base routine baseDealloc().
    subroutine Grid3D_SG_dtor( self )
        implicit none
        !
        type( Grid3D_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor Grid3D_SG_t", self%nx, self%ny, self%nz
        !
        call self%baseDealloc
        !
        call self%deallocateDim
        !
        self%is_allocated = .FALSE.
        !
    end subroutine Grid3D_SG_dtor
    !
    !> No subroutine briefing
    !
    subroutine allocateDim_Grid3D_SG( self )
        implicit none
        !
        class( Grid3D_SG_t ), intent( inout ) :: self
        !
        if( self%is_allocated ) call self%baseDealloc
        !
        call self%deallocateDim
        !
        allocate( self%dx( self%nx ) )
        allocate( self%dy( self%ny ) )
        allocate( self%dz( self%nz ) )
        !
        allocate( self%dx_inv( self%nx ) )
        allocate( self%dy_inv( self%ny ) )
        allocate( self%dz_inv( self%nz ) )
        !
        allocate( self%del_x( self%nx + 1 ) )
        allocate( self%del_y( self%ny + 1 ) )
        allocate( self%del_z( self%nz + 1 ) )
        !
        allocate( self%del_x_inv( self%nx + 1 ) )
        allocate( self%del_y_inv( self%ny + 1 ) )
        allocate( self%del_z_inv( self%nz + 1 ) )
        !
        allocate( self%x_center( self%nx ) )
        allocate( self%y_center( self%ny ) )
        allocate( self%z_center( self%nz ) )
        !
        allocate( self%x_edge( self%nx + 1 ) )
        allocate( self%y_edge( self%ny + 1 ) )
        allocate( self%z_edge( self%nz + 1 ) )
        !
        self%is_allocated = .TRUE.
        !
    end subroutine allocateDim_Grid3D_SG
    !
    !> Deconstructor routine:
    !>     Calls the base routine baseDealloc().
    subroutine deallocateDim_Grid3D_SG( self )
        implicit none
        !
        class( Grid3D_SG_t ), intent( inout ) :: self
        !
        if( allocated( self%dx_inv ) ) deallocate( self%dx_inv )
        if( allocated( self%dy_inv ) ) deallocate( self%dy_inv )
        if( allocated( self%dz_inv ) ) deallocate( self%dz_inv )
        !
        if( allocated( self%del_x ) ) deallocate( self%del_x )
        if( allocated( self%del_y ) ) deallocate( self%del_y )
        if( allocated( self%del_z ) ) deallocate( self%del_z )
        !
        if( allocated( self%del_x_inv ) ) deallocate( self%del_x_inv )
        if( allocated( self%del_y_inv ) ) deallocate( self%del_y_inv )
        if( allocated( self%del_z_inv ) ) deallocate( self%del_z_inv )
        !
        if( allocated( self%x_center ) ) deallocate( self%x_center )
        if( allocated( self%y_center ) ) deallocate( self%y_center )
        if( allocated( self%z_center ) ) deallocate( self%z_center )
        !
        if( allocated( self%x_edge ) ) deallocate( self%x_edge )
        if( allocated( self%y_edge ) ) deallocate( self%y_edge )
        if( allocated( self%z_edge ) ) deallocate( self%z_edge )
        !
    end subroutine deallocateDim_Grid3D_SG
    !
    !> setup does calculations for grid geometry, which cannot be done
    !> until dx, dy, dz, and the origin are set.
    subroutine setup_Grid3D_SG( self, origin )
        implicit none
        !
        class( Grid3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ), optional :: origin(3)
        !
        integer :: ix, iy, iz, i, j, nzAir
        real( kind=prec ) :: xCum, yCum, zCum
        real( kind=prec ) :: ox, oy, oz
        !
        self%dx_inv = 1 / self%dx
        self%dy_inv = 1 / self%dy
        self%dz_inv = 1 / self%dz
        !
        call self%getOrigin( ox, oy, oz )
        !
        if( present( origin ) ) then
            !
            ox = origin(1)
            oy = origin(2)
            oz = origin(3)
            !
            call self%setOrigin(ox, oy, oz)
            !
        endif
        !
        self%x_edge(1) = ox
        self%y_edge(1) = oy
        self%z_edge(1) = 0.0 !> ALWAYS BE ZERO ???? BEFORE WAS oz
        !
        xCum = R_ZERO
        yCum = R_ZERO
        zCum = R_ZERO
        !
        do ix = 1, self%nx
            xCum = xCum + self%dx(ix)
            self%x_edge(ix+1) = xCum + ox
        enddo
        !
        do iy = 1, self%ny
            yCum = yCum + self%dy(iy)
            self%y_edge(iy + 1) = yCum + oy
        enddo
        !
        !> NOTE: adjust for origin later to get airthickness, 
        !> reference to origin at Earth"s surface correct!
        do iz = 1, self%nz
            zCum = zCum + self%dz(iz)
            self%z_edge(iz + 1) = zCum !> CHECK - WHY SHIFT X and Y BUT NOT Z ????
        enddo
        !
        nzAir = self%nzAir
        self%zAirThick = self%z_edge(nzAir + 1)
        !
        !> Distance between center of the selfs
        self%del_x(1) = self%dx(1)
        do ix = 2, self%nx
            self%del_x(ix) = self%dx(ix - 1) + self%dx(ix)
        enddo
        !
        self%del_x(self%nx + 1) = self%dx(self%nx)
        self%del_x = self%del_x / 2.0
        !
        self%del_y(1) = self%dy(1)
        do iy = 2, self%ny
            self%del_y(iy) = self%dy(iy - 1) + self%dy(iy)
        enddo
        !
        self%del_y(self%ny + 1) = self%dy(self%ny)
        self%del_y = self%del_y / 2.0
        !
        self%del_z(1) = self%dz(1)
        do iz = 2, self%nz
            self%del_z(iz) = self%dz(iz - 1) + self%dz(iz)
        enddo
        !
        self%del_z(self%nz + 1) = self%dz(self%nz)
        self%del_z = self%del_z / 2.0
        !
        self%del_x_inv = 1 / self%del_x
        self%del_y_inv = 1 / self%del_y
        self%del_z_inv = 1 / self%del_z
        !
        !> Cumulative distance between the centers, adjusted to model origin
        xCum = R_ZERO
        yCum = R_ZERO
        zCum = R_ZERO
        do ix = 1, self%nx
            xCum = xCum + self%del_x(ix)
            self%x_center(ix) = xCum + ox
        enddo
        do iy = 1, self%ny
            yCum = yCum + self%del_y(iy)
            self%y_center(iy) = yCum + oy
        enddo
        do iz = 1, self%nz
            zCum = zCum + self%del_z(iz)
            self%z_center(iz) = zCum
        enddo
        !
        !> Need to be careful here ... grid origin is given
        !> at Earth"s surface, not top of model domain!
        do iz = 1, self%nz
            self%z_center(iz) = self%z_center(iz) - self%zAirThick + oz
            self%z_edge(iz) = self%z_edge(iz) - self%zAirThick + oz
        enddo
        !
        self%z_edge(self%nz + 1) = self%z_edge(self%nz + 1) - self%zAirThick + oz
        !
    end subroutine setup_Grid3D_SG
    !
    !> No subroutine briefing
    subroutine numberOfEdges_Grid3D_SG( self, n_xedge, n_yedge, n_zedge )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        integer, intent( out ) :: n_xedge, n_yedge, n_zedge
        !
        integer :: nx, ny, nz
        !
        call self%setLimits( XEDGE, nx, ny, nz )
        n_xedge = nx * ny * nz
        !
        call self%setLimits( YEDGE, nx, ny, nz )
        n_yedge  = nx * ny * nz
        !
        call self%setLimits( ZEDGE, nx, ny, nz )
        n_zedge = nx * ny * nz
        !
    end subroutine numberOfEdges_Grid3D_SG
    !
    !> No subroutine briefing
    subroutine numberOfFaces_Grid3D_SG( self, n_xface, n_yface, n_zface ) 
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        integer, intent( out ) :: n_xface, n_yface, n_zface
        !
        integer :: nx, ny, nz
        !
        call self%setLimits( XFACE, nx, ny, nz )
        n_xface = nx * ny * nz
        !
        call self%setLimits( YFACE, nx, ny, nz )
        n_yface = nx * ny * nz
        !
        call self%setLimits( ZFACE, nx, ny, nz )
        n_zface = nx * ny * nz
        !
    end subroutine numberOfFaces_Grid3D_SG
    !
    !> No subroutine briefing
    !
    function numberOfNodes_Grid3D_SG( self ) result(n)
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        !
        integer :: n
        !
        n =(self%nx + 1) * (self%ny + 1) * (self%nz + 1)
        !
    end function numberOfNodes_Grid3D_SG
    !
    !> No function briefing
    !
    function numberOfCells_Grid3D_SG( self ) result( n )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        !
        integer :: n
        !
        call errStop( "numberOfCells_Grid3D_SG > numberOfCells_Grid3D_SG not implemented" )
        !
    end function numberOfCells_Grid3D_SG
    !
    !> No subroutine briefing
    !
    function length_Grid3D_SG( self ) result(n)
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        integer :: n
        !
        n = self%nx * self%ny * self%nz
        !
    end function length_Grid3D_SG
    !
    !> No subroutine briefing
    !
    function slice1D_Grid3D_SG( self ) result( g1D )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        type( Grid1D_t ) :: g1D
        !
        g1D = Grid1D_t( self%nzAir, self%nzEarth, self%dz )
        !
    end function slice1D_Grid3D_SG
    !
    !> No subroutine briefing
    !
    function slice2D_Grid3D_SG( self ) result( g2D )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        type( Grid2D_t ) :: g2D
        !
        !> Should be different for the polarization
        g2D = Grid2D_t( self%ny, self%nzAir, self%nzEarth, self%dy, self%dz )
        !
    end function slice2D_Grid3D_SG
    !
    !
    !
    subroutine write_Grid3D_SG( self )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        !
        write( *, * ) "Grid3D_SG:"
        !
        write( *, * ) "    n_grids: ", self%n_grids
        !
        write( *, * ) "    is_allocated: ", self%is_allocated
        !
        write( *, * ) "    nx, ny, nz: ", self%nx, self%ny, self%nz
        write( *, * ) "    nzAir: ", self%nzAir   !> Number of air layers
        write( *, * ) "    nzEarth: ", self%nzEarth !> Number of earth layers
        !
        write( *, * ) "    geometry: ", self%geometry
        !
        write( *, * ) "    ox, oy, oz: ", self%ox, self%oy, self%oz
        !
        write( *, * ) "    rotDeg: ", self%rotDeg
        !
        write( *, * ) "    zAirThick: ", self%zAirThick
        !
        write( *, * ) "    EDGEb, FACEb, NODEb: ", size( self%EDGEb ), size( self%FACEb ), size( self%NODEb )
        write( *, * ) "    EDGEi, FACEi, NODEi: ", size( self%EDGEi ), size( self%FACEi ), size( self%NODEi )
        !
        write( *, * ) "    dx, dy, dz: ", size( self%dx ), size( self%dy ), size( self%dz )
        write( *, * ) "    dx_inv, dy_inv, dz_inv: ", size( self%dx_inv ), size( self%dy_inv ), size( self%dz_inv )
        !
        write( *, * ) "    del_x, del_y, del_z: ", size( self%del_x ), size( self%del_y ), size( self%del_z )
        write( *, * ) "    del_x_inv, del_y_inv, del_z_inv: ", size( self%del_x_inv ), size( self%del_y_inv ), size( self%del_z_inv )
        !
        write( *, * ) "    x_edge, y_edge, z_edge: ", size( self%x_edge ), size( self%y_edge ), size( self%z_edge )
        !
        write( *, * ) "    x_center, y_center, z_center: ", size( self%x_center ), size( self%y_center ), size( self%z_center )
        !
    end subroutine write_Grid3D_SG
    !
end module Grid3D_SG
