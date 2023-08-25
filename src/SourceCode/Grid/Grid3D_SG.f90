!
!> Implementation of standard 3D Cartesian grid.
!
module Grid3D_SG
    !
    use Grid
    !
    type, extends( Grid_t ) :: Grid3D_SG_t
        !
        !> No derived properties
        !
        contains
            !
            final :: Grid3D_SG_dtor
            !
            procedure, public :: numberOfEdges => numberOfEdges_Grid3D_SG
            procedure, public :: numberOfFaces => numberOfFaces_Grid3D_SG
            procedure, public :: numberOfNodes => numberOfNodes_Grid3D_SG
            procedure, public :: numberOfCells => numberOfCells_Grid3D_SG
            procedure, public :: gridIndex => gridIndex_Grid3D_SG
            procedure, public :: vectorIndex => vectorIndex_Grid3D_SG
            procedure, public :: setLimits => setLimits_Grid3D_SG
            procedure, public :: length => length_Grid3D_SG
            !
            procedure, public :: setup => setup_Grid3D_SG
            !
            procedure, public :: slice1D => slice1D_Grid3D_SG
            procedure, public :: slice2D => slice2D_Grid3D_SG
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
    end subroutine Grid3D_SG_dtor
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
        call self%GetOrigin( ox, oy, oz )
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
        self%z_edge(1) = oz
        !
        xCum = R_ZERO
        yCum = R_ZERO
        zCum = R_ZERO
        !
        do ix = 1, self%nx
            xCum = xCum + self%dx(ix)
            self%x_edge(ix+1) = xCum + ox
        enddo
        do iy = 1, self%ny
            yCum = yCum + self%dy(iy)
            self%y_edge(iy + 1) = yCum + oy
        enddo
        !
        !> NOTE: adjust for origin later to get airthickness, 
        !> reference to origin at Earth"s surface correct!
        do iz = 1, self%nz
            zCum = zCum + self%dz(iz)
            self%z_edge(iz + 1) = zCum
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
    !> Based on matlab method of same name in class Grid_t3D
    !> IndVec is the index within the list of nodes of a fixed type
    !> e.g., among the list of y-Faces.     An offset needs to be
    !> added to get index in list of all faces(for example).
    !
    subroutine gridIndex_Grid3D_SG( self, node_type, ind_vec, i, j, k )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        character(*), intent( in ) :: node_type
        integer, intent( in ) :: ind_vec(:)
        integer, intent( out ) :: i(:), j(:), k(:)
        !
        integer :: nx, ny, nz, nVec, ii
        real(4) :: rNxy, rNx
        !
        call self%setLimits(node_type, nx, ny, nz)
        nVec = size(ind_vec)
        !
        if( nVec .NE. size(i) ) then
            call errStop( "gridIndex_Grid3D_SG > Size of 'ind_vec' and 'i' do not agree." )
        endif
        !
        if( nVec .NE. size(j) ) then
            call errStop( "gridIndex_Grid3D_SG > Size of 'ind_vec' and 'j' do not agree." )
        endif
        !
        if( nVec .NE. size(k) ) then
            call errStop( "gridIndex_Grid3D_SG > Size of 'ind_vec' and 'k' do not agree." )
        endif
        !
        rNxy = float(nx*ny)
        rNx = float(nx)
        !
        do ii = 1, nVec
            i(ii) = mod(ind_vec(ii), nx)
            j(ii) = mod(ceiling(float(ind_vec(ii) )/rNx), ny)
            k(ii) = ceiling(float(ind_vec(ii) )/rNxy)
        enddo
        !
        where( i .EQ. 0 ) i = nx
        where( j .EQ. 0 ) j = ny
        where( k .EQ. 0 ) k = nz
        !
    end subroutine gridIndex_Grid3D_SG
    !
    !> vectorIndex
    !
    !> Based on matlab method of same name in class Grid_t3D
    !> returned array IndVec gives numbering of nodes within
    !> the list for node_type; need to add an offset for position
    !> in full list of all faces or edges(not nodes and cells).
    subroutine vectorIndex_Grid3D_SG( self, node_type, i, j, k, ind_vec )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        character(*), intent( in ) :: node_type
        integer, intent( in ) :: i(:), j(:), k(:)
        integer, intent( out ) :: ind_vec(:)
        !
        integer :: nx, ny, nz, nxy, nVec, ii
        !
        call self%setLimits(node_type, nx, ny, nz)
        !
        nVec = size(ind_vec)
        !
        if( nVec .NE. size(i) ) then
            call errStop( "vectorIndex_Grid3D_SG > Size of 'ind_vec' and 'i' do not agree." )
        endif
        !
        if( nVec .NE. size(J) ) then
            call errStop( "vectorIndex_Grid3D_SG > Size of 'ind_vec' and 'j' do not agree." )
        endif
        !
        if( nVec .NE. size(K) ) then
            call errStop( "vectorIndex_Grid3D_SG > Size of 'ind_cec' and 'k' do not agree." )
        endif
        !
        nxy = nx*ny
        do ii = 1, nVec
            ind_vec(ii) =(K(ii) - 1) * nxy +(j(ii) - 1) * nx + i(ii)
        enddo
        
    end subroutine vectorIndex_Grid3D_SG
    !
    !> No subroutine briefing
    subroutine setLimits_Grid3D_SG( self, node_type, nx, ny, nz )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        character(*), intent( in ) :: node_type
        integer, intent( out ) :: nx, ny, nz
        !
        select case( node_type )
        !
            case( CELL, CELL_EARTH )
                 nx = self%nx
                 ny = self%ny
                 nz = self%nz
            case( NODE )
                 nx = self%nx + 1
                 ny = self%ny + 1
                 nz = self%nz + 1
            case( XEDGE )
                nx = self%nx
                ny = self%ny + 1
                nz = self%nz + 1
                !
            case( XFACE )
                 nx = self%nx + 1
                 ny = self%ny
                 nz = self%nz
            case( YEDGE )
                 nx = self%nx + 1
                 ny = self%ny
                 nz = self%nz + 1
            case( YFACE )
                 nx = self%nx
                 ny = self%ny + 1
                 nz = self%nz
            case( ZEDGE )
                 nx = self%nx + 1
                 ny = self%ny + 1
                 nz = self%nz
            case( ZFACE )
                 nx = self%nx
                 ny = self%ny
                 nz = self%nz + 1
                !
            case default
                !
                call errStop( "setLimits_Grid3D_SG > Undefined node_type ["//node_type//"]" )
                !
        end select
        !
    end subroutine setLimits_Grid3D_SG
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
end module Grid3D_SG
