!
!> Abstract Base class to define a Grid
!> Definition of the base class for grids. 
!> Also defines some auxiliary constants and data types.
!
module Grid
    !
    use Utilities
    use Grid1D
    use Grid2D
    !
    !> Global Grid related properties
    integer :: model_n_air_layer
    !
    real( kind=prec ) :: model_max_height
    !
    character(:), allocatable :: grid_format
    character( len=12 ), parameter :: GRID_SG = "StandardGrid"
    character( len=19 ), parameter :: GRID_MR = "MultiresolutionGrid"
    !
    character(:), allocatable :: model_method
    character( len=12 ), parameter :: MM_METHOD_FIXED_H = "fixed height"
    character( len=6 ), parameter :: MM_METHOD_MIRROR = "mirror"
    !
    !> Base Grid class
    type, abstract :: Grid_t
        !
        integer :: n_grids
        !
        logical :: is_allocated
        !
        integer, allocatable, dimension(:) :: EDGEb, FACEb, NODEb
        integer, allocatable, dimension(:) :: EDGEi, FACEi, NODEi
        integer, allocatable, dimension(:) :: EDGEa, FACEa, NODEa
        !
        character( len=80 ) :: geometry
        !
        integer :: nx, ny, nz    !> Think is usefull have this for SG and MR
        !
        !> USED IN MR setSubGrid !!!!
        real( kind=prec ), allocatable, dimension(:) :: dx, dy, dz
        !
        !> USED IN MR setSubGrid !!!!
        integer :: nzAir
        !
        real( kind=prec ) :: ox, oy, oz
        !
        real( kind=prec ) :: rotDeg
        !
        integer :: nzEarth       ! (JUST SG ????)
        !
        contains
            !
            !> Grid abstract interfaces
            procedure( interface_number_of_edges ), deferred, public :: numberOfEdges
            procedure( interface_number_of_faces ), deferred, public :: numberOfFaces
            procedure( interface_number_of_nodes ), deferred, public :: numberOfNodes
            procedure( interface_number_of_cells ), deferred, public :: numberOfCells
            !
            procedure( interface_allocate_dim_grid ), deferred, public :: allocateDim
            !
            procedure( interface_setup_grid ), deferred, public :: setup
            !
            procedure( interface_slice_1d_grid ), deferred, public :: slice1D
            procedure( interface_slice_2d_grid ), deferred, public :: slice2D
            !
            procedure( interface_write_grid ), deferred, public :: write
            !
            !> Base Grid methods
            procedure, public :: baseInit => initialize_Grid
            procedure, public :: baseDealloc => deallocate_Grid
            !
            procedure, public :: getNGrids => getNGrids_Grid
            !
            procedure, public :: getDimensions => getDimensions_Grid
            !
            procedure, public :: setOrigin => setOrigin_Grid
            procedure, public :: getOrigin => getOrigin_Grid
            !
            procedure, public :: setRotation => setRotation_Grid
            procedure, public :: getRotation => getRotation_Grid
            !
            procedure, public :: setGeometry => setGeometry_Grid
            procedure, public :: getGeometry => getGeometry_Grid
            !
            procedure, public :: setCellSizes => setCellSizes_Grid
            procedure, public :: getCellSizes => getCellSizes_Grid
            !
            procedure, public :: create => create_Grid
            !
            procedure, public :: setupAirLayers => setupAirLayers_Grid
            procedure, public :: updateAirLayers => updateAirLayers_Grid
            !
            procedure, public :: nEdges => nEdges_Grid
            !
            procedure, public :: nFaces => nFaces_Grid
            !
            procedure, public :: setLimits => setLimits_Grid
            procedure, public :: gridIndex => gridIndex_Grid
            procedure, public :: vectorIndex => vectorIndex_Grid
            !
    end type Grid_t
    !
    !> Public Global Generic Grid object
    class( Grid_t ), allocatable, target :: main_grid
    !
    !> Details needed to unambiguosly compute and/or store the air layers;
    !> method options are: mirror; fixed height; read from file
    !> for backwards compatibility, all of the defaults are set to what
    !> was previously hard coded (AK; May 19, 2017)
    !> For backwards compatibility, default is "mirror 10 3. 30."
    !>    GDE 12/17/21 : new method "undefined" is default -- force
    !>        explicit setting of air layers (can still hard code a default
    !>         in the driver program!)
    !
    type :: TAirLayers
        !
        character(:), allocatable :: method
        integer :: nz
        real( kind=prec ) :: maxHeight, minTopDz, alpha
        real( kind=prec ), allocatable, dimension(:) :: dz
        logical :: is_allocated
        !
    end type TAirLayers
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_number_of_edges( self, n_xedge, n_yedge, n_zedge )
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            integer, intent( out ) :: n_xedge, n_yedge, n_zedge
        end subroutine interface_number_of_edges
        !
        !> No interface subroutine briefing
        subroutine interface_number_of_faces( self, n_xface, n_yface, n_zface ) 
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            integer, intent( out ) :: n_xface, n_yface, n_zface
        end subroutine interface_number_of_faces
        !
        !> No interface function briefing
        function interface_number_of_nodes( self ) result(n)
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            integer :: n
        end function interface_number_of_nodes
        !
        !> No interface function briefing
        function interface_number_of_cells( self ) result(n)
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            integer :: n
        end function interface_number_of_cells
        !
        !> No interface subroutine briefing
        !
        subroutine interface_allocate_dim_grid( self )
            import :: Grid_t
            class( Grid_t ), intent( inout ) :: self
        end subroutine interface_allocate_dim_grid
        !
        !> No interface subroutine briefing
        !
        subroutine interface_setup_grid( self, origin )
            import :: Grid_t, prec
            class( Grid_t ), intent( inout ) :: self
            real( kind=prec ), intent( in ), optional :: origin(3)
        end subroutine interface_setup_grid
        !
        !> No interface function briefing
        !
        function interface_slice_1d_grid( self ) result( g1D )
            import :: Grid_t, Grid1D_t
            class( Grid_t ), intent( in ) :: self
            type( Grid1D_t ) :: g1D
        end function interface_slice_1d_grid
        !
        !> No interface function briefing
        !
        function interface_slice_2d_grid( self ) result( g2D )
            import :: Grid_t, Grid2D_t
            class( Grid_t ), intent( in ) :: self
            type( Grid2D_t ) :: g2D
        end function interface_slice_2d_grid
        !
        !> No interface function briefing
        !
        subroutine interface_write_grid( self )
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
        end subroutine interface_write_grid
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    !
    subroutine initialize_Grid( self )
        implicit none
        !
        class( Grid_t ), intent( inout ) :: self
        !
        self%n_grids = 0
        !
        self%geometry = REGION
        !
        self%ox = R_ZERO
        self%oy = R_ZERO
        self%oz = R_ZERO
        !
        self%rotDeg = R_ZERO 
        !
        self%nx = 0
        self%ny = 0
        self%nz = 0
        !
        self%nzAir = 0
        self%nzEarth = 0
        !
        self%is_allocated = .FALSE.
        !
    end subroutine initialize_Grid
    !
    !> No subroutine briefing
    !
    subroutine deallocate_Grid( self )
        implicit none
        !
        class( Grid_t ), intent( inout ) :: self
        !
        self%nx = 0
        self%ny = 0
        self%nz = 0
        !
        if( allocated( self%dx ) ) deallocate( self%dx )
        if( allocated( self%dy ) ) deallocate( self%dy )
        if( allocated( self%dz ) ) deallocate( self%dz )
        !
        self%is_allocated = .FALSE.
        !
    end subroutine deallocate_Grid
    !
    !> No subroutine briefing
    !
    function getNGrids_Grid( self ) result( n_grids )
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
        !
        integer :: n_grids
        !
        n_grids = self%n_grids
        !
    end function getNGrids_Grid
    !
    !> No subroutine briefing
    !
    subroutine getDimensions_Grid( self, nx, ny, nz, nzAir )
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
        integer, intent( out ) :: nx, ny, nz, nzAir
        !
        nx = self%nx
        ny = self%ny
        nz = self%nz
        nzAir = self%nzAir
        !
    end subroutine getDimensions_Grid
    !
    !> No subroutine briefing
    !
    subroutine setOrigin_Grid( self, ox, oy, oz )
        implicit none
        !
        class( Grid_t ), intent(inout) :: self
        real( kind=prec ), intent( in ) :: ox, oy, oz
        !
        self%ox = ox
        self%oy = oy
        self%oz = oz
        !
    end subroutine setOrigin_Grid
    !
    !> No subroutine briefing
    !
    subroutine getOrigin_Grid( self, ox, oy, oz )
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
        real( kind=prec ), intent( out ) :: ox, oy, oz
        !
        ox = self%ox
        oy = self%oy
        oz = self%oz
        !
    end subroutine getOrigin_Grid
    !
    !> No subroutine briefing
    !
    subroutine setRotation_Grid( self, rotDeg )
        implicit none
        !
        class( Grid_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rotDeg
        !
        self%rotDeg = rotDeg
        !
    end subroutine setRotation_Grid
    !
    !> No subroutine briefing
    !
    function getRotation_Grid( self ) result(rotDeg)
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
        !
        real( kind=prec ) :: rotDeg
        !
        rotDeg = self%rotDeg
        !
    end function getRotation_Grid
    !
    !> No subroutine briefing
    !
    subroutine setGeometry_Grid( self, s )
        implicit none
        !
        class( Grid_t ), intent( inout ) :: self
        character(*), intent( in ) :: s
        !
        self%geometry = s
        !
    end subroutine setGeometry_Grid
    !
    !> No function briefing
    !
    function getGeometry_Grid( self ) result(s)
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
        !
        character(80) :: s
        !
        s = self%geometry
        !
    end function getGeometry_Grid
    !
    !> No subroutine briefing
    !
    subroutine setCellSizes_Grid( self, dx, dy, dz )
        implicit none
        !
        class( Grid_t ), intent( inout ) :: self
        real( kind=prec ), dimension(:), intent( in ) :: dx, dy, dz
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setCellSizes_Grid > Grid not allocated." )
        endif
        !
        !> Check dimensions
        if( ( size( dx ) .NE. size( self%dx ) ) .OR. &
            ( size( dy ) .NE. size( self%dy ) ) .OR. &
            ( size( dz ) .NE. size( self%dz ) ) ) then
            !
            call errStop( "setCellSizes_Grid > Incompatible sizes for cell arrays." )
            !
        endif
        !
        self%dx = dx
        self%dy = dy
        self%dz = dz
        !
    end subroutine setCellSizes_Grid
    !
    !> No subroutine briefing
    subroutine getCellSizes_Grid( self, dx, dy, dz )
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
        real( kind=prec ), intent( out ) :: dx(:), dy(:), dz(:)
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "getCellSizes_Grid > Grid not allocated." )
        endif
        !
        !> Check dimensions
        if( ( size( dx ) .NE. size( self%dx ) ) .OR. &
            ( size( dy ) .NE. size( self%dy ) ) .OR. &
            ( size( dz ) .NE. size( self%dz ) ) ) then
            !
            call errStop( "getCellSizes_Grid > Incompatible sizes for cell arrays." )
            !
        endif
        !
        dx = self%dx
        dy = self%dy
        dz = self%dz
        !
    end subroutine getCellSizes_Grid
    !
    !> No subroutine briefing
    subroutine create_Grid( self, nx, ny, nzAir, nzEarth )
        implicit none
        !
        class( Grid_t ), intent( inout ) :: self
        integer, intent( in ) :: nx, ny, nzAir, nzEarth
        !
        integer :: nz
        !
        nz = nzEarth + nzAir
        !
        self%nzAir = nzAir
        self%nzEarth = nzEarth
        !
        self%nx = nx
        self%ny = ny
        self%nz = nz
        !
        call self%allocateDim
        !
    end subroutine create_Grid
    !
    !> setupAirLayers computes the Dz in the airlayers structure
    !> using the grid to get the top layers Dz;
    !> all values expected in km on input
    !> For backwards compatibility, default is "mirror 10 3. 30."
    !> but the use of "fixed height 12 1000" is recommended
    !
    subroutine setupAirLayers_Grid( self, airLayers, method, &
                                    nzAir, maxHeight, minTopDz, &
                                    alpha, dzAir )
        implicit none
        !
        class( Grid_t ), intent( inout ) :: self
        type( TAirLayers ), intent( inout ) :: airlayers
        character(:), allocatable, intent( in ), optional :: method
        integer, intent( in ), optional :: nzAir
        real( kind=prec ), intent( in ), optional :: maxHeight, minTopDz, alpha
        real( kind=prec ), intent( in ), optional, pointer :: dzAir
        !
        integer :: ix, iy, iz, i, j
        integer :: status
        real( kind=prec ) :: z1_log, dlogz, z_log, height1, height2
        !
        airLayers%is_allocated = .FALSE.
        !
        if( present( method ) ) then
            airlayers%method = method
        else
            airlayers%method = "fixed height"
        endif
        !
        if( present( nzAir ) ) then
            airlayers%nz = nzAir
        else
            airlayers%nz = model_n_air_layer
        endif
        !
        if( .NOT. ( index( airLayers%method, "read from file" ) > 0 ) ) then
            if( airLayers%is_allocated ) then
                deallocate( airlayers%dz )
            endif
            allocate( airLayers%dz( airLayers%nz ) )
            airLayers%is_allocated = .TRUE.
        endif
        !
        if( present( maxHeight ) ) then
            airLayers%maxHeight = 1000. * maxHeight
        else
            airLayers%maxHeight = model_max_height
        endif
        
        if( present( minTopDz ) ) then
            airLayers%minTopDz = 1000. * minTopDz
        else
            airLayers%minTopDz = 100.0
        endif
        
        if( present( alpha ) ) then
            airLayers%alpha = alpha
        else
            airLayers%alpha = 3.
        endif
        !
        if( index( airLayers%method, "mirror" ) > 0 ) then
            !
            !> Following is Kush"s approach to setting air layers:
            !> mirror imaging the dz values in the air layer with respect to
            !> earth layer as far as we can using the following formulation
            !> air layer(bottom:top) = (alpha)^(j-1) * earth layer(top:bottom)
            !
            do iz = airLayers%nz, 1, -1
                j = airLayers%nz - iz + 1
                airLayers%dz(iz) = ((airLayers%alpha)**(j - 1) ) * self%dz(self%nzAir + j)
            enddo
            !
            !> The topmost air layer has to be at least 30 km
            if(airLayers%dz(1).lt.airLayers%minTopDz) then
                airLayers%dz(1) = airLayers%minTopDz
            endif

        elseif( index(airLayers%method, "fixed height") > 0 ) then 
            !
            !> ON IMPLEMENTATION
            z1_log = log10( self%Dz( self%NzAir + 1 ) )
            dlogz = ( log10( airlayers%maxHeight ) - z1_log ) / ( airlayers%Nz )
            !
            z_log = z1_log
            do iz = airlayers%Nz, 1, -1
                airlayers%Dz(iz) = 10.**(z_log+dlogz) - 10.**(z_log)
                z_log = z_log + dlogz
            enddo
            !
            !> OTHER IMPLEMENTATION
            !
            ! z1_log = log10(self%dz(self%nzAir + 1) )
            ! dlogz = (log10(airLayers%maxHeight) - z1_log)/(airLayers%nz-1)
            ! z_log = z1_log
            ! height1 = 10.**z1log
            ! airLayers%dz(airLayers%Nz) = height1
            ! do iz = airLayers%Nz-1, 1, -1
                ! z_log = z_log + dlogz
                ! height2 = 10.**z_log 
                ! airLayers%dz(iz) = height2-height1
                ! height1 = height2
            ! enddo
            !
        elseif( index( airLayers%method, "read from file" ) > 0 ) then
            !
            !> Air layers have been read from file and are
            !> already stored in Dz, so only need to reallocate
            !> if passing a new array to it.
            !
            if( present( dzAir) ) then
                !
                if( airLayers%is_allocated) then
                    deallocate( airLayers%dz )
                endif
                !
                allocate( airLayers%dz(airLayers%nz) )
                airLayers%dz = dzAir
                !
            endif
            !
        endif
        !
    end subroutine setupAirLayers_Grid
    !
    !> Procedure updateAirLayers_Grid
    !> Assumes that the grid is already defined, and merely
    !> includes the new air layers in the grid.
    !
    subroutine updateAirLayers_Grid( self, nzAir, dzAir )
        implicit none
        !
        class( Grid_t ), intent( inout ) :: self
        integer, intent( in ) :: nzAir
        real( kind=prec ), intent( in ) :: dzAir(:)
        !
        integer :: nzAir_old, nzEarth_old
        integer :: nx_old, ny_old, nz_old
        real( kind=prec ), allocatable :: dx_old(:), dy_old(:), dz_old(:)
        real( kind=prec ) :: ox_old, oy_old, oz_old
        real( kind=prec ) :: rotDeg_old
        character( len=80 ) :: geometry_old
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "updateAirLayers_Grid > Grid not allocated." )
        endif
        !
        nx_old = self%nx
        ny_old = self%ny
        nz_old = self%nz
        nzAir_old = self%nzair
        nzEarth_old = self%nzEarth
        dx_old = self%dx
        dy_old = self%dy
        dz_old = self%dz
        geometry_old = self%getGeometry()
        !
        ox_old = self%ox
        oy_old = self%oy
        oz_old = self%oz
        !
        rotdeg_old = self%rotdeg
        !
        call self%baseDealloc
        call self%create( nx_old, ny_old, nzAir, nzEarth_old )
        !
        !> Set air layers to dzAir values and copy the rest
        self%dz(1:nzAir) = dzAir
        self%dz(nzAir+1:self%nz) = dz_old( nzAir_old+1:nz_old )
        self%dy = dy_old
        self%dx = dx_old
        !
        call self%setOrigin( ox_old, oy_old, oz_old )
        call self%setRotation( rotdeg_old )
        call self%setGeometry( geometry_old )
        !
        !> setup the rest of the grid from scratch
        call self%setup
        !
    end subroutine updateAirLayers_Grid
    !
    !> No subroutine briefing
    !
    subroutine nEdges_Grid( self, n_xedge, n_yedge, n_zedge )
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
        integer, intent( out )  :: n_xedge, n_yedge, n_zedge
        integer :: nx, ny, nz

        call self%setLimits( XEDGE, nx, ny, nz )
        n_xedge = nx*ny*nz
        !
        call self%setLimits( YEDGE, nx, ny, nz )
        n_yedge = nx*ny*nz
        !
        call self%setLimits( ZEDGE, nx, ny, nz )
        n_zedge = nx*ny*nz
        !
    end subroutine nEdges_Grid
    !
    !> No subroutine briefing
    !
    subroutine nFaces_Grid( self, n_xface, n_yface, n_zface )
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
        integer, intent(out) :: n_xface, n_yface, n_zface
        integer :: nx, ny, nz
        !
        call self%setLimits( XFACE, nx, ny, nz )
        n_xface = nx*ny*nz
        !
        call self%setLimits( YFACE, nx, ny, nz )
        n_yface = nx*ny*nz
        !
        call self%setLimits( ZFACE, nx, ny, nz )
        n_zface = nx*ny*nz
        !
    end subroutine nFaces_Grid
    !
    !> No subroutine briefing
    subroutine setLimits_Grid( self, node_type, nx, ny, nz )
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
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
                call errStop( "setLimits_Grid > Undefined node_type ["//node_type//"]" )
                !
        end select
        !
    end subroutine setLimits_Grid
    !
    !> Based on matlab method of same name in class Grid_t3D
    !> IndVec is the index within the list of nodes of a fixed type
    !> e.g., among the list of y-Faces.     An offset needs to be
    !> added to get index in list of all faces(for example).
    !
    subroutine gridIndex_Grid( self, node_type, ind_vec, i, j, k )
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
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
            call errStop( "gridIndex_Grid > Size of 'ind_vec' and 'i' do not agree." )
        endif
        !
        if( nVec .NE. size(j) ) then
            call errStop( "gridIndex_Grid > Size of 'ind_vec' and 'j' do not agree." )
        endif
        !
        if( nVec .NE. size(k) ) then
            call errStop( "gridIndex_Grid > Size of 'ind_vec' and 'k' do not agree." )
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
    end subroutine gridIndex_Grid
    !
    !> vectorIndex
    !
    !> Based on matlab method of same name in class Grid_t3D
    !> returned array IndVec gives numbering of nodes within
    !> the list for node_type; need to add an offset for position
    !> in full list of all faces or edges(not nodes and cells).
    !
    subroutine vectorIndex_Grid( self, node_type, i, j, k, ind_vec )
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
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
            call errStop( "vectorIndex_Grid > Size of 'ind_vec' and 'i' do not agree." )
        endif
        !
        if( nVec .NE. size(J) ) then
            call errStop( "vectorIndex_Grid > Size of 'ind_vec' and 'j' do not agree." )
        endif
        !
        if( nVec .NE. size(K) ) then
            call errStop( "vectorIndex_Grid > Size of 'ind_cec' and 'k' do not agree." )
        endif
        !
        nxy = nx*ny
        do ii = 1, nVec
            ind_vec(ii) =(K(ii) - 1) * nxy +(j(ii) - 1) * nx + i(ii)
        enddo
        !
    end subroutine vectorIndex_Grid
    !
end module Grid
