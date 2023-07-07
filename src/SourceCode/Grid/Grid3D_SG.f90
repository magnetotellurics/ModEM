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
            procedure, public :: limits => limits_Grid3D_SG
            procedure, public :: length => length_Grid3D_SG
            !
            !
            procedure, public :: create => create_Grid3D_SG
            procedure, public :: allocateDim => allocateDim_Grid3D_SG
            procedure, public :: setup => setup_Grid3D_SG
            procedure, public :: setupAirLayers => setupAirLayers_Grid3D_SG
            procedure, public :: updateAirLayers => updateAirLayers_Grid3D_SG
            !
            procedure, public :: setCellSizes => setCellSizes_Grid3D_SG
            procedure, public :: getCellSizes => getCellSizes_Grid3D_SG
            !
            procedure, public :: slice1D => slice1D_Grid3D_SG
            procedure, public :: slice2D => slice2D_Grid3D_SG
            !
    end type Grid3D_SG_t
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
    interface Grid3D_SG_t
        module procedure Grid3D_SG_t_ctor
    end interface Grid3D_SG_t
    !
contains
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
    !> Class constructor for simple tensor product grid
    !> Usage obj = Grid_t3D(Dx,Dy,Dz,Nza)
    !> Dx, Dy, Dz are cell dimensions for x, y, z direction
    !> Nza is number of air layers to allow (included in Dz)
    !
    function Grid3D_SG_t_ctor( nx, ny, nzAir, nzEarth, dx, dy, dz ) result( self )
        implicit none
        !
        integer, intent( in ) :: nx, ny, nzAir, nzEarth
        real( kind=prec ), dimension(:), intent( in ) :: dx, dy, dz
        !
        type( Grid3D_SG_t ) :: self
        !
        !write( *, * ) "Constructor Grid3D_SG_t"
        !
        call self%init
        !
        call self%create( nx, ny, nzAir, nzEarth )
        !
        call self%setCellSizes( dx, dy, dz )
        !
        call self%setup
        !
    end function Grid3D_SG_t_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine dealloc().
    subroutine Grid3D_SG_dtor( self )
        implicit none
        !
        type( Grid3D_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor Grid3D_SG_t", self%nx, self%ny, self%nz
        !
        call self%dealloc
        !
    end subroutine Grid3D_SG_dtor
    !
    !> No subroutine briefing
    subroutine create_Grid3D_SG( self, nx, ny, nzAir, nzEarth )
        implicit none
        !
        class( Grid3D_SG_t ), intent( inout ) :: self
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
    end subroutine create_Grid3D_SG
    !
    !> No subroutine briefing
    subroutine allocateDim_Grid3D_SG( self )
        implicit none
        !
        class( Grid3D_SG_t ), intent(inout) :: self
        !
        integer :: nx, ny, nz
        !
        if( self%is_allocated ) call self%dealloc
        !
        nx = self%nx
        ny = self%ny
        nz = self%nz
        !
        allocate( self%dx(nx) )
        allocate( self%dy(ny) )
        allocate( self%dz(nz) )
        !
        !> dx_inv = 1/ dx and similarly for dy_inv and dz_inv
        allocate( self%dx_inv(nx) )
        allocate( self%dy_inv(ny) )
        allocate( self%dz_inv(nz) )
        !
        !> del_x, del_y, and del_z are the distances between
        !> the electrical field defined on the center of the
        !> edges in x, y, and z axes, respectively.
        allocate( self%del_x(nx + 1) )
        allocate( self%del_y(ny + 1) )
        allocate( self%del_z(nz + 1) )
        !
        allocate( self%del_x_inv(nx + 1) )
        allocate( self%del_y_inv(ny + 1) )
        allocate( self%del_z_inv(nz + 1) )
        !
        !> x_edge is the array for cumulative distance of the edge
        !> for each grid (starting from the coordinate axes) with
        !> dimensions nx + 1.
        !> x_center is the array for cumulative distance of the center
        !> for each grid (starting from the coordinate axes) with
        !> dimension n.
        !> y_edge, y_center, z_edge, z_center are analagous arrays for
        !> other directions.
        !
        allocate( self%x_edge(nx + 1) )
        allocate( self%y_edge(ny + 1) )
        allocate( self%z_edge(nz + 1) )
        allocate( self%x_center(nx) )
        allocate( self%y_center(ny) )
        allocate( self%z_center(nz) )
        !
        self%is_allocated = .TRUE.
        !
    end subroutine allocateDim_Grid3D_SG
    !
    !> setup does calculations for grid geometry, which cannot be done
    !> until dx, dy, dz, and the origin are set.
    subroutine setup_Grid3D_SG(self, origin)
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
        if(present(origin) ) then
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
        self%del_x = self%del_x/2.0
        !
        self%del_y(1) = self%dy(1)
        do iy = 2, self%ny
            self%del_y(iy) = self%dy(iy - 1) + self%dy(iy)
        enddo
        !
        self%del_y(self%ny + 1) = self%dy(self%ny)
        self%del_y = self%del_y/2.0
        !
        self%del_z(1) = self%dz(1)
        do iz = 2, self%nz
            self%del_z(iz) = self%dz(iz - 1) + self%dz(iz)
        enddo
        !
        self%del_z(self%nz + 1) = self%dz(self%nz)
        self%del_z = self%del_z/2.0
        !
        self%del_x_inv = 1/self%del_x
        self%del_y_inv = 1/self%del_y
        self%del_z_inv = 1/self%del_z
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
    !> setupAirLayers computes the Dz in the airlayers structure
    !> using the grid to get the top layers Dz;
    !> all values expected in km on input
    !> For backwards compatibility, default is "mirror 10 3. 30."
    !> but the use of "fixed height 12 1000" is recommended
    !
    subroutine setupAirLayers_Grid3D_SG( self, airLayers, method, &
                                        nzAir, maxHeight, minTopDz, &
                                        alpha, dzAir )
        implicit none
        !
        class( Grid3D_SG_t ), intent( inout ) :: self
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
        if(index(airLayers%method, "mirror") > 0) then
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

        else if(index(airLayers%method, "fixed height") > 0) then 
            !
            !> ON IMPLEMENTATION
            z1_log = log10( self%Dz( self%NzAir + 1 ) )
            dlogz = ( log10( airlayers%maxHeight ) - z1_log ) / ( airlayers%Nz )

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
        else if(index(airLayers%method, "read from file") > 0) then
            !
            !> Air layers have been read from file and are
            !> already stored in Dz, so only need to reallocate
            !> if passing a new array to it.
            !
            if( present(dzAir) ) then
                if( airLayers%is_allocated) then
                     deallocate( airLayers%dz )
                endif
                allocate( airLayers%dz(airLayers%nz) )
                airLayers%dz = dzAir
            endif
        endif
        !
    end subroutine setupAirLayers_Grid3D_SG
    !
    !> No subroutine briefing
    subroutine DeallocateAirLayers_Grid3D_SG(self, airLayers)
        implicit none
        !
        class( Grid3D_SG_t ), intent( inout ) :: self
        type( TAirLayers ), intent( inout ) :: airLayers
        !
        deallocate( airLayers%dz )
        !
    end subroutine DeallocateAirLayers_Grid3D_SG
    !
    !> Procedure updateAirLayers_Grid3D_SG
    !> Assumes that the grid is already defined, and merely
    !> includes the new air layers in the grid.
    subroutine updateAirLayers_Grid3D_SG( self, nzAir, dzAir )
        implicit none
        !
        class( Grid3D_SG_t ), intent(inout) :: self
        integer, intent( in ) :: nzAir
        real( kind=prec ), intent( in ) :: dzAir(:)
        !
        integer :: nzAir_old, nzEarth_old
        integer :: nx_old, ny_old, nz_old
        real( kind=prec ), allocatable :: dx_old(:), dy_old(:), dz_old(:)
        real( kind=prec ) :: ox_old, oy_old, oz_old
        real( kind=prec ) :: rotDeg_old
        character(len = 80) :: geometry_old
        !
        if(.NOT.self%is_allocated) then
             stop "Error: updateAirLayers_Grid3D_SG > Grid not allocated."
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
        call self%dealloc
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
    end subroutine updateAirLayers_Grid3D_SG
    !
    !> No subroutine briefing
    subroutine setCellSizes_Grid3D_SG(self, dx, dy, dz)
        implicit none
        !
        class( Grid3D_SG_t ), intent(inout) :: self
        real( kind=prec ), dimension(:), intent( in ) :: dx, dy, dz
        !
        if( .NOT. self%is_allocated ) then
             stop "Error: setCellSizes_Grid3D_SG > Grid not allocated."
        endif
        !
        !> Check dimensions
        if((size(dx).NE.size(self%dx) ) .OR. &
                 (size(dy).NE.size(self%dy) ) .OR. &
                 (size(dz).NE.size(self%dz) )) then
             write( *, * ) "Error:Grid3D_SG_t:setCellSizes:"
             stop "    Incompatible sizes for cell arrays."
        endif
        !
        self%dx = dx
        self%dy = dy
        self%dz = dz
        !
    end subroutine setCellSizes_Grid3D_SG
    !
    !> No subroutine briefing
    subroutine getCellSizes_Grid3D_SG(self, dx, dy, dz)
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        real( kind=prec ), intent( out ) :: dx(:), dy(:), dz(:)
        !
        if( .NOT. self%is_allocated ) then
             stop "Error: getCellSizes_Grid3D_SG > Grid not allocated."
        endif
        !
        !> Check dimensions
        if((size(dx).NE.size(self%dx) ) .OR. &
                 (size(dy).NE.size(self%dy) ) .OR. &
                 (size(dz).NE.size(self%dz) )) then
             stop "Error: getCellSizes_Grid3D_SG > Incompatible sizes for cell arrays."
        endif
        !
        dx = self%dx
        dy = self%dy
        dz = self%dz
        !
    end subroutine getCellSizes_Grid3D_SG
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
        call self%limits( XEDGE, nx, ny, nz )
        n_xedge = nx * ny * nz
        !
        call self%limits( YEDGE, nx, ny, nz )
        n_yedge  = nx * ny * nz
        !
        call self%limits( ZEDGE, nx, ny, nz )
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
        
        call self%limits( XFACE, nx, ny, nz )
        n_xface = nx * ny * nz
        
        call self%limits( YFACE, nx, ny, nz )
        n_yface = nx * ny * nz
        
        call self%limits( ZFACE, nx, ny, nz )
        n_zface = nx * ny * nz
        
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
        
        n = (self%nx + 1)*(self%ny + 1)*(self%nz + 1)

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
        stop "Error: numberOfCells_Grid3D_SG not implemented"
        !
    end function numberOfCells_Grid3D_SG
    !
    !> Based on matlab method of same name in class Grid_t3D
    !> IndVec is the index within the list of nodes of a fixed type
    !> e.g., among the list of y-Faces.     An offset needs to be
    !> added to get index in list of all faces (for example).
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

        call self%limits(node_type, nx, ny, nz)
        nVec = size(ind_vec)
        
        if(nVec.NE.size(i) ) then
             stop "Size of 'ind_vec' and 'i' do not agree."
        endif
        
        if(nVec.NE.size(j) ) then
             stop "Size of 'ind_vec' and 'j' do not agree."
        endif
        
        if(nVec.NE.size(k) ) then
             stop "Size of 'ind_vec' and 'k' do not agree."
        endif
        
        rNxy = float(nx*ny)
        rNx = float(nx)
        
        do ii = 1, nVec
            i(ii) = mod(ind_vec(ii), nx)
            j(ii) = mod(ceiling(float(ind_vec(ii) )/rNx), ny)
            k(ii) = ceiling(float(ind_vec(ii) )/rNxy)
        enddo
        
        where(i.EQ.0) i = nx
        where(j.EQ.0) j = ny
        where(k.EQ.0) k = nz

    end subroutine gridIndex_Grid3D_SG
    !
    !> vectorIndex
    !
    !> Based on matlab method of same name in class Grid_t3D
    !> returned array IndVec gives numbering of nodes within
    !> the list for node_type; need to add an offset for position
    !> in full list of all faces or edges (not nodes and cells).
    subroutine vectorIndex_Grid3D_SG( self, node_type, i, j, k, ind_vec )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        character(*), intent( in ) :: node_type
        integer, intent( in ) :: i(:), j(:), k(:)
        integer, intent( out ) :: ind_vec(:)
        !
        integer :: nx, ny, nz, nxy, nVec, ii
        
        call self%limits(node_type, nx, ny, nz)
        
        nVec = size(ind_vec)
        
        if(nVec.NE.size (i) ) then
             stop "Size of 'ind_vec' and 'i' do not agree."
        endif
        
        if(nVec.NE.size (J) ) then
             stop "Size of 'ind_vec' and 'j' do not agree."
        endif
        
        if(nVec.NE.size (K) ) then
             stop "Size of 'ind_cec' and 'k' do not agree."
        endif
        
        nxy = nx*ny
        do ii = 1, nVec
             ind_vec(ii) = (K(ii) - 1) * nxy + (j(ii) - 1) * nx + i(ii)
        enddo
        
    end subroutine vectorIndex_Grid3D_SG
    !
    !> No subroutine briefing
    subroutine limits_Grid3D_SG( self, node_type, nx, ny, nz )
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
                write( *, * ) "Error: limits_Grid3D_SG > Undefined node_type [", node_type, "]"
                stop
                !
        end select
        !
    end subroutine limits_Grid3D_SG
    !
    !> No subroutine briefing
    !
    function length_Grid3D_SG( self ) result(n)
        class( Grid3D_SG_t ), intent( in ) :: self
        integer :: n
        !
        n = self%nx * self%ny * self%nz
        !
    end function length_Grid3D_SG
    !
end module Grid3D_SG
