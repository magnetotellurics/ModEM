!
!> This file is part of the ModEM modeling and inversion package.
!>
!> LICENSING information
!
!> Copyright(C) 2020 ModEM research group.
!> Contact: http://
!>
!> GNU General Public License Usage
!> This file may be used under the terms of the GNU
!> General Public License version 3.0 as published by the Free Software
!> Foundation and appearing in the file LICENSE.GPL included in the
!> packaging of this file.  Please review the following information to
!> ensure the GNU General Public License version 3.0 requirements will be
!> met: http://www.gnu.org/copyleft/gpl.html.
!>
!> SUMMARY
!>
!> No module briefing
!
module Grid3D_MR
    !
    use Utilities
    use Grid3D_SG
    !
    type, extends( Grid_t ) :: Grid3D_MR_t
        !
        !> MR properties
        integer, allocatable, dimension(:,:) :: coarseness, z_limits
        !
        integer, allocatable, dimension(:,:) :: n_active_edge, n_active_face
        !
        integer, allocatable, dimension(:) :: n_active_node, n_active_cell, cs
        !
        type( Grid3D_SG_t ), allocatable, dimension(:) :: sub_grid
        !
        !real( kind=prec ), allocatable, dimension(n_grids) :: z_top
        !
        logical :: is_initialized
        !
        contains
            !
            procedure, public :: numberOfEdges => numberOfEdges_Grid3D_MR
            procedure, public :: numberOfFaces => numberOfFaces_Grid3D_MR
            procedure, public :: numberOfNodes => numberOfNodes_Grid3D_MR
            procedure, public :: numberOfCells => numberOfCells_Grid3D_MR
            !
            procedure, public :: reduceActive => reduceActive_Grid3D_MR
            !
            procedure, public :: setActivelimits => setActiveLimits_Grid3D_MR
            !
            procedure, public :: nActive => nActive_Grid3D_MR
            !
            procedure, public :: active => active_Grid3D_MR
            !
            procedure, public :: allocateDim => allocateDim_Grid3D_MR
            !
            procedure, public :: setup => setup_Grid3D_MR
            !
            procedure, public :: slice1D => slice1D_Grid3D_MR
            procedure, public :: slice2D => slice2D_Grid3D_MR
            !
            procedure, public :: setSubGrid, setupMR
            !
            procedure, public :: updateAirLayers => updateAirLayers_Grid3D_MR
            !
            procedure, public :: write => write_Grid3D_MR
            !
    end type Grid3D_MR_t
    !
    public :: repMat
    !
    interface Grid3D_MR_t
        module procedure Grid3D_MR_t_ctor
    end interface Grid3D_MR_t
    !
contains
    !
    !> Creates and initialize a MR Grid from a standard Grid.
    !
    function Grid3D_MR_t_ctor( nx, ny, nzAir, nzEarth, dx, dy, dz, layers ) result( self )
        implicit none
        !
        integer, intent( in ) :: nx, ny, nzAir, nzEarth
        real( kind=prec ), dimension(:), intent( in ) :: dx, dy, dz
        integer, allocatable, dimension(:), intent( in ) :: layers
        !
        type( Grid3D_MR_t ) :: self
        !  
        integer :: i
        !
        !write( *, * ) "Constructor Grid3D_MR_t", size(dx), size(dy), size(dz)
        !
        call self%baseInit
        !
        call self%create( nx, ny, nzAir, nzEarth )
        !
        call self%setCellSizes( dx, dy, dz )
        !
        self%n_grids = size( layers ) / 2
        !
        self%cs = layers
        !
        if( self%n_grids > 0 ) then
            !
            allocate( self%n_active_edge( 3, self%n_grids ) )
            !
            allocate( self%n_active_face( 3, self%n_grids ) )
            !
            allocate( self%n_active_node( self%n_grids ) )
            !
            allocate( self%n_active_cell( self%n_grids ) )
            !
            allocate( self%z_limits( 2, self%n_grids ) )
            !
            allocate( self%coarseness( self%n_grids, 4 ) )
            !
            allocate( self%sub_grid( self%n_grids ) )
            !
            self%coarseness = 0
            !
            do i = 0, self%n_grids - 1
                !
                self%coarseness( i + 1, 1 ) = layers( 2 * i + 1 )
                !
                self%coarseness( i + 1, 2 ) = layers( 2 * i + 2 )
                !
            enddo
            !
            call self%setup
            !
            call self%setupMR
            !
            self%is_initialized = .TRUE.
            !
            !write( *, * ) "CoarseMAT Row1 (Coarse): [", self%coarseness(:,1), "]"
            !write( *, * ) "CoarseMAT Row2 (Layers): [", self%coarseness(:,2), "]"
            !write( *, * ) "CoarseMAT Row3 (iStart): [", self%coarseness(:,3), "]"
            !write( *, * ) "CoarseMAT Row4 ( iEnd ): [", self%coarseness(:,4), "]"
            !
        else
            self%is_initialized = .FALSE.
        endif
        !
    end function Grid3D_MR_t_ctor
    !
    !> No subroutine briefing
    subroutine allocateDim_Grid3D_MR( self )
        implicit none
        !
        class( Grid3D_MR_t ), intent( inout ) :: self
        !
        if( self%is_allocated ) call self%baseDealloc
        !
        allocate( self%dx( self%nx ) )
        allocate( self%dy( self%ny ) )
        allocate( self%dz( self%nz ) )
        !
        self%is_allocated = .TRUE.
        !
    end subroutine allocateDim_Grid3D_MR
    !
    !> setup does calculations for grid geometry, which cannot be done
    !> until dx, dy, dz, and the origin are set.
    subroutine setup_Grid3D_MR( self, origin )
        implicit none
        !
        class( Grid3D_MR_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ), optional :: origin(3)
        !
        real( kind=prec ) :: ox, oy, oz
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
    end subroutine setup_Grid3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setupMR( self, origin )
        implicit none
        !
        class( Grid3D_MR_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ), optional :: origin(3)
        !
        integer :: i, i1, i2, last
        real( kind=prec ) :: ddz_interface
        !
        ! Check if coarseness parameters are consistent with
        if( sum( self%coarseness(:, 2) ) /= size( self%dz ) ) then
            call errStop( "setupMR > Inconsistent grid coarseness parameter!" )
        endif
        !
        do i = 1, self%n_grids
            !
            i1 = sum( self%coarseness( 1 : i, 2 ) ) - self%coarseness( i, 2 ) + 1
            i2 = sum( self%coarseness( 1 : i, 2 ) )
            !
            call self%setSubGrid( i, i1, i2 )
            !
        enddo
        !
        ! Correct dual edge lengths at interfaces in each subgrid
        do i = 2, self%n_grids
            !
            last = size( self%sub_grid(i - 1)%dz )
            ddz_interface =(self%sub_grid(i - 1)%dz(last) + &
            self%sub_grid(i)%dz(1)) / 2
            !
            last = size( self%sub_grid(i - 1)%del_z)
            self%sub_grid(i - 1)%del_z(last) = ddz_interface
            self%sub_grid(i)%del_z(1) = ddz_interface
            !
        enddo
        !
    end subroutine setupMR
    !
    !> No subroutine briefing
    !
    subroutine setSubGrid( self, k, i1, i2 )
        implicit none
        !
        class( Grid3D_MR_t ), intent( inout ) :: self
        integer, intent( in ) :: k, i1, i2
        !
        integer :: cs
        real( kind=prec ), dimension(:), allocatable :: dx
        real( kind=prec ), dimension(:), allocatable :: dy
        real( kind=prec ), dimension(:), allocatable :: dz
        integer :: nx, ny, nz_air, nz_earth
        !
        cs = 2 ** self%coarseness( k, 1 )
        !
        call redimension_cells( cs, self%dx, dx )
        call redimension_cells( cs, self%dy, dy )
        dz = self%Dz(i1:i2)
        !
        self%coarseness( k, 3 ) = i1
        self%coarseness( k, 4 ) = i2
        !
        nx = size( dx )
        ny = size( dy )
        !
        nz_air = min( max( 0, self%nzAir - i1 + 1 ), self%coarseness( k, 2 ) )
        nz_earth = size(Dz) - nz_air
        !
        self%sub_grid(k) = Grid3D_SG_t( nx, ny, nz_air, nz_earth )
        !
        self%sub_grid(k)%dx = Dx
        self%sub_grid(k)%dy = Dy
        self%sub_grid(k)%dz = Dz
        !
        ! Origin correction for spherical grid
        if(i1 < nz_air) then
            ! Whole sub-grid in the air
            self%sub_grid(k)%oz = sum(self%dz(i2 + 1 : self%nzAir))
            !
        elseif(i1 > self%nzAir + 1) then
            self%sub_grid(k)%oz = -1.0*&
            sum(self%dz(self%nzAir + 1 : i1 - 1))
        else
            self%sub_grid(k)%oz = 0.0
        endif
        !
        call self%sub_grid(k)%setup
        !
        contains
            !
            subroutine Redimension_cells( cs_in, D_in, D_out )
                implicit none
                !
                integer, intent(in) :: cs_in
                real( kind=prec ), intent( in ) :: D_in(:)
                real( kind=prec ), allocatable :: D_out(:)
                !
                integer :: n, iout, k, count
                real( kind=prec ), allocatable :: tmp(:)
                !
                n = size(D_in)
                allocate(tmp(n))
                !
                do iout = 1, n
                    tmp(iout) = 0.0
                    do k = iout, iout -(cs_in - 1), -1
                        if(k < 1) exit
                        tmp(iout) = tmp(iout) + D_in(k)
                    enddo
                enddo
                !
                count = 0
                do k = cs, n, cs
                    count = count + 1
                enddo
                !
                allocate( D_out( count ) )
                !
                count = 1
                do k = cs, n, cs
                    D_out(count) = tmp(k)
                    count = count + 1
                enddo
                !
            end subroutine Redimension_cells
            !
    end subroutine setSubGrid
    !
    !> No subroutine briefing
    !
    subroutine numberOfEdges_Grid3D_MR( self, n_xedge, n_yedge, n_zedge )
        implicit none
        !
        class( Grid3D_MR_t ), intent( in ) :: self
        integer, intent( out ) :: n_xedge, n_yedge, n_zedge
        !
        integer :: k
        integer :: nx, ny, nz
        !
        n_xedge = 0
        n_yedge = 0
        n_zedge = 0
        !
        do k = 1, self%n_grids
            !
            call self%sub_grid(k)%numberOfEdges( nx, ny, nz )
            !
            n_xedge = n_xedge + nx
            n_yedge = n_yedge + ny
            n_zedge = n_zedge + nz
            !
        enddo
        !
    end subroutine numberOfEdges_Grid3D_MR
    !
    !> No subroutine briefing
    !
    subroutine numberOfFaces_Grid3D_MR( self, n_xface, n_yface, n_zface ) 
        implicit none
        !
        class( Grid3D_MR_t ), intent( in ) :: self
        integer, intent( out ) :: n_xface, n_yface, n_zface
        !
        integer :: k
        integer :: nx, ny, nz
        !
        n_xface = 0
        n_yface = 0
        n_zface = 0
        !
        do k = 1, self%n_grids
            !
            call self%sub_grid(k)%numberOfFaces(nx, ny, nz)
            !
            n_xface = n_xface + nx
            n_yface = n_yface + ny
            n_zface = n_zface + nz
            !
        enddo
        !
    end subroutine numberOfFaces_Grid3D_MR
    !
    !> No function briefing
    !
    function numberOfNodes_Grid3D_MR(self) result(n)
        implicit none
        !
        class( Grid3D_MR_t ), intent( in ) :: self
        !
        integer :: i, n
        !
        n = 0
        !
        do i = 1, self%n_grids
            n = n + self%sub_grid(i)%numberOfNodes()
        enddo
        !
    end function numberOfNodes_Grid3D_MR
    !
    !> No function briefing
    !
    function numberOfCells_Grid3D_MR( self ) result( n )
        implicit none
        !
        class( Grid3D_MR_t ), intent( in ) :: self
        !
        integer :: n
        !
        call errStop( "numberOfCells_Grid3D_MR not implemented" )
        !
    end function numberOfCells_Grid3D_MR
    !
    !> For sub-grid iGrid, returns logical array indicating if top/bottom boundaries
    !> are active--TopBottom(1) is for top, TopBottom(2) is for bottom
    !> exactly same function needed for all MR fields!
    !
    function active_Grid3D_MR( self, iGrid ) result( TopBottom ) 
        implicit none
        !
        class( Grid3D_MR_t ), intent( in ) :: self
        integer, intent( in ) :: iGrid
        !
        logical :: TopBottom(2)
        !
        !  perhaps some error checking
        !   e.g., 0<iGrid<nGrid ????
        if( iGrid .EQ. 1 ) then
            TopBottom(1) = .TRUE.
        else
            TopBottom(1) = self%cs( iGrid ) .LT. self%cs( iGrid - 1 )
        endif
        !
        if( iGrid .EQ. self%n_grids ) then
            TopBottom(2) = .TRUE.
        else
            TopBottom(2) = self%cs( iGrid ) .LT. self%cs( iGrid + 1 )
        endif
        !
    end function active_Grid3D_MR
    !
    !> Function to calculate total number of active elements for vectors/scalars of some type
    !> In each sub grid
    !
    function nActive_Grid3D_MR( self, grid_type ) result( n_active )
        implicit none
        !
        class( Grid3D_MR_t ), intent( in ) :: self
        character( len=4 ), intent( in ) :: grid_type
        !
        integer :: iGrid
        integer, dimension(:), allocatable :: n_active
        !
        allocate( n_active( self%n_grids ) )
        !
        do iGrid = 1, self%n_grids
            !
            select case( grid_type )
                !
                case( EDGE )
                    n_active(iGrid) = &
                    self%n_active_edge( 1, iGrid ) + &
                    self%n_active_edge( 2, iGrid ) + &
                    self%n_active_edge( 3, iGrid )
                case( FACE )
                    n_active(iGrid) = &
                    self%n_active_face( 1, iGrid ) + &
                    self%n_active_face( 2, iGrid ) + &
                    self%n_active_face( 3, iGrid )
                case( NODE )
                    n_active(iGrid) = self%n_active_node( iGrid )
                case( CELL )
                    n_active(iGrid) = self%n_active_cell( iGrid )
                case default
                    call errStop( "nActive_Grid3D_MR > Unknown grid_type" )
                !
            end select
            !
        enddo
        !
    end function nActive_Grid3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setActiveLimits_Grid3D_MR( self )
        implicit none
        !
        class( Grid3D_MR_t ), intent( inout ) :: self
        !
        integer :: iGrid
        logical :: TopBottom(2)
        integer :: nX, nY, nZ
        !  just the algorithm here
        !self%Nactive = 0
        !
        do iGrid = 1, self%n_grids
            !  this sets limits(in vertical) of active faces, edges, nodes
            !  First set total--active + inactive for sub-grid
            !
            call self%sub_grid( iGrid )%numberOfEdges( nX, nY, nZ )
            self%n_active_edge( 1, iGrid ) = nX
            self%n_active_edge( 2, iGrid ) = nY
            self%n_active_edge( 3, iGrid ) = nZ
            !
            call self%sub_grid( iGrid )%numberOfFaces( nX, nY, nZ )
            self%n_active_face( 1, iGrid ) = nX
            self%n_active_face( 2, iGrid ) = nY
            self%n_active_face( 3, iGrid ) = nZ
            !
            self%n_active_node( iGrid ) = self%sub_grid( iGrid )%numberOfNodes()
            !  not sure this function is defined -- value his nx*ny*nz
            self%n_active_cell( iGrid ) = self%sub_grid( iGrid )%numberOfCells()
            !
            TopBottom = self%active( iGrid )
            !
            ! IF WHAT ????
            if( TopBottom(1) ) then
                self%z_limits( 1, iGrid ) = 1
            else
                self%z_limits( 1, iGrid ) = 2
                call self%reduceActive( iGrid )
            endif
            !
            if( TopBottom(2) ) then
                self%z_limits( 2, iGrid ) = self%sub_grid( iGrid )%nz + 1
            else
                self%z_limits( 2, iGrid ) = self%sub_grid( iGrid )%nz
                call self%reduceActive( iGrid )
            endif
        enddo
        !
    end subroutine setActiveLimits_Grid3D_MR
    !
    !> just algorithm. -- reduce number of active edges/faces/nodes
    !> For one vertical layer in the sub_grid
    !
    subroutine reduceActive_Grid3D_MR( self, iGrid )
        implicit none
        !
        class( Grid3D_MR_t ), intent( inout ) :: self
        integer, intent( in ) :: iGrid
        !
        self%n_active_edge(1,iGrid) = self%n_active_edge(1,iGrid) -  &
        self%sub_grid(iGrid)%nx *(self%sub_grid(iGrid)%ny+1)
        !
        self%n_active_edge(2,iGrid) = self%n_active_edge(2,iGrid) -  &
        (self%sub_grid(iGrid)%nx+1) * self%sub_grid(iGrid)%ny
        !
        self%n_active_face(3,iGrid) = self%n_active_face(3,iGrid) -  &
        self%sub_grid(iGrid)%nx * self%sub_grid(iGrid)%ny
        !
        self%n_active_node(iGrid) = self%n_active_node(iGrid) -  &
        (self%sub_grid(iGrid)%ny+1) *(self%sub_grid(iGrid)%nx+1)
        !
    end subroutine reduceActive_Grid3D_MR
    !
    !> No function briefing
    !
    function slice1D_Grid3D_MR( self ) result( g1D )
        implicit none
        !
        class( Grid3D_MR_t ), intent( in ) :: self
        type( Grid1D_t ) :: g1D
        !
        g1D = Grid1D_t( self%nzAir, self%nzEarth, self%dz )
        !
    end function slice1D_Grid3D_MR
    !
    !> No function briefing
    !
    function slice2D_Grid3D_MR( self ) result( g2D )
        implicit none
        !
        class( Grid3D_MR_t ), intent( in ) :: self
        type( Grid2D_t ) :: g2D
        !
        call errStop( "slice2D_Grid3D_MR not implemented!" )
        !
    end function slice2D_Grid3D_MR
    !
    function repMat( m_in, nx, ny, nz, transp ) result( m_out )
        implicit none
        !
        real( kind=prec ), dimension(:), intent( in ) :: m_in
        integer, intent( in ) :: nx, ny, nz
        logical, intent( in ) :: transp
        !
        real( kind=prec ), dimension(:,:,:), allocatable :: m_out
        integer :: i, i1, i2, n_in
        !
        n_in = size( m_in )
        !
        if( transp ) then
            !
            allocate( m_out( nx, n_in * ny, nz ) )
            !
            ! Copy along 2nd dimension.
            i1 = 1; i2 = n_in
            do i = 1, ny
                m_out(1, i1:i2, 1) = m_in
                i1 = i2 + 1
                i2 = i2 + n_in
            enddo
            !
            ! Copy along 1st dimension.
            do i = 1, nx
                m_out(i, :, 1) = m_out(1, :, 1)
            enddo
            !
            ! Copy along 3rd dimension
            do i = 1, nz
                m_out(:, :, i) = m_out(:, :, 1)
            enddo
            !
        else
            !
            allocate(m_out(n_in*nx, ny, nz))
            !
            ! Copy along 1st dimension.
            i1 = 1; i2 = n_in
            do i = 1, nx
                m_out(i1:i2, 1, 1) = m_in
                i1 = i2 + 1
                i2 = i2 + n_in
            enddo
            !
            ! Copy along 2nd dimension.
            do i = 1, ny
                m_out(:, i, 1) = m_out(:, 1, 1)
            enddo
            !
            ! Copy along 3rd dimension
            do i = 1, nz
                m_out(:, :, i) = m_out(:, :, 1)
            enddo
            !
        endif
        !
    end function repMat
    !
    !> Procedure updateAirLayers_Grid
    !> Assumes that the grid is already defined, and merely
    !> includes the new air layers in the grid.
    !
    subroutine updateAirLayers_Grid3D_MR( self, nzAir, dzAir )
        implicit none
        !
        class( Grid3D_MR_t ), intent( inout ) :: self
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
        call self%sub_grid(1)%updateAirLayers( nzAir, dzAir )
        !
        call self%setup
        !
    end subroutine updateAirLayers_Grid3D_MR
    !
    !
    !
    subroutine write_Grid3D_MR( self )
        implicit none
        !
        class( Grid3D_MR_t ), intent( in ) :: self
        !
        integer :: i
        !
        write( *, * ) "GridMR:"
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
        write( *, * ) "    EDGEb, FACEb, NODEb: ", size( self%EDGEb ), size( self%FACEb ), size( self%NODEb )
        write( *, * ) "    EDGEi, FACEi, NODEi: ", size( self%EDGEi ), size( self%FACEi ), size( self%NODEi )
        write( *, * ) "    EDGEa, FACEa, NODEa: ", size( self%EDGEa ), size( self%FACEa ), size( self%NODEa )
        !
        write( *, * ) "    dx, dy, dz: ", size( self%dx ), size( self%dy ), size( self%dz )
        !
        do i = 1, self%n_grids
            call self%sub_grid(i)%write
        enddo
        !
    end subroutine write_Grid3D_MR
    !
end module Grid3D_MR
