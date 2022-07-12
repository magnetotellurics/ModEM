!**
! Implementation of standard 3D cartesian grid.
!
!*
module Grid3D_SG
    !
    use Constants
    use Grid
    use Grid1D
    use Grid2D
    use ModEMControlFile
    !
    type, extends( Grid_t ) :: Grid3D_SG_t
        !
        !
        !
        contains
            !
            final :: Grid3D_SG_dtor
            !**
            ! Overriden methods
            !*
            procedure, public :: NumberOfEdges => NumberOfEdgesGrid3D_SG
            procedure, public :: NumberOfFaces => NumberOfFacesGrid3D_SG
            procedure, public :: NumberOfNodes => NumberOfNodesGrid3D_SG
            procedure, public :: GridIndex => GridIndexGrid3D_SG
            procedure, public :: VectorIndex => VectorIndexGrid3D_SG
            procedure, public :: Limits => LimitsGrid3D_SG
            procedure, public :: IsAllocated => IsAllocatedGrid3D_SG
            procedure, public :: Length => LengthGrid3D_SG
            !**
            !
            procedure, public :: Create => CreateGrid3D_SG
            procedure, public :: allocateDim => allocateDimGrid3D_SG
            procedure, public :: Setup => SetupGrid3D_SG
            procedure, public :: SetupAirLayers => SetupAirLayersGrid3D_SG
            procedure, public :: UpdateAirLayers => UpdateAirLayersGrid3D_SG
            ! 
            procedure, public :: GetDimensions => GetDimensionsGrid3D_SG
            procedure, public :: SetCellSizes => SetCellSizesGrid3D_SG
            procedure, public :: GetCellSizes => GetCellSizesGrid3D_SG

            procedure, public :: Copy_from => Copy_fromGrid3D_SG

            procedure, public :: Slice1D => Slice1DGrid3D_SG
            procedure, public :: Slice2D => Slice2DGrid3D_SG
        !
    end type Grid3D_SG_t

    !**
    ! Details needed to unambiguosly compute and/or store the air layers;
    ! method options are: mirror; fixed height; read from file
    ! for backwards compatibility, all of the defaults are set to what
    ! was previously hard coded (AK; May 19, 2017)
    ! For backwards compatibility, default is "mirror 10 3. 30."
    !    GDE 12/17/21 : new method "undefined" is default -- force
    !        explicit setting of air layers (can still hard code a default
    !         in the driver program!)
    !*
    type :: TAirLayers
         character(:), allocatable :: method
         integer :: nz
         real( kind=prec ) :: maxHeight, minTopDz, alpha
         real( kind=prec ), allocatable, dimension(:) :: dz
         logical :: is_allocated
    end type TAirLayers
    !
    interface Grid3D_SG_t
         module procedure Grid3D_SG_t_ctor
    end interface Grid3D_SG_t
    !
contains
    function Slice1DGrid3D_SG(self) result( g1D )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        type( Grid1D_t ) :: g1D
        !
        g1D = Grid1D_t( self%nzAir, self%nzEarth, self%dz )
        !
    end function Slice1DGrid3D_SG
    !
    function Slice2DGrid3D_SG(self) result( g2D )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        type( Grid2D_t ) :: g2D
        !
        ! Should be different for the polarization
        g2D = Grid2D_t( self%ny, self%nzAir, self%nzEarth, self%dy, self%dz )
        !
    end function Slice2DGrid3D_SG
    !**
    ! Class constructor for simple tensor product grid
    ! Usage obj = Grid_t3D(Dx,Dy,Dz,Nza)
    ! Dx, Dy, Dz are cell dimensions for x, y, z direction
    ! Nza is number of air layers to allow (included in Dz)
    !*
    function Grid3D_SG_t_ctor( nx, ny, nzAir, nzEarth, dx, dy, dz ) result(self)
        implicit none
        !
        integer, intent( in ) :: nx, ny, nzAir, nzEarth
        real( kind=prec ), dimension(:), intent( in ) :: dx, dy, dz
        !
        type( Grid3D_SG_t ) :: self
        !
        !write(*,*) "Constructor Grid3D_SG_t"
        !
        call self%init()
        !
        call self%Create( nx, ny, nzAir, nzEarth )
        !
        call self%SetCellSizes( dx, dy, dz )
        !
        call self%Setup()
        !
    end function Grid3D_SG_t_ctor
    !
    ! ModelOperator_MF destructor
    subroutine Grid3D_SG_dtor( self )
        implicit none
        !
        type( Grid3D_SG_t ), intent( inout ) :: self
        !
        !write(*,*) "Destructor Grid3D_SG_t", self%nx, self%ny, self%nz
        !
        call self%dealloc()
        !
    end subroutine Grid3D_SG_dtor
    !
    subroutine CreateGrid3D_SG( self, nx, ny, nzAir, nzEarth )
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
        call self%allocateDim()

    end subroutine CreateGrid3D_SG
    
    subroutine allocateDimGrid3D_SG(self)
        implicit none
        !
        class( Grid3D_SG_t ), intent(inout) :: self
        !
        integer :: nx, ny, nz
        !
        if( self%is_allocated ) call self%dealloc()
        !
        nx = self%nx
        ny = self%ny
        nz = self%nz
        !
        allocate( real( kind=prec ) :: self%dx(nx) )
        allocate( real( kind=prec ) :: self%dy(ny) )
        allocate( real( kind=prec ) :: self%dz(nz) )
        !
        ! dxinv    = 1/ dx and similarly for dyinv and dzinv
        allocate( real( kind=prec ) :: self%dxInv(nx) )
        allocate( real( kind=prec ) :: self%dyInv(ny) )
        allocate( real( kind=prec ) :: self%dzInv(nz) )
        !
        ! delX, delY, and delZ are the distances between
        ! the electrical field defined on the center of the
        ! edges in x, y, and z axes, respectively.
        allocate( real( kind=prec ) :: self%delX(nx + 1) )
        allocate( real( kind=prec ) :: self%delY(ny + 1) )
        allocate( real( kind=prec ) :: self%delZ(nz + 1) )
        !
        allocate( real( kind=prec ) :: self%delXInv(nx + 1) )
        allocate( real( kind=prec ) :: self%delYInv(ny + 1) )
        allocate( real( kind=prec ) :: self%delZInv(nz + 1) )
        !
        ! xEdge is the array for cumulative distance of the edge
        ! for each grid (starting from the coordinate axes) with
        ! dimensions nx + 1.
        ! xCenter is the array for cumulative distance of the center
        ! for each grid (starting from the coordinate axes) with
        ! dimension n.
        ! yEdge, yCenter, zEdge, zCenter are analagous arrays for
        ! other directions.
        allocate( real( kind=prec ) :: self%xEdge(nx + 1) )
        allocate( real( kind=prec ) :: self%yEdge(ny + 1) )
        allocate( real( kind=prec ) :: self%zEdge(nz + 1) )
        allocate( real( kind=prec ) :: self%xCenter(nx) )
        allocate( real( kind=prec ) :: self%yCenter(ny) )
        allocate( real( kind=prec ) :: self%zCenter(nz) )
        !
        self%is_allocated = .TRUE.
        !
    end subroutine allocateDimGrid3D_SG

    subroutine DeallocateDimGrid3D_SG(self)
        implicit none
        !
        class( Grid3D_SG_t ), intent(inout) :: self
        !
        if( .NOT. self%is_allocated ) return

        deallocate(self%dx)
        deallocate(self%dy)
        deallocate(self%dz)
        
        deallocate(self%dxInv)
        deallocate(self%dyInv)
        deallocate(self%dzInv)
        
        deallocate(self%delX)
        deallocate(self%delY)
        deallocate(self%delZ)
        
        deallocate(self%delXInv)
        deallocate(self%delYInv)
        deallocate(self%delZInv)
        
        deallocate(self%xEdge)
        deallocate(self%yEdge)
        deallocate(self%zEdge)
        deallocate(self%xCenter)
        deallocate(self%yCenter)
        deallocate(self%zCenter)
        
        self%is_allocated = .FALSE.
        
    end subroutine DeallocateDimGrid3D_SG
    !**
    ! Setup does calculations for grid geometry, which cannot be done
    ! until dx, dy, dz, and the origin are set.
    !
    !*
    subroutine SetupGrid3D_SG(self, origin)
        implicit none
        !
        class( Grid3D_SG_t ), intent( inout ) :: self
        real( kind=prec ) , intent( in ), optional :: origin(3)
        !
        integer :: ix, iy, iz, i, j, nzAir
        real( kind=prec ) :: xCum, yCum, zCum
        real( kind=prec ) :: ox, oy, oz
        !
        self%dxInv = 1 / self%dx
        self%dyInv = 1 / self%dy
        self%dzInv = 1 / self%dz
        !
        call self%GetOrigin( ox, oy, oz )
        !
        if(present(origin) ) then
            !
            ox = origin(1)
            oy = origin(2)
            oz = origin(3)
            !
            call self%SetOrigin(ox, oy, oz)
            !
        end if
        !
        self%xEdge(1) = ox
        self%yEdge(1) = oy
        self%zEdge(1) = oz
        
        xCum = 0.0
        yCum = 0.0
        zCum = 0.0
        do ix = 1, self%nx
            xCum = xCum + self%dx(ix)
            self%xEdge(ix+1) = xCum + ox
        end do
        do iy = 1, self%ny
            yCum = yCum + self%dy(iy)
            self%yEdge(iy + 1) = yCum + oy
        end do
        
        ! NOTE: adjust for origin later to get airthickness, 
        ! reference to origin at Earth"s surface correct!
        do iz = 1, self%nz
            zCum = zCum + self%dz(iz)
            self%zEdge(iz + 1) = zCum
        end do

        nzAir = self%nzAir
        self%zAirThick = self%zEdge(nzAir + 1)

        ! Distance between center of the selfs
        self%delX(1) = self%dx(1)
        do ix = 2, self%nx
            self%delX(ix) = self%dx(ix - 1) + self%dx(ix)
        end do
        self%delX(self%nx + 1) = self%dx(self%nx)
        self%delX = self%delX/2.0
        
        self%delY(1) = self%dy(1)
        do iy = 2, self%ny
            self%delY(iy) = self%dy(iy - 1) + self%dy(iy)
        end do
        self%delY(self%ny + 1) = self%dy(self%ny)
        self%delY = self%delY/2.0
        
        self%delZ(1) = self%dz(1)
        do iz = 2, self%nz
            self%delZ(iz) = self%dz(iz - 1) + self%dz(iz)
        end do
        self%delZ(self%nz + 1) = self%dz(self%nz)
        self%delZ = self%delZ/2.0
        
        self%delXInv = 1/self%delX
        self%delYInv = 1/self%delY
        self%delZInv = 1/self%delZ

        ! Cumulative distance between the centers, adjusted to model origin
        xCum = 0.0
        yCum = 0.0
        zCum = 0.0
        do ix = 1, self%nx
            xCum = xCum + self%delX(ix)
            self%xCenter(ix) = xCum + ox
        end do
        do iy = 1, self%ny
            yCum = yCum + self%delY(iy)
            self%yCenter(iy) = yCum + oy
        end do
        do iz = 1, self%nz
            zCum = zCum + self%delZ(iz)
            self%zCenter(iz) = zCum
        end do
        ! Need to be careful here ... grid origin is given
        ! at Earth"s surface, not top of model domain!
        do iz = 1, self%nz
            self%zCenter(iz) = self%zCenter(iz) - self%zAirThick + oz
            self%zEdge(iz) = self%zEdge(iz) - self%zAirThick + oz
        end do
        !
        self%zEdge(self%nz + 1) = self%zEdge(self%nz + 1) - self%zAirThick + oz
        !
    end subroutine SetupGrid3D_SG

    !**
    ! SetupAirLayers computes the Dz in the airlayers structure
    ! using the grid to get the top layers Dz;
    ! all values expected in km on input
    ! For backwards compatibility, default is "mirror 10 3. 30."
    ! but the use of "fixed height 12 1000" is recommended
    !*
    subroutine SetupAirLayersGrid3D_SG( self, airLayers, method, &
                                        nzAir, maxHeight, minTopDz, &
                                        alpha, dzAir )
        implicit none
        !
        class( Grid3D_SG_t ), intent( inout ) :: self
        type( TAirLayers ), intent( inout )   :: airlayers
        character(:), allocatable, intent( in ), optional :: method
        integer, intent( in ), optional :: nzAir
        real( kind=prec ), intent( in ), optional :: maxHeight, minTopDz, alpha
        real( kind=prec ), intent( in ), optional, pointer :: dzAir
        !
        integer :: ix, iy, iz, i, j
        integer :: status
        real( kind=prec ) :: z1Log, dlogz, zLog, height1, height2
        !
        airLayers%is_allocated = .FALSE.
        !
        if( present( method ) ) then
            airlayers%method = method
        else
            airlayers%method = "fixed height"
        end if
        !
        if( present( nzAir ) ) then
            airlayers%nz = nzAir
        else
            airlayers%nz = model_n_air_layer
        end if
        !
        if( .NOT. ( index( airLayers%method, "read from file" ) > 0 ) ) then
            if( airLayers%is_allocated ) then
                deallocate( airlayers%dz )
            end if
            allocate( airLayers%dz( airLayers%nz ) )
            airLayers%is_allocated = .TRUE.
        end if
        !
        if( present( maxHeight ) ) then
            airLayers%maxHeight = 1000. * maxHeight
        else
            airLayers%maxHeight = model_max_height
        end if
        
        if( present( minTopDz ) ) then
            airLayers%minTopDz = 1000. * minTopDz
        else
            airLayers%minTopDz = 100.0
        end if
        
        if( present( alpha ) ) then
            airLayers%alpha = alpha
        else
            airLayers%alpha = 3.
        end if
        
        if(index(airLayers%method, "mirror") > 0) then
			!**
			! Following is Kush"s approach to setting air layers:
			! mirror imaging the dz values in the air layer with respect to
			! earth layer as far as we can using the following formulation
			! air layer(bottom:top) = (alpha)^(j-1) * earth layer(top:bottom)
			!*
			do iz = airLayers%nz, 1, -1
				j = airLayers%nz - iz + 1
				airLayers%dz(iz) = ((airLayers%alpha)**(j - 1) ) * self%dz(self%nzAir + j)
			end do
			!
			! The topmost air layer has to be at least 30 km
			if(airLayers%dz(1).lt.airLayers%minTopDz) then
				airLayers%dz(1) = airLayers%minTopDz
			end if

        else if(index(airLayers%method, "fixed height") > 0) then 
			!
			z1Log = log10(self%dz(self%nzAir + 1) )
			dlogz = (log10(airLayers%maxHeight) - z1Log)/(airLayers%nz-1)
			zLog = z1Log
			height1 = 10.**z1log
			airLayers%dz(airLayers%Nz) = height1
			do iz = airLayers%Nz-1, 1, -1
				zLog = zLog + dlogz
				height2 = 10.**zLog 
				airLayers%dz(iz) = height2-height1
				height1 = height2
			end do
			!
        else if(index(airLayers%method, "read from file") > 0) then
            !**
            ! Air layers have been read from file and are
            ! already stored in Dz, so only need to reallocate
            ! if passing a new array to it.
            !*
            if( present(dzAir) ) then
                if( airLayers%is_allocated) then
                     deallocate( airLayers%dz )
                end if
                allocate( airLayers%dz(airLayers%nz) )
                airLayers%dz = dzAir
            end if
        end if
        !
    end subroutine SetupAirLayersGrid3D_SG
    !**
    !
    !*
    subroutine DeallocateAirLayersGrid3D_SG(self, airLayers)
        implicit none
        !
        class( Grid3D_SG_t ), intent( inout ) :: self
        type( TAirLayers ), intent( inout ) :: airLayers
        !
        deallocate( airLayers%dz )
        !
    end subroutine DeallocateAirLayersGrid3D_SG
    !**
    ! Assumes that the grid is already defined, and merely
    ! includes the new air layers in the grid.
    !*
    subroutine UpdateAirLayersGrid3D_SG(self, nzAir, dzAir)
        implicit none
        !
        class( Grid3D_SG_t ), intent(inout) :: self
        integer                     , intent( in )        :: nzAir
        real( kind=prec ) , intent( in )        :: dzAir(:)
        ! Local variables
        integer :: nzAir_old, nzEarth_old
        integer :: nx_old, ny_old, nz_old
        real( kind=prec ), allocatable :: dx_old(:), dy_old(:), dz_old(:)
        real( kind=prec ) :: ox_old, oy_old, oz_old
        real( kind=prec ) :: rotDeg_old
        character(len = 80) :: geometry_old
        
        if(.NOT.self%is_allocated) then
             write( *, * ) "ERROR:Grid3D_SG_t:UpdateAirLayers"
             stop "    Grid not allocated."
        end if

        nx_old = self%nx
        ny_old = self%ny
        nz_old = self%nz
        nzAir_old = self%nzair
        nzEarth_old = self%nzEarth
        dx_old = self%dx
        dy_old = self%dy
        dz_old = self%dz
        geometry_old = self%GetGridGeometry()
    
        ox_old = self%ox
        oy_old = self%oy
        oz_old = self%oz
    
        rotdeg_old = self%rotdeg
        
        call self%dealloc()
        call self%Create(nx_old, ny_old, nzAir, nzEarth_old)
        
        ! Set air layers to dzAir values and copy the rest
        self%dz(1:nzAir) = dzAir
        self%dz(nzAir+1:self%nz) = dz_old(nzAir_old+1:nz_old)
        self%dy = dy_old
        self%dx = dx_old


        call self%SetOrigin(ox_old, oy_old, oz_old)
        call self%SetGridRotation(rotdeg_old)
        call self%SetGridGeometry(geometry_old)

        ! Setup the rest of the grid from scratch
        call self%Setup()
        
    end subroutine UpdateAirLayersGrid3D_SG

    subroutine GetDimensionsGrid3D_SG(self, nx, ny, nz, nzAir)
        implicit none
        !
        class( Grid3D_SG_t ), intent( in )  :: self
        integer, intent( out )             :: nx, ny, nz, nzAir
        !
        nx = self%nx
        ny = self%ny
        nz = self%nz
        nzAir = self%nzAir
        !
    end subroutine GetDimensionsGrid3D_SG

    subroutine SetCellSizesGrid3D_SG(self, dx, dy, dz)
        implicit none
        !
        class( Grid3D_SG_t ), intent(inout) :: self
        real( kind=prec ) , dimension(:), intent( in ) :: dx, dy, dz

        if(.NOT.self%IsAllocated() ) then
             write( *, * ) "ERROR:Grid3D_SG_t:SetCellSizes:"
             stop "    Grid not allocated."
        end if

        ! Check dimensions
        if((size(dx).ne.size(self%dx) ).or.&
                 (size(dy).ne.size(self%dy) ).or.&
                 (size(dz).ne.size(self%dz) )) then
             write( *, * ) "ERROR:Grid3D_SG_t:SetCellSizes:"
             stop "    Incompatible sizes for cell arrays."
        end if

        self%dx = dx
        self%dy = dy
        self%dz = dz

    end subroutine SetCellSizesGrid3D_SG
    
    subroutine GetCellSizesGrid3D_SG(self, dx, dy, dz)
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        real( kind=prec ) , intent( out )  :: dx(:), dy(:), dz(:)
        !
        if(.NOT.self%IsAllocated() ) then
             write( *, * ) "ERROR:Grid3D_SG_t:GetCellSizes:"
             stop "    Grid not allocated."
        end if
        !
        ! Check dimensions
        if((size(dx).ne.size(self%dx) ).or.&
                 (size(dy).ne.size(self%dy) ).or.&
                 (size(dz).ne.size(self%dz) )) then
             write( *, * ) "ERROR:Grid3D_SG_t:GetCellSizes:"
             stop "    Incompatible sizes for cell arrays."
        end if
        !
        dx = self%dx
        dy = self%dy
        dz = self%dz
        !
    end subroutine GetCellSizesGrid3D_SG
    !**
    ! NumberOfEdges
    !
    !*
    subroutine NumberOfEdgesGrid3D_SG(self, nXedge, nYedge, nZedge)
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        integer, intent( out )           :: nXedge, nYedge, nZedge
        ! Local variables
        integer :: nx, ny, nz
        
        call self%Limits( "XEDGE", nx, ny, nz )
        nXedge = nx*ny*nz
        
        call self%Limits( "YEDGE", nx, ny, nz )
        nYedge = nx*ny*nz

        call self%Limits( "ZEDGE", nx, ny, nz )
        nZedge = nx*ny*nz
        
    end subroutine NumberOfEdgesGrid3D_SG
    !**
    ! NumberOfFaces
    !
    !*
    subroutine NumberOfFacesGrid3D_SG( self, nXface, nYface, nZface ) 
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        integer, intent( out )             :: nXface, nYface, nZface
        !
        integer :: nx, ny, nz
        
        call self%Limits("XFACE", nx, ny, nz)
        nXface = nx*ny*nz
        
        call self%Limits("YFACE", nx, ny, nz)
        nYface = nx*ny*nz
        
        call self%Limits("ZFACE", nx, ny, nz)
        nZface = nx*ny*nz
        
    end subroutine NumberOfFacesGrid3D_SG
    
    !**
    ! NumberOfNodes
    !
    !*
    function NumberOfNodesGrid3D_SG(self) result(n)
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        ! Local variables
        integer :: n
        
        n = (self%nx + 1)*(self%ny + 1)*(self%nz + 1)

    end function NumberOfNodesGrid3D_SG
    
    !**
    ! GridIndex
    !
    ! Based on matlab method of same name in class Grid_t3D
    ! IndVec is the index within the list of nodes of a fixed type
    ! e.g., among the list of y-Faces.     An offset needs to be
    ! added to get index in list of all faces (for example).
    !*
    subroutine GridIndexGrid3D_SG( self, nodeType, indVec, i, j, k )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        character(*), intent( in )       :: nodeType
        integer, intent( in )            :: indVec(:)
        integer, intent( out )           :: i(:), j(:), k(:)
        !
        integer :: nx, ny, nz, nVec, ii
        real(4) :: rNxy, rNx

        call self%Limits(nodeType, nx, ny, nz)
        nVec = size(indVec)
        
        if(nVec.ne.size(i) ) then
             stop "Size of 'ind_vec' and 'i' do not agree."
        end if
        
        if(nVec.ne.size(j) ) then
             stop "Size of 'ind_vec' and 'j' do not agree."
        end if
        
        if(nVec.ne.size(k) ) then
             stop "Size of 'ind_vec' and 'k' do not agree."
        end if
        
        rNxy = float(nx*ny)
        rNx    = float(nx)
        
        do ii = 1, nVec
            i(ii) = mod(indVec(ii), nx)
            j(ii) = mod(ceiling(float(indVec(ii) )/rNx), ny)
            k(ii) = ceiling(float(indVec(ii) )/rNxy)
        end do
        
        where(i.eq.0) i = nx
        where(j.eq.0) j = ny
        where(k.eq.0) k = nz

    end subroutine GridIndexGrid3D_SG

    !**
    ! VectorIndex
    !
    ! Based on matlab method of same name in class Grid_t3D
    ! returned array IndVec gives numbering of nodes within
    ! the list for nodeType; need to add an offset for position
    ! in full list of all faces or edges (not nodes and cells).
    !*
    subroutine VectorIndexGrid3D_SG( self, nodeType, i, j, k, indVec )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        character(*), intent( in )       :: nodeType
        integer, intent( in )            :: i(:), j(:), k(:)
        integer, intent( out )           :: indVec(:)
        ! Local variables
        integer :: nx, ny, nz, nxy, nVec, ii
        
        call self%Limits(nodeType, nx, ny, nz)
        
        nVec = size(indVec)
        
        if(nVec.ne.size (i) ) then
             stop "Size of 'ind_vec' and 'i' do not agree."
        end if
        
        if(nVec.ne.size (J) ) then
             stop "Size of 'ind_vec' and 'j' do not agree."
        end if
        
        if(nVec.ne.size (K) ) then
             stop "Size of 'ind_cec' and 'k' do not agree."
        end if
        
        nxy = nx*ny
        do ii = 1, nVec
             indVec(ii) = (K(ii) - 1) * nxy + (j(ii) - 1) * nx + i(ii)
        end do
        
    end subroutine VectorIndexGrid3D_SG
    !**
    ! Limits
    !
    !*
    subroutine LimitsGrid3D_SG( self, nodeType, nx, ny, nz )
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        character(*), intent( in )       :: nodeType
        integer, intent( out )           :: nx, ny, nz
        !
        select case(nodeType)
        case(CENTER)
             nx = self%nx
             ny = self%ny
             nz = self%nz

        case(CORNER)
             nx = self%nx + 1
             ny = self%ny + 1
             nz = self%nz + 1

        case(XEDGE)
             nx = self%nx
             ny = self%ny + 1
             nz = self%nz + 1

        case(XFACE)
             nx = self%nx + 1
             ny = self%ny
             nz = self%nz

        case(YEDGE)
             nx = self%nx + 1
             ny = self%ny
             nz = self%nz + 1

        case(YFACE)
             nx = self%nx
             ny = self%ny + 1
             nz = self%nz
             
        case(ZEDGE)
             nx = self%nx + 1
             ny = self%ny + 1
             nz = self%nz
             
        case(ZFACE)
             nx = self%nx
             ny = self%ny
             nz = self%nz + 1
             
        end select
        
    end subroutine LimitsGrid3D_SG

    function IsAllocatedGrid3D_SG(self) result(f)
        implicit none
        !
        class( Grid3D_SG_t ), intent( in ) :: self
        !
        logical :: f
        !
        f = self%is_allocated
    end function IsAllocatedGrid3D_SG

    function LengthGrid3D_SG(self) result(n)
        class( Grid3D_SG_t ), intent( in ) :: self
        integer :: n
        !
        n = self%nx * self%ny * self%nz
        !
    end function LengthGrid3D_SG

    subroutine Copy_fromGrid3D_SG(self, g)
        implicit none
        !
        class( Grid3D_SG_t ), intent(inout) :: self
        class(Grid_t)         , intent( in )        :: g
    end subroutine Copy_fromGrid3D_SG
    
end module Grid3D_SG
