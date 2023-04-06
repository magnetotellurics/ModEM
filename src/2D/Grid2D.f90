!
!> Implementation of standard 2D cartesian grid.
!>        there are too many things that are hard-coded for 3D in the
!>        abstract Grid class to make this an extension of that
!>        Might make an extension of Grid3D, but I see some issues there also
!
!
module Grid2D
    !
    use Constants
    use Grid1D
    !
    type :: Grid2D_t
         !
         !> Grid Dimensions:
         !>    x is along strike; y and z are coordinates of 2D model
         !> ny is grid dimension (number of cells) in the y-direction
         !> nzEarth is number of earth layers used in the grid modeling
         !> nz is grid dimension (number of cells) in the z-direction:
         !> nz = nzAir + nzEarth
         !
         integer :: ny, nz    !> Number of grid cells in y,z
                                                        !> directions

         integer :: nzAir             !> Number of air layers
         integer :: nzEarth         !> Number of earth layers

         !
         !> Grid geometry:
         !> Dy, DDy are arrays of grid spacings in x-direction
         !> Dy denotes spacings betgween cell edges: dimension: Dy(ny)
         !> DDy denotes spacings between cell centers: dimension: DDy(ny+1)
         !> why dimensions: DDy(ny+1) (DDy(2) is distance between centers of cells
         !> 2 and 1)
         !
         real( kind=prec ), allocatable, dimension(:) :: dy, dz
         real( kind=prec ), allocatable, dimension(:) :: dyInv, dzInv

         real( kind=prec ), allocatable, dimension(:) :: delY, delZ
         real( kind=prec ), allocatable, dimension(:) :: delYInv, delZInv
         !
         !> Book-keeping on cumulative distances
         !
         real( kind=prec ), allocatable, dimension(:) :: yEdge
         real( kind=prec ), allocatable, dimension(:) :: zEdge
         real( kind=prec ), allocatable, dimension(:) :: yCenter
         real( kind=prec ), allocatable, dimension(:) :: zCenter
     
         !
         !> Total thickness of the air above
         !
         real( kind=prec ) :: zAirThick = R_ZERO
 
         logical :: is_allocated = .FALSE. 

     contains
         !>     are all of these needed?
         procedure, public :: NumberOfEdges
         procedure, public :: NumberOfNodes
         procedure, public :: GridIndex
         procedure, public :: VectorIndex
         procedure, public :: Limits
         procedure, public :: Length
         procedure, public :: Slice1D
         
         !
         !
         procedure, public :: create => createGrid2D
         procedure, public :: allocate => allocateGrid2D
         procedure, public :: deallocate => deallocateGrid2D
         procedure, public :: setup =>  setupGrid2D
        
         !> 
         procedure, public :: getDimensions
         procedure, public :: SetCellSizes
         procedure, public :: GetCellSizes

    end type Grid2D_t

    interface Grid2D_t
         module procedure Grid2D_t_ctor
    end interface Grid2D_t
    
contains

    !
    !> Class constructor for simple 2D tensor product grid
    !> Dy, Dz are cell dimensions for y, z direction
    !> Nza is number of air layers to allow (included in Dz)
    !
    !> No subroutine briefing
    !
    function Grid2D_t_ctor(ny, nzAir, nzEarth, dy, dz) result(grid)
        !
        integer, intent( in ) :: ny, nzAir, nzEarth
        real( kind=prec ), dimension(:), intent( in ) :: dy, dz
        !
        type( Grid2D_t ) :: grid
        call grid%create( ny, nzAir, nzEarth )
        call grid%SetCellSizes( dy, dz )
        !
        call grid%setup
        
    end function Grid2D_t_ctor
    !
    !> No subroutine briefing
    subroutine createGrid2D(self, ny, nzAir, nzEarth)
        !
        class(Grid2D_t), intent(inout) :: self
        integer, intent( in ) :: ny, nzAir, nzEarth
        !
        integer :: nz
        
        nz = nzEarth + nzAir
        
        self%nzAir = nzAir
        self%nzEarth = nzEarth
        
        self%ny = ny        
        self%nz = nz
        
        call self%allocate

    end subroutine createGrid2D
    !
    !> No subroutine briefing
    subroutine allocateGrid2D( self )
        !
        class(Grid2D_t), intent(inout) :: self
        !
        integer :: ny, nz
        if( self%is_allocated ) call self%deallocate()

        ny = self%ny; nz = self%nz
        
        allocate(self%dy(ny))
        allocate(self%dz(nz))
        
        allocate(self%dyInv(ny))
        allocate(self%dzInv(nz))
        
        !> delY, and delZ are the distances between
        !> the electrical field defined on the center of the
        !> edges in y, and z axes, respectively.
        allocate(self%delY(ny + 1))
        allocate(self%delZ(nz + 1))
        
        allocate(self%delYInv(ny + 1))
        allocate(self%delZInv(nz + 1))
        
        allocate(self%yEdge(ny + 1))
        allocate(self%zEdge(nz + 1))
        allocate(self%yCenter(ny))
        allocate(self%zCenter(nz))
        
        self%is_allocated = .TRUE.
        
    end subroutine allocateGrid2D
    !
    !> No subroutine briefing
    subroutine deallocateGrid2D( self )
        !
        class(Grid2D_t), intent(inout) :: self
        
        if( .NOT.self%is_allocated ) return

        deallocate(self%dy)
        deallocate(self%dz)
        
        deallocate(self%dyInv)
        deallocate(self%dzInv)
        
        deallocate(self%delY)
        deallocate(self%delZ)
        
        deallocate(self%delYInv)
        deallocate(self%delZInv)
        
        deallocate(self%yEdge)
        deallocate(self%zEdge)
        deallocate(self%yCenter)
        deallocate(self%zCenter)
        
        self%is_allocated = .FALSE.
        
    end subroutine deallocateGrid2D
    !
    !> Procedure setupGrid2D
    !> setupGrid2D does calculations for grid geometry, which cannot be done
    !> until dy, dz, are set -- for 2D I am omitting origin -- can't see the use!
    subroutine setupGrid2D( self )
        implicit none
        !
        class(Grid2D_t), intent(inout) :: self
        !
        integer :: iy, iz, i, j
        integer :: nzAir
        integer :: aStatus
        real( kind=prec ) :: yCum, zCum
        
        self%dyInv = 1/self%dy
        self%dzInv = 1/self%dz

        yCum = R_ZERO
        zCum = R_ZERO
        do iy = 1, self%ny
             yCum = yCum + self%dy(iy)
             self%yEdge(iy + 1) = yCum
        enddo
        
        do iz = 1, self%nz
             zCum = zCum + self%dz(iz)
             self%zEdge(iz + 1) = zCum
        enddo

        nzAir = self%nzAir
        self%zAirThick = self%zEdge(nzAir + 1)

        !> Distance between center of the selfs
        self%delY(1) = self%dy(1)
        do iy = 2, self%ny
             self%delY(iy) = self%dy(iy - 1) + self%dy(iy)
        enddo
        self%delY(self%ny + 1) = self%dy(self%ny)
        self%delY = self%delY/2.0
        
        self%delZ(1) = self%dz(1)
        do iz = 2, self%nz
             self%delZ(iz) = self%dz(iz - 1) + self%dz(iz)
        enddo
        self%delZ(self%nz + 1) = self%dz(self%nz)
        self%delZ = self%delZ/2.0
        
        self%delYInv = 1/self%delY
        self%delZInv = 1/self%delZ

        !> Cumulative distance between the centers
        yCum = R_ZERO
        zCum = R_ZERO
        do iy = 1, self%ny
             yCum = yCum + self%delY(iy)
             self%yCenter(iy) = yCum
        enddo
        do iz = 1, self%nz
             zCum = zCum + self%delZ(iz)
             self%zCenter(iz) = zCum
        enddo
        
        do iz = 1, self%nz
             self%zCenter(iz) = self%zCenter(iz) - self%zAirThick
             self%zEdge(iz) = self%zEdge(iz) - self%zAirThick
        enddo
        self%zEdge(self%nz + 1) = self%zEdge(self%nz + 1) - &
                 self%zAirThick
        
    end subroutine setupGrid2D
    !
    !> No subroutine briefing
    subroutine getDimensions(self, ny, nz, nzAir)
        !
        class(Grid2D_t), intent( in ) :: self
        integer, intent( out ) :: ny, nz, nzAir

        ny = self%ny
        nz = self%nz
        nzAir = self%nzAir
    end subroutine getDimensions
    !
    !> No subroutine briefing
    subroutine SetCellSizes(self, dy, dz)
        !
        class(Grid2D_t), intent(inout) :: self
        real( kind=prec ), dimension(:), intent( in ) :: dy, dz

        if(.NOT.self%is_allocated) then
             write(*, *) 'ERROR:Grid2D_t:SetCellSizes:'
             write(*, *) '    Grid not allocated.'

             stop
        endif

        !> Check dimensions
        if((size(dy).NE.size(self%dy)) .OR. &
                 (size(dz).NE.size(self%dz))) then
             write(*, *) 'ERROR:Grid2D_t:SetCellSizes:'
             write(*, *) '    Incompatible sizes for cell arrays.'

             stop
        endif

        self%dy = dy
        self%dz = dz

    end subroutine SetCellSizes
    !
    !> No subroutine briefing
    subroutine GetCellSizes(self, dy, dz)
        !
        class(Grid2D_t), intent( in ) :: self
        real( kind=prec ), intent( out ) :: dy(:), dz(:)

        if(.NOT.self%is_allocated) then
             write(*, *) 'ERROR:Grid2D_t:SetCellSizes:'
             write(*, *) '    Grid not allocated.'

             stop
        endif

        !> Check dimensions
        if( (size(dy).NE.size(self%dy)) .OR. &
                 (size(dz).NE.size(self%dz))) then
             write(*, *) 'ERROR:Grid2D_t:SetCellSizes:'
             write(*, *) '    Incompatible sizes for cell arrays.'
             
             stop
        endif
        
        dy = self%dy
        dz = self%dz
        
    end subroutine GetCellSizes
    !
    !> No subroutine briefing
    subroutine NumberOfEdges(self, nYedge, nZedge)
        !
        class(Grid2D_t), intent( in ) :: self
        integer, intent( out ) :: nYedge, nZedge
        !
        integer :: ny, nz
        
        call self%Limits('YEDGE', ny, nz)
        nYedge = ny*nz

        call self%Limits('ZEDGE', ny, nz)
        nZedge = ny*nz
        
    end subroutine NumberOfEdges
    !
    !> No subroutine briefing
    !
    function NumberOfNodes( self ) result(n)
        !
        class(Grid2D_t), intent( in ) :: self
        !
        integer :: n
        
        n = (self%ny + 1)*(self%nz + 1)

    end function NumberOfNodes
    !
    !> Procedure GridIndex
    !>
    !> Based on matlab method of same name in class Grid_t3D
    !> IndVec is the index within the list of nodes of a fixed type
    !> e.g., among the list of y-Faces.     An offset needs to be
    !> added to get index in list of all faces (for example).
    subroutine GridIndex(self, nodeType, indVec, j, k)
        !
        class(Grid2D_t), intent( in ) :: self
        character(*), intent( in ) :: nodeType
        integer, intent( in ) :: indVec(:)
        integer, intent( out ) :: j(:), k(:)
        !
        integer :: ny, nz, nVec, ii
        real(4) :: rNy

        call self%Limits(nodeType, ny, nz)
        nVec = size(indVec)
        
        if(nVec.NE.size(j)) then
             print *, 'Size of "ind_vec" and "j" do not agree.'
             stop
        endif
        
        if(nVec.NE.size(k)) then
             print *, 'Size of "ind_vec" and "k" do not agree.'
             stop
        endif
        
        rNy = float(ny)
        
        do ii = 1, nVec
             j(ii) = mod(indVec(ii), ny)
             k(ii) = ceiling(float(indVec(ii))/rNy)
        enddo
        
        where(j.EQ.0) j = ny
        where(k.EQ.0) k = nz

    end subroutine GridIndex
    !
    !> Procedure VectorIndex
    !
    !> Based on matlab method of same name in class Grid_t3D
    !> returned array IndVec gives numbering of nodes within
    !> the list for nodeType; need to add an offset for position
    !> in full list of all faces or edges (not nodes and cells).
    subroutine VectorIndex(self, nodeType, j, k, indVec)
        !
        class(Grid2D_t), intent( in ) :: self
        character(*), intent( in ) :: nodeType
        integer, intent( in ) :: j(:), k(:)
        integer, intent( out ) :: indVec(:)
        !
        integer :: ny, nz, nVec, ii
        
        call self%Limits(nodeType, ny, nz)
        
        nVec = size(indVec)
        
        if(nVec.NE.size (J)) then
             print *, 'Size of "ind_vec" and "j" do not agree.'
             stop
        endif
        
        if(nVec.NE.size (K)) then
             print *, 'Size of "ind_vec" and "k" do not agree.'
             stop
        endif
        
        do ii = 1, nVec
             indVec(ii) = (K(ii) - 1) * ny + j(ii)
        enddo
        
    end subroutine VectorIndex
    !
    !> No subroutine briefing
    subroutine Limits(self, nodeType, ny, nz)
        !
        class(Grid2D_t), intent( in ) :: self
        character(*), intent( in ) :: nodeType
        integer, intent( out ) :: ny, nz

        !>     need to either make node types explicit, or include the parameter definitions
        !>         from abstract grid class
        select case(nodeType)
             case(CELL, CELL_EARTH)
                    ny = self%ny
                    nz = self%nz

             case(NODE)
                    ny = self%ny + 1
                    nz = self%nz + 1

             case(YEDGE)
                    ny = self%ny
                    nz = self%nz + 1

             case(ZEDGE)
                    ny = self%ny + 1
                    nz = self%nz
             
             case DEFAULT
                    print*,'requested nodeType not defined for 2D grids'
                    stop
        end select
        
    end subroutine Limits
    !
    !> No subroutine briefing
    !
    function IsallocateGrid2Dd( self ) result(f)
        !
        class(Grid2D_t), intent( in ) :: self
        !
        logical :: f

        f = self%is_allocated
    end function IsallocateGrid2Dd
    !
    !> No subroutine briefing
    !
    function Length( self ) result(n)
        class(Grid2D_t), intent( in ) :: self
        integer :: n

        n = self%ny*self%nz
    end function Length
    !
    !> No subroutine briefing
    !
    function Slice1D( self ) result(g1D)
        implicit none
        class(Grid2D_t), intent( in ) :: self
        type(Grid1D_t) :: g1d
        !
        !>     not much to this -- just createGrid2D Grid1D object
        g1D = Grid1d_t( self%nzAir, self%nzEarth, self%dz )
        !
    end function Slice1D

end module Grid2D
