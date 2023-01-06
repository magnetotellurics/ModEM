!
!> Abstract Base class to define a Grid
!> Definition of the base class for grids. 
!> Also defines some auxiliary constants and data types.
!
module Grid
    !
    use Constants
    use Grid1D
    use Grid2D
    !
    character(:), allocatable :: grid_type
    character( len=12 ), parameter :: GRID_SG = "StandardGrid"
    character( len=19 ), parameter :: GRID_MR = "MultiresolutionGrid"
    !
    character(:), allocatable :: model_method
    character( len=12 ), parameter :: MM_METHOD_FIXED_H = "fixed height"
    character( len=6 ), parameter :: MM_METHOD_MIRROR  = "mirror"
    !
    integer :: model_n_air_layer
    real( kind=prec ) :: model_max_height
    !
    type, abstract :: Grid_t
        !
        character( len=80 ) :: geometry
        !
        real( kind=prec ) :: ox, oy, oz
        !
        real( kind=prec ) :: rotDeg
        !
        integer :: nx, ny, nz
        !
        integer :: nzAir   !> Number of air layers
        integer :: nzEarth !> Number of earth layers
        !
        real( kind=prec ), allocatable, dimension(:) :: dx, dy, dz
        real( kind=prec ), allocatable, dimension(:) :: dxInv, dyInv, dzInv
        !
        real( kind=prec ), allocatable, dimension(:) :: delX, delY, delZ
        real( kind=prec ), allocatable, dimension(:) :: delXInv, delYInv, delZInv
        !
        real( kind=prec ), allocatable, dimension(:) :: xEdge
        real( kind=prec ), allocatable, dimension(:) :: yEdge
        real( kind=prec ), allocatable, dimension(:) :: zEdge
        real( kind=prec ), allocatable, dimension(:) :: xCenter
        real( kind=prec ), allocatable, dimension(:) :: yCenter
        real( kind=prec ), allocatable, dimension(:) :: zCenter
        !
        real( kind=prec ) :: zAirThick
        !
        logical :: is_allocated
        !
        contains
            !
            procedure, public :: init => initializeGrid
            procedure, public :: dealloc => deallocateGrid
            !
            procedure(interface_NumberOfEdges) , deferred, public :: NumberOfEdges
            procedure(interface_NumberOfFaces) , deferred, public :: NumberOfFaces
            procedure(interface_NumberOfNodes) , deferred, public :: NumberOfNodes
            procedure(interface_GridIndex)    , deferred, public :: GridIndex
            procedure(interface_VectorIndex), deferred, public :: VectorIndex
            procedure(interface_Limits), deferred, public :: Limits
            procedure(interface_IsAllocated), deferred, public :: IsAllocated

            procedure(interface_Copy_from), deferred, public :: Copy_from

            procedure(interface_SetCellSizes), deferred, public :: SetCellSizes

            procedure(interface_slice_1d_grid), deferred, public :: Slice1D
            procedure(interface_slice_2d_grid), deferred, public :: Slice2D
            
            procedure, public :: GetDimensions => GetDimensionsGrid
            
            procedure, public :: SetOrigin
            procedure, public :: GetOrigin

            procedure, public :: SetGridRotation
            procedure, public :: GetGridRotation

            procedure, public :: SetGridGeometry
            procedure, public :: GetGridGeometry

    end type Grid_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_NumberOfEdges( self, nXedge, nYedge, nZedge )
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            integer, intent( out ) :: nXedge, nYedge, nZedge
        end subroutine interface_NumberOfEdges
        !
        !> No interface subroutine briefing
        subroutine interface_NumberOfFaces(self, nXface, nYface, nZface) 
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            integer, intent( out ) :: nXface, nYface, nZface
        end subroutine interface_NumberOfFaces
        !
        !> No interface function briefing
        function interface_NumberOfNodes(self) result(n)
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            integer :: n
        end function interface_NumberOfNodes
        !
        !> GridIndex
        !
        !> Based on matlab method of same name in class Grid_t
        !> IndVec is the index within the list of nodes of a fixed type
        !> e.g., among the list of y-Faces.     An offset needs to be
        !> added to get index in list of all faces (for example).
        subroutine interface_GridIndex(self, nodeType, indVec, i, j, k)
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            character(*), intent( in ) :: nodeType
            integer, dimension (:), intent( in ) :: indVec
            integer, dimension (:), intent( out ) :: i, j, k
        end subroutine interface_GridIndex
        !
        !> VectorIndex
        !
        !> Based on matlab method of same name in class Grid_t
        !> returned array IndVec gives numbering of nodes within
        !> the list for nodeType; need to add an offset for position
        !> in full list of all faces or edges (not nodes and cells).
        subroutine interface_VectorIndex(self, nodeType, i, j, k, indVec)
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            character(*), intent( in ) :: nodeType
            integer, dimension (:), intent( in ) :: i, j, k
            integer, dimension (:), intent( out ) :: indVec
        end subroutine interface_VectorIndex
        !
        !> No interface subroutine briefing
        subroutine interface_Limits(self, nodeType, nx, ny, nz)
            import :: Grid_t
            class( Grid_t ) , intent( in ) :: self
            character(*) , intent( in ) :: nodeType
            integer     , intent( out ) :: nx, ny, nz
        end subroutine interface_Limits
        !
        !> No interface function briefing
        function interface_IsAllocated(self) result(f)
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            logical :: f
        end function interface_IsAllocated
        !
        !> No interface subroutine briefing
        subroutine interface_Copy_from(self, g)
            import :: Grid_t
            class( Grid_t ), intent(inout) :: self
            class( Grid_t ), intent( in ) :: g
        end subroutine interface_Copy_from
        !
        !> No interface subroutine briefing
        subroutine interface_SetCellSizes(self, dx, dy, dz)
            import :: Grid_t, prec
            class( Grid_t ), intent(inout) :: self
            real( kind=prec ) , dimension(:), intent( in ) :: dx, dy, dz
        end subroutine interface_SetCellSizes
        !
        !> No interface function briefing
        function interface_slice_1d_grid(self) result( g1D )
            import :: Grid_t, Grid1D_t
            class( Grid_t ), intent( in ) :: self
            type(Grid1D_t) :: g1D
        end function interface_slice_1d_grid
        !
        !> No interface function briefing
        function interface_slice_2d_grid(self) result( g2D )
            import :: Grid_t, Grid2D_t
            class( Grid_t ), intent( in ) :: self
            type(Grid2D_t) :: g2D
        end function interface_slice_2d_grid
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine GetDimensionsGrid( self, nx, ny, nz, nzAir )
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
    end subroutine GetDimensionsGrid
    !
    !> No subroutine briefing
    subroutine SetOrigin(self, ox, oy, oz)
        implicit none
        !
        class( Grid_t ), intent(inout) :: self
        real( kind=prec ), intent( in ) :: ox, oy, oz
        
        self%ox = ox
        self%oy = oy
        self%oz = oz
    end subroutine SetOrigin
    !
    !> No subroutine briefing
    subroutine GetOrigin(self, ox, oy, oz)
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
        real( kind=prec ), intent( out ) :: ox, oy, oz
        
        ox = self%ox
        oy = self%oy
        oz = self%oz
    end subroutine GetOrigin
    !
    !> No subroutine briefing
    subroutine SetGridRotation(self, rotDeg)
        implicit none
        !
        class( Grid_t ), intent(inout) :: self
        real( kind=prec ), intent( in ) :: rotDeg
        
        self%rotDeg = rotDeg
    end subroutine SetGridRotation
    !
    !> No function briefing
    function GetGridRotation(self) result(rotDeg)
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
        !
        real( kind=prec ) :: rotDeg
        
        rotDeg = self%rotDeg
    end function GetGridRotation
    !
    !> No subroutine briefing
    subroutine SetGridGeometry(self, s)
        implicit none
        !
        class( Grid_t ), intent(inout) :: self
        character(*) , intent( in ) :: s
        !
        self%geometry = s
        !
    end subroutine SetGridGeometry
    !
    !> No function briefing
    function GetGridGeometry(self) result(s)
        implicit none
        !
        class( Grid_t ), intent( in ) :: self
        !
        character(80) :: s
        !
        s = self%geometry
        !
    end function GetGridGeometry
    !
    !> No subroutine briefing
    subroutine initializeGrid( self )
        implicit none
        !
        class( Grid_t ), intent( inout ) :: self
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
        self%zAirThick = 0
        !
        self%is_allocated = .FALSE.
        !
    end subroutine initializeGrid
    !
    !> No subroutine briefing
    subroutine deallocateGrid( self )
        implicit none
        !
        class( Grid_t ), intent( inout ) :: self
        !
        self%nx = 0
        self%ny = 0
        self%nz = 0
        !
        deallocate( self%dx )
        deallocate( self%dy )
        deallocate( self%dz )
        !
        deallocate( self%dxInv )
        deallocate( self%dyInv )
        deallocate( self%dzInv )
        !
        deallocate( self%delX )
        deallocate( self%delY )
        deallocate( self%delZ )
        !
        deallocate( self%delXInv )
        deallocate( self%delYInv )
        deallocate( self%delZInv )
        !
        deallocate( self%xEdge )
        deallocate( self%yEdge )
        deallocate( self%zEdge )
        deallocate( self%xCenter )
        deallocate( self%yCenter )
        deallocate( self%zCenter )
        !
        self%is_allocated = .FALSE.
        !
    end subroutine deallocateGrid
    !
end module Grid
