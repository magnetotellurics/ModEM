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
    !> Global Grid related properties
    integer :: model_n_air_layer
    !
    real( kind=prec ) :: model_max_height
    !
    character(:), allocatable :: grid_type
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
        logical :: is_allocated
        !
        integer :: nx, ny, nz
        integer :: nzAir   !> Number of air layers
        integer :: nzEarth !> Number of earth layers
        !
        character( len=80 ) :: geometry
        !
        real( kind=prec ) :: ox, oy, oz
        !
        real( kind=prec ) :: rotDeg
        !
        real( kind=prec ) :: zAirThick
        !
        integer, dimension(:), allocatable :: ind_interior_edges, ind_interior_nodes
        integer, dimension(:), allocatable :: ind_boundaries_edges, ind_boundaries_nodes
        !
        real( kind=prec ), allocatable, dimension(:) :: dx, dy, dz
        real( kind=prec ), allocatable, dimension(:) :: dx_inv, dy_inv, dz_inv
        !
        real( kind=prec ), allocatable, dimension(:) :: del_x, del_y, del_z
        real( kind=prec ), allocatable, dimension(:) :: del_x_inv, del_y_inv, del_z_inv
        !
        real( kind=prec ), allocatable, dimension(:) :: x_edge, y_edge, z_edge
        !
        real( kind=prec ), allocatable, dimension(:) :: x_center, y_center, z_center
        !
        contains
            !
            !> Grid abstract interfaces
            procedure( interface_number_of_edges ), deferred, public :: numberOfEdges
            procedure( interface_number_of_faces ), deferred, public :: numberOfFaces
            procedure( interface_number_of_nodes ), deferred, public :: numberOfNodes
            procedure( interface_number_of_cells ), deferred, public :: numberOfCells
            procedure( interface_grid_index ), deferred, public :: gridIndex
            procedure( interface_vector_index ), deferred, public :: vectorIndex
            procedure( interface_limits ), deferred, public :: limits
            !
            procedure( interface_set_cell_sizes ), deferred, public :: setCellSizes
            !
            procedure( interface_slice_1d_grid ), deferred, public :: slice1D
            procedure( interface_slice_2d_grid ), deferred, public :: slice2D
            !
            !> Base Grid methods
            procedure, public :: init => initialize_Grid
            procedure, public :: dealloc => deallocate_Grid
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
    end type Grid_t
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
        !> Based on matlab method of same name in class Grid_t
        !> IndVec is the index within the list of nodes of a fixed type
        !> e.g., among the list of y-Faces.     An offset needs to be
        !> added to get index in list of all faces (for example).
        !
        subroutine interface_grid_index( self, node_type, ind_vec, i, j, k )
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            character(*), intent( in ) :: node_type
            integer, dimension (:), intent( in ) :: ind_vec
            integer, dimension (:), intent( out ) :: i, j, k
        end subroutine interface_grid_index
        !
        !> Based on matlab method of same name in class Grid_t
        !> returned array IndVec gives numbering of nodes within
        !> the list for node_type; need to add an offset for position
        !> in full list of all faces or edges (not nodes and cells).
        !
        subroutine interface_vector_index(self, node_type, i, j, k, ind_vec)
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            character(*), intent( in ) :: node_type
            integer, dimension (:), intent( in ) :: i, j, k
            integer, dimension (:), intent( out ) :: ind_vec
        end subroutine interface_vector_index
        !
        !> No interface subroutine briefing
        !
        subroutine interface_limits(self, node_type, nx, ny, nz)
            import :: Grid_t
            class( Grid_t ), intent( in ) :: self
            character(*), intent( in ) :: node_type
            integer, intent( out ) :: nx, ny, nz
        end subroutine interface_limits
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_cell_sizes( self, dx, dy, dz )
            import :: Grid_t, prec
            class( Grid_t ), intent( inout ) :: self
            real( kind=prec ), dimension(:), intent( in ) :: dx, dy, dz
        end subroutine interface_set_cell_sizes
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
        deallocate( self%dx )
        deallocate( self%dy )
        deallocate( self%dz )
        !
        deallocate( self%dx_inv )
        deallocate( self%dy_inv )
        deallocate( self%dz_inv )
        !
        deallocate( self%del_x )
        deallocate( self%del_y )
        deallocate( self%del_z )
        !
        deallocate( self%del_x_inv )
        deallocate( self%del_y_inv )
        deallocate( self%del_z_inv )
        !
        deallocate( self%x_edge )
        deallocate( self%y_edge )
        deallocate( self%z_edge )
        deallocate( self%x_center )
        deallocate( self%y_center )
        deallocate( self%z_center )
        !
        self%is_allocated = .FALSE.
        !
    end subroutine deallocate_Grid
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
end module Grid
