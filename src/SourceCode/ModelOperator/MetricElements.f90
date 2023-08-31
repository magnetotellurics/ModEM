!
!> Abstract Base class to define MetricElements
!
module MetricElements
    !
    use cVector3D_MR
    use iScalar3D_SG
    !
    type, abstract :: MetricElements_t
        !
        !> Generic grid
        class( Grid_t ), pointer :: grid
        !
        class( Vector_t ), allocatable :: edge_length, dual_edge_length
        !
        class( Vector_t ), allocatable :: face_area, dual_face_area
        !
        class( Vector_t ), allocatable :: v_edge
        !
        class( Scalar_t ), allocatable :: v_node, v_cell
        !
     contains
        !
        procedure( interface_set_edge_length_metric_elements ), deferred, public :: setEdgeLength
        procedure( interface_set_face_area_metric_elements ), deferred, public :: setFaceArea
        procedure( interface_set_dual_edge_length_metric_elements ), deferred, public :: setDualEdgeLength
        procedure( interface_set_dual_face_area_metric_elements ), deferred, public :: setDualFaceArea
        procedure( interface_set_cell_volume_metric_elements ), deferred, public :: setCellVolume
        procedure( interface_set_edge_volume_metric_elements ), deferred, public :: setEdgeVolume
        procedure( interface_set_node_volume_metric_elements ), deferred, public :: setNodeVolume
        !
        procedure( interface_set_grid_index_arrays_metric_elements ), deferred, public :: setGridIndexArrays
        !
        !procedure( interface_boundary_index_metric_elements ), deferred, public :: boundaryIndex
        !
        procedure( interface_create_scalar_metric_elements ), deferred, public :: createScalar
        !
        procedure( interface_create_vector_metric_elements ), deferred, public :: createVector
        !
        procedure, public :: baseDealloc => deallocate_MetricElements
        !
        procedure, public :: alloc => allocate_MetricElements
        !
        procedure, public :: setup => setup_MetricElements
        !
    end type MetricElements_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_edge_length_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_edge_length_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_face_area_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_face_area_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_dual_edge_length_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_dual_edge_length_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_dual_face_area_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_dual_face_area_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_cell_volume_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_cell_volume_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_edge_volume_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_edge_volume_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_node_volume_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_node_volume_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_grid_index_arrays_metric_elements( self, grid )
            import :: MetricElements_t, Grid_t
            !
            class( MetricElements_t ), intent( in) :: self
            class( Grid_t ), intent( inout ) :: grid
            !
        end subroutine interface_set_grid_index_arrays_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_limits_metric_elements( self, node_type, nx, ny, nz )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( in ) :: self
            character(*), intent( in ) :: node_type
            integer, intent( out ) :: nx, ny, nz
            !
        end subroutine interface_set_limits_metric_elements
        ! 
        !subroutine interface_boundary_index_metric_elements( self, grid_type, INDb, INDi )
            ! import :: MetricElements_t
            ! !
            ! class( MetricElements_t ), intent( in ) :: self
            ! character(*), intent( in ) :: grid_type
            ! integer, allocatable, dimension(:), intent( inout ) :: INDb, INDi
            ! !
        ! end subroutine interface_boundary_index_metric_elements
        !
        subroutine interface_create_scalar_metric_elements( self, scalar_type, grid_type, scalar )
            import :: MetricElements_t, Scalar_t
            !
            class( MetricElements_t ), intent( in ) :: self
            integer, intent( in ) :: scalar_type
            character( len=4 ), intent( in ) :: grid_type
            class( Scalar_t ), allocatable, intent( out ) :: scalar
            !
        end subroutine interface_create_scalar_metric_elements
        !
        !> Create proper vector from the Grid
        !
        subroutine interface_create_vector_metric_elements( self, vector_type, grid_type, vector )
            import :: MetricElements_t, Vector_t
            !
            class( MetricElements_t ), intent( in ) :: self
            integer, intent( in ) :: vector_type
            character( len=4 ), intent( in ) :: grid_type
            class( Vector_t ), allocatable, intent( out ) :: vector
            !
        end subroutine interface_create_vector_metric_elements
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    !
    subroutine setup_MetricElements( self )
        implicit none
        !
        class( MetricElements_t ), intent( inout ) :: self
        !
        call self%setEdgeLength
        call self%setFaceArea
        call self%setDualEdgeLength
        call self%setDualFaceArea
        call self%setEdgeVolume
        !
        call self%setCellVolume
        call self%setNodeVolume
        !
    end subroutine setup_MetricElements
    !
    !> No subroutine briefing
    !
    subroutine allocate_MetricElements( self )
        implicit none
        !
        class( MetricElements_t ), intent( inout ) :: self
        !
        call self%createVector( real_t, EDGE, self%edge_length )
        call self%createVector( real_t, FACE, self%dual_edge_length )
        call self%createVector( real_t, FACE, self%face_area )
        call self%createVector( real_t, EDGE, self%dual_face_area )
        !
        call self%createVector( real_t, EDGE, self%v_edge )
        !
        call self%createScalar( real_t, NODE, self%v_node )
        call self%createScalar( real_t, CELL, self%v_cell )
        !
    end subroutine allocate_MetricElements
    !
    !> No subroutine briefing
    !
    subroutine deallocate_MetricElements( self )
        implicit none
        !
        class( MetricElements_t ), intent( inout ) :: self
        !
        if( allocated( self%edge_length ) ) deallocate( self%edge_length )
        if( allocated( self%face_area ) ) deallocate( self%face_area )
        if( allocated( self%dual_face_area ) ) deallocate( self%dual_face_area )
        if( allocated( self%dual_edge_length ) ) deallocate( self%dual_edge_length )
        if( allocated( self%v_edge ) ) deallocate( self%v_edge )
        !
        if( allocated( self%v_node ) ) deallocate( self%v_node )
        if( allocated( self%v_cell ) ) deallocate( self%v_cell )
        !
    end subroutine deallocate_MetricElements
    !
end module MetricElements
