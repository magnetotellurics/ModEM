!
!> Abstract Base class to define MetricElements
!
module MetricElements
    !
    use Grid
    use Vector
    use Scalar
    !
    type, abstract :: MetricElements_t
        !
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
        procedure, public :: baseDealloc => deallocateMetricElements
        !
        procedure, public :: setMetricElements
        !
    end type MetricElements_t
    !
    abstract interface
        !
        !> No function briefing
        function interface_get_scalar_metric_elements( self, grid_type, model_operator_type ) result( scalar )
            import :: MetricElements_t, Scalar_t
            class( MetricElements_t ), intent( in ) :: self
            character( len=4 ), intent( in ), optional :: grid_type
            character( len=10 ), intent( in ), optional :: model_operator_type
            class( Scalar_t ), allocatable :: scalar
        end function interface_get_scalar_metric_elements
        !
        !> No function briefing
        function interface_get_vector_metric_elements( self, grid_type, model_operator_type ) result( vector )
            import :: MetricElements_t, Vector_t
            class( MetricElements_t ), intent( in ) :: self
            character( len=4 ), intent( in ), optional :: grid_type
            character( len=10 ), intent( in ), optional :: model_operator_type
            class( Vector_t ), allocatable :: vector
        end function interface_get_vector_metric_elements
        !
        !> No function briefing
        function interface_create_scalar_metric_elements( self, grid_type, model_operator_type ) result( scalar )
            import :: MetricElements_t, Scalar_t
            class( MetricElements_t ), intent( in ) :: self
            character( len=4 ), intent( in ), optional :: grid_type
            character( len=10 ), intent( in ), optional :: model_operator_type
            class( Scalar_t ), allocatable :: scalar
        end function interface_create_scalar_metric_elements
        !
        !> No function briefing
        function interface_create_vector_metric_elements( self, grid_type, model_operator_type ) result( vector )
            import :: MetricElements_t, Vector_t
            class( MetricElements_t ), intent( in ) :: self
            character( len=4 ), intent( in ), optional :: grid_type
            character( len=10 ), intent( in ), optional :: model_operator_type
            class( Vector_t ), allocatable :: vector
        end function interface_create_vector_metric_elements
        !
        !> No function briefing
        function interface_create_field_metric_elements( self, v_type, grid_type, model_operator_type ) result( field )
            import :: MetricElements_t, Field_t
            class( MetricElements_t ), intent( in ) :: self
            character( len=6 ), intent( in ), optional :: v_type
            character( len=4 ), intent( in ), optional :: grid_type
            character( len=10 ), intent( in ), optional :: model_operator_type
            class( Field_t ), allocatable :: field
        end function interface_create_field_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_edge_length_metric_elements( self )
            import :: MetricElements_t
            class( MetricElements_t ), intent( inout ) :: self
        end subroutine interface_set_edge_length_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_face_area_metric_elements( self )
            import :: MetricElements_t
            class( MetricElements_t ), intent( inout ) :: self
        end subroutine interface_set_face_area_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_dual_edge_length_metric_elements( self )
            import :: MetricElements_t
            class( MetricElements_t ), intent( inout ) :: self
        end subroutine interface_set_dual_edge_length_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_dual_face_area_metric_elements( self )
            import :: MetricElements_t
            class( MetricElements_t ), intent( inout ) :: self
        end subroutine interface_set_dual_face_area_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_cell_volume_metric_elements( self )
            import :: MetricElements_t
            class( MetricElements_t ), intent( inout ) :: self
        end subroutine interface_set_cell_volume_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_edge_volume_metric_elements( self )
            import :: MetricElements_t
            class( MetricElements_t ), intent( inout ) :: self
        end subroutine interface_set_edge_volume_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_node_volume_metric_elements( self )
            import :: MetricElements_t
            class( MetricElements_t ), intent( inout ) :: self
        end subroutine interface_set_node_volume_metric_elements
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine setMetricElements( self )
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
    end subroutine setMetricElements
    !
    !> No subroutine briefing
    subroutine deallocateMetricElements( self )
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
    end subroutine deallocateMetricElements
    !
end module MetricElements
