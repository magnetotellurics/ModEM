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
        class( Vector_t ), allocatable :: EdgeLength
        class( Vector_t ), allocatable :: FaceArea
        class( Vector_t ), allocatable :: DualFaceArea
        class( Vector_t ), allocatable :: DualEdgeLength
        class( Vector_t ), allocatable :: Vedge
        !
        class( Scalar_t ), allocatable :: Vnode
        class( Scalar_t ), allocatable :: Vcell
        !
     contains
        !
        procedure, public :: dealloc => deallocateMetricElements
        !
        procedure, public :: setMetricElements
        !
        procedure( interface_set_edge_length_metric_elements ), deferred, public :: SetEdgeLength
        procedure( interface_set_face_area_metric_elements ), deferred, public :: SetFaceArea
        procedure( interface_set_dual_edge_length_metric_elements ), deferred, public :: SetDualEdgeLength
        procedure( interface_set_dual_face_area_metric_elements ), deferred, public :: SetDualFaceArea
        procedure( interface_set_cell_volume_metric_elements ), deferred, public :: SetCellVolume
        procedure( interface_set_edge_volume_metric_elements ), deferred, public :: SetEdgeVolume
        procedure( interface_set_node_volume_metric_elements ), deferred, public :: SetNodeVolume
        !
    end type MetricElements_t
    !
    abstract interface
        !
        !> No function briefing
        function interface_get_scalar_metric_elements( self, grid_type, field_type ) result( scalar )
            import :: MetricElements_t, Scalar_t
            class( MetricElements_t ), intent( in ) :: self
            character( len=4 ), intent( in ), optional :: grid_type
            character( len=10 ), intent( in ), optional :: field_type
            class( Scalar_t ), allocatable :: scalar
        end function interface_get_scalar_metric_elements
        !
        !> No function briefing
        function interface_get_vector_metric_elements( self, grid_type, field_type ) result( vector )
            import :: MetricElements_t, Vector_t
            class( MetricElements_t ), intent( in ) :: self
            character( len=4 ), intent( in ), optional :: grid_type
            character( len=10 ), intent( in ), optional :: field_type
            class( Vector_t ), allocatable :: vector
        end function interface_get_vector_metric_elements
        !
        !> No function briefing
        function interface_create_scalar_metric_elements( self, grid_type, field_type ) result( scalar )
            import :: MetricElements_t, Scalar_t
            class( MetricElements_t ), intent( in ) :: self
            character( len=4 ), intent( in ), optional :: grid_type
            character( len=10 ), intent( in ), optional :: field_type
            class( Scalar_t ), allocatable :: scalar
        end function interface_create_scalar_metric_elements
        !
        !> No function briefing
        function interface_create_vector_metric_elements( self, grid_type, field_type ) result( vector )
            import :: MetricElements_t, Vector_t
            class( MetricElements_t ), intent( in ) :: self
            character( len=4 ), intent( in ), optional :: grid_type
            character( len=10 ), intent( in ), optional :: field_type
            class( Vector_t ), allocatable :: vector
        end function interface_create_vector_metric_elements
        !
        !> No function briefing
        function interface_create_field_metric_elements( self, v_type, grid_type, field_type ) result( field )
            import :: MetricElements_t, Field_t
            class( MetricElements_t ), intent( in ) :: self
            character( len=6 ), intent( in ), optional :: v_type
            character( len=4 ), intent( in ), optional :: grid_type
            character( len=10 ), intent( in ), optional :: field_type
            class( Field_t ), allocatable :: field
        end function interface_create_field_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_edge_length_metric_elements( self )
            import :: MetricElements_t
            class(MetricElements_t), intent( inout ) :: self
        end subroutine interface_set_edge_length_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_face_area_metric_elements( self )
            import :: MetricElements_t
            class(MetricElements_t), intent( inout ) :: self
        end subroutine interface_set_face_area_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_dual_edge_length_metric_elements( self )
            import :: MetricElements_t
            class(MetricElements_t), intent( inout ) :: self
        end subroutine interface_set_dual_edge_length_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_dual_face_area_metric_elements( self )
            import :: MetricElements_t
            class(MetricElements_t), intent( inout ) :: self
        end subroutine interface_set_dual_face_area_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_cell_volume_metric_elements( self )
            import :: MetricElements_t
            class(MetricElements_t), intent( inout ) :: self
        end subroutine interface_set_cell_volume_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_edge_volume_metric_elements( self )
            import :: MetricElements_t
            class(MetricElements_t), intent( inout ) :: self
        end subroutine interface_set_edge_volume_metric_elements
        !
        !> No interface subroutine briefing
        subroutine interface_set_node_volume_metric_elements( self )
            import :: MetricElements_t
            class(MetricElements_t), intent( inout ) :: self
        end subroutine interface_set_node_volume_metric_elements
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine SetMetricElements( self )
        implicit none
        !
        class( MetricElements_t ), intent( inout ) :: self
        !
        call self%SetEdgeLength()
        call self%SetFaceArea()
        call self%SetDualEdgeLength()
        call self%SetDualFaceArea()
        call self%SetEdgeVolume()
        !
        call self%SetCellVolume()
        call self%SetNodeVolume()
        !
    end subroutine SetMetricElements
    !
    !> No subroutine briefing
    subroutine deallocateMetricElements( self )
        implicit none
        !
        class( MetricElements_t ), intent( inout ) :: self
        !
        if( allocated( self%EdgeLength ) ) deallocate( self%EdgeLength )
        if( allocated( self%FaceArea ) ) deallocate( self%FaceArea )
        if( allocated( self%DualFaceArea ) ) deallocate( self%DualFaceArea )
        if( allocated( self%DualEdgeLength ) ) deallocate( self%DualEdgeLength )
        if( allocated( self%Vedge ) ) deallocate( self%Vedge )
        !
        if( allocated( self%Vnode ) ) deallocate( self%Vnode )
        if( allocated( self%Vcell ) ) deallocate( self%Vcell )
        !
    end subroutine deallocateMetricElements
    !
end module MetricElements
