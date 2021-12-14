module MetricElements
   !
   use Grid
   use rVector
   use rScalar
   !
   type, abstract :: MetricElements_t
       !
       class( Grid_t ), pointer :: grid
       !
       class( rVector_t ), allocatable :: EdgeLength
       class( rVector_t ), allocatable :: FaceArea
       class( rVector_t ), allocatable :: DualFaceArea
       class( rVector_t ), allocatable :: DualEdgeLength
       class( rScalar_t ), allocatable :: Vnode
       class( rScalar_t ), allocatable :: Vcell
       class( rVector_t ), allocatable :: Vedge
       !
    contains
       !
       procedure( interface_set_edge_length_metric_elements ), deferred, public      :: SetEdgeLength
       procedure( interface_set_face_area_metric_elements ), deferred, public        :: SetFaceArea
       procedure( interface_set_dual_edge_length_metric_elements ), deferred, public :: SetDualEdgeLength
       procedure( interface_set_dual_face_area_metric_elements ), deferred, public   :: SetDualFaceArea
       procedure( interface_set_cell_volume_metric_elements ), deferred, public      :: SetCellVolume
       procedure( interface_set_edge_volume_metric_elements ), deferred, public      :: SetEdgeVolume
       procedure( interface_set_node_volume_metric_elements ), deferred, public      :: SetNodeVolume
       !
       procedure, public :: setMetricElements
       !
   end type MetricElements_t
   !
   abstract interface
       !
       subroutine interface_set_edge_length_metric_elements(self)
          import :: MetricElements_t
          class(MetricElements_t), intent( inout ) :: self
       end subroutine interface_set_edge_length_metric_elements
       !
       subroutine interface_set_face_area_metric_elements(self)
          import :: MetricElements_t
          class(MetricElements_t), intent( inout ) :: self
       end subroutine interface_set_face_area_metric_elements
       !
       subroutine interface_set_dual_edge_length_metric_elements(self)
          import :: MetricElements_t
          class(MetricElements_t), intent( inout ) :: self
       end subroutine interface_set_dual_edge_length_metric_elements
       !
       subroutine interface_set_dual_face_area_metric_elements(self)
          import :: MetricElements_t
          class(MetricElements_t), intent( inout ) :: self
       end subroutine interface_set_dual_face_area_metric_elements
       !
       subroutine interface_set_cell_volume_metric_elements(self)
          import :: MetricElements_t
          class(MetricElements_t), intent( inout ) :: self
       end subroutine interface_set_cell_volume_metric_elements
       !
       subroutine interface_set_edge_volume_metric_elements(self)
          import :: MetricElements_t
          class(MetricElements_t), intent( inout ) :: self
       end subroutine interface_set_edge_volume_metric_elements
       !
       subroutine interface_set_node_volume_metric_elements(self)
          import :: MetricElements_t
          class(MetricElements_t), intent( inout ) :: self
       end subroutine interface_set_node_volume_metric_elements
       !
   end interface
   !
contains
   !
   subroutine SetMetricElements( self )
      implicit none
      !
      class( MetricElements_t ), intent( inout ) :: self
      !
      call self%SetEdgeLength()
      call self%SetFaceArea()
      call self%SetDualEdgeLength()
      call self%SetDualFaceArea()
      call self%SetCellVolume()
      call self%SetNodeVolume()
      call self%SetEdgeVolume()
      !
   end subroutine SetMetricElements
   !
end module MetricElements
