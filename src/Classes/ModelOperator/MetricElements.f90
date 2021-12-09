module MetricElements
  use Grid
  use rVector
  use rScalar
  
  type, abstract :: MetricElements_t
     !
	 class(Grid_t), pointer :: grid
	 !
     !   these have now been added to abstract class -- for specific extensions
     !     different types of rVector/rScalar may be needed    
     class(rVector_t), allocatable :: EdgeLength
     class(rVector_t), allocatable :: FaceArea
     class(rVector_t), allocatable :: DualFaceArea
     class(rVector_t), allocatable :: DualEdgeLength
     class(rScalar_t), allocatable :: Vnode
     class(rScalar_t), allocatable :: Vcell
     class(rVector_t), allocatable :: Vedge
   contains
     
     procedure(iface_SetEdgeLength)    , deferred, public :: SetEdgeLength
     procedure(iface_SetFaceArea)      , deferred, public :: SetFaceArea
     procedure(iface_SetDualEdgeLength), deferred, public :: SetDualEdgeLength
     procedure(iface_SetDualFaceArea)  , deferred, public :: SetDualFaceArea
     procedure(iface_SetCellVolume)    , deferred, public :: SetCellVolume
     procedure(iface_SetEdgeVolume)    , deferred, public :: SetEdgeVolume
     procedure(iface_SetNodeVolume)    , deferred, public :: SetNodeVolume
     
     procedure, public :: setMetricElements
     
  end type MetricElements_t
  
  abstract interface
     
     subroutine iface_SetEdgeLength(self)
       import :: MetricElements_t
       class(MetricElements_t), intent(inout) :: self
     end subroutine iface_SetEdgeLength
     
     subroutine iface_SetFaceArea(self)
       import :: MetricElements_t
       class(MetricElements_t), intent(inout) :: self
     end subroutine iface_SetFaceArea
     
     subroutine iface_SetDualEdgeLength(self)
       import :: MetricElements_t
       class(MetricElements_t), intent(inout) :: self
     end subroutine iface_SetDualEdgeLength
     
     subroutine iface_SetDualFaceArea(self)
       import :: MetricElements_t
       class(MetricElements_t), intent(inout) :: self
     end subroutine iface_SetDualFaceArea
     
     subroutine iface_SetCellVolume(self)
       import :: MetricElements_t
       class(MetricElements_t), intent(inout) :: self
     end subroutine iface_SetCellVolume
     
     subroutine iface_SetEdgeVolume(self)
       import :: MetricElements_t
       class(MetricElements_t), intent(inout) :: self
     end subroutine iface_SetEdgeVolume
     
     subroutine iface_SetNodeVolume(self)
       import :: MetricElements_t
       class(MetricElements_t), intent(inout) :: self
     end subroutine iface_SetNodeVolume
     
  end interface
  
contains
  
  subroutine SetMetricElements(self)
    implicit none
    !
    class(MetricElements_t), intent(inout) :: self
    
    call self%SetEdgeLength()  
    call self%SetFaceArea()  
    call self%SetDualEdgeLength()  
    call self%SetDualFaceArea()  
    call self%SetCellVolume()  
    call self%SetNodeVolume()  
    call self%SetEdgeVolume()
    
  end subroutine SetMetricElements
  
end module MetricElements
