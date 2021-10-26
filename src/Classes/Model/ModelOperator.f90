module ModelOperator
  !
  use Constants
  use Grid
  use cVector
  use cScalar
  use MetricElements
  use ModelParameter
  !
  type, abstract :: ModelOperator_t
     !
     class( Grid_t ), pointer  :: grid
     real( kind = prec )       :: omega
     logical                   :: is_allocated = .false.
     !
     contains 
     !
     procedure( iface_UpdateFrequency ), deferred, public :: UpdateFrequency
     procedure( iface_SetEquations ), deferred, public    :: SetEquations
     procedure( iface_SetCond )     , deferred, public    :: SetCond
     procedure( iface_AMult )       , deferred, public    :: Amult
     procedure( iface_MultAib )     , deferred, public    :: MultAib
     procedure( iface_MultCurlT )   , deferred, public    :: MultCurlT
     
  end type ModelOperator_t
  
  abstract interface

     subroutine iface_UpdateFrequency(self, omega)
       import :: ModelOperator_t, prec
       class(ModelOperator_t) :: self 
       real(kind = prec), intent(in) :: omega
     end subroutine Iface_UpdateFrequency
     
     !**
     ! SetEquations
     !*
     subroutine iface_SetEquations(self)
       import :: ModelOperator_t
       class(ModelOperator_t), intent(inout) :: self
     end subroutine iface_SetEquations

     !**
     ! SetCond
     subroutine iface_SetCond(self, CondParam) 
       import :: ModelOperator_t, ModelParameter_t
       class(ModelOperator_t) , intent(inout) :: self
       class(ModelParameter_t), intent(inout)    :: CondParam
     end subroutine iface_SetCond
     
     !**
     ! MultAib
     !*
     function iface_MultAib(self, bdry) result(outE)
       import :: ModelOperator_t, cVector_t
       class(ModelOperator_t), intent(in) :: self
       class(cVector_t)      , intent(in)  :: bdry
       class(cVector_t), allocatable :: outE
     end function iface_multAib
     
     !**
     ! MultCurlT
     !*
     subroutine iface_MultCurlT(self, inH, outE)
       import :: ModelOperator_t, cVector_t
       class(ModelOperator_t) , intent(in) :: self
       class(cVector_t)       , intent(in) :: inH
       class(cVector_t)       , intent(out), allocatable :: outE       
     end subroutine iface_MultCurlT
     
     function iface_Amult(self, x, p_adjt) result(y)
       import :: ModelOperator_t, cVector_t, prec
       class(ModelOperator_t), intent(in)  :: self
       class(cVector_t)      , intent(in)  :: x
       logical               , intent(in), optional :: p_adjt    
       class(cVector_t), allocatable :: y
     end function iface_Amult
          
  end interface
  
end module ModelOperator
