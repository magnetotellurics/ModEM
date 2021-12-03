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
     class( Grid_t ), pointer :: grid
     !
     class( MetricElements_t ), allocatable :: metric
     !
     logical                   :: is_allocated = .false.
     !
     contains 
     !
     !procedure( iface_UpdateFrequency ), deferred, public :: UpdateFrequency
     procedure( iface_SetEquations ), deferred, public    :: SetEquations
     procedure( iface_SetCond )     , deferred, public    :: SetCond
     procedure( iface_AMult )       , deferred, public    :: Amult
     procedure( iface_MultAib )     , deferred, public    :: MultAib
     procedure( iface_MultCurlT )   , deferred, public    :: MultCurlT
     !   following procedures are generally used for divergence correction
     !   and might in some cases (e.g., model operator for "SP2" case)
     !   only be implemented as dummy procedures
     procedure( iface_DivCgrad )    , deferred, public    :: DivCgrad
     procedure( iface_DivC )        , deferred, public    :: DivC
     procedure( iface_Grad )        , deferred, public    :: Grad
     procedure( iface_Div )         , deferred, public    :: Div
     !   these will be coded to return cScalar/cVector of type appropriate
     !     for specific ModelOperator implementation
     !    I see no need for real versions -- but we can add if needed!
     !
	 procedure( iface_createScalar )   , deferred, public    :: createScalar
     procedure( iface_createVector )   , deferred, public    :: createVector
     
  end type ModelOperator_t
  
  abstract interface

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
     subroutine iface_MultAib(self, bdry,outE)
       import :: ModelOperator_t, cVector_t
       class(ModelOperator_t), intent(in) :: self
       class(cVector_t)      , intent(in)  :: bdry
       class(cVector_t)      , intent(inout) :: outE
     end subroutine iface_MultAib
     
     !**
     ! MultCurlT
     !*
     subroutine iface_MultCurlT(self, inH, outE)
       import :: ModelOperator_t, cVector_t
       class(ModelOperator_t) , intent(in)    :: self
       class(cVector_t)       , intent(inout) :: inH
       class(cVector_t), allocatable, intent(inout) :: outE
     end subroutine iface_MultCurlT
     
     subroutine iface_Amult(self, omega, x, y, p_adjt)
       import :: ModelOperator_t, cVector_t, prec
       class(ModelOperator_t), intent(in)           :: self
       real( kind = prec ), intent(in), optional    :: omega
       class(cVector_t)      , intent(in)           :: x
       class(cVector_t)      , intent(inout)        :: y
       logical               , intent(in), optional :: p_adjt
     end subroutine iface_Amult
!
!     these are for divergence correction, might be dummies in some cases
     subroutine iface_DivCgrad(self, inPhi, outPhi)
       import :: ModelOperator_t, cScalar_t
       class(ModelOperator_t) , intent(in) :: self
       class(cScalar_t)       , intent(in) :: inPhi
       class(cScalar_t)       , intent(inout) :: outPhi       
     end subroutine iface_DivCgrad

     subroutine iface_DivC(self, inE, outPhi)
       import :: ModelOperator_t, cVector_t, cScalar_t
       class(ModelOperator_t) , intent(in) :: self
       class(cVector_t)       , intent(in) :: inE
       class(cScalar_t)       , intent(inout) :: outPhi       
     end subroutine iface_DivC

     subroutine iface_Grad(self, inPhi, outE)
       import :: ModelOperator_t, cVector_t, cScalar_t
       class(ModelOperator_t) , intent(in) :: self
       class(cScalar_t)       , intent(in) :: inPhi
       class(cVector_t)       , intent(inout) :: outE       
     end subroutine iface_Grad

     subroutine iface_Div(self, inE, outPhi)
       import :: ModelOperator_t, cVector_t, cScalar_t
       class(ModelOperator_t) , intent(in) :: self
       class(cVector_t)       , intent(in) :: inE
       class(cScalar_t)       , intent(inout) :: outPhi       
     end subroutine iface_Div
     !
     function iface_createVector( self, gridType ) result(cVec)
       !
       import :: ModelOperator_t, cVector_t
       class(ModelOperator_t) , intent(in) :: self
       character(len=80), intent(in), optional :: gridType
       class(cVector_t), allocatable :: cVec
     end function iface_createVector

     function iface_createScalar( self, gridType ) result(cSclr)
       !
       import :: ModelOperator_t, cScalar_t
       class(ModelOperator_t) , intent(in) :: self
       character(len=80), intent(in), optional :: gridType
       class(cScalar_t), allocatable :: cSclr
     end function iface_createScalar

  end interface
  
end module ModelOperator
