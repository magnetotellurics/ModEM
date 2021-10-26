!**
! Base MT Source class
! This just has the very basics needed for MT; need to build
! on this with more specialized extensions which load and
! interpolate from a nested model, compute 1D or 2D solutions
! for boundary data, or compute forcing for secondary field
! formulation.
!*
module Source_MT
  use Constants
  use cVector
  use Source
  use ModelOperator
  
  type, extends( Source_t ) :: Source_MT_t
     
     real(kind = prec) :: omega = 0.0         ! Source frequency
     character         :: polarization = 'X'  ! 'X' or 'Y'
     
     contains
        procedure, public :: setSourceParams
        procedure, public :: setRHS
        !
  end type Source_MT_T
  !
  interface Source_MT_t
    module procedure Source_MT_ctor
  end interface Source_MT_t
  !
contains
  !
  !
  function Source_MT_ctor( model_operator, omega, pol ) result( self )
    !
    type( Source_MT_t )  :: self
    class( ModelOperator_t ), target, intent(in) :: model_operator
    real( kind=prec ), intent(in)                 :: omega
    character, intent(in)                         :: pol
    !class( cVector_t ), intent(in), optional     :: E
    !
    write(*,*) "Constructor Source_MT_t"
    !
	!allocate( Source_MT_t :: self )
	!
    call self%init()
    !
    self%model_operator => model_operator
    !
    self%non_zero_source = .false.
    self%adjt = .false.
    !
    call self%setSourceParams( omega, pol )
    !
  end function Source_MT_ctor
  !
  subroutine setSourceParams( self, omega, pol )
    implicit none
    class(Source_MT_t), intent(inout)      :: self
    real(kind = prec), intent(in)          :: omega
    character, intent(in)                  :: pol
    !class(cVector_t), intent(in), optional :: E
    
    self%omega = omega
    self%polarization = pol
    
    !if ( present( E ) ) then
       !self%E = E
    !end if
    
  end subroutine setSourceParams
  
  !**
  ! setRHS
  !*
  subroutine setRHS(self)
    class(Source_MT_t), intent(inout) :: self
    ! Local variables
    class(cVector_t), allocatable:: bdry
    
    bdry = self%E%Boundary()
    
    self%rhs = self%model_operator%MultAib(bdry)
    self%rhs = C_MinusOne * self%rhs
    self%e0 = self%E%Interior()
  end subroutine setRHS
  
end module Source_MT
