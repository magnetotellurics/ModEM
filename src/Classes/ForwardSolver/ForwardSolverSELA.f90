! *************
! 
! Derived class to define a MT ForwardSolver
! 
! Last modified at 26/10/2021 by Paulo Werdt
! 
! *************
! 
module ForwardSolverSELA
  !
  use ForwardSolver
  use cVector3D_SG
  !
  type, extends( ForwardSolver_t ), public :: ForwardSolverSELA_t
    !
    ! PROPERTIES HERE
    !
    contains
      !
      final :: ForwardSolverSELA_dtor
      !
      procedure, public :: getESolution => getESolutionForwardSolverSELA
      !
  end type ForwardSolverSELA_t
  !
  interface ForwardSolverSELA_t
    module procedure ForwardSolverSELA_ctor
  end interface ForwardSolverSELA_t
  !
contains
  !
  !
  function ForwardSolverSELA_ctor() result( self )
    !
    class( ForwardSolverSELA_t ), pointer :: self
    !
    !write(*,*) "Constructor ForwardSolverSELA_t"
    !
    allocate( ForwardSolverSELA_t :: self )
    !
    call self%init()
    !
  end function ForwardSolverSELA_ctor
  !
  ! Destructor
  subroutine ForwardSolverSELA_dtor( self )
    implicit none
    !
    type( ForwardSolverSELA_t ), intent( in out ) :: self
    !
    ! write(*,*) "Destructor ForwardSolverSELA_t"
    !
    call self%dealloc()
    !
  end subroutine ForwardSolverSELA_dtor
  !
  !
  function getESolutionForwardSolverSELA( self, period, imode, source ) result( e_solution )
    implicit none
    !
    class( ForwardSolverSELA_t ), intent( inout ) :: self
    real( kind=prec ), intent(in)                 :: period
    integer, intent(in)                           :: imode
    class( Source_t ), allocatable, intent( in )  :: source
    !
    class( cVector_t ), allocatable :: e_solution
    !
    write(*,*) "Implementing getESolutionForwardSolverSELA"
    !
  end function getESolutionForwardSolverSELA
  !
  !
  subroutine defineSource( self )
    !
    class( ForwardSolverSELA_t ), intent(in)  :: self
    !
    write(*,*) "defineSource ForwardSolverSELA: "
    !
  end subroutine defineSource
  !
end module ForwardSolverSELA
