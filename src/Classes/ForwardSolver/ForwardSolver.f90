! *************
! 
! Base class to define a ForwardSolver
! 
! Last modified at 16/08/2021 by Paulo Werdt
! 
! *************
! 
module ForwardSolver
  !
  use Constants
  use cVector
  use Source
  use ModelOperator
  !
  type, abstract :: ForwardSolver_t
    !
	class( ModelOperator_t ), allocatable :: model_operator
    class( cVector_t ), allocatable :: e_solution
	!
	integer :: nIterTotal
	logical :: failed
    !
  contains
    !
    procedure, public :: init    => initializeFWD
    procedure, public :: dealloc => deallocateFWD
    !
    procedure( interface_get_e_solution_fwd ), deferred, public :: getESolution
    !
  end type ForwardSolver_t
  !
  abstract interface
    !
    function interface_get_e_solution_fwd( self, source ) result( e_solution )
      import :: ForwardSolver_t, prec, cVector_t, Source_t
      !
      class( ForwardSolver_t ), intent( inout )    :: self
      class( Source_t ), allocatable, intent( in ) :: source
      !
      class( cVector_t ), allocatable :: e_solution
      !
    end function interface_get_e_solution_fwd
    !
  end interface
  !
  contains
  !
  subroutine initializeFWD( self )
    class( ForwardSolver_t ), intent( inout ) :: self
    !
  end subroutine initializeFWD
  !
  subroutine deallocateFWD( self )
    class( ForwardSolver_t ), intent( inout )  :: self
    !
    if( allocated( self%e_solution ) ) deallocate( self%e_solution )
    !
  end subroutine deallocateFWD
  !
end module ForwardSolver
