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
  !
  type, abstract :: ForwardSolver_t
    !
    class( cVector_t ), allocatable :: e_solution
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
    function interface_get_e_solution_fwd( self, period, imode, source ) result( e_solution )
      import :: ForwardSolver_t, prec, cVector_t, Source_t
      !
      class( ForwardSolver_t ), intent( inout )  :: self
      real( kind=prec ), intent( in )            :: period
      integer, intent( in )                      :: imode
      class( Source_t ), allocatable, intent( in )   :: source
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
