module Solver
   !
   use Constants
   use cVector
   use ModelOperator
   use PreConditioner
   !
   type, abstract :: Solver_t
      !
      integer                        :: max_iter, n_iter
      real( kind=prec )              :: omega, tolerance
      real( kind=prec ), allocatable :: relErr(:) ! relative error at each iteration
      !
      logical :: failed = .false., converged = .false.
      !
      class( ModelOperator_t ), pointer  :: model_operator
	  !
	  ! PreConditioner as a property of this
      class( PreConditioner_t ), allocatable :: preconditioner
      !
      contains
         !
         procedure, public :: init    => initializeSolver
         procedure, public :: dealloc => deallocateSolver
         !
         procedure, public :: setParameters
         procedure, public :: zeroDiagnostics
         !
   end type Solver_t
   !
contains
    !
    subroutine setParameters( self, max_iter, tolerance )
        !   actually, this set routine will be the same for all
        !    extensions -- should this be abstract
        ! import :: Solver_t  is this needed here????
        class( Solver_t ), intent(inout) :: self
        integer, intent( in )            :: max_iter
        real( kind=prec ), intent(in)    :: tolerance
        integer :: status

        self%max_iter = max_iter
        self%tolerance = tolerance
        !   perhaps check if relErr is already allocated; if so deallocate
        allocate(self%relErr(max_iter),STAT=status)
        !  if we are not going to check "status" of allocate, why 
        !    return this?

    end subroutine setParameters
    !************************************************
    subroutine zeroDiagnostics( self )
        !   zeros diagnostics for solver object
        class( Solver_t ), intent( inout ) :: self
        self%n_iter = 0
        self%relErr = R_ZERO
    end subroutine zeroDiagnostics 
    !
    subroutine initializeSolver( self )
      implicit none
      !
	  class( Solver_t ), intent( inout ) :: self
	  !
	  self%model_operator => null()
	  !
      self%max_iter = 0
	  self%n_iter = 0
      self%omega = 0.0
	  self%tolerance = 0.0
      !
	  self%failed = .false.
	  self%converged = .false.
	  !
   end subroutine initializeSolver
   !
   subroutine deallocateSolver( self )
      implicit none
      !
      class( Solver_t ), intent( inout ) :: self
      !
      !if( allocated( self%preconditioner ) ) deallocate( self%preconditioner )
      !
   end subroutine deallocateSolver
   !
end module Solver
