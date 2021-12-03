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
   use cVector3D_SG
   use Source
   use Solver
   use Grid3D_SG
   !
   type, abstract :: ForwardSolver_t
      !
      real( kind=prec ) :: period
      !
	  ! Pointers to previously instantiated objects
      class( Solver_t ), pointer      :: solver
      !
      class( cVector_t ), allocatable :: e_solution
      !
      integer :: n_iter_total = 0
      logical :: failed = .false.
      !
   contains
      !
      procedure, public :: init    => initializeFWD
      procedure, public :: dealloc => deallocateFWD
      !
      procedure, public :: setSolver         => setSolverFWD
      !
	  procedure( interface_set_period_fwd ), deferred, public     :: setPeriod
      procedure( interface_get_e_solution_fwd ), deferred, public :: getESolution
      !
   end type ForwardSolver_t
   !
   abstract interface
      !
      subroutine interface_set_period_fwd( self, period )
         import :: ForwardSolver_t, prec
         !
         class( ForwardSolver_t ), intent( inout ) :: self
         real( kind=prec ), intent( in )           :: period
         !
      end subroutine interface_set_period_fwd
      !
      function interface_get_e_solution_fwd( self, source, polarization ) result( e_solution )
         import :: ForwardSolver_t, prec, cVector_t, Source_t
         !
         class( ForwardSolver_t ), intent( inout ) :: self
         class( Source_t ), intent( in )           :: source
		 integer, intent( in )                     :: polarization
         !
         class( cVector_t ), allocatable           :: e_solution
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
      self%solver         => null()
	  !
   end subroutine initializeFWD
   !
   subroutine deallocateFWD( self )
      class( ForwardSolver_t ), intent( inout )   :: self
      !
      if( allocated( self%e_solution ) ) deallocate( self%e_solution )
      !
   end subroutine deallocateFWD
   !
   subroutine setSolverFWD( self, solver )
      !
      class( ForwardSolver_t ), intent( inout ) :: self
      class( Solver_t ), target, intent( in )   :: solver
      !
      self%solver => solver
      !
      ! Allocate e_solution based on the grid
      select type( grid => solver%model_operator%grid )
         class is( Grid3D_SG_t )
             allocate( self%e_solution, source=cVector3D_SG_t( grid, EDGE ) )
      end select
      !
   end subroutine setSolverFWD
   !
end module ForwardSolver
