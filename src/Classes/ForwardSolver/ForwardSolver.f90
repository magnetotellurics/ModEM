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
   use Solver
   !
   type, abstract :: ForwardSolver_t
      !
      real( kind=prec ) :: period
      !
      ! Solver as a property of this
      class( Solver_t ), allocatable :: solver
      !
      integer :: n_iter_total = 0
      logical :: failed = .false.
      !
   contains
      !
      procedure, public :: init    => initializeFWD
      procedure, public :: dealloc => deallocateFWD
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
      !    I have eliminated polarization here -- this information could be carried
      !      in source object, if needed (as for Forward_File)
      subroutine interface_get_e_solution_fwd( self, source , e_solution )
         import :: ForwardSolver_t, prec, cVector_t, Source_t
         !
         class( ForwardSolver_t ), intent( inout ) :: self
         !   why should source be inout???
         class( Source_t ), intent( inout )        :: source
         class( cVector_t ), intent( inout)           :: e_solution
         !
      end subroutine interface_get_e_solution_fwd
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
      class( ForwardSolver_t ), intent( inout )   :: self
      !
   end subroutine deallocateFWD
   !
end module ForwardSolver
