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
   use ModelParameter
   use Solver
   !
   type, abstract :: ForwardSolver_t
      !
      real( kind=prec ) :: period = 0.0
      !
      ! Solver as a property of this
      class( Solver_t ), allocatable :: solver
      !
      !   I think these iteration control/diagnositc variables  would always be used
      !     (at least for an iterative solver!)
      !
      !    control parameters: these need to be set before running
      !     set with procedures in instantiable class
      real( kind=prec ) :: tolerance = 0.0  !  target relative residuals
      integer :: max_iter_total = 0     !  limit on number of iteration
      !    iterative solver diagnostics -- these will be set with procedures in
      !      instantiable class, 
      integer :: n_iter_actual = 0           !  actual total number of iterations
      real( kind=prec ) :: relResFinal = 0.0     !  achieved relative residual
      real( kind=prec ), allocatable, dimension(:)  :: relResVec 
                             !  relative residual as a function of iteration
      logical :: failed = .false.   !  flag set to true if target relRes is not achieved
                                    !   maybe this should be an integer "status"?
      !
   contains
      !
      procedure, public :: init    => initializeFWD
      procedure, public :: dealloc => deallocateFWD
      !
      procedure( interface_set_period_fwd ), deferred, public     :: setPeriod
      procedure( interface_set_cond_fwd ), deferred, public       :: setCond
      procedure( interface_set_iter_fwd ), deferred, public       :: setIterControl
      procedure( interface_init_diag_fwd), deferred, public       :: initDiagnostics
      procedure( interface_zero_diag_fwd), deferred, public       :: zeroDiagnostics
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
      subroutine interface_set_cond_fwd( self, modPar )
         import :: ForwardSolver_t, ModelParameter_t
         !
         class( ForwardSolver_t ), intent( inout ) :: self
         class( ModelParameter_t ), intent( in ) :: modPar
         !
      end subroutine interface_set_cond_fwd
      !
      subroutine interface_set_iter_fwd( self, maxit, tol )
         import :: ForwardSolver_t, prec
         class( ForwardSolver_t ), intent( inout )  :: self
         real( kind=prec ), intent(in)                :: tol
         integer, intent(in)                        ::  maxit
         !
      end subroutine interface_set_iter_fwd
      !
      subroutine interface_init_diag_fwd( self )
         import :: ForwardSolver_t
         class( ForwardSolver_t ), intent( inout ) :: self
         !
      end subroutine interface_init_diag_fwd
      !
      subroutine interface_zero_diag_fwd( self )
         import :: ForwardSolver_t
         class( ForwardSolver_t ), intent( inout ) :: self
         !
      end subroutine interface_zero_diag_fwd
      !
      !    I have eliminated polarization here -- this information could be carried
      !      in source object, if needed (as for Forward_File)
      subroutine interface_get_e_solution_fwd( self, source , e_solution )
         import :: ForwardSolver_t, prec, cVector_t, Source_t
         !
         class( ForwardSolver_t ), intent( inout ) :: self
         !   why should source be inout???
         class( Source_t ), intent( inout )        :: source
         class( cVector_t ), intent( inout )       :: e_solution
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
