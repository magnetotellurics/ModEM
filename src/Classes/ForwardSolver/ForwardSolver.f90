! *************
! 
! Base class to define a ForwardSolver
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
    character(:), allocatable :: forward_solver_type
    character ( len=21 ), parameter :: FWD_FILE  = "ForwardSolverFromFile"
    character ( len=15 ), parameter :: FWD_IT    = "ForwardSolverIT"
    character ( len=18 ), parameter :: FWD_IT_DC = "ForwardSolverIT_DC"
    !
    type, abstract :: ForwardSolver_t
        !
        class( Solver_t ), allocatable :: solver
        !
        real( kind=prec ) :: tolerance, relResFinal
        !
        integer           :: max_iter_total, n_iter_actual
        !
        real( kind=prec ), allocatable, dimension(:) :: relResVec
        !
        logical :: failed
        !
        contains
            !
            procedure, public :: init    => initializeForwardSolver
            procedure, public :: dealloc => deallocateForwardSolver
            !
            procedure( interface_set_frequency_fwd ), deferred, public  :: setFrequency
            procedure( interface_set_cond_fwd ), deferred, public       :: setCond
            procedure( interface_set_iter_fwd ), deferred, public       :: setIterControl
            procedure( interface_init_diag_fwd ), deferred, public      :: initDiagnostics
            procedure( interface_zero_diag_fwd ), deferred, public      :: zeroDiagnostics
            procedure( interface_get_e_solution_fwd ), deferred, public :: getESolution
            !
    end type ForwardSolver_t
    !
    abstract interface
        !
        subroutine interface_set_frequency_fwd( self, model_parameter, period )
            import :: ForwardSolver_t, ModelParameter_t, prec
            !
            class( ForwardSolver_t ), intent( inout ) :: self
			class( ModelParameter_t ), intent( in )   :: model_parameter
            real( kind=prec ), intent( in )           :: period
            !
        end subroutine interface_set_frequency_fwd
        !
        subroutine interface_set_cond_fwd( self, model_parameter )
            import :: ForwardSolver_t, ModelParameter_t
            !
            class( ForwardSolver_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in )   :: model_parameter
            !
        end subroutine interface_set_cond_fwd
        !
        subroutine interface_set_iter_fwd( self )
            import :: ForwardSolver_t
            class( ForwardSolver_t ), intent( inout ) :: self
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
        subroutine interface_get_e_solution_fwd( self, source, e_solution )
            import :: ForwardSolver_t, Source_t, cVector_t
            !
            class( ForwardSolver_t ), intent( inout ) :: self
            class( Source_t ), intent( in )           :: source
            class( cVector_t ), intent( inout )       :: e_solution
            !
        end subroutine interface_get_e_solution_fwd
        !
    end interface
    !
    contains
        !
        subroutine initializeForwardSolver( self )
            implicit none
            !
            class( ForwardSolver_t ), intent( inout ) :: self
            !
            self%tolerance = 0.0
            self%max_iter_total = 0
            self%n_iter_actual = 0
            self%relResFinal = 0.0
            !
            self%failed = .FALSE.
            !
        end subroutine initializeForwardSolver
        !
        subroutine deallocateForwardSolver( self )
            implicit none
            !
            class( ForwardSolver_t ), intent( inout ) :: self
            !
            !
            deallocate( self%solver )
            !
            deallocate( self%relResVec )
            !
        end subroutine deallocateForwardSolver
        !
end module ForwardSolver
