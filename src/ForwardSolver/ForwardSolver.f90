!
!> Abstract Base class to define a ForwardSolver
!
module ForwardSolver
    !
    use Constants
    use Vector
    use Source
    use ModelParameter
    use Solver
    !
    character(:), allocatable :: forward_solver_type
    character( len=21 ), parameter :: FWD_FILE = "ForwardSolverFromFile"
    character( len=15 ), parameter :: FWD_IT = "ForwardSolverIT"
    character( len=18 ), parameter :: FWD_IT_DC = "ForwardSolverIT_DC"
    !
    type, abstract :: ForwardSolver_t
        !
        class( Solver_t ), allocatable :: solver
        !
        integer :: max_iter_total, n_iter_actual
        !
        real( kind=prec ) :: tolerance, relResFinal
        !
        real( kind=prec ), allocatable, dimension(:) :: relResVec
        !
        logical :: failed
        !
        contains
            !
            procedure, public :: init => initializeForwardSolver
            procedure, public :: dealloc => deallocateForwardSolver
            !
            procedure( interface_set_frequency_foward_solver ), deferred, public :: setFrequency
            !
            procedure( interface_set_iter_foward_solver ), deferred, public :: setIterControl
            !
            procedure( interface_init_diag_foward_solver ), deferred, public :: initDiagnostics
            !
            procedure( interface_zero_diag_foward_solver ), deferred, public :: zeroDiagnostics
            !
            procedure( interface_create_e_solution_foward_solver ), deferred, public :: createESolution
			!
            procedure( interface_copy_from_foward_solver ), deferred, public :: copyFrom
            generic :: assignment(=) => copyFrom
            !
    end type ForwardSolver_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_set_frequency_foward_solver( self, sigma, period )
            import :: ForwardSolver_t, ModelParameter_t, prec
            !
            class( ForwardSolver_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in ) :: sigma
            real( kind=prec ), intent( in ) :: period
            !
        end subroutine interface_set_frequency_foward_solver
        !
        !> No interface subroutine briefing
        subroutine interface_set_iter_foward_solver( self )
            !
            import :: ForwardSolver_t
            !
            class( ForwardSolver_t ), intent( inout ) :: self
            !
        end subroutine interface_set_iter_foward_solver
        !
        !> No interface subroutine briefing
        subroutine interface_init_diag_foward_solver( self )
            import :: ForwardSolver_t
            !
            class( ForwardSolver_t ), intent( inout ) :: self
            !
        end subroutine interface_init_diag_foward_solver
        !
        !> No interface subroutine briefing
        subroutine interface_zero_diag_foward_solver( self )
            import :: ForwardSolver_t
            !
            class( ForwardSolver_t ), intent( inout ) :: self
            !
        end subroutine interface_zero_diag_foward_solver
        !
        !> No interface function briefing
        subroutine interface_create_e_solution_foward_solver( self, pol, source, e_solution )
            import :: ForwardSolver_t, Source_t, Vector_t
            !
            class( ForwardSolver_t ), intent( inout ) :: self
            integer, intent( in ) :: pol
            class( Source_t ), intent( in ) :: source
            class( Vector_t ), intent( inout ) :: e_solution
            !
        end subroutine interface_create_e_solution_foward_solver
        !
        !> No interface subroutine briefing
        subroutine interface_copy_from_foward_solver( self, rhs )
            import :: ForwardSolver_t            
            class( ForwardSolver_t ), intent( inout ) :: self
            class( ForwardSolver_t ), intent( in ) :: rhs
        end subroutine interface_copy_from_foward_solver
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine initializeForwardSolver( self )
        implicit none
        !
        class( ForwardSolver_t ), intent( inout ) :: self
        !
        self%max_iter_total = 0
        !
        self%n_iter_actual = 0
        !
        self%tolerance = R_ZERO
        !
        self%relResFinal = R_ZERO
        !
        self%failed = .FALSE.
        !
    end subroutine initializeForwardSolver
    !
    !> No subroutine briefing
    subroutine deallocateForwardSolver( self )
        implicit none
        !
        class( ForwardSolver_t ), intent( inout ) :: self
        !
        deallocate( self%relResVec )
        !
        deallocate( self%solver )
        !
    end subroutine deallocateForwardSolver
    !
end module ForwardSolver
