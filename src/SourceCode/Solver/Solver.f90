!
!> Abstract Base class to define a Solver
!
module Solver
    !
    use PreConditioner
    !
    !> Solver parameters
    integer :: max_solver_iters, max_divcor_calls, max_divcor_iters
    !
    real( kind=prec ) :: tolerance_divcor, tolerance_solver
    !
    character(:), allocatable :: solver_type
    character( len=3 ), parameter :: QMR = "QMR"
    character( len=3 ), parameter :: PCG = "PCG"
    character( len=4 ), parameter :: BiCG = "BiCG"
    !
    !> Solver Base Type
    type, abstract :: Solver_t
        !
        integer :: max_iters, n_iter
        real( kind=prec ) :: omega, tolerance
        real( kind=prec ), allocatable :: relErr(:) !> relative error at each iteration
        !
        logical :: failed, converged
        !
        class( PreConditioner_t ), allocatable :: preconditioner
        !
        contains
            !
            procedure( interface_set_solver_defaults ), deferred, public :: setDefaults
            !
            procedure, public :: baseInit => initializeSolver
            procedure, public :: baseDealloc => deallocateSolver
            !
            procedure, public :: setParameters => setParametersSolver
            procedure, public :: zeroDiagnostics => zeroDiagnosticsSolver
            !
    end type Solver_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_set_solver_defaults( self )
           import :: Solver_t
           class( Solver_t ), intent( inout ) :: self
        end subroutine interface_set_solver_defaults
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine setParametersSolver( self, max_iters, tolerance )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        integer, intent( in ) :: max_iters
        real( kind=prec ), intent( in ) :: tolerance
        !
        self%max_iters = max_iters
        !
        allocate( self%relErr( max_iters ) )
        !
        self%relErr = R_ZERO
        !
        self%tolerance = tolerance
        !
    end subroutine setParametersSolver
    !
    !> No subroutine briefing
    subroutine zeroDiagnosticsSolver( self )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        !
        self%n_iter = 0
        self%relErr = R_ZERO
        !
    end subroutine zeroDiagnosticsSolver 
    !
    !> No subroutine briefing
    subroutine initializeSolver( self )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        !
        self%max_iters = 0
        self%n_iter = 0
        self%omega = R_ZERO
        self%tolerance = R_ZERO
        !
        self%failed = .FALSE.
        self%converged = .FALSE.
        !
    end subroutine initializeSolver
    !
    !> No subroutine briefing
    subroutine deallocateSolver( self )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        !
        if( allocated( self%preconditioner ) ) deallocate( self%preconditioner )
        !
        if( allocated( self%relErr ) ) deallocate( self%relErr )
        !
    end subroutine deallocateSolver
    !
end module Solver
