!
!> Abstract Base class to define a Solver
!
module Solver
    !
    use Utilities
    use PreConditioner
    !
    !> Solver parameters
    integer :: max_solver_iters, max_solver_calls, max_divcor_iters
    !
    real( kind=prec ) :: tolerance_solver, tolerance_divcor
    !
    character(:), allocatable :: solver_type
    character( len=3 ), parameter :: SLV_QMR = "QMR"
    character( len=4 ), parameter :: SLV_BICG = "BICG"
    !
    !> Solver Base Type
    type, abstract :: Solver_t
        !
        integer :: iter, n_iter, max_iters
        real( kind=prec ) :: omega, tolerance
        real( kind=prec ), allocatable :: relErr(:) !> relative error at each iteration
        !
        logical :: converged
        !
        class( PreConditioner_t ), allocatable :: preconditioner
        !
        contains
            !
            procedure, public :: set => set_Solver
            !
            procedure, public :: baseInit => initialize_Solve
            procedure, public :: baseDealloc => deallocate_Solver
            !
            procedure, public :: zeroDiagnostics => zeroDiagnostics_Solver
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
    subroutine set_Solver( self, max_iters, tolerance )
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
    end subroutine set_Solver
    !
    !> No subroutine briefing
    subroutine zeroDiagnostics_Solver( self )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        !
        self%n_iter = 0
        self%relErr = R_ZERO
        !
    end subroutine zeroDiagnostics_Solver 
    !
    !> No subroutine briefing
    subroutine initialize_Solve( self )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        !
        self%max_iters = 0
        self%n_iter = 0
        self%omega = R_ZERO
        self%tolerance = R_ZERO
        !
        self%converged = .FALSE.
        !
    end subroutine initialize_Solve
    !
    !> No subroutine briefing
    subroutine deallocate_Solver( self )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        !
        if( allocated( self%preconditioner ) ) deallocate( self%preconditioner )
        !
        if( allocated( self%relErr ) ) deallocate( self%relErr )
        !
    end subroutine deallocate_Solver
    !
end module Solver
