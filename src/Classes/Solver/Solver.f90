module Solver
    !
    use Constants
    use cVector
    use PreConditioner
    !
    ! Solver defaults
    integer :: max_iter
    real( kind=prec ) :: tolerance
    !
    character(:), allocatable :: solver_type
    character ( len=3 ), parameter :: QMR  = "QMR"
    character ( len=3 ), parameter :: PCG  = "PCG"
    character ( len=4 ), parameter :: BiCG = "BiCG"
    !
    ! Solver Base Type
    type, abstract :: Solver_t
        !
        integer                        :: max_iter, n_iter
        real( kind=prec )              :: omega, tolerance
        real( kind=prec ), allocatable :: relErr(:) ! relative error at each iteration
        !
        logical :: failed, converged
        !
        !
        class( PreConditioner_t ), allocatable :: preconditioner
        !
        contains
           !    deferred (abstract class) proedures
           procedure( interface_set_solver_defaults), deferred, public :: setDefaults
           !
           procedure, public :: init    => initializeSolver
           procedure, public :: dealloc => deallocateSolver
           !
           procedure, public :: setParameters => setParametersSolver
           procedure, public :: zeroDiagnostics => zeroDiagnosticsSolver
           !
    end type Solver_t
    !
    abstract interface
        !
        subroutine interface_set_solver_defaults(self)
           import :: Solver_t
           class( Solver_t ), intent( inout ) :: self
        end subroutine interface_set_solver_defaults
        !
    end interface
    !
contains
    !
    subroutine setParametersSolver( self, max_iter, tolerance )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        integer, intent( in )              :: max_iter
        real( kind=prec ), intent( in )    :: tolerance
        !
        self%max_iter = max_iter
        !
        if( allocated( self%relErr ) ) deallocate( self%relErr )
        allocate( self%relErr( max_iter ) )
        !
        self%tolerance = tolerance
        !
    end subroutine setParametersSolver
    !
    !********
    !
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
    !********
    !
    subroutine initializeSolver( self )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        !
        self%n_iter = 0
        self%omega = 0.0
        !
        self%failed = .FALSE.
        self%converged = .FALSE.
        !
    end subroutine initializeSolver
    !
    subroutine deallocateSolver( self )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        !
        deallocate( self%preconditioner )
        !
        deallocate( self%relErr )
        !
    end subroutine deallocateSolver
    !
end module Solver
