module Solver
    !
    use Constants
    use cVector
    use PreConditioner
    !
    character(:), allocatable :: solver_type
    character ( len=3 ), parameter :: QMR  = "QMR"
    character ( len=3 ), parameter :: PCG  = "PCG"
    character ( len=4 ), parameter :: BiCG = "BiCG"
    !
    ! SOLVER DEFAULTS
    integer :: maxIter = 20
    real( kind=prec ) :: tolerance = 0.0000000001
    !
    type, abstract :: Solver_t
        !
        integer                        :: max_iter, n_iter
        real( kind=prec )              :: omega, tolerance
        real( kind=prec ), allocatable :: relErr(:) ! relative error at each iteration
        !
        logical :: failed, converged
        !
        ! PreConditioner as a property of this
        class( PreConditioner_t ), allocatable :: preconditioner
        !
        contains
           !    deferred (abstract class) proedures
           procedure( interface_set_solver_defaults), deferred, public :: SetDefaults
           !
           procedure, public :: init    => initializeSolver
           procedure, public :: dealloc => deallocateSolver
           !
           procedure, public :: setParameters
           procedure, public :: zeroDiagnostics
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
    subroutine setParameters( self, max_iter, tolerance )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        integer, intent( in )              :: max_iter
        real( kind=prec ), intent( in )    :: tolerance
        !
        self%max_iter = max_iter
        self%tolerance = tolerance
        !
        if( allocated( self%relErr ) ) deallocate( self%relErr )
        allocate( self%relErr( max_iter ) )
        !
    end subroutine setParameters
    !
    !********
    !
    subroutine zeroDiagnostics( self )
        implicit none
        !
        class( Solver_t ), intent( inout ) :: self
        !
        self%n_iter = 0
        self%relErr = R_ZERO
        !
    end subroutine zeroDiagnostics 
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
