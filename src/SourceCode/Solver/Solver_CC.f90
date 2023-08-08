!
!> Abstract Base class to define a Solver
!
module Solver_CC
    !
    use Solver
    !
    !> Solver_CC Type
    type, abstract, extends( Solver_t ) :: Solver_CC_t
        !
        ! No derived properties
        !
        contains
            !
            procedure( interface_solve_solver ), deferred, public :: solve
            !
    end type Solver_CC_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_solve_solver( self, b, x )
            import :: Solver_CC_t, Vector_t
            !
            class( Solver_CC_t ), intent( inout ) :: self
            class( Vector_t ), intent( in ) :: b
            class( Vector_t ), intent( inout ) :: x
            !
        end subroutine interface_solve_solver
        !
    end interface
    !
!contains
    !
end module Solver_CC
