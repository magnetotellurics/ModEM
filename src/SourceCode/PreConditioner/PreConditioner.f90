!
!> Abstract Base class to define a PreConditioner
!>
!> Properties declared in abstract class will be: ModOp, omega (in all
!> cases I can now think of, preconditioner depends on the full operator,
!> which depends on omega);    setPreconditioner(self,omega) should be a procedure
!> in the abstract class, since we will assume that this is always defined.
!
module PreConditioner
    !
    use ModelOperator
    !
    type, abstract :: PreConditioner_t
        !
        real( kind=prec ) :: omega
        !
        class( ModelOperator_t ), pointer :: model_operator
        !
        contains
            !
            procedure( interface_set_preconditioner ), deferred, public :: setPreconditioner
            !
            procedure( interface_ltsolve_preconditioner ), deferred, public :: LTsolve
            procedure( interface_utsolve_preconditioner ), deferred, public :: UTsolve
            procedure( interface_lusolve_preconditioner ), deferred, public :: LUsolve
            !
    end type PreConditioner_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_preconditioner( self, omega )
            import :: PreConditioner_t, prec
            !
            class( PreConditioner_t ), intent( inout ) :: self
            real( kind=prec ), intent( in ) :: omega
            !
        end subroutine interface_set_preconditioner
        !
        !> No interface subroutine briefing
        !
        subroutine interface_ltsolve_preconditioner( self, in_e, out_e, adjoint )
            import :: PreConditioner_t, Vector_t
            !
            class( PreConditioner_t ), intent( inout ) :: self
            class( Vector_t ), intent( in ) :: in_e
            class( Vector_t ), intent( inout ) :: out_e
            logical, intent( in ) :: adjoint
            !
        end subroutine interface_ltsolve_preconditioner
        !
        !> No interface subroutine briefing
        !
        subroutine interface_utsolve_preconditioner( self, in_e, out_e, adjoint )
            import :: PreConditioner_t, Vector_t
            !
            class( PreConditioner_t ), intent( inout ) :: self
            class( Vector_t ), intent( in ) :: in_e
            class( Vector_t ), intent( inout ) :: out_e
            logical, intent( in ) :: adjoint
            !
        end subroutine interface_utsolve_preconditioner
        !
        !> No interface subroutine briefing
        !
        subroutine interface_lusolve_preconditioner( self, in_phi, out_phi )
            import :: PreConditioner_t, Scalar_t
            !
            class( PreConditioner_t ), intent( inout ) :: self
            class( Scalar_t ), intent( in ) :: in_phi
            class( Scalar_t ), intent( inout ) :: out_phi
            !
        end subroutine interface_lusolve_preconditioner
        !
    end interface
    !
end module PreConditioner
!