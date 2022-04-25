!**
! Properties declared in abstract class will be: ModOp, omega (in all
! cases I can now think of, preconditioner depends on the full operator,
! which depends on omega);    setPreconditioner(self,omega) should be a procedure
! in the abstract class, since we will assume that this is always defined.
!*
module PreConditioner
    !
    use cVector
    use cScalar
    use ModelOperator
    !
    type, abstract :: PreConditioner_t
        !
        class( ModelOperator_t ), pointer :: model_operator
        !
        contains
            !
            procedure( iface_set_preconditioner ), deferred, public :: setPreconditioner
            !
            procedure( iface_ltsolve_preconditioner ), deferred, public :: LTsolve
            procedure( iface_utsolve_preconditioner ), deferred, public :: UTsolve
            procedure( iface_lusolve_preconditioner ), deferred, public :: LUsolve
            !
    end type PreConditioner_t
    !
    abstract interface
        !
        subroutine iface_set_preconditioner( self, omega )
            import :: PreConditioner_t, prec
            !
            class( PreConditioner_t ), intent( inout ) :: self
            real( kind=prec ), intent( in )            :: omega
        end subroutine iface_set_preconditioner
        !
        subroutine iface_ltsolve_preconditioner( self, inE, outE, adjt )
            import :: PreConditioner_t, cVector_t
            !
            class( PreConditioner_t ), intent( inout ) :: self
            class( cVector_t ), intent( in )           :: inE
            class( cVector_t ), intent( inout )        :: outE
            logical, intent( in )                      :: adjt
        end subroutine iface_ltsolve_preconditioner
        !
        subroutine iface_utsolve_preconditioner( self, inE, outE, adjt )
            import :: PreConditioner_t, cVector_t
            !
            class( PreConditioner_t ), intent( inout ) :: self
            class( cVector_t ), intent( in )           :: inE
            class( cVector_t ), intent( inout )        :: outE
            logical, intent( in )                      :: adjt
        end subroutine iface_utsolve_preconditioner
        !
        subroutine iface_lusolve_preconditioner( self, inPhi, outPhi )
            import :: PreConditioner_t, cScalar_t
            !
            class( PreConditioner_t ), intent( inout ) :: self
            class( cScalar_t ), intent( in )           :: inPhi
            class( cScalar_t ), intent( inout )        :: outPhi
        end subroutine iface_lusolve_preconditioner
        !
    end interface
    !
end module PreConditioner
