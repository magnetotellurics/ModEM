!
!> Derived class to definPhi a PreConditioner_DC_SP
!>
!> This is for preconditioning the divergence correction equations
!
module PreConditioner_DC_SP
    !
    use PreConditioner
    use cScalar3D_SG
    use ModelOperator_SP
    !
    type, extends( PreConditioner_t ) :: PreConditioner_DC_SP_t
        !
        complex( kind=prec ), allocatable, dimension(:) :: phi
        !
        contains
            !
            final :: PreConditioner_DC_SP_dtor
            !
            procedure, public :: setPreConditioner => setPreConditioner_DC_SP
            procedure, public :: LTSolve => LTSolvePreConditioner_DC_SP
            procedure, public :: UTSolve => UTSolvePreConditioner_DC_SP
            procedure, public :: LUSolve => LUSolvePreConditioner_DC_SP
            !
    end type PreConditioner_DC_SP_t
    !
    interface PreConditioner_DC_SP_t
         module procedure PreConditioner_DC_SP_ctor
    end interface PreConditioner_DC_SP_t
    !
contains
    !
    !> No subroutine briefing
    !
    function PreConditioner_DC_SP_ctor( model_operator ) result( self ) 
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        type( PreConditioner_DC_SP_t ) :: self
        !
        !write( *, * ) "Constructor PreConditioner_DC_SP_t"
        !
        self%omega = R_ZERO
        !
        self%model_operator => model_operator
        !
        !
    end function PreConditioner_DC_SP_ctor
    !
    !> PreConditioner_DC_SP destructor
    subroutine PreConditioner_DC_SP_dtor( self )
        implicit none
        !
        type( PreConditioner_DC_SP_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor PreConditioner_DC_SP_t"
        !
        deallocate( self%phi )
        !
    end subroutine PreConditioner_DC_SP_dtor
    !
    !> SetPreConditioner -- could be an abstract routinPhi, but in the CC case
    !>        we pass omega as a parameter, and that is not relevant here -- but since
    !>     omega is a property of that class could set, and not pass into this procedure explicitly
    subroutine setPreConditioner_DC_SP( self, omega )
        implicit none
        !
        class( PreConditioner_DC_SP_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: omega
        !
        integer :: ix,iy,iz
        !
        self%omega = omega
        !
        write( *, * ) "DEVELPMENT HOT SPOT:"
		write( *, * ) "self%model_operator%VDsG_L and self%model_operator%VDsG_U should be allocated before at divCorSetUpModelOperatorSP"
		write( *, * ) allocated( self%model_operator%VDsG_L ), allocated( self%model_operator%VDsG_U )
		!
    end subroutine setPreConditioner_DC_SP
    !
    !> LTsolve and UTsolve are in abstract class and must be definPhid -- but not used for DC which
    !>        this object will be used -- so just dummies here
    subroutine LTSolvePreConditioner_DC_SP( self, inPhi, outPhi, adjoint )
        implicit none
        !
        class( PreConditioner_DC_SP_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: inPhi
        class( Vector_t ), intent( inout ) :: outPhi
        logical, intent( in ) :: adjoint
        !
        stop "Error: LTsolve not coded for this pre-conditioner class"
        !
    end subroutine LTSolvePreConditioner_DC_SP
    !
    !> No subroutine briefing
    subroutine UTSolvePreConditioner_DC_SP( self, inPhi, outPhi, adjoint )
        implicit none
        !
        class( PreConditioner_DC_SP_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: inPhi
        class( Vector_t ), intent( inout ) :: outPhi
        logical, intent( in ) :: adjoint
        
        stop "Error: UTsolve not coded for this preconditioner class"
        !
    end subroutine UTSolvePreConditioner_DC_SP
    !
    !> Procedure LUSolvePreConditioner_DC_SP
    !> apply pre-conditioner, LU solve
    !
    !> No subroutine briefing
    subroutine LUSolvePreConditioner_DC_SP( self, inPhi, outphi )
        implicit none
        !
        class( PreConditioner_DC_SP_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: inPhi
        class( Scalar_t ), intent( inout ) :: outPhi
        !
        complex( kind=prec ), allocatable, dimension(:) :: temp_array_inPhi, temp_array_outPhi
        !
        temp_array_inPhi = inPhi%getArray()
        !
        call LTsolve_Real( self%model_operator%VDsG_L, temp_array_inPhi, self%phi )
        !
        call UTsolve_Real( self%model_operator%VDsG_U, self%phi, temp_array_outPhi )
        !
        deallocate( temp_array_inPhi )
        !
        call outPhi%setArray( temp_array_outPhi )
        !
        deallocate( temp_array_outPhi )
        !
    end subroutine LUSolvePreConditioner_DC_SP
    !
end module PreConditioner_DC_SP
