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
        if( allocated( self%phi ) ) deallocate( self%phi )
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
        !write( *, * ) "setPreConditioner_DC_SP"
        !
        self%omega = omega
        !
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_SP_t )
                !
                if( allocated( self%phi ) ) deallocate( self%phi )
                !
                allocate( self%phi( size( model_operator%NODEi ) ) )
                !
            class default
                call errStop( "setPreConditioner_DC_SP > Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine setPreConditioner_DC_SP
    !
    !> LTsolve and UTsolve are in abstract class and must be definPhid -- but not used for DC which
    !>        this object will be used -- so just dummies here
    subroutine LTSolvePreConditioner_DC_SP( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_DC_SP_t ), intent( inout ) :: self
        class( Vector_t ), intent( inout ) :: in_e, out_e
        logical, intent( in ) :: adjoint
        !
        call errStop( "LTSolvePreConditioner_DC_SP not implemented" )
        !
    end subroutine LTSolvePreConditioner_DC_SP
    !
    !> No subroutine briefing
    subroutine UTSolvePreConditioner_DC_SP( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_DC_SP_t ), intent( inout ) :: self
        class( Vector_t ), intent( inout ) :: in_e, out_e
        logical, intent( in ) :: adjoint
        !
        call errStop( "UTSolvePreConditioner_DC_SP not implemented" )
        !
    end subroutine UTSolvePreConditioner_DC_SP
    !
    !> Procedure LUSolvePreConditioner_DC_SP
    !> apply pre-conditioner, LU solve
    !
    !> No subroutine briefing
    subroutine LUSolvePreConditioner_DC_SP( self, in_phi, out_phi )
        implicit none
        !
        class( PreConditioner_DC_SP_t ), intent( inout ) :: self
        class( Scalar_t ), intent( inout ) :: in_phi, out_phi
        !
        complex( kind=prec ), allocatable, dimension(:) :: in_phi_v, out_phi_v
        !
        !write(*,*) "LUSolvePreConditioner_DC_SP: ", in_phi%length(), out_phi%length()
        !
        if( .NOT. in_phi%is_allocated ) then
            call errStop( "LUSolvePreConditioner_DC_SP > in_phi not allocated yet" )
        endif
        !
        if( .NOT. out_phi%is_allocated ) then
            call errStop( "LUSolvePreConditioner_DC_SP > out_phi not allocated yet" )
        endif
        !
        in_phi_v = in_phi%getArray()
        !
        out_phi_v = out_phi%getArray()
        !
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_SP_t )
                !
                !write(*,*) "LTsolve_Real: ", model_operator%VDsG_L%nCol, model_operator%VDsG_L%nRow, size( in_phi_v ), size( self%phi )
                !
                call LTsolve_Real( model_operator%VDsG_L, in_phi_v, self%phi )
                !
                !write(*,*) "UTsolve_Real: ", model_operator%VDsG_U%nCol, model_operator%VDsG_U%nRow, size( self%phi ), size( out_phi_v )
                !
                call UTsolve_Real( model_operator%VDsG_U, self%phi, out_phi_v )
                !
                call out_phi%setArray( out_phi_v )
                !
            class default
                call errStop( "LUSolvePreConditioner_DC_SP > Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine LUSolvePreConditioner_DC_SP
    !
end module PreConditioner_DC_SP
