!
!> Derived class to define a Divergence Correction PreConditioner
!> Using Matrix Free System
!
module PreConditioner_DC_MF_SG
    !
    use PreConditioner
    use ModelOperator_MF_SG
    use cScalar3D_SG
    !
    type, extends( PreConditioner_t ) :: PreConditioner_DC_MF_SG_t
        !
        type( cScalar3D_SG_t ) :: d
        !
        contains
            !
            procedure, public :: setPreConditioner => setPreConditioner_DC_MF_SG
            procedure, public :: LTSolve => LTSolve_PreConditioner_DC_MF_SG
            procedure, public :: UTSolve => UTSolve_PreConditioner_DC_MF_SG
            procedure, public :: LUSolve => LUSolve_PreConditioner_DC_MF_SG
            !
    end type PreConditioner_DC_MF_SG_t
    !
    interface PreConditioner_DC_MF_SG_t
         module procedure PreConditioner_DC_MF_SG_ctor
    end interface PreConditioner_DC_MF_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function PreConditioner_DC_MF_SG_ctor( model_operator ) result( self ) 
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        type( PreConditioner_DC_MF_SG_t ) :: self
        !
        !write( *, * ) "Constructor PreConditioner_DC_MF_SG_t"
        !
        self%omega = R_ZERO
        !
        self%model_operator => model_operator
        !
        self%d = cScalar3D_SG_t( self%model_operator%metric%grid, NODE )
        !
        call self%d%zeros
        !
    end function PreConditioner_DC_MF_SG_ctor
    !
    !> setPreConditioner: could be an abstract routine, but in the CC case
    !>     we pass omega as a parameter, and that is not relevant here -- but since
    !>     omega is a property of that class could set, and not pass into this procedure explicitly
    !
    subroutine setPreConditioner_DC_MF_SG( self, omega )
        implicit none
        !
        class( PreConditioner_DC_MF_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: omega
        !
        integer :: ix, iy, iz
        !
        !write( *, * ) "setPreConditioner_DC_MF_SG"
        !
        self%omega = omega
        !
        !> Compute inverse diagonal elements for D-ILU (interior nodes only)
        !> set top nodes to 1.0
        self%d%v(1,:,:) = 1.0
        self%d%v(:,1,:) = 1.0
        self%d%v(:,:,1) = 1.0
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                do iz = 2, model_operator%metric%grid%nz
                    do iy = 2, model_operator%metric%grid%ny
                        do ix = 2, model_operator%metric%grid%nx
                            !
                            self%d%v(ix, iy, iz) = model_operator%c%v(ix, iy, iz) - &
                            model_operator%db1%x(ix,iy,iz)*model_operator%db2%x(ix-1,iy,iz) * &
                            self%d%v(ix-1,iy,iz)- &
                            model_operator%db1%y(ix,iy,iz)*model_operator%db2%y(ix,iy-1,iz) * &
                            self%d%v(ix,iy-1,iz)- &
                            model_operator%db1%z(ix,iy,iz)*model_operator%db2%z(ix,iy,iz-1) * &
                            self%d%v(ix,iy,iz-1)
                            !
                            self%d%v(ix, iy, iz) = 1.0/ self%d%v(ix, iy, iz)
                            !
                        enddo
                    enddo
                enddo
                !
            class default
                call errStop( "setPreConditioner_DC_MF_SG > Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine setPreConditioner_DC_MF_SG
    !
    !> LTsolve and UTsolve are in abstract class and must be defined -- but not used for DC which
    !>     this object will be used -- so just dummies here
    !
    subroutine LTSolve_PreConditioner_DC_MF_SG( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_DC_MF_SG_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        logical, intent( in ) :: adjoint
        !
        call errStop( "LTSolve_PreConditioner_DC_MF_SG not implemented yet" )
        !
    end subroutine LTSolve_PreConditioner_DC_MF_SG
    !
    !> No subroutine briefing
    !
    subroutine UTSolve_PreConditioner_DC_MF_SG( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_DC_MF_SG_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        logical, intent( in ) :: adjoint
        !
        call errStop( "UTSolve_PreConditioner_DC_MF_SG not implemented yet" )
        !
    end subroutine UTSolve_PreConditioner_DC_MF_SG
    !
    !> Procedure LUSolve_PreConditioner_DC_MF_SG
    !> apply pre-conditioner, LU solve
    !
    !> No subroutine briefing
    !
    subroutine LUSolve_PreConditioner_DC_MF_SG( self, in_phi, out_phi )
        implicit none
        !
        class( PreConditioner_DC_MF_SG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: in_phi
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        type( cScalar3D_SG_t ) :: out_phi_copy
        integer :: ix, iy, iz
        !
        if( .NOT. in_phi%is_allocated ) then
            call errStop( "LUSolve_PreConditioner_DC_MF_SG > in_phi not allocated yet" )
        endif
        !
        if( .NOT. out_phi%is_allocated ) then
            call errStop( "LUSolve_PreConditioner_DC_MF_SG > out_phi not allocated yet" )
        endif
        !
        out_phi_copy = out_phi
        call out_phi_copy%zeros
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                !> Instantiate the ModelOperator object
                select type( in_phi )
                    !
                    class is( cScalar3D_SG_t )
                        !
                        !> forward substitution (Solve lower triangular system)
                        !> the coefficients are only for the interior nodes
                        do iz = 2, in_phi%nz
                            do iy = 2, in_phi%ny
                                do ix = 2, in_phi%nx
                                    out_phi_copy%v(ix, iy, iz) = in_phi%v(ix, iy, iz) &
                                    - out_phi_copy%v(ix-1,iy,iz)*model_operator%db1%x(ix,iy,iz)&
                                    *self%d%v(ix-1,iy,iz) &
                                    - out_phi_copy%v(ix,iy-1,iz)*model_operator%db1%y(ix,iy,iz)&
                                    *self%d%v(ix,iy-1,iz) &
                                    - out_phi_copy%v(ix,iy,iz-1)*model_operator%db1%z(ix,iy,iz)&
                                    *self%d%v(ix,iy,iz-1)
                                enddo
                            enddo
                        enddo
                        !
                        !> backward substitution (Solve upper triangular system)
                        !> the coefficients are only for the interior nodes
                        do iz = in_phi%nz,2,-1
                            do iy = in_phi%ny,2,-1
                                do ix = in_phi%nx,2,-1
                                    out_phi_copy%v(ix, iy, iz) = (out_phi_copy%v(ix, iy, iz) &
                                    - out_phi_copy%v(ix+1, iy, iz)*model_operator%db2%x(ix, iy, iz) &
                                    - out_phi_copy%v(ix, iy+1, iz)*model_operator%db2%y(ix, iy, iz) &
                                    - out_phi_copy%v(ix, iy, iz+1)*model_operator%db2%z(ix, iy, iz)) &
                                    *self%d%v(ix, iy, iz)
                                enddo
                            enddo
                        enddo
                        !
                        out_phi = out_phi_copy
                        !
                    class default
                        call errStop( "UTSolvePreConditioner_CC_MF_SG > Unclassified in_phi" )
                    !
                end select
                !
            class default
                call errStop( "UTSolvePreConditioner_CC_MF_SG > Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine LUSolve_PreConditioner_DC_MF_SG
    !
end module PreConditioner_DC_MF_SG
