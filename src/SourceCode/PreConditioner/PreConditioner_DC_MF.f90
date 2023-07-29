!
!> Derived class to define a PreConditioner_DC_MF
!>
!> This is for preconditioning the divergence correction equations
!
module PreConditioner_DC_MF
    !
    use PreConditioner
    use ModelOperator_MF_SG
    !
    type, extends( PreConditioner_t ) :: PreConditioner_DC_MF_t
        !
        class( Scalar_t ), allocatable :: d
        !
        contains
            !
            procedure, public :: setPreConditioner => setPreConditioner_DC_MF
            procedure, public :: LTSolve => LTSolve_PreConditioner_DC_MF
            procedure, public :: UTSolve => UTSolve_PreConditioner_DC_MF
            procedure, public :: LUSolve => LUSolve_PreConditioner_DC_MF
            !
    end type PreConditioner_DC_MF_t
    !
    interface PreConditioner_DC_MF_t
         module procedure PreConditioner_DC_MF_ctor
    end interface PreConditioner_DC_MF_t
    !
contains
    !
    !> No subroutine briefing
    !
    function PreConditioner_DC_MF_ctor( model_operator ) result( self ) 
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        type( PreConditioner_DC_MF_t ) :: self
        !
        !write( *, * ) "Constructor PreConditioner_DC_MF_t"
        !
        self%omega = R_ZERO
        !
        self%model_operator => model_operator
        !
        call self%model_operator%metric%createScalar( complex_t, NODE, self%d )
        !
        call self%d%zeros
        !
    end function PreConditioner_DC_MF_ctor
    !
    !> SetPreConditioner -- could be an abstract routine, but in the CC case
    !>        we pass omega as a parameter, and that is not relevant here -- but since
    !>     omega is a property of that class could set, and not pass into this procedure explicitly
    subroutine setPreConditioner_DC_MF( self, omega )
        implicit none
        !
        class( PreConditioner_DC_MF_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: omega
        !
        integer :: ix,iy,iz
        complex( kind=prec ), allocatable, dimension(:, :, :) :: d_v, c_v
        !
        !write( *, * ) "setPreConditioner_DC_MF"
        !
        d_v = self%d%getV()
        !
        self%omega = omega
        !
        !> Compute inverse diagonal elements for D-ILU (interior nodes only)
        !> set top nodes to 1.0
        d_v(1,:,:) = 1.0
        d_v(:,1,:) = 1.0
        d_v(:,:,1) = 1.0
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                c_v = model_operator%c%getV()
                !
                do iz = 2, model_operator%metric%grid%nz
                    do iy = 2, model_operator%metric%grid%ny
                        do ix = 2, model_operator%metric%grid%nx
                            d_v(ix, iy, iz) = c_v(ix, iy, iz) - &
                            model_operator%db1%x(ix,iy,iz)*model_operator%db2%x(ix-1,iy,iz) * &
                            d_v(ix-1,iy,iz)- &
                            model_operator%db1%y(ix,iy,iz)*model_operator%db2%y(ix,iy-1,iz) * &
                            d_v(ix,iy-1,iz)- &
                            model_operator%db1%z(ix,iy,iz)*model_operator%db2%z(ix,iy,iz-1) * &
                            d_v(ix,iy,iz-1)
                            !
                            d_v(ix, iy, iz) = 1.0/ d_v(ix, iy, iz)
                        enddo
                    enddo
                enddo
                !
                call self%d%setV( d_v )
                !
            class default
                call errStop( "setPreConditioner_DC_MF > Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine setPreConditioner_DC_MF
    !
    !> LTsolve and UTsolve are in abstract class and must be defined -- but not used for DC which
    !>        this object will be used -- so just dummies here
    subroutine LTSolve_PreConditioner_DC_MF( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_DC_MF_t ), intent( inout ) :: self
        class( Vector_t ), intent( inout ) :: in_e, out_e
        logical, intent( in ) :: adjoint
        !
        call errStop( "LTSolve_PreConditioner_DC_MF not implemented yet" )
        !
    end subroutine LTSolve_PreConditioner_DC_MF
    !
    !> No subroutine briefing
    subroutine UTSolve_PreConditioner_DC_MF( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_DC_MF_t ), intent( inout ) :: self
        class( Vector_t ), intent( inout ) :: in_e, out_e
        logical, intent( in ) :: adjoint
        !
        call errStop( "UTSolve_PreConditioner_DC_MF not implemented yet" )
        !
    end subroutine UTSolve_PreConditioner_DC_MF
    !
    !> Procedure LUSolve_PreConditioner_DC_MF
    !> apply pre-conditioner, LU solve
    !
    !> No subroutine briefing
    subroutine LUSolve_PreConditioner_DC_MF( self, in_phi, out_phi )
        implicit none
        !
        class( PreConditioner_DC_MF_t ), intent( inout ) :: self
        class( Scalar_t ), intent( inout ) :: in_phi, out_phi
        !
        integer :: ix, iy, iz
        complex( kind=prec ), allocatable, dimension(:, :, :) :: in_phi_v, out_phi_v, d_v
        !
        !write(*,*) "LUSolve_PreConditioner_DC_MF: ", in_phi%length(), out_phi%length()
        !
        if( .NOT. in_phi%is_allocated ) then
            call errStop( "LUSolve_PreConditioner_DC_MF > in_phi not allocated yet" )
        endif
        !
        if( .NOT. out_phi%is_allocated ) then
            call errStop( "LUSolve_PreConditioner_DC_MF > out_phi not allocated yet" )
        endif
        !
        in_phi_v = in_phi%getV()
        !
        d_v = self%d%getV()
        !
        call out_phi%zeros
        !
        out_phi_v = out_phi%getV()
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                !> forward substitution (Solve lower triangular system)
                !> the coefficients are only for the interior nodes
                do iz = 2, in_phi%nz
                    do iy = 2, in_phi%ny
                        do ix = 2, in_phi%nx
                            out_phi_v(ix, iy, iz) = in_phi_v(ix, iy, iz) &
                            - out_phi_v(ix-1,iy,iz)*model_operator%db1%x(ix,iy,iz)&
                            *d_v(ix-1,iy,iz) &
                            - out_phi_v(ix,iy-1,iz)*model_operator%db1%y(ix,iy,iz)&
                            *d_v(ix,iy-1,iz) &
                            - out_phi_v(ix,iy,iz-1)*model_operator%db1%z(ix,iy,iz)&
                            *d_v(ix,iy,iz-1)
                        enddo
                    enddo
                enddo
                !
                !write(*,*) "LTsolve_Real: ", size( d_v ), size( in_phi_v ), size( out_phi_v )
                !
                !> backward substitution (Solve upper triangular system)
                !> the coefficients are only for the interior nodes
                do iz = in_phi%nz,2,-1
                    do iy = in_phi%ny,2,-1
                        do ix = in_phi%nx,2,-1
                            out_phi_v(ix, iy, iz) = (out_phi_v(ix, iy, iz) &
                            - out_phi_v(ix+1, iy, iz)*model_operator%db2%x(ix, iy, iz) &
                            - out_phi_v(ix, iy+1, iz)*model_operator%db2%y(ix, iy, iz) &
                            - out_phi_v(ix, iy, iz+1)*model_operator%db2%z(ix, iy, iz)) &
                            *d_v(ix, iy, iz)
                        enddo
                    enddo
                enddo
                !
                !write(*,*) "UTsolve_Real: ", size( d_v ), size( in_phi_v ), size( out_phi_v )
                !
                call out_phi%setV( out_phi_v )
                !
            class default
                call errStop( "LUSolve_PreConditioner_DC_MF > Unclassified ModelOperator" )
        end select
        !
    end subroutine LUSolve_PreConditioner_DC_MF
    !
end module PreConditioner_DC_MF
