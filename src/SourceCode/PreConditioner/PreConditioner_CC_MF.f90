!
!> Derived class to define a PreConditioner_CC_MF
!>
!> This specific version will only be used with matrix-free,
!> which is only implemented for CSG.
!
module PreConditioner_CC_MF
    !
    use PreConditioner
    use ModelOperator_MF_SG
    !
    type, extends( PreConditioner_t ) :: PreConditioner_CC_MF_t
        !
        class( Vector_t ), allocatable :: Dilu
        !
        contains
            !
            procedure, public :: setPreConditioner => setPreConditioner_CC_MF !> This needs to be called by Solver    object
            !
            procedure, public :: LTSolve => LTSolvePreConditioner_CC_MF !> These are left (M1) and right (M2)
            procedure, public :: UTSolve => UTSolvePreConditioner_CC_MF !> preconditioning matrices for curl-curl equation.
            procedure, public :: LUSolve => LUSolvePreConditioner_CC_MF !> preconditoner for symmetric divCGrad operator
            !
    end type PreConditioner_CC_MF_t
    !
    interface PreConditioner_CC_MF_t
        module procedure PreConditioner_CC_MF_ctor
    end interface PreConditioner_CC_MF_t
    !
contains
    !
    !> No subroutine briefing
    !
    function PreConditioner_CC_MF_ctor( model_operator ) result( self ) 
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        type( PreConditioner_CC_MF_t ) :: self
        !
        !write( *, * ) "Constructor PreConditioner_CC_MF_t"
        !
        self%omega = R_ZERO
        !
        self%model_operator => model_operator
        !
        call self%model_operator%metric%createVector( complex_t, EDGE, self%Dilu )
        !
        call self%Dilu%zeros
        !
    end function PreConditioner_CC_MF_ctor
    !
    !> SetPreConditioner
    !
    subroutine setPreConditioner_CC_MF( self, omega )
        implicit none
        !
        class( PreConditioner_CC_MF_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: omega
        !
        integer :: status, ix, iy, iz
        complex( kind=prec ) :: c_factor
        complex( kind=prec ), allocatable, dimension(:, :, :) :: dilu_x, dilu_y, dilu_z
        complex( kind=prec ), allocatable, dimension(:, :, :) :: sigma_e_x, sigma_e_y, sigma_e_z
        !
        dilu_x = self%Dilu%getX()
        dilu_y = self%Dilu%getY()
        dilu_z = self%Dilu%getZ()
        !
        !> Save omega in object, to record
        self%omega = omega
        !
        c_factor = ONE_I * omega * isign * mu_0
        !
        !> Initialize the non-interior values
        !> only the interior edge values are really used
        dilu_x(:,1,:) = C_ONE
        dilu_x(:,:,1) = C_ONE
        dilu_y(1,:,:) = C_ONE
        dilu_y(:,:,1) = C_ONE
        dilu_z(1,:,:) = C_ONE
        dilu_z(:,1,:) = C_ONE
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                sigma_e_x = model_operator%sigma_e%getX()
                sigma_e_y = model_operator%sigma_e%getY()
                sigma_e_z = model_operator%sigma_e%getZ()
                !
                !> Now set interior values
                do ix = 1, model_operator%metric%grid%nx
                    do iy = 2, model_operator%metric%grid%ny
                        do iz = 2, model_operator%metric%grid%nz
                            dilu_x(ix, iy, iz) = model_operator%xXO(iy,iz) + &
                            c_factor * sigma_e_x(ix, iy, iz)    &
                            - model_operator%xXY(iy, 1)*model_operator%xXY(iy-1, 2) &
                            *dilu_x(ix,iy-1,iz) &
                            - model_operator%xXZ(iz, 1)*model_operator%xXZ(iz-1, 2) &
                            *dilu_x(ix,iy,iz-1)
                            dilu_x(ix, iy, iz) = C_ONE/dilu_x(ix, iy, iz)
                        enddo
                    enddo
                enddo
                !
                !> The coefficients for y are only for the interior nodes
                !> but need to initialize edges for recursive algorithm.
                do iy = 1, model_operator%metric%grid%ny
                    do iz = 2, model_operator%metric%grid%nz
                        do ix = 2, model_operator%metric%grid%nx
                            dilu_y(ix, iy, iz) = model_operator%yYO(ix,iz) + &
                            c_factor * sigma_e_y(ix, iy, iz) &
                            - model_operator%yYZ(iz, 1)*model_operator%yYZ(iz-1, 2) &
                            *dilu_y(ix, iy, iz-1) &
                            - model_operator%yYX(ix, 1)*model_operator%yYX(ix-1, 2) &
                            *dilu_y(ix-1, iy, iz)
                            dilu_y(ix, iy, iz) = C_ONE/dilu_y(ix, iy, iz)
                        enddo
                    enddo
                enddo
                !
                !> The coefficients for z are only for the interior nodes
                !> but need to initialize edges for recursive algorithm.
                do iz = 1, model_operator%metric%grid%nz
                    do ix = 2, model_operator%metric%grid%nx
                        do iy = 2, model_operator%metric%grid%ny
                            dilu_z(ix, iy, iz) = model_operator%zZO(ix,iy) + &
                            c_factor * sigma_e_z(ix, iy, iz) &
                            - model_operator%zZX(ix, 1)*model_operator%zZX(ix-1, 2)*    &
                            dilu_z(ix-1, iy, iz) &
                            - model_operator%zZY(iy, 1)*model_operator%zZY(iy-1, 2) &
                            *dilu_z(ix, iy-1, iz)
                            dilu_z(ix, iy, iz) = C_ONE/dilu_z(ix, iy, iz)
                        enddo
                    enddo
                enddo
                !
                call self%Dilu%setX( dilu_x )
                call self%Dilu%setY( dilu_y )
                call self%Dilu%setZ( dilu_z )
                !
            class default
                call errStop( "setPreConditioner_CC_MF > Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine setPreConditioner_CC_MF
    !
    !> Procedure LTSolvePreConditioner_CC_MF
    !> Purpose: to solve the lower triangular system (or it"s adjoint);
    !> for the d-ilu pre-conditioner.
    !
    subroutine LTSolvePreConditioner_CC_MF( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_CC_MF_t ), intent( inout ) :: self
        class( Vector_t ), intent( inout ) :: in_e, out_e
        logical, intent( in ) :: adjoint
        !
        integer :: ix, iy, iz
        complex( kind=prec ), allocatable, dimension(:, :, :) :: in_e_x, in_e_y, in_e_z
        complex( kind=prec ), allocatable, dimension(:, :, :) :: out_e_x, out_e_y, out_e_z
        complex( kind=prec ), allocatable, dimension(:, :, :) :: dilu_x, dilu_y, dilu_z
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "LTSolvePreConditioner_CC_MF > in_e not allocated yet" )
        endif
        !
        in_e_x = in_e%getX()
        in_e_y = in_e%getY()
        in_e_z = in_e%getZ()
        !
        dilu_x = self%Dilu%getX()
        dilu_y = self%Dilu%getY()
        dilu_z = self%Dilu%getZ()
        !
        !> as usual I am cutting some of the error checking, which is not
        !> consistent with new classes
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "LTSolvePreConditioner_CC_MF > out_e not allocated yet" )
        endif
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                if( .NOT. adjoint ) then
                    !
                    out_e = in_e
                    !
                    call out_e%div( model_operator%Metric%v_edge )
                    !
                    out_e_x = out_e%getX()
                    out_e_y = out_e%getY()
                    out_e_z = out_e%getZ()
                    !
                    do ix = 1, in_e%nx
                        do iz = 2, in_e%nz
                            do iy = 2, in_e%ny
                                out_e_x(ix, iy, iz) = (out_e_x(ix, iy, iz) - &
                                out_e_x(ix, iy-1, iz)*model_operator%xXY(iy, 1) - &
                                out_e_x(ix, iy, iz-1)*model_operator%xXZ(iz, 1))* &
                                dilu_x(ix, iy, iz)
                            enddo
                        enddo
                    enddo
                    !
                    do iy = 1, in_e%ny
                        do iz = 2, in_e%nz
                            do ix = 2, in_e%nx
                                out_e_y(ix, iy, iz) = (out_e_y(ix, iy, iz) - &
                                out_e_y(ix, iy, iz-1)*model_operator%yYZ(iz, 1) - &
                                out_e_y(ix-1, iy, iz)*model_operator%yYX(ix, 1))* &
                                dilu_y(ix, iy, iz)
                            enddo
                        enddo
                    enddo
                    !
                    do iz = 1, in_e%nz
                        do iy = 2, in_e%ny
                                do ix = 2, in_e%nx
                                out_e_z(ix, iy, iz) = (out_e_z(ix, iy, iz) - &
                                out_e_z(ix-1, iy, iz)*model_operator%zZX(ix, 1) - &
                                out_e_z(ix, iy-1, iz)*model_operator%zZY(iy, 1))* &
                                dilu_z(ix, iy, iz)
                            enddo
                        enddo
                    enddo
                    !
                    call out_e%setX( out_e_x )
                    call out_e%setY( out_e_y )
                    call out_e%setZ( out_e_z )
                    !
                else
                    !> adjoint = .TRUE. -- reverse mapping in to out
                    !>     need to make sure that out_e is zero on boundaries initially -- this is not
                    !>        done explicitly in ModEM stable!
                    call out_e%zeros
                    !
                    out_e_x = out_e%getX()
                    out_e_y = out_e%getY()
                    out_e_z = out_e%getZ()
                    !
                    do ix = 1, in_e%nx
                        do iy = in_e%ny, 2, -1
                            do iz = in_e%nz, 2, -1
                                out_e_x(ix, iy, iz) = (in_e_x(ix, iy, iz) - &
                                out_e_x(ix, iy+1, iz)*model_operator%xXY(iy+1, 1) - &
                                out_e_x(ix, iy, iz+1)*model_operator%xXZ(iz+1, 1))* &
                                conjg(dilu_x(ix, iy, iz))
                            enddo
                        enddo
                    enddo
                    !
                    !> The coefficients for y are only for the interior nodes
                    do iy = 1, in_e%ny
                        do ix = in_e%nx, 2, -1
                            do iz = in_e%nz, 2, -1
                                out_e_y(ix, iy, iz) = (in_e_y(ix, iy, iz) - &
                                out_e_y(ix, iy, iz+1)*model_operator%yYZ(iz+1, 1) - &
                                out_e_y(ix+1, iy, iz)*model_operator%yYX(ix+1, 1))* &
                                conjg(dilu_y(ix, iy, iz))
                            enddo
                        enddo
                    enddo
                    !
                    do iz = 1, in_e%nz
                        do ix = in_e%nx, 2, -1
                            do iy = in_e%ny, 2, -1
                                out_e_z(ix, iy, iz) = (in_e_z(ix, iy, iz) - &
                                out_e_z(ix+1, iy, iz)*model_operator%zZX(ix+1, 1) - &
                                out_e_z(ix, iy+1, iz)*model_operator%zZY(iy+1, 1))* &
                                conjg(dilu_z(ix, iy, iz))
                            enddo
                        enddo
                    enddo
                    !
                    call out_e%setX( out_e_x )
                    call out_e%setY( out_e_y )
                    call out_e%setZ( out_e_z )
                    !
                    !>     for adjoint to the division by volume elements last
                    call out_e%div( model_operator%metric%v_edge )
                    !
                endif
                !
            class default
                call errStop( "LTSolvePreConditioner_CC_MF > Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine LTSolvePreConditioner_CC_MF
    !
    !> Procedure UTSolvePreConditioner_CC_MF
    !> Purpose: to solve the upper triangular system (or it"s adjoint);
    !> for the d-ilu pre-condtioner
    !
    subroutine UTSolvePreConditioner_CC_MF( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_CC_MF_t ), intent( inout ) :: self
        class( Vector_t ), intent( inout ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        logical, intent( in ) :: adjoint
        !
        integer :: ix, iy, iz
        complex( kind=prec ), allocatable, dimension(:, :, :) :: in_e_x, in_e_y, in_e_z
        complex( kind=prec ), allocatable, dimension(:, :, :) :: out_e_x, out_e_y, out_e_z
        complex( kind=prec ), allocatable, dimension(:, :, :) :: dilu_x, dilu_y, dilu_z
        !
        in_e_x = in_e%getX()
        in_e_y = in_e%getY()
        in_e_z = in_e%getZ()
        !
        dilu_x = self%Dilu%getX()
        dilu_y = self%Dilu%getY()
        dilu_z = self%Dilu%getZ()
        !
        !>    to be safe, zero out outR
        call out_e%zeros
        !
        out_e_x = out_e%getX()
        out_e_y = out_e%getY()
        out_e_z = out_e%getZ()
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                if( .NOT. adjoint ) then
                    !> for standard upper triangular solution
                    !
                    do ix = 1, in_e%nx
                        do iz = in_e%nz, 2, -1
                            do iy = in_e%ny, 2, -1
                                out_e_x(ix, iy, iz) = in_e_x(ix, iy, iz) - &
                                ( out_e_x(ix, iy+1, iz)*model_operator%xXY(iy, 2) &
                                + out_e_x(ix, iy, iz+1)*model_operator%xXZ(iz, 2))* &
                                dilu_x(ix, iy, iz)
                            enddo
                        enddo
                    enddo
                    !
                    do iy = 1, in_e%ny
                        do iz = in_e%nz, 2, -1
                            do ix = in_e%nx, 2, -1
                                out_e_y(ix, iy, iz) = in_e_y(ix, iy, iz) - &
                                ( out_e_y(ix, iy, iz+1)*model_operator%yYZ(iz, 2) &
                                + out_e_y(ix+1, iy, iz)*model_operator%yYX(ix, 2))* &
                                dilu_y(ix, iy, iz)
                            enddo
                        enddo
                    enddo
                    !
                    do iz = 1, in_e%nz
                        do iy = in_e%ny, 2, -1
                            do ix = in_e%nx, 2, -1
                                out_e_z(ix, iy, iz) = in_e_z(ix, iy, iz) - &
                                ( out_e_z(ix+1, iy, iz)*model_operator%zZX(ix, 2) &
                                + out_e_z(ix, iy+1, iz)*model_operator%zZY(iy, 2))* &
                                dilu_z(ix, iy, iz)
                            enddo
                        enddo
                    enddo
                else
                    !> adjoint = .TRUE.
                    do ix = 1, in_e%nx
                        do iz = 2, in_e%nz
                            do iy = 2, in_e%ny
                                out_e_x(ix, iy, iz) = in_e_x(ix, iy, iz) &
                                - out_e_x(ix, iy-1, iz)*model_operator%xXY(iy-1, 2) &
                                * conjg(dilu_x(ix,iy-1,iz))     &
                                - out_e_x(ix, iy, iz-1)*model_operator%xXZ(iz-1, 2) &
                                * conjg(dilu_x(ix, iy, iz-1))
                            enddo
                        enddo
                    enddo
                    !
                    do iy = 1, in_e%ny
                        do iz = 2, in_e%nz
                            do ix = 2, in_e%nx
                                out_e_y(ix, iy, iz) = in_e_y(ix, iy, iz) &
                                - out_e_y(ix, iy, iz-1)*model_operator%yYZ(iz-1, 2) &
                                * conjg(dilu_y(ix,iy,iz-1)) &
                                - out_e_y(ix-1, iy, iz)*model_operator%yYX(ix-1, 2) &
                                * conjg(dilu_y(ix-1, iy, iz))
                            enddo
                        enddo
                    enddo
                    !
                    do iz = 1, in_e%nz
                        do iy = 2, in_e%ny
                            do ix = 2, in_e%nx
                                out_e_z(ix, iy, iz) = in_e_z(ix, iy, iz) &
                                - out_e_z(ix-1, iy, iz)*model_operator%zZX(ix-1, 2) &
                                * conjg(dilu_z(ix-1,iy,iz)) &
                                - out_e_z(ix, iy-1, iz)*model_operator%zZY(iy-1, 2) &
                                * conjg(dilu_z(ix, iy-1, iz))
                            enddo
                        enddo
                    enddo
                    !
                endif
                !
                call out_e%setX( out_e_x )
                call out_e%setY( out_e_y )
                call out_e%setZ( out_e_z )
                !
            class default
                call errStop( "UTSolvePreConditioner_CC_MF > Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine UTSolvePreConditioner_CC_MF
    !
    !> Procedure LUSolvePreConditioner_CC_MF
    !> this is dummy routine required by abstract preconditioner class
    subroutine LUSolvePreConditioner_CC_MF( self, in_phi, out_phi )
        implicit none
        !
        class( PreConditioner_CC_MF_t ), intent( inout ) :: self
        class( Scalar_t ), intent( inout ) :: in_phi
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        call errStop( "LUSolvePreConditioner_CC_MF not implemented" )
        !
    end subroutine LUSolvePreConditioner_CC_MF
    !
end module PreConditioner_CC_MF

