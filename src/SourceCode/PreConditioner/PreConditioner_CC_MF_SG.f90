!
!> Derived class to define a Curl-Curl PreConditioner
!> Using Matrix Free System
!
module PreConditioner_CC_MF_SG
    !
    use PreConditioner
    use ModelOperator_MF_SG
    use cVector3D_SG
    !
    type, extends( PreConditioner_t ) :: PreConditioner_CC_MF_SG_t
        !
        type( cVector3D_SG_t ) :: Dilu
        !
        contains
            !
            procedure, public :: setPreConditioner => setPreConditioner_CC_MF_SG !> This needs to be called by Solver    object
            !
            procedure, public :: LTSolve => LTSolvePreConditioner_CC_MF_SG !> These are left (M1) and right (M2)
            procedure, public :: UTSolve => UTSolvePreConditioner_CC_MF_SG !> preconditioning matrices for curl-curl equation.
            procedure, public :: LUSolve => LUSolvePreConditioner_CC_MF_SG !> preconditoner for symmetric divCGrad operator
            !
    end type PreConditioner_CC_MF_SG_t
    !
    interface PreConditioner_CC_MF_SG_t
        module procedure PreConditioner_CC_MF_SG_ctor
    end interface PreConditioner_CC_MF_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function PreConditioner_CC_MF_SG_ctor( model_operator ) result( self ) 
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        type( PreConditioner_CC_MF_SG_t ) :: self
        !
        !write( *, * ) "Constructor PreConditioner_CC_MF_SG_t"
        !
        self%omega = R_ZERO
        !
        self%model_operator => model_operator
        !
        self%Dilu = cVector3D_SG_t( self%model_operator%metric%grid, EDGE )
        !
    end function PreConditioner_CC_MF_SG_ctor
    !
    !> SetPreConditioner
    !
    subroutine setPreConditioner_CC_MF_SG( self, omega )
        implicit none
        !
        class( PreConditioner_CC_MF_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: omega
        !
        integer :: status, ix, iy, iz
        complex( kind=prec ) :: c_factor
        !
        !> Save omega in object, to record
        self%omega = omega
        !
        c_factor = ONE_I * omega * isign * mu_0
        !
        !> Initialize the non-interior values
        !> only the interior edge values are really used
        self%Dilu%x(:,1,:) = C_ONE
        self%Dilu%x(:,:,1) = C_ONE
        self%Dilu%y(1,:,:) = C_ONE
        self%Dilu%y(:,:,1) = C_ONE
        self%Dilu%z(1,:,:) = C_ONE
        self%Dilu%z(:,1,:) = C_ONE
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                !> Now set interior values
                do ix = 1, model_operator%metric%grid%nx
                    do iy = 2, model_operator%metric%grid%ny
                        do iz = 2, model_operator%metric%grid%nz
                            self%Dilu%x(ix, iy, iz) = model_operator%xXO(iy,iz) + &
                            c_factor * model_operator%sigma_e%x(ix, iy, iz)    &
                            - model_operator%xXY(iy, 1)*model_operator%xXY(iy-1, 2) &
                            *self%Dilu%x(ix,iy-1,iz) &
                            - model_operator%xXZ(iz, 1)*model_operator%xXZ(iz-1, 2) &
                            *self%Dilu%x(ix,iy,iz-1)
                            self%Dilu%x(ix, iy, iz) = C_ONE/self%Dilu%x(ix, iy, iz)
                        enddo
                    enddo
                enddo
                !
                !> The coefficients for y are only for the interior nodes
                !> but need to initialize edges for recursive algorithm.
                do iy = 1, model_operator%metric%grid%ny
                    do iz = 2, model_operator%metric%grid%nz
                        do ix = 2, model_operator%metric%grid%nx
                            self%Dilu%y(ix, iy, iz) = model_operator%yYO(ix,iz) + &
                            c_factor * model_operator%sigma_e%y(ix, iy, iz) &
                            - model_operator%yYZ(iz, 1)*model_operator%yYZ(iz-1, 2) &
                            *self%Dilu%y(ix, iy, iz-1) &
                            - model_operator%yYX(ix, 1)*model_operator%yYX(ix-1, 2) &
                            *self%Dilu%y(ix-1, iy, iz)
                            self%Dilu%y(ix, iy, iz) = C_ONE/self%Dilu%y(ix, iy, iz)
                        enddo
                    enddo
                enddo
                !
                !> The coefficients for z are only for the interior nodes
                !> but need to initialize edges for recursive algorithm.
                do iz = 1, model_operator%metric%grid%nz
                    do ix = 2, model_operator%metric%grid%nx
                        do iy = 2, model_operator%metric%grid%ny
                            self%Dilu%z(ix, iy, iz) = model_operator%zZO(ix,iy) + &
                            c_factor * model_operator%sigma_e%z(ix, iy, iz) &
                            - model_operator%zZX(ix, 1)*model_operator%zZX(ix-1, 2)*    &
                            self%Dilu%z(ix-1, iy, iz) &
                            - model_operator%zZY(iy, 1)*model_operator%zZY(iy-1, 2) &
                            *self%Dilu%z(ix, iy-1, iz)
                            self%Dilu%z(ix, iy, iz) = C_ONE/self%Dilu%z(ix, iy, iz)
                        enddo
                    enddo
                enddo
                !
            class default
                call errStop( "setPreConditioner_CC_MF_SG > Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine setPreConditioner_CC_MF_SG
    !
    !> Procedure LTSolvePreConditioner_CC_MF_SG
    !> Purpose: to solve the lower triangular system (or it"s adjoint);
    !> for the d-ilu pre-conditioner.
    !
    subroutine LTSolvePreConditioner_CC_MF_SG( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_CC_MF_SG_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        logical, intent( in ) :: adjoint
        !
        type( cVector3D_SG_t ) :: out_e_copy
        integer :: ix, iy, iz
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "LTSolvePreConditioner_CC_MF_SG > in_e not allocated yet" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "LTSolvePreConditioner_CC_MF_SG > out_e not allocated yet" )
        endif
        !
        out_e_copy = out_e
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                !> Instantiate the ModelOperator object
                select type( in_e )
                    !
                    class is( cVector3D_SG_t )
                        !
                        if( .NOT. adjoint ) then
                            !
                            out_e_copy = in_e
                            !
                            call out_e_copy%div( model_operator%Metric%v_edge )
                            !
                            do ix = 1, in_e%nx
                                do iz = 2, in_e%nz
                                    do iy = 2, in_e%ny
                                        out_e_copy%x(ix, iy, iz) = (out_e_copy%x(ix, iy, iz) - &
                                        out_e_copy%x(ix, iy-1, iz)*model_operator%xXY(iy, 1) - &
                                        out_e_copy%x(ix, iy, iz-1)*model_operator%xXZ(iz, 1))* &
                                        self%Dilu%x(ix, iy, iz)
                                    enddo
                                enddo
                            enddo
                            !
                            do iy = 1, in_e%ny
                                do iz = 2, in_e%nz
                                    do ix = 2, in_e%nx
                                        out_e_copy%y(ix, iy, iz) = (out_e_copy%y(ix, iy, iz) - &
                                        out_e_copy%y(ix, iy, iz-1)*model_operator%yYZ(iz, 1) - &
                                        out_e_copy%y(ix-1, iy, iz)*model_operator%yYX(ix, 1))* &
                                        self%Dilu%y(ix, iy, iz)
                                    enddo
                                enddo
                            enddo
                            !
                            do iz = 1, in_e%nz
                                do iy = 2, in_e%ny
                                        do ix = 2, in_e%nx
                                        out_e_copy%z(ix, iy, iz) = (out_e_copy%z(ix, iy, iz) - &
                                        out_e_copy%z(ix-1, iy, iz)*model_operator%zZX(ix, 1) - &
                                        out_e_copy%z(ix, iy-1, iz)*model_operator%zZY(iy, 1))* &
                                        self%Dilu%z(ix, iy, iz)
                                    enddo
                                enddo
                            enddo
                            !
                        else
                            !> adjoint = .TRUE. -- reverse mapping in to out
                            !>    need to make sure that out_e_copy is zero on boundaries initially -- this is not
                            !>    done explicitly in ModEM stable!
                            call out_e_copy%zeros
                            !
                            do ix = 1, in_e%nx
                                do iy = in_e%ny, 2, -1
                                    do iz = in_e%nz, 2, -1
                                        out_e_copy%x(ix, iy, iz) = (in_e%x(ix, iy, iz) - &
                                        out_e_copy%x(ix, iy+1, iz)*model_operator%xXY(iy+1, 1) - &
                                        out_e_copy%x(ix, iy, iz+1)*model_operator%xXZ(iz+1, 1))* &
                                        conjg(self%Dilu%x(ix, iy, iz))
                                    enddo
                                enddo
                            enddo
                            !
                            !> The coefficients for y are only for the interior nodes
                            do iy = 1, in_e%ny
                                do ix = in_e%nx, 2, -1
                                    do iz = in_e%nz, 2, -1
                                        out_e_copy%y(ix, iy, iz) = (in_e%y(ix, iy, iz) - &
                                        out_e_copy%y(ix, iy, iz+1)*model_operator%yYZ(iz+1, 1) - &
                                        out_e_copy%y(ix+1, iy, iz)*model_operator%yYX(ix+1, 1))* &
                                        conjg(self%Dilu%y(ix, iy, iz))
                                    enddo
                                enddo
                            enddo
                            !
                            do iz = 1, in_e%nz
                                do ix = in_e%nx, 2, -1
                                    do iy = in_e%ny, 2, -1
                                        out_e_copy%z(ix, iy, iz) = (in_e%z(ix, iy, iz) - &
                                        out_e_copy%z(ix+1, iy, iz)*model_operator%zZX(ix+1, 1) - &
                                        out_e_copy%z(ix, iy+1, iz)*model_operator%zZY(iy+1, 1))* &
                                        conjg(self%Dilu%z(ix, iy, iz))
                                    enddo
                                enddo
                            enddo
                            !
                            !>     for adjoint to the division by volume elements last
                            call out_e_copy%div( model_operator%metric%v_edge )
                            !
                        endif
                        !
                        out_e = out_e_copy
                        !
                    class default
                        call errStop( "LTSolvePreConditioner_CC_MF_SG > Unclassified in_e" )
                    !
                end select
                !
            class default
                call errStop( "LTSolvePreConditioner_CC_MF_SG > Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine LTSolvePreConditioner_CC_MF_SG
    !
    !> Procedure UTSolvePreConditioner_CC_MF_SG
    !> Purpose: to solve the upper triangular system (or it"s adjoint);
    !> for the d-ilu pre-condtioner
    !
    subroutine UTSolvePreConditioner_CC_MF_SG( self, in_e, out_e, adjoint )
        implicit none
        !
        class( PreConditioner_CC_MF_SG_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        logical, intent( in ) :: adjoint
        !
        type( cVector3D_SG_t ) :: out_e_copy
        integer :: ix, iy, iz
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "UTSolvePreConditioner_CC_MF_SG > in_e not allocated yet" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "UTSolvePreConditioner_CC_MF_SG > out_e not allocated yet" )
        endif
        !
        !> to be safe, zero out outR
        out_e_copy = out_e
        call out_e_copy%zeros
        !
        !> Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                !> Instantiate the ModelOperator object
                select type( in_e )
                    !
                    class is( cVector3D_SG_t )
                        !
                        if( .NOT. adjoint ) then
                            !> for standard upper triangular solution
                            !
                            do ix = 1, in_e%nx
                                do iz = in_e%nz, 2, -1
                                    do iy = in_e%ny, 2, -1
                                        out_e_copy%x(ix, iy, iz) = in_e%x(ix, iy, iz) - &
                                        ( out_e_copy%x(ix, iy+1, iz)*model_operator%xXY(iy, 2) &
                                        + out_e_copy%x(ix, iy, iz+1)*model_operator%xXZ(iz, 2))* &
                                        self%Dilu%x(ix, iy, iz)
                                    enddo
                                enddo
                            enddo
                            !
                            do iy = 1, in_e%ny
                                do iz = in_e%nz, 2, -1
                                    do ix = in_e%nx, 2, -1
                                        out_e_copy%y(ix, iy, iz) = in_e%y(ix, iy, iz) - &
                                        ( out_e_copy%y(ix, iy, iz+1)*model_operator%yYZ(iz, 2) &
                                        + out_e_copy%y(ix+1, iy, iz)*model_operator%yYX(ix, 2))* &
                                        self%Dilu%y(ix, iy, iz)
                                    enddo
                                enddo
                            enddo
                            !
                            do iz = 1, in_e%nz
                                do iy = in_e%ny, 2, -1
                                    do ix = in_e%nx, 2, -1
                                        out_e_copy%z(ix, iy, iz) = in_e%z(ix, iy, iz) - &
                                        ( out_e_copy%z(ix+1, iy, iz)*model_operator%zZX(ix, 2) &
                                        + out_e_copy%z(ix, iy+1, iz)*model_operator%zZY(iy, 2))* &
                                        self%Dilu%z(ix, iy, iz)
                                    enddo
                                enddo
                            enddo
                        else
                            !> adjoint = .TRUE.
                            do ix = 1, in_e%nx
                                do iz = 2, in_e%nz
                                    do iy = 2, in_e%ny
                                        out_e_copy%x(ix, iy, iz) = in_e%x(ix, iy, iz) &
                                        - out_e_copy%x(ix, iy-1, iz)*model_operator%xXY(iy-1, 2) &
                                        * conjg(self%Dilu%x(ix,iy-1,iz))     &
                                        - out_e_copy%x(ix, iy, iz-1)*model_operator%xXZ(iz-1, 2) &
                                        * conjg(self%Dilu%x(ix, iy, iz-1))
                                    enddo
                                enddo
                            enddo
                            !
                            do iy = 1, in_e%ny
                                do iz = 2, in_e%nz
                                    do ix = 2, in_e%nx
                                        out_e_copy%y(ix, iy, iz) = in_e%y(ix, iy, iz) &
                                        - out_e_copy%y(ix, iy, iz-1)*model_operator%yYZ(iz-1, 2) &
                                        * conjg(self%Dilu%y(ix,iy,iz-1)) &
                                        - out_e_copy%y(ix-1, iy, iz)*model_operator%yYX(ix-1, 2) &
                                        * conjg(self%Dilu%y(ix-1, iy, iz))
                                    enddo
                                enddo
                            enddo
                            !
                            do iz = 1, in_e%nz
                                do iy = 2, in_e%ny
                                    do ix = 2, in_e%nx
                                        out_e_copy%z(ix, iy, iz) = in_e%z(ix, iy, iz) &
                                        - out_e_copy%z(ix-1, iy, iz)*model_operator%zZX(ix-1, 2) &
                                        * conjg(self%Dilu%z(ix-1,iy,iz)) &
                                        - out_e_copy%z(ix, iy-1, iz)*model_operator%zZY(iy-1, 2) &
                                        * conjg(self%Dilu%z(ix, iy-1, iz))
                                    enddo
                                enddo
                            enddo
                            !
                        endif
                        !
                        out_e = out_e_copy
                        !
                    class default
                        call errStop( "UTSolvePreConditioner_CC_MF_SG > Unclassified in_e" )
                    !
                end select
                !
            class default
                call errStop( "UTSolvePreConditioner_CC_MF_SG > Unclassified ModelOperator" )
            !
        end select
        !
    end subroutine UTSolvePreConditioner_CC_MF_SG
    !
    !> Procedure LUSolvePreConditioner_CC_MF_SG
    !> this is dummy routine required by abstract preconditioner class
    !
    subroutine LUSolvePreConditioner_CC_MF_SG( self, in_phi, out_phi )
        implicit none
        !
        class( PreConditioner_CC_MF_SG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: in_phi
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        call errStop( "LUSolvePreConditioner_CC_MF_SG not implemented" )
        !
    end subroutine LUSolvePreConditioner_CC_MF_SG
    !
end module PreConditioner_CC_MF_SG

