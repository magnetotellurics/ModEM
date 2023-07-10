!
!> Derived class to define a PreConditioner_CC_MF
!>
!> This specific version will only be used with matrix-free,
!> which is only implemented for CSG.
!
module PreConditioner_CC_MF
    !
    use PreConditioner
    use ModelOperator_MF
    !
    type, extends( PreConditioner_t ) :: PreConditioner_CC_MF_t
        !
        type( cVector3D_SG_t ) :: Dilu
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
        self%Dilu = cVector3D_SG_t( self%model_operator%metric%grid, EDGE )
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
            class is( ModelOperator_MF_t )
                !
                !> Now set interior values
                do ix = 1, model_operator%metric%grid%nx
                    do iy = 2, model_operator%metric%grid%ny
                        do iz = 2, model_operator%metric%grid%nz
                            self%Dilu%x(ix, iy, iz) = model_operator%xXO(iy,iz) + &
                            c_factor*model_operator%Sigma_E%x(ix, iy, iz)    &
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
                            c_factor*model_operator%Sigma_E%y(ix, iy, iz) &
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
                            c_factor*model_operator%Sigma_E%z(ix, iy, iz) &
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
                stop "Error: setPreConditioner_CC_MF > Unclassified ModelOperator"
            !
        end select
        !
    end subroutine setPreConditioner_CC_MF
    !
    !> Procedure LTSolvePreConditioner_CC_MF
    !> Purpose: to solve the lower triangular system (or it"s adjoint);
    !> for the d-ilu pre-conditioner.
    !
    subroutine LTSolvePreConditioner_CC_MF( self, inE, outE, adjoint )
        implicit none
        !
        class( PreConditioner_CC_MF_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: inE
        class( Vector_t ), intent( inout ) :: outE
        logical, intent( in ) :: adjoint
        !
        integer :: ix, iy, iz
        !
        !>     as usual I am cutting some of the error checking, which is not
        !>     consistent with new classes
        !
        select type( inE )
            !
            class is( cVector3D_SG_t )
                !
                select type( outE )
                    !
                    class is( cVector3D_SG_t )
                        !
                        if( .NOT. outE%is_allocated ) then
                            stop "Error: LTSolvePreConditioner_CC_MF > outE not allocated yet"
                        endif
                        !
                        !> Instantiate the ModelOperator object
                        select type( model_operator => self%model_operator )
                            !
                            class is( ModelOperator_MF_t )
                                !
                                if( .NOT. adjoint ) then
                                    !
                                    outE = inE
                                    !
                                    call outE%div( model_operator%Metric%Vedge )
                                    !
                                    do ix = 1, inE%nx
                                        do iz = 2, inE%nz
                                            do iy = 2, inE%ny
                                                outE%x(ix, iy, iz) = (outE%x(ix, iy, iz) - &
                                                outE%x(ix, iy-1, iz)*model_operator%xXY(iy, 1) - &
                                                outE%x(ix, iy, iz-1)*model_operator%xXZ(iz, 1))* &
                                                self%Dilu%x(ix, iy, iz)
                                            enddo
                                        enddo
                                    enddo

                                    do iy = 1, inE%ny
                                        do iz = 2, inE%nz
                                            do ix = 2, inE%nx
                                                outE%y(ix, iy, iz) = (outE%y(ix, iy, iz) - &
                                                outE%y(ix, iy, iz-1)*model_operator%yYZ(iz, 1) - &
                                                outE%y(ix-1, iy, iz)*model_operator%yYX(ix, 1))* &
                                                self%Dilu%y(ix, iy, iz)
                                            enddo
                                        enddo
                                    enddo

                                    do iz = 1, inE%nz
                                        do iy = 2, inE%ny
                                                do ix = 2, inE%nx
                                                outE%z(ix, iy, iz) = (outE%z(ix, iy, iz) - &
                                                outE%z(ix-1, iy, iz)*model_operator%zZX(ix, 1) - &
                                                outE%z(ix, iy-1, iz)*model_operator%zZY(iy, 1))* &
                                                self%Dilu%z(ix, iy, iz)
                                            enddo
                                        enddo
                                    enddo
                                else
                                    !> adjoint = .TRUE. -- reverse mapping in to out
                                    !>     need to make sure that outE is zero on boundaries initially -- this is not
                                    !>        done explicitly in ModEM stable!
                                    call outE%zeros
                                    !
                                    do ix = 1, inE%nx
                                        do iy = inE%ny, 2, -1
                                            do iz = inE%nz, 2, -1
                                                outE%x(ix, iy, iz) = (inE%x(ix, iy, iz) - &
                                                outE%x(ix, iy+1, iz)*model_operator%xXY(iy+1, 1) - &
                                                outE%x(ix, iy, iz+1)*model_operator%xXZ(iz+1, 1))* &
                                                conjg(self%Dilu%x(ix, iy, iz))
                                            enddo
                                        enddo
                                    enddo
                                    !
                                    !> The coefficients for y are only for the interior nodes
                                    do iy = 1, inE%ny
                                        do ix = inE%nx, 2, -1
                                            do iz = inE%nz, 2, -1
                                                outE%y(ix, iy, iz) = (inE%y(ix, iy, iz) - &
                                                outE%y(ix, iy, iz+1)*model_operator%yYZ(iz+1, 1) - &
                                                outE%y(ix+1, iy, iz)*model_operator%yYX(ix+1, 1))* &
                                                conjg(self%Dilu%y(ix, iy, iz))
                                            enddo
                                        enddo
                                    enddo
                                    !
                                    do iz = 1, inE%nz
                                        do ix = inE%nx, 2, -1
                                            do iy = inE%ny, 2, -1
                                                outE%z(ix, iy, iz) = (inE%z(ix, iy, iz) - &
                                                outE%z(ix+1, iy, iz)*model_operator%zZX(ix+1, 1) - &
                                                outE%z(ix, iy+1, iz)*model_operator%zZY(iy+1, 1))* &
                                                conjg(self%Dilu%z(ix, iy, iz))
                                            enddo
                                        enddo
                                    enddo
                                    !
                                    !>     for adjoint to the division by volume elements last
                                    call outE%div( model_operator%Metric%Vedge )
                                    !
                                endif
                                !
                            class default
                                stop "Error: LTSolvePreConditioner_CC_MF > Unclassified ModelOperator"
                        !
                        end select
                        !
                    class default
                        stop "Error: LTSolvePreConditioner_CC_MF > Unclassified OutE"
                    !
                end select
                !
            class default
                stop "Error: LTSolvePreConditioner_CC_MF > Unclassified InE"
            !
        end select
        !
    end subroutine LTSolvePreConditioner_CC_MF
    !
    !> Procedure UTSolvePreConditioner_CC_MF
    !> Purpose: to solve the upper triangular system (or it"s adjoint);
    !> for the d-ilu pre-condtioner
    !
    subroutine UTSolvePreConditioner_CC_MF( self, inE, outE, adjoint )
        implicit none
        !
        class( PreConditioner_CC_MF_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: inE
        class( Vector_t ), intent( inout ) :: outE
        logical, intent( in ) :: adjoint
        !
        integer :: ix, iy, iz
        !
        select type( inE )
            !
            class is( cVector3D_SG_t )
                !
                select type( outE )
                    !
                    class is( cVector3D_SG_t )
                        !
                        !>    to be safe, zero out outR
                        call outE%zeros
                        !
                        !> Instantiate the ModelOperator object
                        select type( model_operator => self%model_operator )
                            !
                            class is( ModelOperator_MF_t )
                                !
                                if( .NOT. adjoint ) then
                                    !> for standard upper triangular solution
                                    !
                                    do ix = 1, inE%nx
                                        do iz = inE%nz, 2, -1
                                            do iy = inE%ny, 2, -1
                                                outE%x(ix, iy, iz) = inE%x(ix, iy, iz) - &
                                                ( outE%x(ix, iy+1, iz)*model_operator%xXY(iy, 2) &
                                                + outE%x(ix, iy, iz+1)*model_operator%xXZ(iz, 2))* &
                                                self%Dilu%x(ix, iy, iz)
                                            enddo
                                        enddo
                                    enddo
                                    !
                                    do iy = 1, inE%ny
                                        do iz = inE%nz, 2, -1
                                            do ix = inE%nx, 2, -1
                                                outE%y(ix, iy, iz) = inE%y(ix, iy, iz) - &
                                                ( outE%y(ix, iy, iz+1)*model_operator%yYZ(iz, 2) &
                                                + outE%y(ix+1, iy, iz)*model_operator%yYX(ix, 2))* &
                                                self%Dilu%y(ix, iy, iz)
                                            enddo
                                        enddo
                                    enddo
                                    !
                                    do iz = 1, inE%nz
                                        do iy = inE%ny, 2, -1
                                            do ix = inE%nx, 2, -1
                                                outE%z(ix, iy, iz) = inE%z(ix, iy, iz) - &
                                                ( outE%z(ix+1, iy, iz)*model_operator%zZX(ix, 2) &
                                                + outE%z(ix, iy+1, iz)*model_operator%zZY(iy, 2))* &
                                                self%Dilu%z(ix, iy, iz)
                                            enddo
                                        enddo
                                    enddo
                                else
                                    !> adjoint = .TRUE.
                                    do ix = 1, inE%nx
                                        do iz = 2, inE%nz
                                            do iy = 2, inE%ny
                                                outE%x(ix, iy, iz) = inE%x(ix, iy, iz) &
                                                - outE%x(ix, iy-1, iz)*model_operator%xXY(iy-1, 2) &
                                                * conjg(self%Dilu%x(ix,iy-1,iz))     &
                                                - outE%x(ix, iy, iz-1)*model_operator%xXZ(iz-1, 2) &
                                                * conjg(self%Dilu%x(ix, iy, iz-1))
                                            enddo
                                        enddo
                                    enddo
                                    !
                                    do iy = 1, inE%ny
                                        do iz = 2, inE%nz
                                            do ix = 2, inE%nx
                                                outE%y(ix, iy, iz) = inE%y(ix, iy, iz) &
                                                - outE%y(ix, iy, iz-1)*model_operator%yYZ(iz-1, 2) &
                                                * conjg(self%Dilu%y(ix,iy,iz-1)) &
                                                - outE%y(ix-1, iy, iz)*model_operator%yYX(ix-1, 2) &
                                                * conjg(self%Dilu%y(ix-1, iy, iz))
                                            enddo
                                        enddo
                                    enddo
                                    !
                                    do iz = 1, inE%nz
                                        do iy = 2, inE%ny
                                            do ix = 2, inE%nx
                                                outE%z(ix, iy, iz) = inE%z(ix, iy, iz) &
                                                - outE%z(ix-1, iy, iz)*model_operator%zZX(ix-1, 2) &
                                                * conjg(self%Dilu%z(ix-1,iy,iz)) &
                                                - outE%z(ix, iy-1, iz)*model_operator%zZY(iy-1, 2) &
                                                * conjg(self%Dilu%z(ix, iy-1, iz))
                                            enddo
                                        enddo
                                    enddo
                                endif
                                !
                            class default
                                stop "Error: UTSolvePreConditioner_CC_MF > Unclassified ModelOperator"
                            !
                        end select
                    !
                    class default
                        stop "Error: UTSolvePreConditioner_CC_MF > Unclassified OutE"
                !
                end select
                !
            class default
                stop "Error: UTSolvePreConditioner_CC_MF > Unclassified InE"
            !
        end select
        !
    end subroutine UTSolvePreConditioner_CC_MF
    !
    !> Procedure LUSolvePreConditioner_CC_MF
    !> this is dummy routine required by abstract preconditioner class
    subroutine LUSolvePreConditioner_CC_MF( self, inPhi, outPhi )
        implicit none
        !
        class( PreConditioner_CC_MF_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: inPhi
        class( Scalar_t ), intent( inout ) :: outPhi
        !
        stop "Error: LUSolvePreConditioner_CC_MF not implemented"
        !
    end subroutine LUSolvePreConditioner_CC_MF
    !
end module PreConditioner_CC_MF

