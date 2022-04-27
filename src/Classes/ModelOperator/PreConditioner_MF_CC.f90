!**
! This specific version will only be used with matrix-free,
! which is only implemented for CSG.
!*
module PreConditioner_MF_CC
    !
    use Constants
    use Grid
    use cVector
    use cVector3D_SG
    use cScalar3D_SG
    use ModelOperator_MF
    use PreConditioner
    !
    type, extends( PreConditioner_t ) :: PreConditioner_MF_CC_t
        !
        type( cVector3D_SG_t ), allocatable :: Dilu
        !
        contains
            !
            final :: PreConditioner_MF_CC_dtor
            !
            procedure, public :: setPreConditioner => setPreConditioner_MF_CC ! This needs to be called by Solver    object
            !
            procedure, public :: LTSolve => LTSolvePreConditioner_MF_CC ! These are left (M1) and right (M2)
            procedure, public :: UTSolve => UTSolvePreConditioner_MF_CC ! preconditioning matrices for curl-curl equation.
            procedure, public :: LUSolve => LUSolvePreConditioner_MF_CC ! preconditoner for symmetric divCgrad operator
            !
    end type PreConditioner_MF_CC_t
    !
    interface PreConditioner_MF_CC_t
        module procedure PreConditioner_MF_CC_ctor
    end interface PreConditioner_MF_CC_t
    !
contains
    !**
    ! Class constructor
    !*
    function PreConditioner_MF_CC_ctor( model_operator ) result( self ) 
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        type( PreConditioner_MF_CC_t ) :: self
        !
        write(*,*) "Constructor PreConditioner_MF_CC_t"
        !
        self%omega = 0.0
        !
        self%model_operator => model_operator
        !
        select type( grid => model_operator%metric%grid )
            class is( Grid3D_SG_t )
                !
                allocate( self%Dilu, source = cVector3D_SG_t( grid, EDGE ) )
                !
        end select
        !
    end function PreConditioner_MF_CC_ctor
    !
    ! Destructor
    subroutine PreConditioner_MF_CC_dtor( self )
      implicit none
      !
      type( PreConditioner_MF_CC_t ), intent( inout ) :: self
      !
      write(*,*) "Destructor PreConditioner_MF_CC"
      !
      deallocate( self%Dilu )
      !
    end subroutine PreConditioner_MF_CC_dtor
    !**
    ! SetPreConditioner
    !*
    subroutine setPreConditioner_MF_CC( self, omega )
        implicit none
        !
        class( PreConditioner_MF_CC_t ), intent( inout ) :: self
        real( kind=prec ), intent( in )                  :: omega
        !
        integer :: status, ix, iy, iz
        complex( kind=prec ) :: cFac
        !
        ! Save omega in object, to record
        self%omega = omega
        cFac = ISIGN*ONE_I*omega*MU_0
        !
        ! Initialize the non-interior values
        ! only the interior edge values are really used
        self%Dilu%x(:,1,:) = C_ONE
        self%Dilu%x(:,:,1) = C_ONE
        self%Dilu%y(1,:,:) = C_ONE
        self%Dilu%y(:,:,1) = C_ONE
        self%Dilu%z(1,:,:) = C_ONE
        self%Dilu%z(:,1,:) = C_ONE
        !
        ! Instantiate the ModelOperator object
        select type( model_operator => self%model_operator )
            !
            class is( ModelOperator_MF_t )
                !
                ! Now set interior values
                do ix = 1, model_operator%metric%grid%nx
                    do iy = 2, model_operator%metric%grid%ny
                        do iz = 2, model_operator%metric%grid%nz
                            self%Dilu%x(ix, iy, iz) = model_operator%xXO(iy,iz) + &
                            cFac*model_operator%Sigma_E%x(ix, iy, iz)    &
                            - model_operator%xXY(iy, 1)*model_operator%xXY(iy-1, 2) &
                            *self%Dilu%x(ix,iy-1,iz) &
                            - model_operator%xXZ(iz, 1)*model_operator%xXZ(iz-1, 2) &
                            *self%Dilu%x(ix,iy,iz-1)
                            self%Dilu%x(ix, iy, iz) = C_ONE/self%Dilu%x(ix, iy, iz)
                        end do
                    end do
                end do
                !
                ! The coefficients for y are only for the interior nodes
                ! but need to initialize edges for recursive algorithm.
                do iy = 1, model_operator%metric%grid%ny
                    do iz = 2, model_operator%metric%grid%nz
                        do ix = 2, model_operator%metric%grid%nx
                            self%Dilu%y(ix, iy, iz) = model_operator%yYO(ix,iz) + &
                            cFac*model_operator%Sigma_E%y(ix, iy, iz) &
                            - model_operator%yYZ(iz, 1)*model_operator%yYZ(iz-1, 2) &
                            *self%Dilu%y(ix, iy, iz-1) &
                            - model_operator%yYX(ix, 1)*model_operator%yYX(ix-1, 2) &
                            *self%Dilu%y(ix-1, iy, iz)
                            self%Dilu%y(ix, iy, iz) = C_ONE/self%Dilu%y(ix, iy, iz)
                        end do
                    end do
                end do
                !
                ! The coefficients for z are only for the interior nodes
                ! but need to initialize edges for recursive algorithm.
                do iz = 1, model_operator%metric%grid%nz
                    do ix = 2, model_operator%metric%grid%nx
                        do iy = 2, model_operator%metric%grid%ny
                            self%Dilu%z(ix, iy, iz) = model_operator%zZO(ix,iy) + &
                            cFac*model_operator%Sigma_E%z(ix, iy, iz) &
                            - model_operator%zZX(ix, 1)*model_operator%zZX(ix-1, 2)*    &
                            self%Dilu%z(ix-1, iy, iz) &
                            - model_operator%zZY(iy, 1)*model_operator%zZY(iy-1, 2) &
                            *self%Dilu%z(ix, iy-1, iz)
                            self%Dilu%z(ix, iy, iz) = C_ONE/self%Dilu%z(ix, iy, iz)
                        end do
                    end do
                end do
                !
                !
            class default
                 stop "setPreConditioner_MF_CC: Unclassified ModelOperator"
            !
        end select
        !
    end subroutine setPreConditioner_MF_CC
    
    !**
    ! Purpose: to solve the lower triangular system (or it"s adjoint);
    ! for the d-ilu pre-condtioner.
    !*
    subroutine LTSolvePreConditioner_MF_CC( self, inE, outE, adjt )
        implicit none
        !
        class( PreConditioner_MF_CC_t ), intent( inout ) :: self
        class( cVector_t ), intent( in )                 :: inE
        class( cVector_t ), intent( inout )              :: outE
        logical, intent( in )                            :: adjt
        !
        integer :: ix, iy, iz
        !
        !     as usual I am cutting some of the error checking, which is not
        !     consistent with new classes
        !
        select type( inE )
        class is( cVector3D_SG_t )
            !
            select type( outE )
            class is( cVector3D_SG_t )
                !
                if ( .NOT. outE%is_allocated ) then
                    STOP "outE in LTsolve not allocated yet"
                endif
                !
                ! Instantiate the ModelOperator object
                select type( model_operator => self%model_operator )
                    !
                    class is( ModelOperator_MF_t )
                        !
                        if ( .NOT. adjt ) then
                            !    we will need element/by element division (rdvide    in matlab)
                            !                        Call diagDiv(inE, V_E, outE)    !     this is ModEM routine
                            !     I am assuming that this TVector function implements inE./Vedge
                            !
                            outE = inE    ! assuming this works as copy?
                            !
                            call outE%divs( model_operator%Metric%Vedge )
                            !
                            do ix = 1, inE%nx
                                do iz = 2, inE%nz
                                    do iy = 2, inE%ny
                                        outE%x(ix, iy, iz) = (outE%x(ix, iy, iz) - &
                                        outE%x(ix, iy-1, iz)*model_operator%xXY(iy, 1) - &
                                        outE%x(ix, iy, iz-1)*model_operator%xXZ(iz, 1))* &
                                        self%Dilu%x(ix, iy, iz)
                                    end do
                                end do
                            end do

                            do iy = 1, inE%ny
                                do iz = 2, inE%nz
                                    do ix = 2, inE%nx
                                        outE%y(ix, iy, iz) = (outE%y(ix, iy, iz) - &
                                        outE%y(ix, iy, iz-1)*model_operator%yYZ(iz, 1) - &
                                        outE%y(ix-1, iy, iz)*model_operator%yYX(ix, 1))* &
                                        self%Dilu%y(ix, iy, iz)
                                    end do
                                end do
                            end do

                            do iz = 1, inE%nz
                                do iy = 2, inE%ny
                                        do ix = 2, inE%nx
                                        outE%z(ix, iy, iz) = (outE%z(ix, iy, iz) - &
                                        outE%z(ix-1, iy, iz)*model_operator%zZX(ix, 1) - &
                                        outE%z(ix, iy-1, iz)*model_operator%zZY(iy, 1))* &
                                        self%Dilu%z(ix, iy, iz)
                                    end do
                                end do
                            end do
                        else
                            ! adjoint = .true. -- reverse mapping in to out
                            !     need to make sure that outE is zero on boundaries initially -- this is not
                            !        done explicitly in ModEM stable!
                            call outE%zeros()     !    let"s do explicitly -- but consider if necessary
                            !
                            do ix = 1, inE%nx
                                do iy = inE%ny, 2, -1
                                    do iz = inE%nz, 2, -1
                                        outE%x(ix, iy, iz) = (inE%x(ix, iy, iz) - &
                                        outE%x(ix, iy+1, iz)*model_operator%xXY(iy+1, 1) - &
                                        outE%x(ix, iy, iz+1)*model_operator%xXZ(iz+1, 1))* &
                                        conjg(self%Dilu%x(ix, iy, iz))
                                    end do
                                end do
                            end do
                            !
                            ! The coefficients for y are only for the interior nodes
                            do iy = 1, inE%ny
                                do ix = inE%nx, 2, -1
                                    do iz = inE%nz, 2, -1
                                        outE%y(ix, iy, iz) = (inE%y(ix, iy, iz) - &
                                        outE%y(ix, iy, iz+1)*model_operator%yYZ(iz+1, 1) - &
                                        outE%y(ix+1, iy, iz)*model_operator%yYX(ix+1, 1))* &
                                        conjg(self%Dilu%y(ix, iy, iz))
                                    end do
                                end do
                            end do

                            do iz = 1, inE%nz
                                do ix = inE%nx, 2, -1
                                    do iy = inE%ny, 2, -1
                                        outE%z(ix, iy, iz) = (inE%z(ix, iy, iz) - &
                                        outE%z(ix+1, iy, iz)*model_operator%zZX(ix+1, 1) - &
                                        outE%z(ix, iy+1, iz)*model_operator%zZY(iy+1, 1))* &
                                        conjg(self%Dilu%z(ix, iy, iz))
                                    end do
                                end do
                            end do

                            !     for adjoint to the division by volume elements last
                            call outE%divs(model_operator%Metric%Vedge)
                            !
                        end if 
                    class default
                         stop "LTSolvePreConditioner_MF_CC: Unclassified ModelOperator"
                    !
                    end select
                    !
            end select
        end select
        !
    end subroutine LTSolvePreConditioner_MF_CC
    
    !**
    ! Purpose: to solve the upper triangular system (or it"s adjoint);
    ! for the d-ilu pre-condtioner
    !*
    subroutine UTSolvePreConditioner_MF_CC( self, inE, outE, adjt )
        implicit none
        !
        class( PreConditioner_MF_CC_t ), intent( inout ) :: self
        class( cVector_t ), intent( in )                 :: inE
        class( cVector_t ), intent( inout )              :: outE
        logical, intent( in )                            :: adjt
        !
        integer :: ix, iy, iz
        !
        select type( inE )
        class is( cVector3D_SG_t )
            !
            select type( outE )
            class is( cVector3D_SG_t )
                !
                !    to be safe, zero out outR
                call outE%zeros()
                !
                ! Instantiate the ModelOperator object
                select type( model_operator => self%model_operator )
                    !
                    class is( ModelOperator_MF_t )
                        !
                        if (.NOT.adjt) then
                            ! for standard upper triangular solution
                            !
                            do ix = 1, inE%nx
                                do iz = inE%nz, 2, -1
                                    do iy = inE%ny, 2, -1
                                        outE%x(ix, iy, iz) = inE%x(ix, iy, iz) - &
                                        ( outE%x(ix, iy+1, iz)*model_operator%xXY(iy, 2) &
                                        + outE%x(ix, iy, iz+1)*model_operator%xXZ(iz, 2))* &
                                        self%Dilu%x(ix, iy, iz)
                                    end do
                                end do
                            end do
                            !
                            do iy = 1, inE%ny
                                do iz = inE%nz, 2, -1
                                    do ix = inE%nx, 2, -1
                                        outE%y(ix, iy, iz) = inE%y(ix, iy, iz) - &
                                        ( outE%y(ix, iy, iz+1)*model_operator%yYZ(iz, 2) &
                                        + outE%y(ix+1, iy, iz)*model_operator%yYX(ix, 2))* &
                                        self%Dilu%y(ix, iy, iz)
                                    end do
                                end do
                            end do
                            !
                            do iz = 1, inE%nz
                                do iy = inE%ny, 2, -1
                                    do ix = inE%nx, 2, -1
                                        outE%z(ix, iy, iz) = inE%z(ix, iy, iz) - &
                                        ( outE%z(ix+1, iy, iz)*model_operator%zZX(ix, 2) &
                                        + outE%z(ix, iy+1, iz)*model_operator%zZY(iy, 2))* &
                                        self%Dilu%z(ix, iy, iz)
                                    end do
                                end do
                            end do
                        else
                            ! adjoint = .true.
                            do ix = 1, inE%nx
                                do iz = 2, inE%nz
                                    do iy = 2, inE%ny
                                        outE%x(ix, iy, iz) = inE%x(ix, iy, iz) &
                                        - outE%x(ix, iy-1, iz)*model_operator%xXY(iy-1, 2) &
                                        * conjg(self%Dilu%x(ix,iy-1,iz))     &
                                        - outE%x(ix, iy, iz-1)*model_operator%xXZ(iz-1, 2) &
                                        * conjg(self%Dilu%x(ix, iy, iz-1))
                                    end do
                                end do
                            end do
                            !
                            do iy = 1, inE%ny
                                do iz = 2, inE%nz
                                    do ix = 2, inE%nx
                                        outE%y(ix, iy, iz) = inE%y(ix, iy, iz) &
                                        - outE%y(ix, iy, iz-1)*model_operator%yYZ(iz-1, 2) &
                                        * conjg(self%Dilu%y(ix,iy,iz-1)) &
                                        - outE%y(ix-1, iy, iz)*model_operator%yYX(ix-1, 2) &
                                        * conjg(self%Dilu%y(ix-1, iy, iz))
                                    end do
                                end do
                            end do
                            !
                            do iz = 1, inE%nz
                                do iy = 2, inE%ny
                                    do ix = 2, inE%nx
                                        outE%z(ix, iy, iz) = inE%z(ix, iy, iz) &
                                        - outE%z(ix-1, iy, iz)*model_operator%zZX(ix-1, 2) &
                                        * conjg(self%Dilu%z(ix-1,iy,iz)) &
                                        - outE%z(ix, iy-1, iz)*model_operator%zZY(iy-1, 2) &
                                        * conjg(self%Dilu%z(ix, iy-1, iz))
                                    end do
                                end do
                            end do
                        end if
                    class default
                         stop "setPreConditioner_MF_DC: Unclassified ModelOperator"
                    !
                end select
                !
            end select
        end select
        !
    end subroutine UTSolvePreConditioner_MF_CC
    !**
    ! this is dummy routine required by abstract preconditioner class
    !*
    subroutine LUSolvePreConditioner_MF_CC( self, inPhi, outPhi )
        implicit none
        !
        class( PreConditioner_MF_CC_t ), intent( inout ) :: self
        class( cScalar_t ), intent( in )                 :: inPhi
        class( cScalar_t ), intent( inout )              :: outPhi
        !
        STOP "ERROR: LUsolve is not coded for this pre-conditioner class"
        !
    end subroutine LUSolvePreConditioner_MF_CC
    !
end module PreConditioner_MF_CC

