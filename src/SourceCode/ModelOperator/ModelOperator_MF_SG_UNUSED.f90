    !
    !> maxwell computes the finite difference equation in complex vectors
    !> for del X del X E +/- i*omega*mu*conductivity*E in unsymmetrical form.
    !> Note that the difference equation is only for interior edges. However,
    !> it does use the contribution from the boundary edges. The coefficients
    !> are  calculated in CurlcurleSetUp. Remember, in the operators that are
    !> used in iterative fashion, output is always initialized outside
    !
    !> Obs.: input and output electrical field as complex vector
    !
    subroutine maxwell_ModelOperator_MF_SG( self, in_e, adjoint, out_e )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: in_e
        logical, intent( in ) :: adjoint
        class( Vector_t ), intent( inout ) :: out_e
        !
        integer :: diag_sign, ix, iy, iz
        complex( kind=prec ), allocatable, dimension(:,:,:) :: in_e_x, in_e_y, in_e_z
        complex( kind=prec ), allocatable, dimension(:,:,:) :: Adiag_e_v, out_e_v
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "maxwell_ModelOperator_MF_SG > in_e not allocated" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "maxwell_ModelOperator_MF_SG > out_e not allocated" )
        endif
        !
        !> Check whether the bounds are the same
        if( ( in_e%nx == out_e%nx ) .AND. &
            ( in_e%ny == out_e%ny ) .AND. &
            ( in_e%nz == out_e%nz ) ) then
            !
            if( in_e%grid_type == out_e%grid_type ) then
                !
                if( adjoint ) then
                    diag_sign = -1 * ISIGN
                else
                    diag_sign = ISIGN
                end if
                !
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ix,iy,iz)
                !
                !> Apply difference equation to compute Ex (only on interior nodes)
                !> the diagonal nodes have the imaginary component added
!!$OMP DO SCHEDULE(STATIC)
                !
                in_e_x = in_e%getX()
                in_e_y = in_e%getY()
                in_e_z = in_e%getZ()
                !
                !> Ex
                Adiag_e_v = self%Adiag%getX()
                out_e_v = out_e%getX()
                !
                do iz = 2, in_e%nz
                    do iy = 2, in_e%ny
                        do ix = 1, in_e%nx
                            out_e_v(ix,iy,iz) = self%xY(ix,iy)*(in_e_y(ix+1,iy,iz)-&
                            in_e_y(ix,iy,iz)-in_e_y(ix+1,iy-1,iz)&
                            +in_e_y(ix,iy-1,iz))+&
                            self%xZ(ix,iz)*(in_e_z(ix+1,iy,iz)-in_e_z(ix,iy,iz)&
                            -in_e_z(ix+1,iy,iz-1)+in_e_z(ix,iy,iz-1))+&
                            self%xXY(iy,2)*in_e_x(ix,iy+1,iz)+&
                            self%xXY(iy,1)*in_e_x(ix,iy-1,iz)+&
                            self%xXZ(iz,2)*in_e_x(ix,iy,iz+1)+&
                            self%xXZ(iz,1)*in_e_x(ix,iy,iz-1)+&
                            (self%xXO(iy,iz) + diag_sign * Adiag_e_v(ix,iy,iz))*in_e_x(ix,iy,iz)
                        enddo
                    enddo
                enddo
                !
!!$OMP END DO NOWAIT
                !
                call out_e%setX( out_e_v )
                !
                !> Ey
                Adiag_e_v = self%Adiag%getY()
                out_e_v = out_e%getY()
                !
                !> Apply difference equation to compute Ey (only on interior nodes)
                !> the diagonal nodes have the imaginary component added
!!$OMP DO SCHEDULE(STATIC)
                !
                do iz = 2, in_e%nz
                    do iy = 1, in_e%ny
                        do ix = 2, in_e%nx
                        out_e_v(ix,iy,iz) = self%yZ(iy,iz)*(in_e_z(ix,iy+1,iz)-&
                        in_e_z(ix,iy,iz)-in_e_z(ix,iy+1,iz-1)+in_e_z(ix,iy,iz-1))&
                        +self%yX(ix,iy)*(in_e_x(ix,iy+1,iz)-in_e_x(ix,iy,iz)&
                        -in_e_x(ix-1,iy+1,iz)+in_e_x(ix-1,iy,iz))+&
                        self%yYZ(iz,2)*in_e_y(ix,iy,iz+1)+&
                        self%yYZ(iz,1)*in_e_y(ix,iy,iz-1)+&
                        self%yYX(ix,2)*in_e_y(ix+1,iy,iz)+&
                        self%yYX(ix,1)*in_e_y(ix-1,iy,iz)+&
                        (self%yYO(ix,iz)+diag_sign*Adiag_e_v(ix,iy,iz))*in_e_y(ix,iy,iz)
                        enddo
                    enddo
                enddo
                !
!!$OMP END DO NOWAIT
                !
                !> Apply difference equation to compute Ey (only on interior nodes)
                !> the diagonal nodes have the imaginary component added
!!$OMP DO SCHEDULE(STATIC)
                !
                call out_e%setY( out_e_v )
                !
                !> Ez
                Adiag_e_v = self%Adiag%getZ()
                out_e_v = out_e%getZ()
                !
                do iz = 1, in_e%nz
                    do iy = 2, in_e%ny
                        do ix = 2, in_e%nx
                            out_e_v(ix,iy,iz) = self%zX(ix,iz)*(in_e_x(ix,iy,iz+1)-&
                            in_e_x(ix,iy,iz)-in_e_x(ix-1,iy,iz+1)+in_e_x(ix-1,iy,iz))&
                            +self%zY(iy,iz)*(in_e_y(ix,iy,iz+1)-in_e_y(ix,iy,iz)&
                            -in_e_y(ix,iy-1,iz+1)+in_e_y(ix,iy-1,iz))+&
                            self%zZX(ix,2)*in_e_z(ix+1,iy,iz)+&
                            self%zZX(ix,1)*in_e_z(ix-1,iy,iz)+&
                            self%zZY(iy,2)*in_e_z(ix,iy+1,iz)+&
                            self%zZY(iy,1)*in_e_z(ix,iy-1,iz)+&
                            (self%zZO(ix,iy)+diag_sign*Adiag_e_v(ix,iy,iz))*in_e_z(ix,iy,iz)
                        enddo
                    enddo
                enddo
                !
                call out_e%setZ( out_e_v )
                !
!!$OMP END DO NOWAIT
                !
!!$OMP END PARALLEL
            !
            else
                call errStop( "maxwell_ModelOperator_MF_SG > not compatible usage for existing data types" )
            end if
        !
        else
            call errStop( "maxwell_ModelOperator_MF_SG > Error-complex vectors are not of same size" )
        end if
        !
    end subroutine maxwell_ModelOperator_MF_SG
    !
    !> Purpose: to solve the lower triangular system (or it's adjoint);
    !> for the d-ilu pre-condtioner.
    !
    subroutine M1solve_ModelOperator_MF_SG( self, in_e, adjoint, out_e )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: in_e
        logical, intent( in ) :: adjoint
        class( Vector_t ), intent( inout ) :: out_e
        !
        integer :: ix, iy, iz
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "M1solve_ModelOperator_MF_SG > in_e not allocated" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "M1solve_ModelOperator_MF_SG > out_e not allocated" )
        endif
        !
        !> Check whether all the vector nodes are of the same size
        if( ( in_e%nx == out_e%nx ) .AND. &
            ( in_e%ny == out_e%ny ) .AND. &
            ( in_e%nz == out_e%nz ) ) then
            !
            if( in_e%grid_type == out_e%grid_type ) then
                !
                if( .NOT. adjoint ) then
                    ! adjoint = .FALSE.
                    call diagDiv( in_e, V_E, out_e )
                    !
                    ! ... note that we only parallelize the outer loops
!!$OMP PARALLEL DEFAULT(SHARED)
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
                    !
                    do ix = 1, in_e%nx
!!$OMP ORDERED
                        do iz = 2, in_e%nz
                            do iy = 2, in_e%ny
                                !
                                out_e%x(ix, iy, iz) = (out_e%x(ix, iy, iz) - &
                                out_e%x(ix, iy-1, iz)*xXY(iy, 1) - &
                                out_e%x(ix, iy, iz-1)*xXZ(iz, 1))* &
                                Dilu%x(ix, iy, iz)
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
!!$OMP END DO
                    !
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
                    do iy = 1, in_e%ny
                        !$OMP ORDERED
                        do iz = 2, in_e%nz
                            do ix = 2, in_e%nx
                                !
                                out_e%y(ix, iy, iz) = (out_e%y(ix, iy, iz) - &
                                out_e%y(ix, iy, iz-1)*yYZ(iz, 1) - &
                                out_e%y(ix-1, iy, iz)*yYX(ix, 1))* &
                                Dilu%y(ix, iy, iz)
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
!!$OMP END DO
                    !
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
                    do iz = 1, in_e%nz
                        !$OMP ORDERED
                        do iy = 2, in_e%ny
                            do ix = 2, in_e%nx
                                !
                                out_e%z(ix, iy, iz) = (out_e%z(ix, iy, iz) - &
                                out_e%z(ix-1, iy, iz)*zZX(ix, 1) - &
                                out_e%z(ix, iy-1, iz)*zZY(iy, 1))* &
                                Dilu%z(ix, iy, iz)
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
!!$OMP END DO
!!$OMP END PARALLEL
                    ! adjoint = .TRUE.
                else
                    !
                    ! ... note that we only parallelize the outer loops
!!$OMP PARALLEL DEFAULT(SHARED)
                    ! the coefficients for x are only for the interior nodes
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
                    do ix = 1, in_e%nx
!!$OMP ORDERED
                        do iy = in_e%ny, 2, -1
                            do iz = in_e%nz, 2, -1
                                !
                                out_e%x(ix, iy, iz) = (in_e%x(ix, iy, iz) - &
                                out_e%x(ix, iy+1, iz)*xXY(iy+1, 1) - &
                                out_e%x(ix, iy, iz+1)*xXZ(iz+1, 1))* &
                                conjg(Dilu%x(ix, iy, iz))
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
!!$OMP END DO
                    !
                    ! the coefficients for y are only for the interior nodes
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
                    do iy = 1, in_e%ny
!!$OMP ORDERED
                        do ix = in_e%nx, 2, -1
                            do iz = in_e%nz, 2, -1
                                !
                                out_e%y(ix, iy, iz) = (in_e%y(ix, iy, iz) - &
                                out_e%y(ix, iy, iz+1)*yYZ(iz+1, 1) - &
                                out_e%y(ix+1, iy, iz)*yYX(ix+1, 1))* &
                                conjg(Dilu%y(ix, iy, iz))
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
!!$OMP END DO
                    !
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
                    do iz = 1, in_e%nz
!!$OMP ORDERED
                        do ix = in_e%nx, 2, -1
                            do iy = in_e%ny, 2, -1
                                !
                                out_e%z(ix, iy, iz) = (in_e%z(ix, iy, iz) - &
                                out_e%z(ix+1, iy, iz)*zZX(ix+1, 1) - &
                                out_e%z(ix, iy+1, iz)*zZY(iy+1, 1))* &
                                conjg(Dilu%z(ix, iy, iz))
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
!!$OMP END DO
!!$OMP END PARALLEL
                    !
                    call diagDiv(out_e, V_E, out_e)
                    !
                endif
                !
            else
                call errStop( "M1solve_ModelOperator_MF_SG > not compatible usage for M1solve" )
            end if
            !
        else
            call errStop( "M1solve_ModelOperator_MF_SG > lower triangular: vectors not same size" )
        end if
        !
    end subroutine M1solve_ModelOperator_MF_SG
    !
    !> Purpose: to solve the upper triangular system (or it's adjoint);
    !> for the d-ilu pre-condtioner
    !
    subroutine M2solve_ModelOperator_MF_SG( self, in_e, adjoint, out_e )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: in_e
        logical, intent( in ) :: adjoint
        class( Vector_t ), intent( inout ) :: out_e
        !
        integer :: ix, iy, iz
        !
        if( .NOT. in_e%is_allocated ) then
            call errStop( "M2solve_ModelOperator_MF_SG > in_e not allocated" )
        endif
        !
        if( .NOT. out_e%is_allocated ) then
            call errStop( "M2solve_ModelOperator_MF_SG > out_e not allocated" )
        endif
        !
        !> Check whether all the vector nodes are of the same size
        if( ( in_e%nx == out_e%nx ) .AND. &
            ( in_e%ny == out_e%ny ) .AND. &
            ( in_e%nz == out_e%nz ) ) then
            !
            if( in_e%grid_type == out_e%grid_type ) then
                !
                ! adjoint = .FALSE.
                if( .NOT. adjoint ) then
                    ! for standard upper triangular solution
                    ! ... note that we only parallelize the outer loops
!!$OMP PARALLEL DEFAULT(SHARED)
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
                    do ix = 1, in_e%nx
!!$OMP ORDERED
                        do iz = in_e%nz, 2, -1
                            do iy = in_e%ny, 2, -1
                                !
                                out_e%x(ix, iy, iz) = in_e%x(ix, iy, iz) - &
                                ( out_e%x(ix, iy+1, iz)*xXY(iy, 2) &
                                + out_e%x(ix, iy, iz+1)*xXZ(iz, 2))* &
                                Dilu%x(ix, iy, iz)
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
!!$OMP END DO
                    !
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
                    do iy = 1, in_e%ny
!!$OMP ORDERED
                        do iz = in_e%nz, 2, -1
                            do ix = in_e%nx, 2, -1
                                !
                                out_e%y(ix, iy, iz) = in_e%y(ix, iy, iz) - &
                                ( out_e%y(ix, iy, iz+1)*yYZ(iz, 2) &
                                + out_e%y(ix+1, iy, iz)*yYX(ix, 2))* &
                                Dilu%y(ix, iy, iz)
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
!!$OMP END DO
                    !
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
                    do iz = 1, in_e%nz
!!$OMP ORDERED
                        do iy = in_e%ny, 2, -1
                            do ix = in_e%nx, 2, -1
                                !
                                out_e%z(ix, iy, iz) = in_e%z(ix, iy, iz) - &
                                ( out_e%z(ix+1, iy, iz)*zZX(ix, 2) &
                                + out_e%z(ix, iy+1, iz)*zZY(iy, 2))* &
                                Dilu%z(ix, iy, iz)
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
                    !
!!$OMP END DO
!!$OMP END PARALLEL
                    ! adjoint = .true.
                    !
                else
                    !
                    ! ... note that we only parallelize the outer loops
!!$OMP PARALLEL DEFAULT(SHARED)
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
                    do ix = 1, in_e%nx
!!$OMP ORDERED
                        do iz = 2, in_e%nz
                            do iy = 2, in_e%ny
                                !
                                out_e%x(ix, iy, iz) = in_e%x(ix, iy, iz) &
                                - out_e%x(ix, iy-1, iz)*xXY(iy-1, 2) &
                                * conjg(Dilu%x(ix,iy-1,iz))   &
                                - out_e%x(ix, iy, iz-1)*xXZ(iz-1, 2) &
                                * conjg(Dilu%x(ix, iy, iz-1))
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
!!$OMP END DO
                    !
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
                    do iy = 1, in_e%ny
!!$OMP ORDERED
                        do iz = 2, in_e%nz
                            do ix = 2, in_e%nx
                                !
                                out_e%y(ix, iy, iz) = in_e%y(ix, iy, iz) &
                                - out_e%y(ix, iy, iz-1)*yYZ(iz-1, 2) &
                                * conjg(Dilu%y(ix,iy,iz-1)) &
                                - out_e%y(ix-1, iy, iz)*yYX(ix-1, 2) &
                                * conjg(Dilu%y(ix-1, iy, iz))
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
!!$OMP END DO
                    !
!!$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
                    do iz = 1, in_e%nz
!!$OMP ORDERED
                        do iy = 2, in_e%ny
                            do ix = 2, in_e%nx
                                !
                                out_e%z(ix, iy, iz) = in_e%z(ix, iy, iz) &
                                - out_e%z(ix-1, iy, iz)*zZX(ix-1, 2) &
                                * conjg(Dilu%z(ix-1,iy,iz)) &
                                - out_e%z(ix, iy-1, iz)*zZY(iy-1, 2) &
                                * conjg(Dilu%z(ix, iy-1, iz))
                                !
                            enddo
                        enddo
!!$OMP END ORDERED
                    enddo
!!$OMP END DO
!!$OMP END PARALLEL
                    !
                endif
                !
            else
                call errStop( "M1solve_ModelOperator_MF_SG > not compatible usage for M1solve" )
            end if
            !
        else
            call errStop( "M1solve_ModelOperator_MF_SG > lower triangular: vectors not same size" )
        end if
        !
    end subroutine M2solve ! M2solve
    !
    !> Adiag sets up the diagonal nodes with the imaginary part added to it
    !> SetUp routines do basic calculations (maybe one time deal or more than one)
    !
    subroutine AdiagSetUp_ModelOperator_MF_SG( self, omega_in )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: omega_in
        !
        complex( kind=prec ), allocatable, dimension(:,:,:) :: Adiag_e_v
        integer :: ix, iy, iz
        !
        if( .NOT. self%Adiag%is_allocated ) then
            call errStop( "AdiagSetUp_ModelOperator_MF_SG > Adiag not allocated" )
        endif
        !
        !> Ex
        Adiag_e_v = self%Adiag%getX()
        !
        do ix = 1, self%metric%grid%nx
            Adiag_e_v(ix,:,:) = ONE_I * omega_in * MU_0 * self%sigma_e%x(ix,:,:)
        enddo
        !
        call self%Adiag%setX( Adiag_e_v )
        !
        !> Ey
        Adiag_e_v = self%Adiag%getY()
        !
        do iy = 1, self%metric%grid%ny
            Adiag_e_v(:,iy,:) = ONE_I * omega_in * MU_0 * self%sigma_e%y(:,iy,:)
        enddo
        !
        call self%Adiag%setY( Adiag_e_v )
        !
        !> Ez
        Adiag_e_v = self%Adiag%getZ()
        !
        do iz = 1, self%metric%grid%nz
            Adiag_e_v(:,:,iz) = ONE_I * omega_in * MU_0 * self%sigma_e%z(:,:,iz)
        enddo
        !
        call self%Adiag%setZ( Adiag_e_v )
        !
    end subroutine AdiagSetUp_ModelOperator_MF_SG
    !