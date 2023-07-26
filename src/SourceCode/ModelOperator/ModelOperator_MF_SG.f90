!
!> Derived class to define a ModelOperator_MF_SG
!
module ModelOperator_MF_SG
    !
    use ModelOperator
    use MetricElements_CSG
    !
    type, extends( ModelOperator_t ) :: ModelOperator_MF_SG_t
        !
        real( kind=prec ), allocatable, dimension(:,:) :: xXY, xXZ
        real( kind=prec ), allocatable, dimension(:,:) :: xY, xZ
        real( kind=prec ), allocatable, dimension(:,:) :: xXO
        real( kind=prec ), allocatable, dimension(:,:) :: yYX, yYZ
        real( kind=prec ), allocatable, dimension(:,:) :: yX, yZ
        real( kind=prec ), allocatable, dimension(:,:) :: yYO
        real( kind=prec ), allocatable, dimension(:,:) :: zZX, zZY
        real( kind=prec ), allocatable, dimension(:,:) :: zX, zY
        real( kind=prec ), allocatable, dimension(:,:) :: zZO
        !
        type( rVector3D_SG_t ) :: sigma_e
        !
        type( rVector3D_SG_t ) :: db1, db2
        !
        type( rScalar3D_SG_t ) :: c
        !
        contains
            !
            final :: ModelOperator_MF_SG_dtor
            !
            !> Setup
            procedure, public :: setEquations => setEquations_ModelOperator_MF_SG
            procedure, public :: setCond => setCond_ModelOperator_MF_SG
            !
            !procedure, public :: divCorInit => divCorInit_ModelOperator_MF_SG
            procedure, public :: divCorSetUp => divCorSetUp_ModelOperator_MF_SG
            !
            !> Operations
            procedure, public :: amult => amult_ModelOperator_MF_SG
            procedure, public :: multAib => multAib_ModelOperator_MF_SG
            !
            procedure, public :: div => div_ModelOperator_MF_SG
            procedure, public :: divC => divC_ModelOperator_MF_SG
            procedure, public :: divCGrad => divCGrad_ModelOperator_MF_SG
            !
            procedure, public :: grad => grad_ModelOperator_MF_SG
            !
            !> Alloc/Dealloc
            procedure :: create => create_ModelOperator_MF_SG 
            procedure :: alloc => allocate_ModelOperator_MF_SG
            procedure :: dealloc => deallocate_ModelOperator_MF_SG
            !
            !> Miscellaneous
            procedure, public :: print => print_ModelOperator_MF_SG
            !
    end type ModelOperator_MF_SG_t
    !
    interface ModelOperator_MF_SG_t
        module procedure ModelOperator_MF_SG_ctor
    end interface ModelOperator_MF_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelOperator_MF_SG_ctor( grid ) result( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        !
        type( ModelOperator_MF_SG_t ) :: self
        !
        !write( *, * ) "Constructor ModelOperator_MF_SG"
        !
        call self%baseInit
        !
        !> Instantiation of the specific object MetricElements
        allocate( self%metric, source = MetricElements_CSG_t( grid ) )
        !
        call self%create( grid )
        !
    end function ModelOperator_MF_SG_ctor
    !
    !> ModelOperator_MF_SG destructor
    subroutine ModelOperator_MF_SG_dtor( self )
        implicit none
        !
        type( ModelOperator_MF_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor ModelOperator_MF_SG_t"
        !
        call self%baseDealloc()
        !
    end subroutine ModelOperator_MF_SG_dtor
    !
    !> No subroutine briefing
    !
    subroutine setEquations_ModelOperator_MF_SG( self )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz 
        !
        do iy = 2, self%metric%grid%ny
            self%xXY(iy, 2) = -1.0 / (self%metric%grid%del_y(iy) * self%metric%grid%dy(iy))
            self%xXY(iy, 1) = -1.0 / (self%metric%grid%del_y(iy) * self%metric%grid%dy(iy-1))
        enddo
        !
        do iz = 2, self%metric%grid%nz
            self%xXZ(iz, 2) = -1.0 / (self%metric%grid%del_z(iz) * self%metric%grid%dz(iz))
            self%xXZ(iz, 1) = -1.0 / (self%metric%grid%del_z(iz) * self%metric%grid%dz(iz-1))
        enddo
        !
        do iy = 2, self%metric%grid%ny
            do iz = 2, self%metric%grid%nz
                self%xXO(iy, iz) = -(self%xXY(iy,1) + self%xXY(iy,2) + &
                self%xXZ(iz,1) + self%xXZ(iz,2))
            enddo
        enddo
        !
        do ix = 1, self%metric%grid%nx
            do iy = 2, self%metric%grid%ny
                self%xY(ix, iy) = 1.0 / (self%metric%grid%del_y(iy)*self%metric%grid%dx(ix))
            enddo
        enddo
        !
        do ix = 1, self%metric%grid%nx
            do iz = 2, self%metric%grid%nz
                self%xZ(ix, iz) = 1.0 / (self%metric%grid%del_z(iz)*self%metric%grid%dx(ix))
            enddo
        enddo
        !
        do iz = 2, self%metric%grid%nz
            self%yYZ(iz, 2) = -1.0 / (self%metric%grid%del_z(iz)*self%metric%grid%dz(iz))
            self%yYZ(iz, 1) = -1.0 / (self%metric%grid%del_z(iz)*self%metric%grid%dz(iz-1))
        enddo
        !
        do ix = 2, self%metric%grid%nx
            self%yYX(ix, 2) = -1.0 / (self%metric%grid%del_x(ix)*self%metric%grid%dx(ix))
            self%yYX(ix, 1) = -1.0 / (self%metric%grid%del_x(ix)*self%metric%grid%dx(ix-1))
        enddo
        !
        do ix = 2, self%metric%grid%nx
            do iz = 2, self%metric%grid%nz
                self%yYO(ix, iz) = -(self%yYX(ix,1) + self%yYX(ix,2) + &
                self%yYZ(iz,1) + self%yYZ(iz,2))
            enddo
        enddo
        !
        do iy = 1, self%metric%grid%ny
            do iz = 2, self%metric%grid%nz
                self%yZ(iy, iz) = 1.0 / (self%metric%grid%del_z(iz)*self%metric%grid%dy(iy))
            enddo
        enddo
        !
        do ix = 2, self%metric%grid%nx
            do iy = 1, self%metric%grid%ny
                self%yX(ix, iy) = 1.0 / (self%metric%grid%del_x(ix)*self%metric%grid%dy(iy))
            enddo
        enddo
        !
        do ix = 2, self%metric%grid%nx
            self%zZX(ix, 2) = -1.0 / (self%metric%grid%del_x(ix)*self%metric%grid%dx(ix))
            self%zZX(ix, 1) = -1.0 / (self%metric%grid%del_x(ix)*self%metric%grid%dx(ix-1))
        enddo
        !
        do iy = 2, self%metric%grid%ny
            self%zZY(iy, 2) = -1.0 / (self%metric%grid%del_y(iy)*self%metric%grid%dy(iy))
            self%zZY(iy, 1) = -1.0 / (self%metric%grid%del_y(iy)*self%metric%grid%dy(iy-1))
        enddo
        !
        do ix = 2, self%metric%grid%nx
            do iy = 2, self%metric%grid%ny
                self%zZO(ix, iy) = -(self%zZX(ix,1) + self%zZX(ix,2) + &
                self%zZY(iy,1) + self%zZY(iy,2))
            enddo
        enddo
        !
        do ix = 2, self%metric%grid%nx
             do iz = 1, self%metric%grid%nz
                 self%zX(ix, iz) = 1.0 / (self%metric%grid%del_x(ix)*self%metric%grid%dz(iz))
             enddo
        enddo
        !
        do iy = 2, self%metric%grid%ny
            do iz = 1, self%metric%grid%nz
                self%zY(iy, iz) = 1.0 / (self%metric%grid%del_y(iy)*self%metric%grid%dz(iz))
            enddo
        enddo
        !
        self%eqset = .TRUE.
        !
    end subroutine setEquations_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    subroutine setCond_ModelOperator_MF_SG( self, sigma, omega )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( inout ) :: sigma
        real( kind=prec ), intent( in ), optional :: omega
        !
        call sigma%PDEmapping( self%sigma_e )
        !
    end subroutine setCond_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    subroutine divCorSetUp_ModelOperator_MF_SG( self )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent(inout) :: self
        !
        integer :: ix, iy, iz
        !
        do iz = 2, self%metric%grid%nz
            do iy = 2, self%metric%grid%ny
                do ix = 2, self%metric%grid%nx
                    !
                    self%db1%x(ix, iy, iz) = self%sigma_e%x(ix - 1, iy, iz)/ &
                    (self%metric%grid%dx(ix - 1)*self%metric%grid%del_x(ix))
                    !
                    self%db2%x(ix, iy, iz) = self%sigma_e%x(ix, iy, iz)/ &
                    (self%metric%grid%dx(ix)*self%metric%grid%del_x(ix))
                    !
                    self%db1%y(ix, iy, iz) = self%sigma_e%y(ix, iy - 1, iz)/ &
                    (self%metric%grid%dy(iy - 1)*self%metric%grid%del_y(iy))
                    !
                    self%db2%y(ix, iy, iz) = self%sigma_e%y(ix, iy, iz)/ &
                    (self%metric%grid%dy(iy)*self%metric%grid%del_y(iy))
                    !
                    self%db1%z(ix, iy, iz) = self%sigma_e%z(ix, iy, iz - 1)/ &
                    (self%metric%grid%dz(iz - 1)*self%metric%grid%del_z(iz))
                    !
                    self%db2%z(ix, iy, iz) = self%sigma_e%z(ix, iy, iz)/ &
                    (self%metric%grid%dz(iz)*self%metric%grid%del_z(iz))
                    !
                    self%c%v(ix, iy, iz) = -(self%db1%x(ix, iy, iz) + &
                    self%db2%x(ix, iy, iz) + &
                    self%db1%y(ix, iy, iz) + &
                    self%db2%y(ix, iy, iz) + &
                    self%db1%z(ix, iy, iz) + &
                    self%db2%z(ix, iy, iz) )
                    !
                enddo
            enddo
        enddo
        !
        select type( v_node => self%metric%v_node )
            !
            class is( rScalar3D_SG_t )
                !
                do iz = 2, self%metric%grid%nz
                    do iy = 2, self%metric%grid%ny
                        do ix = 2, self%metric%grid%nx
                            !
                            self%db1%x(ix, iy, iz) = self%db1%x(ix, iy, iz) * &
                            v_node%v(ix,iy,iz)
                            !
                            self%db1%y(ix, iy, iz) = self%db1%y(ix, iy, iz) * &
                            v_node%v(ix,iy,iz)
                            !
                            self%db1%z(ix, iy, iz) = self%db1%z(ix, iy, iz) * &
                            v_node%v(ix,iy,iz)
                            !
                            self%db2%x(ix, iy, iz) = self%db2%x(ix, iy, iz) * &
                            v_node%v(ix,iy,iz)
                            !
                            self%db2%y(ix, iy, iz) = self%db2%y(ix, iy, iz) * &
                            v_node%v(ix,iy,iz)
                            !
                            self%db2%z(ix, iy, iz) = self%db2%z(ix, iy, iz) * &
                            v_node%v(ix,iy,iz)
                            !
                        enddo
                    enddo
                enddo
                !
            class default
                call errStop( "divCorSetUp_ModelOperator_MF_SG > undefined grid" )
                !
        end select
        !
        call self%c%mult( self%metric%v_node )
        !
        self%db1%x(2, :, :) = R_ZERO
        self%db1%y(:, 2, :) = R_ZERO
        self%db1%z(:, :, 2) = R_ZERO
        !
        self%db2%x(self%metric%grid%nx, :, :) = R_ZERO
        self%db2%y(:, self%metric%grid%ny, :) = R_ZERO
        self%db2%z(:, :, self%metric%grid%nz) = R_ZERO
        !
    end subroutine divCorSetUp_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    !
    subroutine amult_ModelOperator_MF_SG( self, omega, in_e, out_e, p_adjoint )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( in ) :: self
        real( kind=prec ), intent( in ), optional :: omega
        class( Vector_t ), intent( inout ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        logical, intent( in ), optional :: p_adjoint
        !
        integer :: ix, iy, iz
        complex( kind=prec ) :: cvalue
        logical :: adjoint
        !
        if( present( p_adjoint ) ) then
            adjoint = p_adjoint
        else
            adjoint = .FALSE.
        endif
        !
        if( adjoint ) then
            cvalue = -ONE_I * omega * isign * mu_0
        else
            cvalue = ONE_I * omega * isign * mu_0
        endif
        !
        select type( in_e )
            !
            class is( cVector3D_SG_t )
                !
                if( .NOT. out_e%is_allocated ) then
                    call errStop( "amult_ModelOperator_MF_SG > output vector out_e not allocated" )
                endif
                !
                select type( out_e )
                    !
                    class is( cVector3D_SG_t )
                        !
                        call out_e%zeros
                        !
                        do iz = 2, in_e%nz
                            do iy = 2, in_e%ny
                                do ix = 1, in_e%nx
                                    out_e%x(ix, iy, iz) = self%xY(ix, iy)*(in_e%y(ix + 1, iy, iz) - &
                                    in_e%y(ix, iy, iz) - in_e%y(ix + 1, iy - 1, iz) + &
                                    in_e%y(ix, iy - 1, iz)) + &
                                    self%xZ(ix, iz) * (in_e%z(ix + 1, iy, iz) - in_e%z(ix, iy, iz) - &
                                    in_e%z(ix + 1, iy, iz - 1) + in_e%z(ix, iy, iz - 1)) + &
                                    self%xXY(iy, 2) * in_e%x(ix, iy + 1, iz) + &
                                    self%xXY(iy, 1) * in_e%x(ix, iy - 1, iz) + &
                                    self%xXZ(iz, 2) * in_e%x(ix, iy, iz + 1) + &
                                    self%xXZ(iz, 1) * in_e%x(ix, iy, iz - 1) + &
                                    (self%xXO(iy, iz)+cvalue*self%sigma_e%x(ix,iy,iz)) * in_e%x(ix, iy, iz)
                                enddo
                            enddo
                        enddo
                        !
                        do iz = 2, in_e%nz
                            do iy = 1, in_e%ny
                                do ix = 2, in_e%nx
                                    out_e%y(ix, iy, iz) = self%yZ(iy, iz) * (in_e%z(ix, iy + 1, iz) - &
                                    in_e%z(ix, iy, iz) - in_e%z(ix, iy + 1, iz - 1) + in_e%z(ix, iy, iz - 1)) + &
                                    self%yX(ix, iy) * (in_e%x(ix, iy + 1, iz) - in_e%x(ix, iy, iz) - &
                                    in_e%x(ix - 1, iy + 1, iz) + in_e%x(ix - 1, iy, iz)) + &
                                    self%yYZ(iz, 2) * in_e%y(ix, iy, iz + 1) + &
                                    self%yYZ(iz, 1) * in_e%y(ix, iy, iz - 1) + &
                                    self%yYX(ix, 2) * in_e%y(ix + 1, iy, iz) + &
                                    self%yYX(ix, 1) * in_e%y(ix - 1, iy, iz) + &
                                    (self%yYO(ix, iz)+cvalue*self%sigma_e%y(ix,iy,iz)) * in_e%y(ix, iy, iz)
                                enddo
                            enddo
                        enddo
                        !
                        do iz = 1, in_e%nz
                            do iy = 2, in_e%ny
                                do ix = 2, in_e%nx
                                    out_e%z(ix, iy, iz) = self%zX(ix, iz) * (in_e%x(ix, iy, iz + 1) - &
                                    in_e%x(ix, iy, iz) - in_e%x(ix - 1, iy, iz + 1) + in_e%x(ix - 1, iy, iz)) + &
                                    self%zY(iy,iz) * (in_e%y(ix, iy, iz + 1) - in_e%y(ix, iy, iz) - &
                                    in_e%y(ix, iy - 1, iz + 1) + in_e%y(ix, iy - 1, iz)) + &
                                    self%zZX(ix, 2) * in_e%z(ix + 1, iy, iz) + &
                                    self%zZX(ix, 1) * in_e%z(ix - 1, iy, iz) + &
                                    self%zZY(iy, 2) * in_e%z(ix, iy + 1, iz) + &
                                    self%zZY(iy, 1) * in_e%z(ix, iy - 1, iz) + &
                                    (self%zZO(ix, iy)+cvalue*self%sigma_e%z(ix,iy,iz)) * in_e%z(ix, iy, iz)
                                enddo
                            enddo
                        enddo
                        !
                        call out_e%mult( self%metric%v_edge )
                        !
                    class default
                        call errStop( "amult_ModelOperator_MF_SG > Undefined out_e." )
                        !
                end select
                !
            class default
                call errStop( "amult_ModelOperator_MF_SG > Undefined in_e." )
                !
        end select
        !
    end subroutine amult_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    !
    subroutine multAib_ModelOperator_MF_SG( self, in_e, out_e )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        !
        real( kind=prec ) :: omega
        !
        if(.NOT. out_e%is_allocated) then
            call errStop( "multAib_ModelOperator_MF_SG > out_e not allocated" )
        endif
        !
        omega = R_ZERO
        !
        call self%amult( omega, in_e, out_e ) 
        !
    end subroutine multAib_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    !
    subroutine div_ModelOperator_MF_SG( self, in_e, out_phi )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: in_e
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        integer :: ix, iy, iz
        !
        select type( out_phi )
            !
            class is( cScalar3D_SG_t )
                !
                if( .NOT. out_phi%is_allocated ) then
                    call errStop( "div_ModelOperator_MF_SG > Output cScalar object not allocated" )
                endif
                !
                select type( in_e )
                    !
                    class is( cVector3D_SG_t )
                        !
                        call out_phi%zeros
                        !
                        do ix = 2, out_phi%nx
                            do iy = 2, out_phi%ny
                                do iz = 2, out_phi%grid%nz
                                    out_phi%v(ix, iy, iz) = &
                                    (in_e%x(ix, iy, iz) - in_e%x(ix - 1, iy, iz)) * &
                                    in_e%grid%del_x_inv(ix) + &
                                    (in_e%y(ix, iy, iz) - in_e%y(ix, iy - 1, iz)) * &
                                    in_e%grid%del_y_inv(iy) + &
                                    (in_e%z(ix, iy, iz) - in_e%z(ix, iy, iz - 1)) * &
                                    in_e%grid%del_z_inv(iz)
                                enddo
                            enddo
                        enddo
                        !
                    class default
                        call errStop( "div_ModelOperator_MF_SG> in_e type unknown" )
                end select
                !
            class default
                call errStop( "div_ModelOperator_MF_SG > out_phi type unknown" )
                !
        end select
        !
    end subroutine div_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    !
    subroutine divC_ModelOperator_MF_SG( self, in_e, out_phi )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: in_e
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        integer :: ix, iy, iz
        !
        select type( out_phi )
            !
            class is( cScalar3D_SG_t )
                !
                if( .NOT. out_phi%is_allocated ) then
                    call errStop( "divC_ModelOperator_MF_SG > out_phi not allocated" )
                endif
                !
                select type( in_e )
                    !
                    class is ( cVector3D_SG_t )
                        !
                        call out_phi%zeros
                        !
                        do ix = 2, out_phi%nx
                            do iy = 2, out_phi%ny
                                !
                                do iz = 2, out_phi%grid%nzAir
                                    out_phi%v(ix, iy, iz) = &
                                    SIGMA_AIR * (in_e%x(ix, iy, iz) - in_e%x(ix - 1, iy, iz)) * &
                                    in_e%grid%del_x_inv(ix) + &
                                    SIGMA_AIR * (in_e%y(ix, iy, iz) - in_e%y(ix, iy - 1, iz)) * &
                                    in_e%grid%del_y_inv(iy) + &
                                    SIGMA_AIR * (in_e%z(ix, iy, iz) - in_e%z(ix, iy, iz - 1)) * &
                                    in_e%grid%del_z_inv(iz)
                                enddo
                                !
                                !> FOR NODES AT THE AIR-EARTH INTERFACE
                                iz = out_phi%grid%nzAir + 1
                                !
                                out_phi%v(ix, iy, iz) = &
                                (self%sigma_e%x(ix, iy, iz) * in_e%x(ix, iy, iz) -         &
                                self%sigma_e%x(ix - 1, iy, iz) * in_e%x(ix - 1, iy, iz)) * &
                                in_e%grid%del_x_inv(ix) + &
                                (self%sigma_e%y(ix, iy, iz) * in_e%y(ix, iy, iz) -         &
                                self%sigma_e%y(ix, iy - 1, iz) * in_e%y(ix, iy - 1, iz)) * &
                                in_e%grid%del_y_inv(iy) + &
                                (self%sigma_e%z(ix, iy, iz) * in_e%z(ix, iy, iz) -         &
                                SIGMA_AIR * in_e%z(ix, iy, iz - 1)) * &
                                in_e%grid%del_z_inv(iz)
                                !
                                !> FOR NODES INSIDE THE EARTH ONLY
                                !> THE TOP MOST EARTH NODE HAS AN INTERFACE WITH
                                !> AIR, THEREFORE THAT ONE IS SKIPPED HERE
                                do iz = out_phi%grid%nzAir + 2, out_phi%nz
                                    out_phi%v(ix, iy, iz) = &
                                    (self%sigma_e%x(ix,iy,iz)*in_e%x(ix, iy, iz) -                 &
                                    self%sigma_e%x(ix - 1,iy,iz)*in_e%x(ix - 1, iy, iz)) * &
                                    in_e%grid%del_x_inv(ix)            &
                                    +    (self%sigma_e%y(ix,iy,iz)*in_e%y(ix, iy, iz) -            &
                                    self%sigma_e%y(ix,iy - 1,iz)*in_e%y(ix, iy - 1, iz)) * &
                                    in_e%grid%del_y_inv(iy)            &
                                    +    (self%sigma_e%z(ix,iy,iz)*in_e%z(ix, iy, iz) -            &
                                    self%sigma_e%z(ix,iy,iz - 1)*in_e%z(ix, iy, iz - 1)) * &
                                    in_e%grid%del_z_inv(iz)
                                enddo
                            enddo
                        enddo
                        !
                    class default
                        call errStop( "divC_ModelOperator_MF_SG > in_e type unknown" )
                end select
                !
            class default
                call errStop( "divC_ModelOperator_MF_SG > out_phi type unknown" )
                !
        end select
        !
    end subroutine divC_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    subroutine divCGrad_ModelOperator_MF_SG( self, in_phi, out_phi )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( in ) :: self
        class( Scalar_t ), intent( inout ) :: in_phi
        class( Scalar_t ), intent( inout ) :: out_phi
        !
        integer :: ix, iy, iz
        !
        select type( out_phi )
            !
            class is( cScalar3D_SG_t )
                !
                if( .NOT. out_phi%is_allocated) then
                    call errStop( "divCGrad_ModelOperator_MF_SG > out_phi not allocated" )
                endif
                !
                select type( in_phi )
                    !
                    class is( cScalar3D_SG_t )
                        !
                        !> zero output (already allocated) to start
                        call out_phi%zeros
                        !
                        !> The coefficients are only for interior nodes
                        do iz = 2, in_phi%nz
                            do iy = 2, in_phi%ny
                                do ix = 2, in_phi%nx
                                    out_phi%v(ix, iy, iz) = &
                                    in_phi%v(ix + 1, iy, iz) * self%db2%x(ix,iy,iz) + &
                                    in_phi%v(ix - 1, iy, iz) * self%db1%x(ix, iy, iz) + &
                                    in_phi%v(ix, iy + 1, iz) * self%db2%y(ix, iy, iz) + &
                                    in_phi%v(ix, iy - 1, iz) * self%db1%y(ix, iy, iz) + &
                                    in_phi%v(ix, iy, iz + 1) * self%db2%z(ix, iy, iz) + &
                                    in_phi%v(ix, iy, iz - 1) * self%db1%z(ix, iy, iz) + &
                                    in_phi%v(ix, iy, iz) * self%c%v(ix, iy, iz)
                                enddo
                            enddo
                        enddo
                        !
                    class default
                        call errStop( "divCGrad_ModelOperator_MF_SG > Incompatible in_phi." )
                        !
                end select
                !
            class default
                call errStop( "divCGrad_ModelOperator_MF_SG > Incompatible out_phi." )
                !
        end select
        !
    end subroutine divCGrad_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    !
    subroutine grad_ModelOperator_MF_SG( self, in_phi, out_e )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( in ) :: self
        class( Scalar_t ), intent( inout ) :: in_phi
        class( Vector_t ), intent( inout ) :: out_e
        !
        integer :: ix, iy, iz
        !
        select type( out_e )
            !
            class is ( cVector3D_SG_t )
                !
                if( .NOT. out_e%is_allocated ) then
                    call errStop( "grad_ModelOperator_MF_SG > out_e not allocated" )
                endif
                !
                select type( in_phi )
                    !
                    class is ( cScalar3D_SG_t )
                        !
                        call out_e%Zeros
                        !
                        do ix = 1, self%metric%grid%nx 
                            do iy = 2, self%metric%grid%ny
                                do iz = 2, self%metric%grid%nz
                                    out_e%x(ix, iy, iz) = (in_phi%v(ix + 1, iy, iz) - &
                                    in_phi%v(ix, iy, iz)) / self%metric%grid%dx(ix)
                                enddo
                            enddo
                        enddo
                        !
                        do ix = 2, self%metric%grid%nx 
                            do iy = 1, self%metric%grid%ny
                                do iz = 2, self%metric%grid%nz
                                    out_e%y(ix, iy, iz) = (in_phi%v(ix, iy + 1, iz) - &
                                    in_phi%v(ix, iy, iz)) / self%metric%grid%dy(iy)
                                enddo
                            enddo
                        enddo
                        !
                        do ix = 2, self%metric%grid%nx 
                            do iy = 2, self%metric%grid%ny
                                do iz = 1, self%metric%grid%nz    
                                    out_e%z(ix, iy, iz) = (in_phi%v(ix, iy, iz + 1) - &
                                    in_phi%v(ix, iy, iz)) / self%metric%grid%dz(iz)
                                enddo
                            enddo
                        enddo
                        !
                    class default
                        call errStop( "grad_ModelOperator_MF_SG > in_phi type unknown" )
                end select
                !
            class default
                call errStop( "grad_ModelOperator_MF_SG > out_e type unknown" )
        end select
        !
    end subroutine grad_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    !
    subroutine create_ModelOperator_MF_SG( self, grid )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( inout ) :: self
        class( Grid_t ), target, intent( in ) :: grid
        !
        self%is_allocated = .FALSE.
        !
        self%metric%grid => grid
        !
        call self%alloc
        !
    end subroutine create_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    !
    subroutine allocate_ModelOperator_MF_SG( self )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( inout ) :: self
        !
        allocate( self%xXY( self%metric%grid%ny + 1, 2 ) )
        allocate( self%xXZ( self%metric%grid%nz + 1, 2 ) )
        allocate( self%xY( self%metric%grid%nx, self%metric%grid%ny + 1 ) )
        allocate( self%xZ( self%metric%grid%nx, self%metric%grid%nz + 1 ) )
        allocate( self%xXO( self%metric%grid%ny, self%metric%grid%nz) )
        !
        allocate( self%yYZ( self%metric%grid%nz + 1, 2) )
        allocate( self%yYX( self%metric%grid%nx + 1, 2) )
        allocate( self%yZ( self%metric%grid%ny, self%metric%grid%nz + 1 ) )
        allocate( self%yX( self%metric%grid%nx + 1, self%metric%grid%ny ) )
        allocate( self%yYO( self%metric%grid%nx, self%metric%grid%nz ) )
        !
        allocate( self%zZX( self%metric%grid%nx + 1, 2 ) )
        allocate( self%zZY( self%metric%grid%ny + 1, 2) )
        allocate( self%zX( self%metric%grid%nx + 1, self%metric%grid%nz ) )
        allocate( self%zY( self%metric%grid%ny + 1, self%metric%grid%nz ) )
        allocate( self%zZO( self%metric%grid%nx, self%metric%grid%ny ) )
        !
        self%xXY = R_ZERO
        self%xXZ = R_ZERO
        self%xY = R_ZERO
        self%xZ = R_ZERO
        self%xXO = R_ZERO
        self%yYX = R_ZERO
        self%yYZ = R_ZERO
        self%yX = R_ZERO
        self%yZ = R_ZERO
        self%zZX = R_ZERO
        self%zZY = R_ZERO
        self%zX = R_ZERO
        self%zY = R_ZERO
        self%zZO = R_ZERO
        !
        self%sigma_e = rVector3D_SG_t( self%metric%grid, EDGE )
        self%db1 = rVector3D_SG_t( self%metric%grid, EDGE )
        self%db2 = rVector3D_SG_t( self%metric%grid, EDGE )
        self%c = rScalar3D_SG_t( self%metric%grid, NODE )
        !
        self%is_allocated = .TRUE.
        !
    end subroutine allocate_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    subroutine deallocate_ModelOperator_MF_SG( self )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( inout ) :: self
        !
        call self%baseDealloc
        !
        deallocate( self%xXY )
        deallocate( self%xXZ )
        deallocate( self%xY )
        deallocate( self%xZ )
        deallocate( self%xXO )
        !
        deallocate( self%yYZ )
        deallocate( self%yYX )
        deallocate( self%yZ )
        deallocate( self%yX )
        deallocate( self%yYO )
        !
        deallocate( self%zZX )
        deallocate( self%zZY )
        deallocate( self%zX )
        deallocate( self%zY )
        deallocate( self%zZO )
        !
        self%is_allocated = .FALSE.
        !
    end subroutine deallocate_ModelOperator_MF_SG
    !
    !> No subroutine briefing
    subroutine print_ModelOperator_MF_SG( self )
        implicit none
        !
        class( ModelOperator_MF_SG_t ), intent( in ) :: self
        !
        call errStop( "print_ModelOperator_MF_SG not implemented yet" )
        !
    end subroutine print_ModelOperator_MF_SG
    !
end module ModelOperator_MF_SG
