!
!> Derived class to define a ModelOperator_MF
!
module ModelOperator_MF
    !
    use Constants
    use MetricElements_CSG
    use ModelParameter
    use ModelOperator
    !
    type, extends( ModelOperator_t ) :: ModelOperator_MF_t
         !
         logical :: eqset 
         !
         integer :: mKey(8)
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
         type( rVector3D_SG_t ) :: Sigma_E
         !
         type( rVector3D_SG_t ) :: db1, db2
         !
         type( rScalar3D_SG_t ) :: c
         !
         contains
              !
              final :: ModelOperator_MF_dtor
              !
              procedure, public :: setEquations => setEquationsModelOperatorMF
              procedure, public :: setCond => setCondModelOperatorMF
              procedure, public :: amult => amultModelOperatorMF
              procedure, public :: multAib => multAibModelOperatorMF
              procedure, public :: multCurlT => multCurlTModelOperatorMF
              procedure, public :: divCorSetUp => divCorsetUpModelOperatorMF
              !
              procedure, public :: AdjtBC => AdjtBCModelOperatorMF
              !
              procedure :: divCgrad => divCgradModelOperatorMF
              procedure :: divC => divCModelOperatorMF
              procedure :: grad => gradModelOperatorMF
              procedure :: div => divModelOperatorMF
              !
              procedure :: create => createModelOperatorMF 
              procedure :: allocate => allocateModelOperatorMF
              procedure :: deallocate => deallocateModelOperatorMF
              !
              procedure, public :: print => printModelOperatorMF
              !
    end type ModelOperator_MF_t
    !
    interface ModelOperator_MF_t
        module procedure ModelOperator_MF_ctor
    end interface ModelOperator_MF_t
    !
contains
    !
    !> No function briefing
    function ModelOperator_MF_ctor( grid ) result( self )
        implicit none
        !
        class( Grid3D_SG_t ), target, intent( in ) :: grid
        !
        type( ModelOperator_MF_t ) :: self
        !
        !write( *, * ) "Constructor ModelOperator_MF"
        !
        call self%init()
        !
        self%eqset = .FALSE.
        !
        call date_and_time( values=self%mKey )
        !
        !> Instantiation of the specific object MetricElements
        allocate( self%metric, source = MetricElements_CSG_t( grid ) )
        !
        call self%create( grid )
        !
    end function ModelOperator_MF_ctor
    !
    !> ModelOperator_MF destructor
    subroutine ModelOperator_MF_dtor( self )
        implicit none
        !
        type( ModelOperator_MF_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor ModelOperator_MF_t"
        !
        call self%dealloc()
        !
        call self%deallocate()
        !
    end subroutine ModelOperator_MF_dtor
    !
    !> No subroutine briefing
    subroutine createModelOperatorMF( self, grid )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( inout ) :: self
        class( Grid_t ), target, intent( in ) :: grid
        !
        self%is_allocated = .FALSE.
        !
        self%metric%grid => grid
        !
        call self%allocate()
        !
    end subroutine createModelOperatorMF
    !
    !> No subroutine briefing
    subroutine allocateModelOperatorMF( self )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( inout ) :: self
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
        self%xY  = R_ZERO
        self%xZ  = R_ZERO
        self%xXO = R_ZERO
        self%yYX = R_ZERO
        self%yYZ = R_ZERO
        self%yX  = R_ZERO
        self%yZ  = R_ZERO
        self%zZX = R_ZERO
        self%zZY = R_ZERO
        self%zX  = R_ZERO
        self%zY  = R_ZERO
        self%zZO = R_ZERO
        !
        self%Sigma_E = rVector3D_SG_t( self%metric%grid, EDGE )
        self%db1 = rVector3D_SG_t( self%metric%grid, EDGE )
        self%db2 = rVector3D_SG_t( self%metric%grid, EDGE )
        self%c = rScalar3D_SG_t( self%metric%grid, NODE )
        !
        self%is_allocated = .TRUE.
        !
    end subroutine allocateModelOperatorMF
    !
    !> No subroutine briefing
    subroutine deallocateModelOperatorMF( self )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( inout ) :: self
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
    end subroutine deallocateModelOperatorMF
    !
    !> No subroutine briefing
    subroutine setEquationsModelOperatorMF( self )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz 
        !
        do iy = 2, self%metric%grid%ny
            self%xXY(iy, 2) = -1.0 / (self%metric%grid%delY(iy) * self%metric%grid%dy(iy))
            self%xXY(iy, 1) = -1.0 / (self%metric%grid%delY(iy) * self%metric%grid%dy(iy-1))
        enddo
        !
        do iz = 2, self%metric%grid%nz
            self%xXZ(iz, 2) = -1.0 / (self%metric%grid%delZ(iz) * self%metric%grid%dz(iz))
            self%xXZ(iz, 1) = -1.0 / (self%metric%grid%delZ(iz) * self%metric%grid%dz(iz-1))
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
                self%xY(ix, iy) = 1.0 / (self%metric%grid%delY(iy)*self%metric%grid%dx(ix))
            enddo
        enddo
        !
        do ix = 1, self%metric%grid%nx
            do iz = 2, self%metric%grid%nz
                self%xZ(ix, iz) = 1.0 / (self%metric%grid%delZ(iz)*self%metric%grid%dx(ix))
            enddo
        enddo
        !
        do iz = 2, self%metric%grid%nz
            self%yYZ(iz, 2) = -1.0 / (self%metric%grid%delZ(iz)*self%metric%grid%dz(iz))
            self%yYZ(iz, 1) = -1.0 / (self%metric%grid%delZ(iz)*self%metric%grid%dz(iz-1))
        enddo
        !
        do ix = 2, self%metric%grid%nx
            self%yYX(ix, 2) = -1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dx(ix))
            self%yYX(ix, 1) = -1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dx(ix-1))
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
                self%yZ(iy, iz) = 1.0 / (self%metric%grid%delZ(iz)*self%metric%grid%dy(iy))
            enddo
        enddo
        !
        do ix = 2, self%metric%grid%nx
            do iy = 1, self%metric%grid%ny
                self%yX(ix, iy) = 1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dy(iy))
            enddo
        enddo
        !
        do ix = 2, self%metric%grid%nx
            self%zZX(ix, 2) = -1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dx(ix))
            self%zZX(ix, 1) = -1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dx(ix-1))
        enddo
        !
        do iy = 2, self%metric%grid%ny
            self%zZY(iy, 2) = -1.0 / (self%metric%grid%delY(iy)*self%metric%grid%dy(iy))
            self%zZY(iy, 1) = -1.0 / (self%metric%grid%delY(iy)*self%metric%grid%dy(iy-1))
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
                    self%zX(ix, iz) = 1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dz(iz))
             enddo
        enddo
        !
        do iy = 2, self%metric%grid%ny
            do iz = 1, self%metric%grid%nz
                self%zY(iy, iz) = 1.0 / (self%metric%grid%delY(iy)*self%metric%grid%dz(iz))
            enddo
        enddo
        !
        self%eqset = .TRUE.
        !
    end subroutine setEquationsModelOperatorMF
    !
    !> No subroutine briefing
    subroutine setCondModelOperatorMF( self, sigma )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: sigma
        !
        call sigma%PDEmapping( self%sigma_E )
        !

    end subroutine setCondModelOperatorMF
    !
    !> No subroutine briefing
    subroutine divCorsetUpModelOperatorMF( self )
        implicit none
        !
        class( ModelOperator_MF_t ), intent(inout) :: self
        !
        integer :: ix, iy, iz
        !
        do iz = 2, self%metric%grid%nz
            do iy = 2, self%metric%grid%ny
                do ix = 2, self%metric%grid%nx
                    self%db1%x(ix, iy, iz) = self%Sigma_E%x(ix - 1, iy, iz)/ &
                    (self%metric%grid%dx(ix - 1)*self%metric%grid%delX(ix))
                    !
                    self%db2%x(ix, iy, iz) = self%Sigma_E%x(ix, iy, iz)/ &
                    (self%metric%grid%dx(ix)*self%metric%grid%delX(ix))
                    !
                    self%db1%y(ix, iy, iz) = self%Sigma_E%y(ix, iy - 1, iz)/ &
                    (self%metric%grid%dy(iy - 1)*self%metric%grid%delY(iy))
                    !
                    self%db2%y(ix, iy, iz) = self%Sigma_E%y(ix, iy, iz)/ &
                    (self%metric%grid%dy(iy)*self%metric%grid%delY(iy))
                    !
                    self%db1%z(ix, iy, iz) = self%Sigma_E%z(ix, iy, iz - 1)/ &
                    (self%metric%grid%dz(iz - 1)*self%metric%grid%delZ(iz))
                    !
                    self%db2%z(ix, iy, iz) = self%Sigma_E%z(ix, iy, iz)/ &
                    (self%metric%grid%dz(iz)*self%metric%grid%delZ(iz))
                    !
                    self%c%v(ix, iy, iz) = - (self%db1%x(ix, iy, iz) + &
                    self%db2%x(ix, iy, iz) + &
                    self%db1%y(ix, iy, iz) + &
                    self%db2%y(ix, iy, iz) + &
                    self%db1%z(ix, iy, iz) + &
                    self%db2%z(ix, iy, iz))
                enddo
            enddo
        enddo
        !
        select type( vnode => self%Metric%Vnode )
            class is( rScalar3D_SG_t )
                !
                do iz = 2, self%metric%grid%nz
                    do iy = 2, self%metric%grid%ny
                        do ix = 2, self%metric%grid%nx
                            self%db1%x(ix, iy, iz) = self%db1%x(ix, iy, iz) * &
                            vnode%v(ix,iy,iz)
                            !
                            self%db1%y(ix, iy, iz) = self%db1%y(ix, iy, iz) * &
                            vnode%v(ix,iy,iz)
                            !
                            self%db1%z(ix, iy, iz) = self%db1%z(ix, iy, iz) * &
                            vnode%v(ix,iy,iz)
                            !
                            self%db2%x(ix, iy, iz) = self%db2%x(ix, iy, iz) * &
                            vnode%v(ix,iy,iz)
                            !
                            self%db2%y(ix, iy, iz) = self%db2%y(ix, iy, iz) * &
                            vnode%v(ix,iy,iz)
                            !
                            self%db2%z(ix, iy, iz) = self%db2%z(ix, iy, iz) * &
                            vnode%v(ix,iy,iz)
                        enddo
                    enddo
                enddo
                    !
            class default
                stop "Error: getFullVectorCSparsevector3D_SG > undefined grid"
                !
        end select
        !
        call self%c%mult( self%Metric%Vnode )
        !
        self%db1%x(2, :, :) = R_ZERO
        self%db1%y(:, 2, :) = R_ZERO
        self%db1%z(:, :, 2) = R_ZERO
        !
        self%db2%x(self%metric%grid%nx, :, :) = R_ZERO
        self%db2%y(:, self%metric%grid%ny, :) = R_ZERO
        self%db2%z(:, :, self%metric%grid%nz) = R_ZERO
        !
    end subroutine divCorsetupModelOperatorMF
    !
    !> No subroutine briefing
    subroutine amultModelOperatorMF( self, omega, inE, outE, p_adjoint )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        real( kind=prec ), intent( in ), optional :: omega
        class( Vector_t ), intent( in ) :: inE
        class( Vector_t ), intent( inout ) :: outE
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
            cvalue = -ONE_I * omega * isign * MU_0
        else
            cvalue = ONE_I * omega * isign * MU_0
        endif
        !
        select type( inE )
            !
            class is( cVector3D_SG_t )
                !
                if( .NOT. outE%is_allocated ) then
                    stop "Error: amultModelOperatorMF > output vector outE not allocated"
                endif
                !
                select type( outE )
                    !
                    class is( cVector3D_SG_t )
                        !
                        call outE%Zeros()
                        !
                        do iz = 2, inE%nz
                            do iy = 2, inE%ny
                                do ix = 1, inE%nx
                                    outE%x(ix, iy, iz) = self%xY(ix, iy)*(inE%y(ix + 1, iy, iz) - &
                                    inE%y(ix, iy, iz) - inE%y(ix + 1, iy - 1, iz) + &
                                    inE%y(ix, iy - 1, iz)) + &
                                    self%xZ(ix, iz) * (inE%z(ix + 1, iy, iz) - inE%z(ix, iy, iz) - &
                                    inE%z(ix + 1, iy, iz - 1) + inE%z(ix, iy, iz - 1)) + &
                                    self%xXY(iy, 2) * inE%x(ix, iy + 1, iz) + &
                                    self%xXY(iy, 1) * inE%x(ix, iy - 1, iz) + &
                                    self%xXZ(iz, 2) * inE%x(ix, iy, iz + 1) + &
                                    self%xXZ(iz, 1) * inE%x(ix, iy, iz - 1) + &
                                    (self%xXO(iy, iz)+cvalue*self%Sigma_E%x(ix,iy,iz)) * inE%x(ix, iy, iz)
                                enddo
                            enddo
                        enddo
                        !
                        do iz = 2, inE%nz
                            do iy = 1, inE%ny
                                do ix = 2, inE%nx
                                    outE%y(ix, iy, iz) = self%yZ(iy, iz) * (inE%z(ix, iy + 1, iz) - &
                                    inE%z(ix, iy, iz) - inE%z(ix, iy + 1, iz - 1) + inE%z(ix, iy, iz - 1)) + &
                                    self%yX(ix, iy) * (inE%x(ix, iy + 1, iz) - inE%x(ix, iy, iz) - &
                                    inE%x(ix - 1, iy + 1, iz) + inE%x(ix - 1, iy, iz)) + &
                                    self%yYZ(iz, 2) * inE%y(ix, iy, iz + 1) + &
                                    self%yYZ(iz, 1) * inE%y(ix, iy, iz - 1) + &
                                    self%yYX(ix, 2) * inE%y(ix + 1, iy, iz) + &
                                    self%yYX(ix, 1) * inE%y(ix - 1, iy, iz) + &
                                    (self%yYO(ix, iz)+cvalue*self%Sigma_E%y(ix,iy,iz)) * inE%y(ix, iy, iz)
                                enddo
                            enddo
                        enddo
                        !
                        do iz = 1, inE%nz
                            do iy = 2, inE%ny
                                do ix = 2, inE%nx
                                    outE%z(ix, iy, iz) = self%zX(ix, iz) * (inE%x(ix, iy, iz + 1) - &
                                    inE%x(ix, iy, iz) - inE%x(ix - 1, iy, iz + 1) + inE%x(ix - 1, iy, iz)) + &
                                    self%zY(iy,iz) * (inE%y(ix, iy, iz + 1) - inE%y(ix, iy, iz) - &
                                    inE%y(ix, iy - 1, iz + 1) + inE%y(ix, iy - 1, iz)) + &
                                    self%zZX(ix, 2) * inE%z(ix + 1, iy, iz) + &
                                    self%zZX(ix, 1) * inE%z(ix - 1, iy, iz) + &
                                    self%zZY(iy, 2) * inE%z(ix, iy + 1, iz) + &
                                    self%zZY(iy, 1) * inE%z(ix, iy - 1, iz) + &
                                    (self%zZO(ix, iy)+cvalue*self%Sigma_E%z(ix,iy,iz)) * inE%z(ix, iy, iz)
                                enddo
                            enddo
                        enddo
                        !
                        call outE%mult( self%metric%Vedge )
                        !
                    class default
                        stop "Error: amultModelOperatorMF > Undefined outE."
                        !
                end select
                !
            class default
                stop "Error: amultModelOperatorMF > Undefined inE."
                !
        end select
        !
    end subroutine amultModelOperatorMF
    !
    !> No subroutine briefing
    subroutine multAibModelOperatorMF( self, bdry, outE )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: bdry
        class( Vector_t ), intent( inout ) :: outE
        !
        real( kind=prec ) :: omega
        !
        if(.NOT. outE%is_allocated) then
            stop "Error: multAibModelOperatorMF > output vector not allocated"
        endif
        !
        omega = R_ZERO
        !
        call self%amult( omega, bdry, outE ) 
        !
    end subroutine multAibModelOperatorMF
    !
    !> No subroutine briefing
    subroutine multCurlTModelOperatorMF( self, inH, outE )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: inH
        class( Vector_t ), allocatable, intent( inout ) :: outE
        !
        integer :: ix, iy, iz
        !
        select type( inH )
            class is( cVector3D_SG_t )
                !
                call inH%div( self%Metric%FaceArea )
                !
                if(.NOT.outE%is_allocated) then
                     write( *, * ) "Error:  multCurlTModelOperatorMF > output vector not allocated"
                endif
                !
                select type( outE )
                    class is( cVector3D_SG_t )
                        !
                        !> Ex
                        do iy = 2, inH%Ny
                            do iz = 2, inH%Nz
                                outE%x(:, iy, iz) =    (inH%z(:, iy, iz) - &
                                inH%z(:, iy - 1, iz)) - &
                                (inH%y(:, iy, iz) - inH%y(:, iy, iz - 1))
                            enddo
                        enddo
                        !
                        !> Ey
                        do iz = 2, inH%Nz
                            do ix = 2, inH%Nx
                                outE%y(ix, :, iz) = (inH%x(ix, :, iz) - &
                                inH%x(ix, :, iz - 1)) - &
                                (inH%z(ix, :, iz) - inH%z(ix - 1, :, iz))
                            enddo
                        enddo
                        !
                        !> Ez
                        do ix = 2, inH%Nx
                            do iy = 2, inH%Ny
                                outE%z(ix,iy,:) = (inH%y(ix, iy, :) - &
                                inH%y(ix - 1, iy, :)) - &
                                (inH%x(ix, iy, :) - inH%x(ix, iy - 1, :))
                            enddo
                        enddo
                        !
                    class default
                        stop "Error: multCurlTModelOperatorMF > Incompatible input [outE]"
                end select
                    !> 
            class default
                stop "Error: multCurlTModelOperatorMF > Incompatible input [inH]"
                !
        end select
        !
        call outE%mult( self%metric%EdgeLength )
        !
    end subroutine multCurlTModelOperatorMF
    !
    !> No subroutine briefing
    subroutine divCgradModelOperatorMF( self, inPhi, outPhi )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        class( Scalar_t ), intent( in ) :: inPhi
        class( Scalar_t ), intent( inout ) :: outPhi
        !
        integer :: ix, iy, iz
        !
        select type( outphi )
            class is( cScalar3D_SG_t )
                !
                if(.NOT.outPhi%is_allocated) then
                    stop "Error: divCgradModelOperatorMF > Output cScalar object not allocated"
                endif
                !
                select type( inPhi )
                class is( cScalar3D_SG_t )
                    !
                    !> zero output (already allocated) to start
                    call outPhi%Zeros()
                    !
                    !> The coefficients are only for interior nodes
                    do iz = 2, inPhi%nz
                        do iy = 2, inPhi%ny
                            do ix = 2, inPhi%nx
                                outPhi%v(ix, iy, iz) = &
                                inPhi%v(ix + 1, iy, iz) * self%db2%x(ix,iy,iz) + &
                                inPhi%v(ix - 1, iy, iz) * self%db1%x(ix, iy, iz) + &
                                inPhi%v(ix, iy + 1, iz) * self%db2%y(ix, iy, iz) + &
                                inPhi%v(ix, iy - 1, iz) * self%db1%y(ix, iy, iz) + &
                                inPhi%v(ix, iy, iz + 1) * self%db2%z(ix, iy, iz) + &
                                inPhi%v(ix, iy, iz - 1) * self%db1%z(ix, iy, iz) + &
                                inPhi%v(ix, iy, iz) * self%c%v(ix, iy, iz)
                            enddo
                        enddo
                    enddo
                    !
                class default
                    stop "Error: divCgradModelOperatorMF > Incompatible input [inPhi]."
                    !
            end select
                    !
            class default
                stop "Error: divCgradModelOperatorMF > Incompatible input [x]."
                !
        end select
        !
    end subroutine divCgradModelOperatorMF
    !
    !> No subroutine briefing
    subroutine divCModelOperatorMF( self, inE, outPhi )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: inE
        class( Scalar_t ), intent( inout ) :: outPhi
        !
        integer :: ix, iy, iz
        !
        select type( outPhi )
            class is( cScalar3D_SG_t )
                !
                if( .NOT. outPhi%is_allocated) then
                    stop "Error: divCModelOperatorMF > Output cScalar object not allocated"
                endif
                !
                select type( inE )
                    class is ( cVector3D_SG_t )
                        !
                        call outPhi%Zeros()
                        !
                        do ix = 2, outPhi%nx
                            do iy = 2, outPhi%ny
                                !
                                do iz = 2, outPhi%grid%nzAir
                                    outPhi%v(ix, iy, iz) = &
                                    SIGMA_AIR * (inE%x(ix, iy, iz) - inE%x(ix - 1, iy, iz)) * &
                                    inE%grid%delXinv(ix) + &
                                    SIGMA_AIR * (inE%y(ix, iy, iz) - inE%y(ix, iy - 1, iz)) * &
                                    inE%grid%delYinv(iy) + &
                                    SIGMA_AIR * (inE%z(ix, iy, iz) - inE%z(ix, iy, iz - 1)) * &
                                    inE%grid%delZinv(iz)
                                enddo
                                
                                !> FOR NODES AT THE AIR-EARTH INTERFACE
                                iz = outPhi%grid%nzAir + 1
                                !
                                outPhi%v(ix, iy, iz) = &
                                (self%Sigma_E%x(ix, iy, iz) * inE%x(ix, iy, iz) -         &
                                self%Sigma_E%x(ix - 1, iy, iz) * inE%x(ix - 1, iy, iz)) * &
                                inE%grid%delXinv(ix) + &
                                (self%Sigma_E%y(ix, iy, iz) * inE%y(ix, iy, iz) -         &
                                self%Sigma_E%y(ix, iy - 1, iz) * inE%y(ix, iy - 1, iz)) * &
                                inE%grid%delYinv(iy) + &
                                (self%Sigma_E%z(ix, iy, iz) * inE%z(ix, iy, iz) -         &
                                SIGMA_AIR * inE%z(ix, iy, iz - 1)) * &
                                inE%grid%delZinv(iz)

                                !> FOR NODES INSIDE THE EARTH ONLY
                                !> THE TOP MOST EARTH NODE HAS AN INTERFACE WITH
                                !> AIR, THEREFORE THAT ONE IS SKIPPED HERE
                                do iz = outPhi%grid%nzAir + 2, outPhi%nz
                                    outPhi%v(ix, iy, iz) = &
                                    (self%Sigma_E%x(ix,iy,iz)*inE%x(ix, iy, iz) -                 &
                                    self%Sigma_E%x(ix - 1,iy,iz)*inE%x(ix - 1, iy, iz)) * &
                                    inE%grid%delXinv(ix)            &
                                    +    (self%Sigma_E%y(ix,iy,iz)*inE%y(ix, iy, iz) -            &
                                    self%Sigma_E%y(ix,iy - 1,iz)*inE%y(ix, iy - 1, iz)) * &
                                    inE%grid%delYinv(iy)            &
                                    +    (self%Sigma_E%z(ix,iy,iz)*inE%z(ix, iy, iz) -            &
                                    self%Sigma_E%z(ix,iy,iz - 1)*inE%z(ix, iy, iz - 1)) * &
                                    inE%grid%delZinv(iz)
                                enddo
                            enddo
                        enddo
                        !
                    class default
                        stop "Error: divCModelOperatorMF > inE type unknown"
                end select
                !
            class default
                stop "Error: divCModelOperatorMF > outPhi type unknown"
                !
        end select
        !
    end subroutine divCModelOperatorMF
    !
    !> No subroutine briefing
    subroutine gradModelOperatorMF( self, inPhi, outE )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        class( Scalar_t ), intent( in ) :: inPhi
        class( Vector_t ), intent( inout ) :: outE
        !
        integer :: ix, iy, iz
        !
        select type( outE )
            class is ( cVector3D_SG_t )
                !
                if(.NOT.outE%is_allocated) then
                    stop "Error: gradModelOperatorMF > Output cVector object not allocated"
                endif
                !
                select type( inPhi )
                    class is ( cScalar3D_SG_t )
                        !
                        call outE%Zeros
                        !
                        do ix = 1, self%metric%grid%nx 
                            do iy = 2, self%metric%grid%ny
                                do iz = 2, self%metric%grid%nz
                                    outE%x(ix, iy, iz) = (inPhi%v(ix + 1, iy, iz) - &
                                    inPhi%v(ix, iy, iz)) / self%metric%grid%dx(ix)
                                enddo
                            enddo
                        enddo
                        !
                        do ix = 2, self%metric%grid%nx 
                            do iy = 1, self%metric%grid%ny
                                do iz = 2, self%metric%grid%nz
                                    outE%y(ix, iy, iz) = (inPhi%v(ix, iy + 1, iz) - &
                                    inPhi%v(ix, iy, iz)) / self%metric%grid%dy(iy)
                                enddo
                            enddo
                        enddo
                        !
                        do ix = 2, self%metric%grid%nx 
                            do iy = 2, self%metric%grid%ny
                                do iz = 1, self%metric%grid%nz    
                                    outE%z(ix, iy, iz) = (inPhi%v(ix, iy, iz + 1) - &
                                    inPhi%v(ix, iy, iz)) / self%metric%grid%dz(iz)
                                enddo
                            enddo
                        enddo
                        !
                    class default
                        stop "Error: gradModelOperatorMF > inPhi type unknown"
                end select
            class default
                stop "Error: gradModelOperatorMF > outE type unknown"
        end select
        !
    end subroutine gradModelOperatorMF
    !
    !> No subroutine briefing
    subroutine divModelOperatorMF( self, inE, outPhi )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: inE
        class( Scalar_t ), intent( inout ) :: outPhi
        !
        integer :: ix, iy, iz
        !
        select type(outPhi)
            class is(cScalar3D_SG_t)
                !
                if(.NOT.outPhi%is_allocated) then
                    stop "Error: divModelOperatorMF > Output cScalar object not allocated"
                endif
                !
                select type( inE )
                    class is(cVector3D_SG_t)
                        !
                        call outPhi%Zeros()
                        !
                        do ix = 2, outPhi%nx
                            do iy = 2, outPhi%ny
                                do iz = 2, outPhi%grid%nz
                                    outPhi%v(ix, iy, iz) = &
                                    (inE%x(ix, iy, iz) - inE%x(ix - 1, iy, iz)) * &
                                    inE%grid%delXinv(ix) + &
                                    (inE%y(ix, iy, iz) - inE%y(ix, iy - 1, iz)) * &
                                    inE%grid%delYinv(iy) + &
                                    (inE%z(ix, iy, iz) - inE%z(ix, iy, iz - 1)) * &
                                    inE%grid%delZinv(iz)
                                enddo
                            enddo
                        enddo
                        class default
                        stop "Error: divModelOperatorMF> inE type unknown"
                end select
            class default
                stop "Error: divModelOperatorMF > outPhi type unknown"
        end select
        !
    end subroutine divModelOperatorMF
    !
    !  subroutine AdjtBC uses (adjoint) interior node solution to compute
    !  boundary node values for adjoint (or transpose) solution
    !   (NOTE: because off-diagonal part of EM operator is real this works
    !  Assuming boundary conditions for forward problem are
    !  specified tangential E fields, adjoint BC are  homogeneous (to solve for
    !   interior nodes), and solution on boundary nodes is determined from
    !   interior solution via:    E_B - adjt(A_IB)*E_I = 0
    !    where E_B is boundary part of adjoint system rhs, and E_I
    !    is the interior node solution of the adjoint system (solved with
    !    homogeneous tangential BC). This operator computes adjt(A_IB)*E_I.
    !   Output is required for calculating sensitivities
    !     of data to errors in specified BC (and potentially for other sorts
    !     of sensitivities which require the boundary nodes of the adjoint or
    !     transpose solution).
    !    NOTE: this routine can be used for both complex conjugate
    !         transpose and transpose cases.

    !   Uses curl_curl coefficients, available to all routines in this module
    !   NOTE: Must call CurlCurlSetup before use of this routine
    !
    subroutine AdjtBCModelOperatorMF( self, eIn, BC )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: eIn
        ! OUTPUT: boundary condition structure: should be allocated
        !   and initialized before call to this routine
        class( Vector_t ), intent( inout) :: BC
        !
        integer                   :: ix,iy,iz,nx,ny,nz
        !
        !  Multiply FD electric field vector defined on interior nodes (eIn) by
        !  adjoint of A_IB, the interior/boundary sub-block of the differential
        !  operator.
        !
        nx = eIn%nx
        ny = eIn%ny
        nz = eIn%nz
        !
        select type( eIn )
            class is( cVector3D_SG_t )
                !
                if( .NOT. eIn%is_allocated ) then
                    stop "Error: AdjtBCModelOperatorMF > eIn not allocated"
                endif
                !
                select type( BC )
                    class is( cVector3D_SG_t )
						!
						BC%y(:,1,1) = C_ZERO
						BC%y(:,ny,1) = C_ZERO
						BC%y(:,1,nz+1) = C_ZERO
						BC%y(:,ny,nz+1) = C_ZERO
						!
                        do ix = 1, nx
                            do iz = 2, nz
                                !
								BC%y(ix,:,iz) = (- self%yX(ix,1)*Ein%y(ix,1,iz)       &
                                + self%yX(ix+1,1)*Ein%y(ix+1,1,iz)   &
                                + self%xXY(2,1)*Ein%x(ix,2,iz))
                                BC%y(ix,:,iz) = (+ self%yX(ix,ny)*Ein%y(ix,ny,iz)     &
                                - self%yX(ix+1,ny)*Ein%y(ix+1,ny,iz) &
                                + self%xXY(ny,2)*Ein%x(ix,ny,iz))
								!
                            enddo
                        enddo
						!
						BC%y(1,1,:) = C_ZERO
						BC%y(1,ny,:) = C_ZERO
						BC%y(nx+1,1,:) = C_ZERO
						BC%y(nx+1,ny,:) = C_ZERO
						!
                        do iz = 1, nz
                            do ix = 2, nx
                                !
								BC%y(ix,:,iz) = (- self%yZ(1,iz)*Ein%y(ix,1,iz)        &
                                + self%yZ(1,iz+1)*Ein%y(ix,1,iz+1)    &
                                + self%zZY(2,1)*Ein%z(ix,2,iz))
                                BC%y(ix, :, iz) = (+ self%yZ(ny,iz)*Ein%y(ix,ny,iz)      &
                                - self%yZ(ny,iz+1)*Ein%y(ix,ny,iz+1)  &
                                + self%zZY(ny,2)*Ein%z(ix,ny,iz))
								!
                            enddo
                        enddo
						!
						BC%x(1,:,1) = C_ZERO
						BC%x(nx,:,1) = C_ZERO
						BC%x(1,:,nz+1) = C_ZERO
						BC%x(nx,:,nz+1) = C_ZERO
						!
                        do iy = 1, ny
                            do iz = 2, nz
                                !
								BC%x(:,iy,iz) = (- self%xY(1,iy)*Ein%x(1,iy,iz)        &
                                + self%xY(1,iy+1)*Ein%x(1,iy+1,iz)    &
                                + self%yYX(2,1)*Ein%y(2,iy,iz))
                                BC%x(:,iy,iz) = (+ self%xY(nx,iy)*Ein%x(nx,iy,iz)      &
                                - self%xY(nx,iy+1)*Ein%x(nx,iy+1,iz)  &
                                + self%yYX(nx,2)*Ein%y(nx,iy,iz))
								!
                            enddo
                        enddo
						!
						BC%x(1,1,:) = C_ZERO
						BC%x(nx,1,:) = C_ZERO
						BC%x(1,ny+1,:) = C_ZERO
						BC%x(nx,ny+1,:) = C_ZERO
						!
                        do iz = 1, nz
                            do iy = 2, ny
                                !
								BC%x(:,iy,iz) = (- self%xZ(1,iz)*Ein%x(1,iy,iz)       &
                                + self%xZ(1,iz+1)*Ein%x(1,iy,iz+1)   &
                                + self%zZX(2,1)*Ein%z(2,iy,iz))
                                BC%x(:,iy,iz) = (+ self%xZ(nx,iz)*Ein%x(nx,iy,iz)     &
                                - self%xZ(nx,iz+1)*Ein%x(nx,iy,iz+1) &
                                + self%zZX(nx,2)*Ein%z(nx,iy,iz))
								!
                            enddo
                        enddo
						!
						BC%z(:,1,1) = C_ZERO
						BC%z(:,1,nz) = C_ZERO
						BC%z(:,ny+1,1) = C_ZERO
						BC%z(:,ny+1,nz) = C_ZERO
						!
                        do ix = 1, nx
                            do iy = 2, ny
                                !
								BC%z(ix,iy,:) = (- self%zX(ix,1)*Ein%z(ix,iy,1)       &
                                + self%zX(ix+1,1)*Ein%z(ix+1,iy,1)   &
                                + self%xXZ(2,1)*Ein%x(ix,iy,2))
                                BC%z(ix,iy,:) = (+ self%zX(ix,nz)*Ein%z(ix,iy,nz)     &
                                - self%zX(ix+1,nz)*Ein%z(ix+1,iy,nz) &
                                + self%xXZ(nz,2)*Ein%x(ix,iy,nz))
								!
                            enddo
                        enddo
						!
						BC%z(1,:,1) = C_ZERO
						BC%z(1,:,nz) = C_ZERO
						BC%z(nx+1,:,1) = C_ZERO
						BC%z(nx+1,:,nz) = C_ZERO
						!
                        do iy = 1, ny
                            do ix = 2, nx
                                !
								BC%z(ix,iy,:) = (- self%zY(iy,1)*Ein%z(ix,iy,1)        &
                                + self%zY(iy+1,1)*Ein%z(ix,iy+1,1)    &
                                + self%yYZ(2,1)*Ein%y(ix,iy,2))
                                BC%z(ix,iy,:) = (+ self%zY(iy,nz)*Ein%z(ix,iy,nz)      &
                                - self%zY(iy+1,nz)*Ein%z(ix,iy+1,nz)  &
                                + self%yYZ(nz,2)*Ein%y(ix,iy,nz) )
								!
                            enddo
                        enddo
						!
                        class default
                        stop "Error: AdjtBCModelOperatorMF> Ein type unknown"
                        !
                end select
                !
            class default
                stop "Error: AdjtBCModelOperatorMF > BC type unknown"
                !
        end select
        !
    end subroutine AdjtBCModelOperatorMF
    !
    !> No subroutine briefing
    subroutine printModelOperatorMF( self )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        !
        stop "Subroutine print not implemented for ModelOperator_MF"
        !
    end subroutine printModelOperatorMF
    !
end module ModelOperator_MF
