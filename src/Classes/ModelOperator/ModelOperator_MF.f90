module ModelOperator_MF
    !
    use Constants
    use Grid3D_SG
    use rScalar3D_SG
    use cScalar3D_SG
    use cVector3D_SG
    use rVector3D_SG
    use MetricElements_CSG
    use ModelParameter
    use ModelOperator
    !
    type, extends( ModelOperator_t ) :: ModelOperator_MF_t
         !
         logical :: eqset      ! set to true after equations,
                               ! (part independent of sigma) is set.
         integer :: mKey(8)    ! For use with FORTRAN DATE_AND_TIME subroutine.
         !
         integer :: nx, ny, nz ! Redundant with what is in grid,
                               ! but makes code easier to read
         !**
         ! These are the arrays used to represent the curl-curl
         ! operator in MF implementation aBC is for the Ea equation
         ! using the b component in the c direction.
         ! The two array elements correspond to coefficients required to form
         ! the derivative using adjacent faces
         ! E.g., xXY --> coefficient for the Ex equation using Ex variables
         ! in the y direction.
         ! xY is the product of grid spacings in X and Y-directions, etc.
         ! The two array elements again correspond to adjacent faces
         ! xXO is the sum of the all the products for the spacing in all the
         ! directions for adjacent faces for Ex term at ix, ij, ik point
         ! (collection of left/ right horizontal and top/ bottom vertical faces),
         ! etc.
         !*
         real( kind=prec ), allocatable, dimension(:,:) :: xXY, xXZ
         real( kind=prec ), allocatable, dimension(:,:) :: xY, xZ
         real( kind=prec ), allocatable, dimension(:,:) :: xXO
         real( kind=prec ), allocatable, dimension(:,:) :: yYX, yYZ
         real( kind=prec ), allocatable, dimension(:,:) :: yX, yZ
         real( kind=prec ), allocatable, dimension(:,:) :: yYO
         real( kind=prec ), allocatable, dimension(:,:) :: zZX, zZY
         real( kind=prec ), allocatable, dimension(:,:) :: zX, zY
         real( kind=prec ), allocatable, dimension(:,:) :: zZO
         
         ! Coefficients of diagonal of (unweighted) A operator.
         ! Here we store this as real - still need to multiply by i omega
         ! Thus, don"t really need to store Adiag -- just store Sigma_E.
         type( rVector3D_SG_t ) :: Sigma_E
         
         ! Operators for divergence correction (DC) will also be included here
         ! implementation will be trhough a separate module, managed through
         ! "solver" object.     Preconditioners will also be managed through the solver
         ! Might consider making a "base" version w/o DC, and have an extension that
         ! includes these arrayso.

         ! db1%x contains coefficients of the stencil for shift -1 of ix index
         ! (%y,%z give corresponding coefficients for iy, iz)
         ! db2    contains coefficients for corresponding shift of +1
         type( rVector3D_SG_t ) :: db1, db2
         
         ! c contains the coefficients for div sigma grad
         ! operator diagonal. Note that divergence and gradient
         ! (also needed for DC) can be implemented
         ! directly using Metric Elements.
         type( rScalar3D_SG_t ) :: c
         !
         contains
              !
              final :: ModelOperator_MF_dtor
              !
              procedure, public :: setEquations => setEquationsModelOperatorMF
              procedure, public :: setCond      => setCondModelOperatorMF
              procedure, public :: amult        => amultModelOperatorMF
              procedure, public :: multAib      => multAibModelOperatorMF
              procedure, public :: multCurlT    => multCurlTModelOperatorMF
              procedure, public :: divCorSetUp  => divCorsetUpModelOperatorMF
              !
              procedure :: divCgrad => divCgradModelOperatorMF
              procedure :: divC     => divCModelOperatorMF
              procedure :: grad     => gradModelOperatorMF
              procedure :: div      => divModelOperatorMF
              !
              !procedure, public :: createScalar => createScalarModelOperatorMF
              !procedure, public :: createVector => createVectorModelOperatorMF
              !
              ! Private methods
              procedure :: create      => createModelOperatorMF 
              procedure :: allocate    => allocateModelOperatorMF
              procedure :: deallocate  => deallocateModelOperatorMF
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
    function ModelOperator_MF_ctor( grid ) result( self )
        implicit none
        !
        class( Grid3D_SG_t ), target, intent( in ) :: grid
        !
        type( ModelOperator_MF_t ) :: self
        !
        !write(*,*) "Constructor ModelOperator_MF"
        !
        call self%init()
        !
        self%eqset = .FALSE.
        !
        call date_and_time( values=self%mKey )
        !
        self%nx = 0
        self%ny = 0
        self%nz = 0
        !
        ! Instantiation of the specific object MetricElements
        allocate( self%metric, source = MetricElements_CSG_t( grid ) )
        !
        call self%create( grid )
        !
    end function ModelOperator_MF_ctor
    !
    ! ModelOperator_MF destructor
    subroutine ModelOperator_MF_dtor( self )
        implicit none
        !
        type( ModelOperator_MF_t ), intent( inout ) :: self
        !
        !write(*,*) "Destructor ModelOperator_MF_t"
        !
        call self%dealloc()
        !
        call self%deallocate()
        !
    end subroutine ModelOperator_MF_dtor
    !
    !**
    ! Create -- since this just calls allocateModelOperatorMF.
    !*
    subroutine createModelOperatorMF( self, inGrid )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( inout ) :: self
        class( Grid_t ), target, intent( in )        :: inGrid
        !
        self%is_allocated = .FALSE.
        !
        self%metric%grid => inGrid
        self%nx = inGrid%nx
        self%ny = inGrid%ny
        self%nz = inGrid%nz
        !
        call self%allocate()
        !
    end subroutine createModelOperatorMF
    
    !**
    ! allocateModelOperatorMF
    !*
    subroutine allocateModelOperatorMF( self )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( inout ) :: self
        !
        ! Allocate memory for del x del operator coefficient arrays
        ! Coefficients for difference equation only uses interior
        ! nodes. however, we need boundary nodes for the adjoint problem
        allocate( self%xXY( self%ny + 1, 2 ) )
        allocate( self%xXZ( self%nz + 1, 2 ) )
        allocate( self%xY( self%nx, self%ny + 1 ) )
        allocate( self%xZ( self%nx, self%nz + 1 ) )
        allocate( self%xXO( self%ny, self%nz) )
        !
        allocate( self%yYZ( self%nz + 1, 2) )
        allocate( self%yYX( self%nx + 1, 2) )
        allocate( self%yZ( self%ny, self%nz + 1 ) )
        allocate( self%yX( self%nx + 1, self%ny ) )
        allocate( self%yYO( self%nx, self%nz ) )

        allocate( self%zZX( self%nx + 1, 2 ) )
        allocate( self%zZY( self%ny + 1, 2) )
        allocate( self%zX( self%nx + 1, self%nz ) )
        allocate( self%zY( self%ny + 1, self%nz ) )
        allocate( self%zZO( self%nx, self%ny ) )
        !
        ! Initalize all coefficients to zero (some remain zero)
        self%xXY = 0.0
        self%xXZ = 0.0
        self%xY  = 0.0
        self%xZ  = 0.0
        self%xXO = 0.0
        self%yYX = 0.0
        self%yYZ = 0.0
        self%yX  = 0.0
        self%yZ  = 0.0
        self%zZX = 0.0
        self%zZY = 0.0
        self%zX  = 0.0
        self%zY  = 0.0
        self%zZO = 0.0
        !
        select type( grid => self%metric%grid )
            class is( Grid3D_SG_t )
                !
                self%Sigma_E = rVector3D_SG_t( grid, EDGE )
                self%db1 = rVector3D_SG_t( grid, EDGE )
                self%db2 = rVector3D_SG_t( grid, EDGE )
                self%c = rScalar3D_SG_t( grid, NODE )
                !
        end select
        !
        self%is_allocated = .TRUE.
        !
    end subroutine allocateModelOperatorMF
    !**
    ! deallocateModelOperatorMF
    !*
    subroutine deallocateModelOperatorMF( self )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( inout ) :: self
        !
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
    !***
    !     createVector
    !***
    !function createVectorModelOperatorMF( self, gridType ) result( cVec )
        !implicit none
        !
        !     this just returns a cVector of correct type
        !         for this model operator -- 
        !class( ModelOperator_MF_t ), intent( in )   :: self
        !character( len=80 ), intent( in ), optional :: gridType
        !class( cVector_t ), allocatable             :: cVec
        !
        ! do we really need select type here?????
        !select type( grid => self%metric%grid )
            !class is( Grid3D_SG_t )
            !
                !if( present( gridType ) ) then
                    !
                    !allocate( cVec, source = cVector3D_SG_t( grid, gridType ) )
                !else
                    !
                    !    EDGE gridType by default
                    !allocate( cVec, source = cVector3D_SG_t( grid, EDGE ) )
                !endif
        !end select
        !
    !end function createVectorModelOperatorMF
    !***
    !     createScalar
    !***
    !function createScalarModelOperatorMF( self, gridType ) result( c_scalar )
        !implicit none
        !
        !     this just returns a cScalar of correct type
        !         for this model operator -- 
        !class( ModelOperator_MF_t ), intent( in )   :: self
        !character( len=80 ), intent( in ), optional :: gridType
        !class( cScalar_t ), allocatable             :: c_scalar
        !
        !     do we really need select type here?????
        !select type( grid => self%metric%grid )
            !class is( Grid3D_SG_t )
                !if( present( gridType ) ) then
                    !allocate( c_scalar, source = cScalar3D_SG_t( grid, gridType ) )
                !else
                    !    NODE gridType by default
                    !allocate( c_scalar, source = cScalar3D_SG_t( grid, NODE ) )
                !endif
        !end select
        !
    !end function createScalarModelOperatorMF
    !**
    ! setEquations
    !*
    subroutine setEquationsModelOperatorMF( self )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( inout ) :: self
        !
        integer ::ix, iy, iz 
        !
        ! Coefficients for curlcurlE (del X del X E)
        ! coefficents for calculating Ex ; only loop over internal edges
        do iy = 2, self%ny
            self%xXY(iy, 2) = -1.0 / (self%metric%grid%delY(iy) * self%metric%grid%dy(iy))
            self%xXY(iy, 1) = -1.0 / (self%metric%grid%delY(iy) * self%metric%grid%dy(iy-1))
        end do
        !
        do iz = 2, self%nz
            self%xXZ(iz, 2) = -1.0 / (self%metric%grid%delZ(iz) * self%metric%grid%dz(iz))
            self%xXZ(iz, 1) = -1.0 / (self%metric%grid%delZ(iz) * self%metric%grid%dz(iz-1))
        end do
        !
        do iy = 2, self%ny
            do iz = 2, self%nz
                self%xXO(iy, iz) = -(self%xXY(iy,1) + self%xXY(iy,2) + &
                self%xXZ(iz,1) + self%xXZ(iz,2))
            end do
        end do
        !
        do ix = 1, self%nx
            do iy = 2, self%ny
                self%xY(ix, iy) = 1.0 / (self%metric%grid%delY(iy)*self%metric%grid%dx(ix))
            end do
        end do
        !
        do ix = 1, self%nx
            do iz = 2, self%nz
                self%xZ(ix, iz) = 1.0 / (self%metric%grid%delZ(iz)*self%metric%grid%dx(ix))
            end do
        end do
        ! End of Ex coefficients
        !
        ! Coefficents for calculating Ey; only loop over internal edges
        do iz = 2, self%nz
            self%yYZ(iz, 2) = -1.0 / (self%metric%grid%delZ(iz)*self%metric%grid%dz(iz))
            self%yYZ(iz, 1) = -1.0 / (self%metric%grid%delZ(iz)*self%metric%grid%dz(iz-1))
        end do
        !
        do ix = 2, self%nx
            self%yYX(ix, 2) = -1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dx(ix))
            self%yYX(ix, 1) = -1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dx(ix-1))
        end do
        !
        do ix = 2, self%nx
            do iz = 2, self%nz
                self%yYO(ix, iz) = -(self%yYX(ix,1) + self%yYX(ix,2) + &
                self%yYZ(iz,1) + self%yYZ(iz,2))
            end do
        end do
        !
        do iy = 1, self%ny
            do iz = 2, self%nz
                self%yZ(iy, iz) = 1.0 / (self%metric%grid%delZ(iz)*self%metric%grid%dy(iy))
            end do
        end do
        !
        do ix = 2, self%nx
            do iy = 1, self%ny
                self%yX(ix, iy) = 1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dy(iy))
            end do
        end do
        ! End of Ey coefficients
        !
        ! Coefficents for calculating Ez; only loop over internal edges
        do ix = 2, self%nx
            self%zZX(ix, 2) = -1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dx(ix))
            self%zZX(ix, 1) = -1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dx(ix-1))
        end do
        !
        do iy = 2, self%ny
            self%zZY(iy, 2) = -1.0 / (self%metric%grid%delY(iy)*self%metric%grid%dy(iy))
            self%zZY(iy, 1) = -1.0 / (self%metric%grid%delY(iy)*self%metric%grid%dy(iy-1))
        end do
        !
        do ix = 2, self%nx
            do iy = 2, self%ny
                self%zZO(ix, iy) = -(self%zZX(ix,1) + self%zZX(ix,2) + &
                self%zZY(iy,1) + self%zZY(iy,2))
            end do
        enddo
        !
        do ix = 2, self%nx
             do iz = 1, self%nz
                    self%zX(ix, iz) = 1.0 / (self%metric%grid%delX(ix)*self%metric%grid%dz(iz))
             end do
        end do
        !
        do iy = 2, self%ny
            do iz = 1, self%nz
                self%zY(iy, iz) = 1.0 / (self%metric%grid%delY(iy)*self%metric%grid%dz(iz))
            end do
        end do
        !
        self%eqset = .TRUE.
        !
        ! Note that divergence correction coefficients depend on conductivity
        ! and need to be reset whenever conductivity changes -- so these are
        ! not set here.
        !
    end subroutine setEquationsModelOperatorMF
    !**
    ! setCond
    !*
    subroutine setCondModelOperatorMF( self, ModPar )
        implicit none
        !
        class( ModelOperator_MF_t), intent( inout ) :: self
        class( ModelParameter_t), intent( in )      :: ModPar
        !
        call ModPar%PDEmapping( self%sigma_E )
        !
    end subroutine setCondModelOperatorMF
    !**
    ! divcorsetup
    !*
    subroutine divCorsetUpModelOperatorMF( self )
        implicit none
        !
        ! Arguments
        class( ModelOperator_MF_t ), intent(inout) :: self
        ! Local variables
        integer :: ix, iy, iz
        !
        ! The coefficients are only for the interior nodes
        ! these coefficients have not been multiplied by volume elements yet
        do iz = 2, self%nz
            do iy = 2, self%ny
                do ix = 2, self%nx
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
                end do
            end do
        end do
        !
        ! multiply by corner volume elements to make operator symmetric
        ! interior nodes only     -- note that 6 edges meet in each node, and
        !        coefficients for the corresponding edge are stored in x, y, z conmponens
        !        of rVector objects (some components are zero/not used near boundaries)
        
        select type( vnode => self%Metric%Vnode )
        class is( rScalar3D_SG_t )
            !
            do iz = 2, self%nz
                do iy = 2, self%ny
                    do ix = 2, self%nx
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
                    end do
                end do
            end do
            !
        end select
        !
        call self%c%mult( self%Metric%Vnode )
        !
        ! To be explicit set coefficients that multiply edges connected to
        ! boundary nodes to zero (this gaurantees that the BC on the potential
        ! is phi = 0):
        self%db1%x(2, :, :) = R_ZERO
        self%db1%y(:, 2, :) = R_ZERO
        self%db1%z(:, :, 2) = R_ZERO
        self%db2%x(self%nx, :, :) = R_ZERO
        self%db2%y(:, self%ny, :) = R_ZERO
        self%db2%z(:, :, self%nz) = R_ZERO
        !
    end subroutine divCorsetupModelOperatorMF
    !**
    ! Amult -- converted to subroutine
    ! This does the matrix-vector multiply (A+iwB)inE = outE
    ! parameterized by input frequency omega
    !*
    subroutine amultModelOperatorMF( self, omega, x, y, p_adjt )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        real( kind=prec ), intent( in ), optional :: omega
        class( cVector_t ), intent( in )          :: x
        class( cVector_t ), intent( inout )       :: y
        logical, intent( in ), optional           :: p_adjt
        !
        integer :: ix, iy, iz
        complex( kind=prec ) :: c
        logical :: adjt
        !
        if ( present( p_adjt ) ) then
            adjt = p_adjt
        else
            adjt = .FALSE.
        end if
        !
        if (adjt) then
            c = -ONE_I * omega * ISIGN * MU_0
        else
            c = ONE_I * omega * ISIGN * MU_0
        end if
        !
        select type(x)
        class is(cVector3D_SG_t)
            !
            if( .NOT. y%is_allocated ) then
                write(*,*) "ERROR: amult in    ModelOperator_MF"
                stop         "output vector y not allocated"
            endif
            !
            select type(y)
            class is(cVector3D_SG_t)
                call y%Zeros()
                ! Apply difference equation to compute Ex
                ! (compute only on interior nodes)
                do iz = 2, x%nz
                    do iy = 2, x%ny
                        do ix = 1, x%nx
                            y%x(ix, iy, iz) = self%xY(ix, iy)*(x%y(ix + 1, iy, iz) - &
                            x%y(ix, iy, iz) - x%y(ix + 1, iy - 1, iz) + &
                            x%y(ix, iy - 1, iz)) + &
                            self%xZ(ix, iz) * (x%z(ix + 1, iy, iz) - x%z(ix, iy, iz) - &
                            x%z(ix + 1, iy, iz - 1) + x%z(ix, iy, iz - 1)) + &
                            self%xXY(iy, 2) * x%x(ix, iy + 1, iz) + &
                            self%xXY(iy, 1) * x%x(ix, iy - 1, iz) + &
                            self%xXZ(iz, 2) * x%x(ix, iy, iz + 1) + &
                            self%xXZ(iz, 1) * x%x(ix, iy, iz - 1) + &
                            (self%xXO(iy, iz)+c*self%Sigma_E%x(ix,iy,iz)) * x%x(ix, iy, iz)
                        end do
                    end do
                end do
                !
                ! Apply difference equation to compute Ey (only on interior nodes)
                do iz = 2, x%nz
                    do iy = 1, x%ny
                        do ix = 2, x%nx
                            y%y(ix, iy, iz) = self%yZ(iy, iz) * (x%z(ix, iy + 1, iz) - &
                            x%z(ix, iy, iz) - x%z(ix, iy + 1, iz - 1) + x%z(ix, iy, iz - 1)) + &
                            self%yX(ix, iy) * (x%x(ix, iy + 1, iz) - x%x(ix, iy, iz) - &
                            x%x(ix - 1, iy + 1, iz) + x%x(ix - 1, iy, iz)) + &
                            self%yYZ(iz, 2) * x%y(ix, iy, iz + 1) + &
                            self%yYZ(iz, 1) * x%y(ix, iy, iz - 1) + &
                            self%yYX(ix, 2) * x%y(ix + 1, iy, iz) + &
                            self%yYX(ix, 1) * x%y(ix - 1, iy, iz) + &
                            (self%yYO(ix, iz)+c*self%Sigma_E%y(ix,iy,iz)) * x%y(ix, iy, iz)
                        end do
                    end do
                end do
                !
                ! Apply difference equation to compute Ez (only on interior nodes)
                do iz = 1, x%nz
                    do iy = 2, x%ny
                        do ix = 2, x%nx
                            y%z(ix, iy, iz) = self%zX(ix, iz) * (x%x(ix, iy, iz + 1) - &
                            x%x(ix, iy, iz) - x%x(ix - 1, iy, iz + 1) + x%x(ix - 1, iy, iz)) + &
                            self%zY(iy,iz) * (x%y(ix, iy, iz + 1) - x%y(ix, iy, iz) - &
                            x%y(ix, iy - 1, iz + 1) + x%y(ix, iy - 1, iz)) + &
                            self%zZX(ix, 2) * x%z(ix + 1, iy, iz) + &
                            self%zZX(ix, 1) * x%z(ix - 1, iy, iz) + &
                            self%zZY(iy, 2) * x%z(ix, iy + 1, iz) + &
                            self%zZY(iy, 1) * x%z(ix, iy - 1, iz) + &
                            (self%zZO(ix, iy)+c*self%Sigma_E%z(ix,iy,iz)) * x%z(ix, iy, iz)
                        end do
                    end do
                end do
                !
                ! Modified to add diagonal part in main loop ...
                ! Finally multiply by VEdge (in place)
                call y%mults3( self%metric%Vedge )
            !
            end select
            !
            class default
                write( *, * ) "ERROR:ModelOperator_MF::amult:"
                stop        "            Incompatible input [x]. Exiting."
                !
        end select

    end subroutine amultModelOperatorMF
    !**
    ! multAib
    ! This multiplies boundary values by the real matrix Aib -- converts
    ! boundary data to interior forcing.
    ! NOTE: Make sure interior edges are zeroed before calling
    ! this function.
    !*
    subroutine multAibModelOperatorMF( self, bdry, outE )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        !     again do these need to be abstract -- and OK if they are?
        class( cVector_t ), intent( in )    :: bdry
        class( cVector_t ), intent( inout ) :: outE
        !
        real( kind=prec ) :: omega
        !
        if(.NOT. outE%is_allocated) then
            write(*,*) "ERROR: multAib in    ModelOperator_MF"
            stop         "output vector not allocated"
        endif
        !
        omega = R_ZERO     ! diagonal part or A does not enter into this
        !
        call self%amult( omega, bdry, outE ) 
        !
    end subroutine multAibModelOperatorMF
    
    !**
    ! This multiplies Vector defined on faces (hence inH) by the transpose of
    ! the curl operator (which normally maps from edges to faces).     Not same
    ! as adjoint curl, on a non-uniform grid -- implemented using metric elements
    ! + matrix free adjoint curl on unit grid.
    ! curl.    This will enable simplified (and easily generalized to sparse
    ! matrix Model Operator implementations) computation of data functionals
    ! for computing magnetic fields.
    !*
    subroutine multCurlTModelOperatorMF( self, inH, outE )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in )        :: self
        class( cVector_t ), intent( inout )              :: inH
        class( cVector_t ), allocatable, intent( inout ) :: outE
        !
        integer :: ix, iy, iz

        !  NOTE: this only computes outputs for interior edges --  fine
        !    as long as observation locations are not in cell adjacent to boundary --
        !      but othewise should also use boundary edeges to compute H!
        !      Also would need to modify to compute Lrows, since don"t want to
        !         force on boundary!
        !
        !     modified so that there are no local cVectors allocated
        !
        !     do we really need this here?
        select type( inH )
        class is( cVector3D_SG_t )
            !
            !     this overwrites input inH    -- need to be aware of this in using!
            call inH%divs(self%Metric%FaceArea)
            !**
            ! Apply adjoint curl on unit grid
            !*
            !  EVERY OTHER ROUTINE IN MODEL OPERATOR ASSUMES THAT the output
            !     is allocated before calling -- going to follow this convention
            !      for this subroutine also!!!!
            ! allocate( outE, source = cVector3D_SG_t( inH%grid, EDGE ) )
            !
            if(.NOT.outE%is_allocated) then
                 write(*,*) "ERROR: multCurlT in    ModelOperator_MF"
                 stop    "output vector not allocated"
            endif

            select type( outE )
            class is( cVector3D_SG_t )
                !
                ! Ex
                do iy = 2, inH%Ny
                    do iz = 2, inH%Nz
                        outE%x(:, iy, iz) =    (inH%z(:, iy, iz) - &
                        inH%z(:, iy - 1, iz)) - &
                        (inH%y(:, iy, iz) - inH%y(:, iy, iz - 1))
                    end do
                end do
                !
                ! Ey
                do iz = 2, inH%Nz
                    do ix = 2, inH%Nx
                        outE%y(ix, :, iz) = (inH%x(ix, :, iz) - &
                        inH%x(ix, :, iz - 1)) - &
                        (inH%z(ix, :, iz) - inH%z(ix - 1, :, iz))
                    end do
                end do
                !
                ! Ez
                do ix = 2, inH%Nx
                    do iy = 2, inH%Ny
                        outE%z(ix,iy,:) = (inH%y(ix, iy, :) - &
                        inH%y(ix - 1, iy, :)) - &
                        (inH%x(ix, iy, :) - inH%x(ix, iy - 1, :))
                    end do
                end do
                !
            class default
                write( *, * ) "ERROR:ModelOperator_MF::multCurlT:"
                stop          "            Incompatible input [outE]"
            end select
            ! 
        class default
            write( *, * ) "ERROR:ModelOperator_MF::multCurlT:"
            stop          "            Incompatible input [inH]"
            !
        end select
        !
        !**
        ! Post multiply by edge length
        !*
        call outE%mults( self%metric%EdgeLength )
        !
    end subroutine multCurlTModelOperatorMF
    !**
    ! multiply by divergence correction operator.
    !    SEE NOTE BELOW ABOUT SIGN
    !*
    subroutine divCgradModelOperatorMF( self, inPhi, outPhi )
        implicit none
        !
        ! Arguments
        class( ModelOperator_MF_t ), intent( in ) :: self
        ! inphi, outphi have to be abstract to use in generic solver??
        class( cScalar_t ), intent( in )    :: inPhi
        class( cScalar_t ), intent( inout ) :: outPhi
        !
        integer :: ix, iy, iz
        !
        select type( outphi )
        class is( cScalar3D_SG_t )
            !
            if(.NOT.outPhi%is_allocated) then
                write( *, * ) "ERROR:ModelOperator_MF::divCgrad"
                stop          "         Output cScalar object not allocated"
            endif
            !
            select type( inPhi )
            class is( cScalar3D_SG_t )
                !
                ! zero output (already allocated) to start
                call outPhi%Zeros()
                !
                ! The coefficients are only for interior nodes
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
                        end do
                    end do
                end do
                ! NOTE: the coefficients db1, db2, c implement:
                !      phiOut = [div sigma grad] phiIn
                !  IN FACT this is a NEGATIVE DEFINITE OPERATOR!
                !    FLIP SIGN OF phiOut before output to make this 
                !     POSITIVE DEFINITE -- note that (a) solving with PCG
                !        seems to work in ModEM stable without this change, but
                !        this probably depends on specific preconditioner used
                !        (b)  preconditioner, and rhs need to also have signs flipped!
                !outPhi%v = outPhi%v
            end select
            !
        class default
            write( *, * ) "ERROR:ModelOperator_MF::amult:"
            stop          "            Incompatible input [x]. Exiting."
            !
        end select
        !
    end subroutine divCgradModelOperatorMF
    !**
    ! multiply by edge conductivity, apply
    ! divergence ==> source for divergence
    ! correction solver.
    !*
    subroutine divCModelOperatorMF( self, inE, outPhi )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        class( cVector_t ), intent( in )          :: inE
        class( cScalar_t ), intent( inout )       :: outPhi
        !
        integer :: ix, iy, iz
        !
        select type(outPhi)
        class is(cScalar3D_SG_t)
            !
            if( .NOT. outPhi%is_allocated) then
                write( *, * ) "ERROR:ModelOperator_MF::divC"
                stop          "         Output cScalar object not allocated"
            endif
            !
            select type( inE )
            class is ( cVector3D_SG_t )
                !
                call outPhi%Zeros()
                !
                ! Computation done only for internal nodes
                do ix = 2, outPhi%nx
                    do iy = 2, outPhi%ny
                        ! FOR NODES IN THE AIR ONLY
                        do iz = 2, outPhi%grid%nzAir
                            outPhi%v(ix, iy, iz) = &
                            SIGMA_AIR * (inE%x(ix, iy, iz) - inE%x(ix - 1, iy, iz)) * &
                            inE%grid%delXinv(ix) + &
                            SIGMA_AIR * (inE%y(ix, iy, iz) - inE%y(ix, iy - 1, iz)) * &
                            inE%grid%delYinv(iy) + &
                            SIGMA_AIR * (inE%z(ix, iy, iz) - inE%z(ix, iy, iz - 1)) * &
                            inE%grid%delZinv(iz)
                        end do
                        
                        ! FOR NODES AT THE AIR-EARTH INTERFACE
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

                        ! FOR NODES INSIDE THE EARTH ONLY
                        ! THE TOP MOST EARTH NODE HAS AN INTERFACE WITH
                        ! AIR, THEREFORE THAT ONE IS SKIPPED HERE
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
                        end do
                    end do
                end do
                !
            class default
                write( *, * ) "ERROR:ModelOperator_MF_t::divC:"
                STOP                "inE type unknow"
            end select
            !
        class default
            write( *, * ) "ERROR:ModelOperator_MF_t::divC:"
            STOP                "outPhi type unknow"
            !
        end select
        !
    end subroutine divCModelOperatorMF
    !**
    ! gradient
    !*
    subroutine gradModelOperatorMF( self, inPhi, outE )
        implicit none
        !
        class( ModelOperator_MF_t ), intent( in ) :: self
        class( cScalar_t ), intent( in )          :: inPhi
        class( cVector_t ), intent( inout )       :: outE
        !
        integer :: ix, iy, iz
        !
        select type( outE )
        class is(cVector3D_SG_t)
            !
            if(.NOT.outE%is_allocated) then
                write( *, * ) "ERROR:ModelOperator_MF::grad"
                stop          "         Output cVector object not allocated"
            endif
            !
            select type( inPhi )
            class is ( cScalar3D_SG_t )
                !
                call outE%Zeros
                ! Outputs only for interior edges
                do ix = 1, self%nx 
                    do iy = 2, self%ny
                        do iz = 2, self%nz
                            outE%x(ix, iy, iz) = (inPhi%v(ix + 1, iy, iz) - &
                            inPhi%v(ix, iy, iz)) / self%metric%grid%dx(ix)
                        end do
                    end do
                end do
                !
                do ix = 2, self%nx 
                    do iy = 1, self%ny
                        do iz = 2, self%nz
                            outE%y(ix, iy, iz) = (inPhi%v(ix, iy + 1, iz) - &
                            inPhi%v(ix, iy, iz)) / self%metric%grid%dy(iy)
                        end do
                    end do
                end do
                !
                do ix = 2, self%nx 
                    do iy = 2, self%ny
                        do iz = 1, self%nz    
                            outE%z(ix, iy, iz) = (inPhi%v(ix, iy, iz + 1) - &
                            inPhi%v(ix, iy, iz)) / self%metric%grid%dz(iz)
                        end do
                    end do
                end do
                !
            class default
                write( *, * ) "ERROR:ModelOperator_MF_t::grad:"
                STOP                "inPhi type unknow"
            end select
        class default
            write( *, * ) "ERROR:ModelOperator_MF_t::divC:"
            STOP                "outE type unknow"
        end select
        !
    end subroutine gradModelOperatorMF
    !**
    ! divergence, without multiplication by edge conductivity
    !*
    subroutine divModelOperatorMF( self, inE, outPhi )
        implicit none
          !
        ! Arguments
        class( ModelOperator_MF_t ), intent( in ) :: self
        class( cVector_t ), intent( in )          :: inE
        class( cScalar_t ), intent( inout )       :: outPhi
        !
        integer :: ix, iy, iz
        !
        select type(outPhi)
        class is(cScalar3D_SG_t)
            !
            if(.NOT.outPhi%is_allocated) then
                write( *, * ) "ERROR:ModelOperator_MF::div"
                stop          "         Output cScalar object not allocated"
            endif
            !
            select type( inE )
            class is(cVector3D_SG_t)
                !
                call outPhi%Zeros()
                !
                ! Computation done only for internal nodes
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
                        end do
                    end do
                end do
                class default
                write( *, * ) "ERROR:ModelOperator_MF_t:div:"
                STOP                "inE type unknow"
            end select
        class default
            write( *, * ) "ERROR:ModelOperator_MF_t:div:"
            STOP                "outPhi type unknow"
        end select
        !
    end subroutine divModelOperatorMF
    !
    subroutine printModelOperatorMF( self )
        implicit none
        !
        class( ModelOperator_MF_t ), intent(in) :: self
        !
        stop "subroutine print not implemented for ModelOperator_MF"
        !
    end subroutine printModelOperatorMF
    !
end module ModelOperator_MF
