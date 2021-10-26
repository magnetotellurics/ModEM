module ModelOperator_MF
  !
  use Constants
  use Grid3D_SG
  use rScalar3D_SG
  use cScalar3D_SG
  use cVector3D_SG
  use rVector
  use rVector3D_SG  
  use MetricElements_CSG
  use ModelParameter
  use ModelOperator
  !
  type, extends(ModelOperator_t) :: ModelOperator_MF_t
     
     !**
     ! Variables:
     !
     ! Only SG grid is supported for MF.
     !**
     class(MetricElements_CSG_t) , pointer :: metric => null() ! Pointer to metric element object.
     logical :: eqSet = .false.  ! Set to true after equations,
                                 ! (part independent of sigma) is set.
     integer(8) :: mKey = 0      ! For use with fortran DATE_AND_TIME subroutine.
     
     integer :: nx = 0, ny = 0, nz = 0  ! Redundant with what is in grid,
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
     real(kind = prec), allocatable, dimension(:,:) :: xXY, xXZ
     real(kind = prec), allocatable, dimension(:,:) :: xY, xZ
     real(kind = prec), allocatable, dimension(:,:) :: xXO
     real(kind = prec), allocatable, dimension(:,:) :: yYX, yYZ
     real(kind = prec), allocatable, dimension(:,:) :: yX, yZ
     real(kind = prec), allocatable, dimension(:,:) :: yYO
     real(kind = prec), allocatable, dimension(:,:) :: zZX, zZY
     real(kind = prec), allocatable, dimension(:,:) :: zX, zY
     real(kind = prec), allocatable, dimension(:,:) :: zZO
     
     ! Coefficients of diagonal of (unweighted) A operator.
     ! Here we store this as real - still need to multiply by i omega
     ! Thus, don't really need to store Adiag -- just store Sigma_E.
     type(rVector3D_SG_t) :: Sigma_E
     
     ! Operators for divergence correction (DC) will also be included here
     ! implementation will be trhough a separate module, managed through
     ! "solver" object.   Preconditioners will also be managed through the solver
     ! Might consider making a "base" version w/o DC, and have an extension that
     ! includes these arrayso.

     ! db1%x contains coefficients of the stencil for shift -1 of ix index
     ! (%y,%z give corresponding coefficients for iy, iz)
     ! db2  contains coefficients for corresponding shift of +1
     type(rVector3D_SG_t) :: db1, db2
     
     ! c contains the coefficients for div sigma grad
     ! operator diagonal. Note that divergence and gradient
     ! (also needed for DC) can be implemented
     ! directly using Metric Elements.
     type(rScalar3D_SG_t) :: c     
     
   contains
     private

     procedure, public :: UpdateFrequency
     procedure, public :: SetEquations
     procedure, public :: SetCond
     procedure, public :: Amult
     procedure, public :: MultAib
     procedure, public :: MultCurlT

     ! Private methods
     procedure :: create     => createModelOperatorMF 
     procedure :: allocate   => allocateModelOperatorMF
     procedure :: deallocate => deallocateModelOperatorMF

     procedure :: DivCorSetUp
     procedure :: divCgrad 
     procedure :: divC
     procedure :: Grad
     
  end type ModelOperator_MF_t
  
  interface ModelOperator_MF_t
     module procedure ModelOperator_MF_ctor
  end interface ModelOperator_MF_t
  
contains
  
  ! class constructor
  function ModelOperator_MF_ctor(inGrid) result(ModOp)
    class(Grid3D_SG_t), target, intent(in) :: inGrid
    type(ModelOperator_MF_t) :: ModOp
    !
    call ModOp%create( inGrid )
    !
  end function ModelOperator_MF_ctor
  
  !**
  ! Create -- since this just calls allocateModelOperatorMF.
  !*
  subroutine createModelOperatorMF(self, inGrid)
    ! Arguments
    class(ModelOperator_MF_t)  , intent(inout) :: self
    class(Grid3D_SG_t), target , intent(in)    :: inGrid
    
    self%grid => inGrid
    self%nx = inGrid%nx
    self%ny = inGrid%ny
    self%nz = inGrid%nz

    call self%allocate()

  end subroutine createModelOperatorMF
  
  !**
  ! allocateModelOperatorMF
  !*
  subroutine allocateModelOperatorMF(self)
    ! Arguments
    class(ModelOperator_MF_t), intent(inout) :: self
    ! Local variables
    integer :: status

    ! Allocate memory for del x del operator coefficient arrays
    ! Coefficients for difference equation only uses interior
    ! nodes. however, we need boundary nodes for the adjoint problem
    allocate(self%xXY(self%ny + 1, 2)     , STAT = status)
    allocate(self%xXZ(self%nz + 1, 2)     , STAT = status)   
    allocate(self%xY(self%nx, self%ny + 1), STAT = status)
    allocate(self%xZ(self%nx, self%nz + 1), STAT = status)
    allocate(self%xXO(self%ny, self%nz)   , STAT = status)

    allocate(self%yYZ(self%nz + 1, 2)     , STAT = status)   
    allocate(self%yYX(self%nx + 1, 2)     , STAT = status)   
    allocate(self%yZ(self%ny, self%nz + 1), STAT = status)
    allocate(self%yX(self%nx + 1, self%ny), STAT = status)
    allocate(self%yYO(self%nx, self%nz)   , STAT = status)

    allocate(self%zZX(self%nx + 1, 2)     , STAT = status)   
    allocate(self%zZY(self%ny + 1, 2)     , STAT = status)   
    allocate(self%zX(self%nx + 1, self%nz), STAT = status)
    allocate(self%zY(self%ny + 1, self%nz), STAT = status)
    allocate(self%zZO(self%nx, self%ny)   , STAT = status)

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

    select type(grid => self%grid)
    class is(Grid3D_SG_t)

      ! Initialize diagonal -- real Vector_t-- defined on edges
      self%Sigma_E = rVector3D_SG_t( grid, EDGE )

      ! Create Vectors, Scalar for divergence correction
      self%db1 = rVector3D_SG_t(grid, EDGE)
      self%db2 = rVector3D_SG_t(grid, EDGE)
      self%c   = rScalar3D_SG_t(grid, NODE)

      ! The metric element type for MF.
      ! Will always be CSG.
      !
      allocate(self%metric, source = MetricElements_CSG_t(grid))
       
    end select

    self%is_allocated = .true.

  end subroutine allocateModelOperatorMF
  
  !**
  ! deallocateModelOperatorMF
  !*
  subroutine deallocateModelOperatorMF(self)
    ! Arguments
    class(ModelOperator_MF_t), intent(inout) :: self
    ! Local variables
    integer :: status
    
    ! deallocateModelOperatorMF memory for del x del operator coefficient arrays
    ! Coefficients for difference equation only uses interior
    ! nodes. however, we need boundary nodes for the adjoint problem
    deallocate(self%xXY, STAT = status)
    deallocate(self%xXZ, STAT = status)
    deallocate(self%xY , STAT = status)
    deallocate(self%xZ , STAT = status)
    deallocate(self%xXO, STAT = status)
    
    deallocate(self%yYZ, STAT = status)
    deallocate(self%yYX, STAT = status)
    deallocate(self%yZ , STAT = status)
    deallocate(self%yX , STAT = status)
    deallocate(self%yYO, STAT = status)
    
    deallocate(self%zZX, STAT = status)
    deallocate(self%zZY, STAT = status)
    deallocate(self%zX , STAT = status)
    deallocate(self%zY , STAT = status)
    deallocate(self%zZO, STAT = status)
    
    self%is_allocated = .false.
  end subroutine deallocateModelOperatorMF

  subroutine updateFrequency(self, omega)
    class(ModelOperator_MF_t) :: self
    real(kind = prec), intent(in) :: omega

    self%omega = omega
  end subroutine updateFrequency
  
  !**
  ! SetEquations
  !*
  subroutine SetEquations(self)
    ! Arguments
    class(ModelOperator_MF_t), intent(inout) :: self
    ! Local variables
    integer :: status
    integer :: ix, iy, iz 
    
    ! Coefficients for curlcurlE (del X del X E)
    ! coefficents for calculating Ex ; only loop over internal edges
    do iy = 2, self%ny
       self%xXY(iy, 2) = -1.0 / (self%grid%delY(iy) * self%grid%dy(iy))
       self%xXY(iy, 1) = -1.0 / (self%grid%delY(iy) * self%grid%dy(iy-1))
    end do
    
    do iz = 2, self%nz
       self%xXZ(iz, 2) = -1.0 / (self%grid%delZ(iz) * self%grid%dz(iz))
       self%xXZ(iz, 1) = -1.0 / (self%grid%delZ(iz) * self%grid%dz(iz-1))
    end do
    
    do iy = 2, self%ny
       do iz = 2, self%nz
          self%xXO(iy, iz) = -(self%xXY(iy,1) + self%xXY(iy,2) + &
               self%xXZ(iz,1) + self%xXZ(iz,2))
       end do
    end do
    
    do ix = 1, self%nx
       do iy = 2, self%ny
          self%xY(ix, iy) = 1.0 / (self%grid%delY(iy)*self%grid%dx(ix))
       end do
    end do
    
    do ix = 1, self%nx
       do iz = 2, self%nz
          self%xZ(ix, iz) = 1.0 / (self%grid%delZ(iz)*self%grid%dx(ix))
       end do
    end do
    ! End of Ex coefficients
    
    ! Coefficents for calculating Ey; only loop over internal edges
    do iz = 2, self%nz
       self%yYZ(iz, 2) = -1.0 / (self%grid%delZ(iz)*self%grid%dz(iz))
       self%yYZ(iz, 1) = -1.0 / (self%grid%delZ(iz)*self%grid%dz(iz-1))
    end do
    
    do ix = 2, self%nx
       self%yYX(ix, 2) = -1.0 / (self%grid%delX(ix)*self%grid%dx(ix))
       self%yYX(ix, 1) = -1.0 / (self%grid%delX(ix)*self%grid%dx(ix-1))
    end do
    
    do ix = 2, self%nx
       do iz = 2, self%nz
          self%yYO(ix, iz) = -(self%yYX(ix,1) + self%yYX(ix,2) + &
               self%yYZ(iz,1) + self%yYZ(iz,2))
       end do
    end do

    do iy = 1, self%ny
       do iz = 2, self%nz
          self%yZ(iy, iz) = 1.0 / (self%grid%delZ(iz)*self%grid%dy(iy))
       end do
    end do
    
    do ix = 2, self%nx
       do iy = 1, self%ny
          self%yX(ix, iy) = 1.0 / (self%grid%delX(ix)*self%grid%dy(iy))
       end do
    end do
    ! End of Ey coefficients
    
    ! Coefficents for calculating Ez; only loop over internal edges
    do ix = 2, self%nx
       self%zZX(ix, 2) = -1.0 / (self%grid%delX(ix)*self%grid%dx(ix))
       self%zZX(ix, 1) = -1.0 / (self%grid%delX(ix)*self%grid%dx(ix-1))
    end do
    
    do iy = 2, self%ny
       self%zZY(iy, 2) = -1.0 / (self%grid%delY(iy)*self%grid%dy(iy))
       self%zZY(iy, 1) = -1.0 / (self%grid%delY(iy)*self%grid%dy(iy-1))
    end do
    
    do ix = 2, self%nx
       do iy = 2, self%ny
          self%zZO(ix, iy) = -(self%zZX(ix,1) + self%zZX(ix,2) + &
               self%zZY(iy,1) + self%zZY(iy,2))
       end do
    enddo

    do ix = 2, self%nx
       do iz = 1, self%nz
          self%zX(ix, iz) = 1.0 / (self%grid%delX(ix)*self%grid%dz(iz))
       end do
    end do
    
    do iy = 2, self%ny
       do iz = 1, self%nz
          self%zY(iy, iz) = 1.0 / (self%grid%delY(iy)*self%grid%dz(iz))
       end do
    end do
    ! End of Ez coefficients

    !**
    ! Equations are set, so change flag to true
    !*
    self%eqSet = .true.
    
    ! Note that divergence correction coefficients depend on conductivity
    ! and need to be reset whenever conductivity changes -- so these are
    ! not set here.
    
  end subroutine SetEquations
  
  !**
  ! SetCond
  !*
  subroutine SetCond(self, condParam)
    ! Arguments
    class(ModelOperator_MF_t), intent(inout) :: self
    class(ModelParameter_t)  , intent(inout)    :: condParam   
    ! Local variables
    complex(kind = prec) :: c
    !
    self%sigma_E = condParam%PDEmapping()
    !
    call self%DivCorSetUp()
    
  end subroutine SetCond
  
  !**
  ! DivcorSetup
  !*
  subroutine DivCorSetUp(self)
    ! Arguments
    class(ModelOperator_MF_t) :: self
    ! Local variables
    integer :: ix, iy, iz
    
    do iz = 1, self%grid%nzAir
       do iy = 1, self%ny
          do ix = 1, self%nx
             self%sigma_E%x(ix, iy, iz) = SIGMA_AIR
             self%sigma_E%y(ix, iy, iz) = SIGMA_AIR
             self%sigma_E%z(ix, iy, iz) = SIGMA_AIR
          end do
       end do
    end do
    
    ! The coefficients are only for the interior nodes
    ! these coefficients have not been multiplied by volume elements yet
    do iz = 2, self%nz
       do iy = 2, self%ny
          do ix = 2, self%nx
             self%db1%x(ix, iy, iz) = self%Sigma_E%x(ix - 1, iy, iz)/ &
                  (self%grid%dx(ix - 1)*self%grid%delX(ix))
             
             self%db2%x(ix, iy, iz) = self%Sigma_E%x(ix, iy, iz)/ &
                  (self%grid%dx(ix)*self%grid%delX(ix))
             
             self%db1%y(ix, iy, iz) = self%Sigma_E%y(ix, iy - 1, iz)/ &
                  (self%grid%dy(iy - 1)*self%grid%delY(iy))
             
             self%db2%y(ix, iy, iz) = self%Sigma_E%y(ix, iy, iz)/ &
                  (self%grid%dy(iy)*self%grid%delY(iy))
             
             self%db1%z(ix, iy, iz) = self%Sigma_E%z(ix, iy, iz - 1)/ &
                  (self%grid%dz(iz - 1)*self%grid%delZ(iz))
             
             self%db2%z(ix, iy, iz) = self%Sigma_E%z(ix, iy, iz)/ &
                  (self%grid%dz(iz)*self%grid%delZ(iz))
             
             self%c%v(ix, iy, iz) = - (self%db1%x(ix, iy, iz) + &
                  self%db2%x(ix, iy, iz) + &
                  self%db1%y(ix, iy, iz) + &
                  self%db2%y(ix, iy, iz) + &
                  self%db1%z(ix, iy, iz) + &
                  self%db2%z(ix, iy, iz))
          end do
       end do
    end do
    
    ! Multiply by corner volume elements to make operator symmetric
    do iz = 2, self%nz
       do iy = 2, self%ny
          do ix = 2, self%nx
             self%db1%x(ix, iy, iz) = self%db1%x(ix, iy, iz) * &
                  self%Metric%Vnode%v(ix,iy,iz)
             
             self%db1%y(ix, iy, iz) = self%db1%y(ix, iy, iz) * &
                  self%Metric%Vnode%v(ix,iy,iz)
             
             self%db1%z(ix, iy, iz) = self%db1%z(ix, iy, iz) * &
                  self%Metric%Vnode%v(ix,iy,iz)
             
             self%db2%x(ix, iy, iz) = self%db2%x(ix, iy, iz) * &
                  self%Metric%Vnode%v(ix,iy,iz)
             
             self%db2%y(ix, iy, iz) = self%db2%y(ix, iy, iz) * &
                  self%Metric%Vnode%v(ix,iy,iz)
             
             self%db2%z(ix, iy, iz) = self%db2%z(ix, iy, iz) * &
                  self%Metric%Vnode%v(ix,iy,iz)
          end do
       end do
    end do

    self%c = self%c * self%Metric%Vnode
    
    ! To be explicit set coefficients that multiply edges connected to
    ! boundary nodes to zero (this gaurantees that the BC on the potential
    ! is phi = 0):
    self%db1%x(2, :, :) = R_ZERO
    self%db1%y(:, 2, :) = R_ZERO
    self%db1%z(:, :, 2) = R_ZERO
    self%db2%x(self%nx, :, :) = R_ZERO
    self%db2%y(:, self%ny, :) = R_ZERO
    self%db2%z(:, :, self%nz) = R_ZERO

  end subroutine DivCorSetup
  
  !**
  ! AMult
  ! This does the matrix-vector multiply (A+iwB)inE = outE
  ! parameterized by input frequency omega
  !*
  function Amult(self, x, p_adjt) result(y)
    ! Arguments
    class(ModelOperator_MF_t), intent(in) :: self
    class(cVector_t)         , intent(in) :: x
    logical                  , intent(in), optional :: p_adjt    
    ! Local variables
    class(cVector_t), allocatable :: y
    ! Local variables
    integer :: ix, iy, iz
    complex(kind = prec) :: c
    class(cVector3D_SG_t), allocatable :: workE1, workE2
    logical :: adjt

    if (present(p_adjt)) then
       adjt = p_adjt
    else
       adjt = .false.
    end if
    
    select type(x)
    class is(cVector3D_SG_t)
       allocate(y, source = x)
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
                        self%xXO(iy, iz) * x%x(ix, iy, iz)
                end do
             end do
          end do
          
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
                        self%yYO(ix, iz) * x%y(ix, iy, iz)
                end do
             end do
          end do
          
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
                        self%zZO(ix, iy) * x%z(ix, iy, iz)
                end do
             end do
          end do
          
          ! Add diagonal part, multiplied by i*omega (and ISIGN MU_0)
          ! mult needs to know how to muliply complex vector by real
          ! first need to create WorkE1, WorkE2.
          select type(grid => self%grid)
          class is(Grid3D_SG_t)
             allocate(WorkE1, source = cVector3D_SG_t(grid, EDGE))
             allocate(WorkE2, source = cVector3D_SG_t(grid, EDGE))
          end select
          
          if (adjt) then
             c = -C_ONE*self%omega*ISIGN*MU_0
          else
             c = C_ONE*self%omega*ISIGN*MU_0
          end if
          
          WorkE1 = x *self%Sigma_E
          WorkE2 = c*WorkE1 + C_ONE*y
          
          ! Finally multiply by VEdge
          y = WorkE2 * self%metric%Vedge
       end select
    class default
       write(*, *) 'ERROR:ModelOperator_MF::Amult:'
       write(*, *) '      Incompatible input [x]. Exiting.'
       
       STOP
    end select
  end function Amult
  
  !**
  ! MultAib
  ! This multiplies boundary values by the real matrix Aib -- converts
  ! boundary data to interior forcing.
  ! NOTE: Make sure interior edges are zeroed before calling
  ! this function.
  !*
  function MultAib(self, bdry) result(outE)
    ! Arguments
    class(ModelOperator_MF_t), intent(in)  :: self
    class(cVector_t)         , intent(in)  :: bdry
    ! Local variables
    class(cVector_t), allocatable :: outE
    ! Local variables
    logical :: adjt    
    
    adjt = .false. 
    outE = self%AMult(bdry, adjt)
    
  end function multAib
  
  !**
  ! This multiplies Vector defined on faces (hence inH) by the transpose of
  ! the curl operator (which normally maps from edges to faces).   Not same
  ! as adjoint curl, on a non-uniform grid -- implemented using metric elements
  ! + matrix free adjoint curl on unit grid.
  ! curl.  This will enable simplified (and easily generalized to sparse
  ! matrix Model Operator implementations) computation of data functionals
  ! for computing magnetic fields.
  !*
  subroutine MultCurlT(self, inH, outE)
    implicit none
    ! Arguments
    class(ModelOperator_MF_t), intent(in) :: self
    class(cVector_t)         , intent(in) :: inH
    class(cVector_t)         , intent(out), allocatable :: outE
    ! Local variables
    integer :: ix, iy, iz
    type(cVector3D_SG_t) :: workH, workE

    select type(grid => self%grid)
    class is(Grid3D_SG_t)
       select type(inH)
       class is(cVector3D_SG_t)
          workE = cVector3D_SG_t(grid, EDGE)
          workH = cVector3D_SG_t(grid, FACE)
          
          workH = inH / self%Metric%FaceArea
          
          !**
          ! Apply adjoint curl on uint grid
          !*
          !
          ! Ex
          do iy = 2, workH%Ny
             do iz = 2, workH%Nz
                workE%x(:, iy, iz) =  (workH%z(:, iy, iz) - &
                     workH%z(:, iy - 1, iz)) - &
                     (workH%y(:, iy, iz) - workH%y(:, iy, iz - 1))
             end do
          end do
          
          ! Ey
          do iz = 2, workH%Nz
             do ix = 2, workH%Nx
                workE%y(ix, :, iz) = (workH%x(ix, :, iz) - &
                     workH%x(ix, :, iz - 1)) - &
                     (workH%z(ix, :, iz) - workH%z(ix - 1, :, iz))
             end do
          end do
          
          ! Ez
          do ix = 2, workH%Nx
             do iy = 2, workH%Ny
                workE%z(ix,iy,:) = (workH%y(ix, iy, :) - &
                     workH%y(ix - 1, iy, :)) - &
                     (workH%x(ix, iy, :) - workH%x(ix, iy - 1, :))
             end do
          end do
          !
          !**
          ! Post multiply by edge length
          !*
          allocate(outE, source = workE)
          
          outE = workE * self%Metric%EdgeLength
       class default
          write(*, *) 'ERROR:ModelOperator_MF::MultCurlT:'
          write(*, *) '      Incompatible input [inH]'
          
          STOP
       end select
    end select
  end subroutine MultCurlT
  
  !**
  ! Multiply by divergence correction operator.
  !*
  function divCgrad(self, inPhi) result(outPhi) 
    ! Arguments
    class(ModelOperator_MF_t), intent(in) :: self
    class(cScalar3D_SG_t)    , intent(in) :: inPhi
    ! Local variables
    class(cScalar3D_SG_t), allocatable :: outPhi    
    integer :: ix, iy, iz

    allocate(outPhi, source = inPhi)
    call outPhi%Zeros()
    
    ! The coefficients are only for interior nodes
    do iz = 2, inPhi%nz
       do iy = 2, inPhi%ny
          do ix = 2, inPhi%nx
             outPhi%v(ix, iy, iz) = inPhi%v(ix + 1, iy, iz) * self%db2%x(ix,iy,iz) + &
                  inPhi%v(ix - 1, iy, iz) * self%db1%x(ix, iy, iz) + &
                  inPhi%v(ix, iy + 1, iz) * self%db2%y(ix, iy, iz) + &
                  inPhi%v(ix, iy - 1, iz) * self%db1%y(ix, iy, iz) + &
                  inPhi%v(ix, iy, iz + 1) * self%db2%z(ix, iy, iz) + &
                  inPhi%v(ix, iy, iz - 1) * self%db1%z(ix, iy, iz) + &
                  inPhi%v(ix, iy, iz) * self%c%v(ix, iy, iz)
          end do
       end do
    end do
  end function divCgrad
  
  !**
  ! Multiply by edge conductivity, apply
  ! divergence ==> source for divergence
  ! correction solver.
  !*
  function divC(self, inE) result(outSc)
    ! Arguments
    class(ModelOperator_MF_t), intent(in) :: self
    class(cVector3D_SG_t)    , intent(in) :: inE
    ! Local variables
    type(cScalar3D_SG_t) :: outSc
    integer :: ix, iy, iz

    outSc = cScalar3D_SG_t(inE%grid, CORNER)
    call outSc%Zeros()
    
    ! Computation done only for internal nodes
    do ix = 2, outSc%nx
       do iy = 2, outSc%ny
          ! FOR NODES IN THE AIR ONLY
          do iz = 2, outSc%grid%nzAir
             outSc%v(ix, iy, iz) = &
                  SIGMA_AIR * (inE%x(ix, iy, iz) - inE%x(ix - 1, iy, iz)) * &
                  inE%grid%delXinv(ix) + &
                  SIGMA_AIR * (inE%y(ix, iy, iz) - inE%y(ix, iy - 1, iz)) * &
                  inE%grid%delYinv(iy) + &
                  SIGMA_AIR * (inE%z(ix, iy, iz) - inE%z(ix, iy, iz - 1)) * &
                  inE%grid%delZinv(iz)
          end do
          
          ! FOR NODES AT THE AIR-EARTH INTERFACE
          iz = outSc%grid%nzAir + 1
          outSc%v(ix, iy, iz) = &
               (self%Sigma_E%x(ix, iy, iz) * inE%x(ix, iy, iz) -     &
               self%Sigma_E%x(ix - 1, iy, iz) * inE%x(ix - 1, iy, iz)) * &
               inE%grid%delXinv(ix) + &
               (self%Sigma_E%y(ix, iy, iz) * inE%y(ix, iy, iz) -     &
               self%Sigma_E%y(ix, iy - 1, iz) * inE%y(ix, iy - 1, iz)) * &
               inE%grid%delYinv(iy) + &
               (self%Sigma_E%z(ix, iy, iz) * inE%z(ix, iy, iz) -     &
               SIGMA_AIR * inE%z(ix, iy, iz - 1)) * &
               inE%grid%delZinv(iz)
          
          ! FOR NODES INSIDE THE EARTH ONLY
          ! THE TOP MOST EARTH NODE HAS AN INTERFACE WITH
          ! AIR, THEREFORE THAT ONE IS SKIPPED HERE
          do iz = outSc%grid%nzAir + 2, outSc%nz
             outSc%v(ix, iy, iz) = &
                  (self%Sigma_E%x(ix,iy,iz)*inE%x(ix, iy, iz) -         &
                  self%Sigma_E%x(ix - 1,iy,iz)*inE%x(ix - 1, iy, iz)) * &
                  inE%grid%delXinv(ix)      &
                  +  (self%Sigma_E%y(ix,iy,iz)*inE%y(ix, iy, iz) -      &
                  self%Sigma_E%y(ix,iy - 1,iz)*inE%y(ix, iy - 1, iz)) * &
                  inE%grid%delYinv(iy)      &
                  +  (self%Sigma_E%z(ix,iy,iz)*inE%z(ix, iy, iz) -      &
                  self%Sigma_E%z(ix,iy,iz - 1)*inE%z(ix, iy, iz - 1)) * &
                  inE%grid%delZinv(iz)
          end do
       end do
    end do
  end function divC
  
  !**
  ! Gradient
  !*
  function Grad(self, inSc) result(outE)
    ! Arguments
    class(ModelOperator_MF_t), intent(in) :: self
    type(cScalar3D_SG_t)     , intent(in) :: inSc
    ! Local variables
    type(cVector3D_SG_t), allocatable :: outE       
    integer :: ix, iy, iz   
    
    ! Outputs only for interior edges
    do ix = 1, self%nx 
       do iy = 2, self%ny
          do iz = 2, self%nz
             outE%x(ix, iy, iz) = (inSc%v(ix + 1, iy, iz) - &
                  inSc%v(ix, iy, iz)) / self%grid%dx(ix)
          end do
       end do
    end do
    
    do ix = 2, self%nx 
       do iy = 1, self%ny
          do iz = 2, self%nz
             outE%y(ix, iy, iz) = (inSc%v(ix, iy + 1, iz) - &
                  inSc%v(ix, iy, iz)) / self%grid%dy(iy)
          end do
       end do
    end do
    
    do ix = 2, self%nx 
       do iy = 2, self%ny
          do iz = 1, self%nz  
             outE%z(ix, iy, iz) = (inSc%v(ix, iy, iz + 1) - &
                  inSc%v(ix, iy, iz)) / self%grid%dz(iz)
          end do
       end do
    end do
  end function Grad
  
end module ModelOperator_MF
