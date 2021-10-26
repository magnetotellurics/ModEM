!**
! This specific version will only be used with matrix-free,
! which is only implemented for CSG.
!*
module Preconditioner_MF
  use Constants
  use Grid
  use cVector
  use cVector3D_SG
  use cScalar3D_SG
  use ModelOperator_MF
  use PreConditioner

  implicit none

  private
  
  public ::  PreConditioner_MF_t
  
  type, extends(PreConditioner_t) :: PreConditioner_MF_t
     
     real(kind = prec) :: omega
     class(cVector3D_SG_t), allocatable :: Dilu  !   diagonal of D-ILU used for 
                                                 !   preconditioner of curl-curl equation
     class(cScalar3D_SG_t), allocatable :: d     !   diagonal of D-ILU used for
                                                 !   preconditioner of divergence correction equation
      
     ! Pointer to model operator defining system
     ! of equations that this preconditoner
     class(ModelOperator_MF_t), pointer	:: ModOp
     
   contains
     ! Main routines used externally
     procedure, public :: Create
     procedure, public :: Allocate
     procedure, public :: DeAllocate
     procedure, public :: SetPreconditioner ! This needs to be called by Solver  object
                                            ! every time any part of the operator changes.
     procedure         :: SetPreconditioner_DC
     procedure, public :: LTSolve           ! These are left (M1) and right (M2)
     procedure, public :: UTSolve           ! preconditioning matrices for curl-curl equation.
     procedure, public :: divCGradILU
     
  end type PreConditioner_MF_t
  
  interface PreConditioner_MF_t
     module procedure PreConditioner_MF_ctor
  end interface PreConditioner_MF_t
  
contains
  
  !**
  ! Class constructor
  !*
  function PreConditioner_MF_ctor(ModOp) result(PreCond) 
    class(ModelOperator_MF_t), pointer :: ModOp
    type(PreConditioner_MF_t) :: Precond
    
    call PreCond%Create(ModOp)
    
  end function PreConditioner_MF_ctor
  
  !**
  ! Create
  !*
  subroutine Create(self, ModOp)
    class(PreConditioner_MF_t) :: self
    class(ModelOperator_MF_t), pointer :: ModOp
    
    self%ModOp => ModOp
    call self%Allocate()
  end subroutine Create
  
  !**
  ! Allocate
  !*
  subroutine Allocate(self)
    class(PreConditioner_MF_t), intent(inout) :: self
    
    if (.not.self%Dilu%isAllocated) then
       allocate(self%Dilu, source = cVector3D_SG_t(self%ModOp%grid, EDGE))
    else
       ! could check if Dilu is consistent with grid, and do nothing if so
       allocate(self%Dilu, source = cVector3D_SG_t(self%ModOp%grid, EDGE))
    end if
  end subroutine Allocate

  !**
  ! DeAllcoate
  !*
  subroutine DeAllocate(self)
    class(PreConditioner_MF_t) :: self
  end subroutine DeAllocate
  
  !**
  ! SetPreconditioner
  !*
  subroutine SetPreconditioner(self, omega)
    class(Preconditioner_MF_t) :: self
    real(kind = prec), intent(in) :: omega
    ! Local variables
    integer :: status
    integer :: ix, iy, iz
    
    ! Save omega in object, to record
    self%omega = omega
    
    ! Initialize the non-interior values
    ! only the interior edge values are really used
    self%Dilu%x(:,1,:) = C_ZERO
    self%Dilu%x(:,:,1) = C_ZERO
    self%Dilu%y(1,:,:) = C_ZERO
    self%Dilu%y(:,:,1) = C_ZERO
    self%Dilu%z(1,:,:) = C_ZERO
    self%Dilu%z(:,1,:) = C_ZERO
    
    ! Now set interior values
    do ix = 1, self%ModOp%grid%nx
       do iy = 2, self%ModOp%grid%ny
          do iz = 2, self%ModOp%grid%nz
             self%Dilu%x(ix, iy, iz) = self%ModOp%xXO(iy,iz) - &
                  C_ONE*omega*self%ModOp%Sigma_E%x(ix, iy, iz)  &
                  - self%ModOp%xXY(iy, 1)*self%ModOp%xXY(iy-1, 2) &
                  *self%Dilu%x(ix,iy-1,iz) &
                  - self%ModOp%xXZ(iz, 1)*self%ModOp%xXZ(iz-1, 2) &
                  *self%Dilu%x(ix,iy,iz-1)
             self%Dilu%x(ix, iy, iz) = ONE/self%Dilu%x(ix, iy, iz)
          end do
       end do
    end do
    
    ! The coefficients for y are only for the interior nodes
    ! but need to initialize edges for recursive algorithm.
    do iy = 1, self%ModOp%grid%ny
       do iz = 2, self%ModOp%grid%nz
          do ix = 2, self%ModOp%grid%nx
             self%Dilu%y(ix, iy, iz) = self%ModOp%yYO(ix,iz) - &
                  C_ONE*omega*self%ModOp%Sigma_E%y(ix, iy, iz) &
                  - self%ModOp%yYZ(iz, 1)*self%ModOp%yYZ(iz-1, 2) &
                  *self%Dilu%y(ix, iy, iz-1) &
                  - self%ModOp%yYX(ix, 1)*self%ModOp%yYX(ix-1, 2) &
                  *self%Dilu%y(ix-1, iy, iz)
             self%Dilu%y(ix, iy, iz) = ONE/self%Dilu%y(ix, iy, iz)
          end do
       end do
    end do
    
    ! The coefficients for z are only for the interior nodes
    ! but need to initialize edges for recursive algorithm.
    do iz = 1, self%ModOp%grid%nz
       do ix = 2, self%ModOp%grid%nx
          do iy = 2, self%ModOp%grid%ny
             self%Dilu%z(ix, iy, iz) = self%ModOp%zZO(ix,iy) - &
                  C_ONE*omega*self%Modop%Sigma_E%z(ix, iy, iz) &
                  - self%ModOp%zZX(ix, 1)*self%ModOp%zZX(ix-1, 2)*  &
                  self%Dilu%z(ix-1, iy, iz) &
                  - self%ModOp%zZY(iy, 1)*self%ModOp%zZY(iy-1, 2) &
                  *self%Dilu%z(ix, iy-1, iz)
             self%Dilu%z(ix, iy, iz) = ONE/self%Dilu%z(ix, iy, iz)
          end do
       end do
    end do
    
    ! This is for MF, with divergence correction -- if we implement
    ! modified equations in MF, will need to implement a new preconditioner
    ! anyway!
    ! BUT -- don't need to reset this when/if frequency changes (in contrast
    ! to Dilu).   So might remove this call, and make an explicit call where 
    ! needed.
    call self%setPreconditioner_DC()
    
  end subroutine SetPreConditioner
  
  !**
  ! Purpose: to solve the lower triangular system (or it's adjoint);
  ! for the d-ilu pre-condtioner.
  !*
  function LTsolve(self, inE, adjt) result(outE)
    class(PreConditioner_MF_t) :: self
    class(cVector_t), intent(in) :: inE
    logical         , intent(in) :: adjt
    ! Local variables
    integer :: ix, iy, iz
    class(cVector_t), allocatable :: outE
    
    if (.not.inE%isAllocated) then
       write(0,*) 'inE in M1solve not allocated yet'
       stop
    end if
    
    if (.not.outE%isAllocated) then
       write(0,*) 'outE in M1solve not allocated yet'
       stop
    end if

    allocate(outE, source = inE)
    call outE%Zeros()
    
    !   as usual I am cutting some of the error checking, which is not
    !   consistent with new classes
    
    if (.not.adjt) then
       ! adjoint = .false.
       !  we will need element/by element division (rdvide  in matlab)
       !            Call diagDiv(inE, V_E, outE)  !   this is ModEM routine
       !   I am assuming that this TVector function implements inE./Vedge
       outE = inE/self%ModOp%Metric%Vedge
       
       do ix = 1, inE%nx
          do iz = 2, inE%nz
             do iy = 2, inE%ny
                outE%x(ix, iy, iz) = (outE%x(ix, iy, iz) - &
                     outE%x(ix, iy-1, iz)*self%ModOp%xXY(iy, 1) - &
                     outE%x(ix, iy, iz-1)*self%ModOp%xXZ(iz, 1))* &
                     self%Dilu%x(ix, iy, iz)
             end do
          end do
       end do

       do iy = 1, inE%ny
          do iz = 2, inE%nz
             do ix = 2, inE%nx
                outE%y(ix, iy, iz) = (outE%y(ix, iy, iz) - &
                     outE%y(ix, iy, iz-1)*self%ModOp%yYZ(iz, 1) - &
                     outE%y(ix-1, iy, iz)*self%ModOp%yYX(ix, 1))* &
                     self%Dilu%y(ix, iy, iz)
             end do
          end do
       end do
       
       do iz = 1, inE%nz
          do iy = 2, inE%ny
             do ix = 2, inE%nx
                outE%z(ix, iy, iz) = (outE%z(ix, iy, iz) - &
                     outE%z(ix-1, iy, iz)*self%ModOp%zZX(ix, 1) - &
                     outE%z(ix, iy-1, iz)*self%ModOp%zZY(iy, 1))* &
                     self%Dilu%z(ix, iy, iz)
             end do
          end do
       end do
    else
       ! adjoint = .true.
       
       do ix = 1, inE%nx
          do iy = inE%ny, 2, -1
             do iz = inE%nz, 2, -1
                outE%x(ix, iy, iz) = (inE%x(ix, iy, iz) - &
                     outE%x(ix, iy+1, iz)*self%ModOp%xXY(iy+1, 1) - &
                     outE%x(ix, iy, iz+1)*self%ModOp%xXZ(iz+1, 1))* &
                     conjg(self%Dilu%x(ix, iy, iz))
             end do
          end do
       end do
       
       ! The coefficients for y are only for the interior nodes
       do iy = 1, inE%ny
          do ix = inE%nx, 2, -1
             do iz = inE%nz, 2, -1
                outE%y(ix, iy, iz) = (inE%y(ix, iy, iz) - &
                     outE%y(ix, iy, iz+1)*self%ModOp%yYZ(iz+1, 1) - &
                     outE%y(ix+1, iy, iz)*self%ModOp%yYX(ix+1, 1))* &
                     conjg(self%Dilu%y(ix, iy, iz))
             end do
          end do
       end do
       
       do iz = 1, inE%nz
          do ix = inE%nx, 2, -1
             do iy = inE%ny, 2, -1
                outE%z(ix, iy, iz) = (inE%z(ix, iy, iz) - &
                     outE%z(ix+1, iy, iz)*self%ModOp%zZX(ix+1, 1) - &
                     outE%z(ix, iy+1, iz)*self%ModOp%zZY(iy+1, 1))* &
                     conjg(self%Dilu%z(ix, iy, iz))
             end do
          end do
       end do

       ! IN ModEM we had allowed overwriting of first argument by output;
       ! IF WE DO NOT DO THIS, will have to make a copy of outE before
       ! dividing
       outE = outE/self%ModOp%Metric%Vedge
       
    end if
    
  end function LTsolve
  
  !**
  ! Purpose: to solve the upper triangular system (or it's adjoint);
  ! for the d-ilu pre-condtioner
  !*
  function UTsolve(self, inE, adjt) result(outE)
    class(PreConditioner_MF_t) :: self
    class(cVector_t), intent(in) :: inE
    logical         , intent(in) :: adjt
    ! Local variables
    integer :: ix, iy, iz
    class(cVector_t), allocatable :: outE

    allocate(outE, source = inE)
    call outE%Zeros()
    
    if (.not.adjt) then
       ! adjoint = .false.
       ! for standard upper triangular solution
       
       do ix = 1, inE%nx
          do iz = inE%nz, 2, -1
             do iy = inE%ny, 2, -1
                outE%x(ix, iy, iz) = inE%x(ix, iy, iz) - &
                     ( outE%x(ix, iy+1, iz)*self%ModOp%xXY(iy, 2) &
                     + outE%x(ix, iy, iz+1)*self%ModOp%xXZ(iz, 2))* &
                     self%Dilu%x(ix, iy, iz)
             end do
          end do
       end do
       
       do iy = 1, inE%ny
          do iz = inE%nz, 2, -1
             do ix = inE%nx, 2, -1
                outE%y(ix, iy, iz) = inE%y(ix, iy, iz) - &
                     ( outE%y(ix, iy, iz+1)*self%ModOp%yYZ(iz, 2) &
                     + outE%y(ix+1, iy, iz)*self%ModOp%yYX(ix, 2))* &
                     self%Dilu%y(ix, iy, iz)
             end do
          end do
       end do

       do iz = 1, inE%nz
          do iy = inE%ny, 2, -1
             do ix = inE%nx, 2, -1
                outE%z(ix, iy, iz) = inE%z(ix, iy, iz) - &
                     ( outE%z(ix+1, iy, iz)*self%ModOp%zZX(ix, 2) &
                     + outE%z(ix, iy+1, iz)*self%ModOp%zZY(iy, 2))* &
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
                     - outE%x(ix, iy-1, iz)*self%ModOp%xXY(iy-1, 2) &
                     * conjg(self%Dilu%x(ix,iy-1,iz))   &
                     - outE%x(ix, iy, iz-1)*self%ModOp%xXZ(iz-1, 2) &
                     * conjg(self%Dilu%x(ix, iy, iz-1))
             end do
          end do
       end do
       
       do iy = 1, inE%ny
          do iz = 2, inE%nz
             do ix = 2, inE%nx
                outE%y(ix, iy, iz) = inE%y(ix, iy, iz) &
                     - outE%y(ix, iy, iz-1)*self%ModOp%yYZ(iz-1, 2) &
                     * conjg(self%Dilu%y(ix,iy,iz-1)) &
                     - outE%y(ix-1, iy, iz)*self%ModOp%yYX(ix-1, 2) &
                     * conjg(self%Dilu%y(ix-1, iy, iz))
             end do
          end do
       end do
       
       do iz = 1, inE%nz
          do iy = 2, inE%ny
             do ix = 2, inE%nx
                outE%z(ix, iy, iz) = inE%z(ix, iy, iz) &
                     - outE%z(ix-1, iy, iz)*self%ModOp%zZX(ix-1, 2) &
                     * conjg(self%Dilu%z(ix-1,iy,iz)) &
                     - outE%z(ix, iy-1, iz)*self%ModOp%zZY(iy-1, 2) &
                     * conjg(self%Dilu%z(ix, iy-1, iz))
             end do
          end do
       end do
    end if
  end function UTsolve
  
  !************************************************************
  !     Divergence Correction Preconditioner Routines  
  !************************************************************
  
  !**
  ! Set up preconditioner for divergence correction equations
  ! div sigma grad phi = s
  ! where s is div sigma E_n   where E_n is trial solution to curl-curl
  ! This has to be called AFTER setting divergence correction operators
  ! with call to ModOp%DivCorSetup  ; this has to be called after changing
  ! conductivity parameter.
  !*
  subroutine SetPreconditioner_DC(self)
    class(Preconditioner_MF_t), intent(inout) :: self
    !  local variables
    integer :: ix,iy,iz
    
    ! Compute inverse diagonal elements for D-ILU (interior nodes only)
    ! set top nodes to 1.0
    self%d%v(1,:,:) = 1.0
    self%d%v(:,1,:) = 1.0
    self%d%v(:,:,1) = 1.0
    do iz = 2, self%ModOp%nz
       do iy = 2, self%ModOp%ny
          do ix = 2, self%ModOp%nx
             self%d%v(ix, iy, iz) = self%ModOp%c%v(ix, iy, iz) - &
                  self%ModOp%db1%x(ix,iy,iz)*self%ModOp%db2%x(ix-1,iy,iz) * &
                  self%d%v(ix-1,iy,iz)- &
                  self%ModOp%db1%y(ix,iy,iz)*self%ModOp%db2%y(ix,iy-1,iz) * &
                  self%d%v(ix,iy-1,iz)- &
                  self%ModOp%db1%z(ix,iy,iz)*self%ModOp%db2%z(ix,iy,iz-1) * &
                  self%d%v(ix,iy,iz-1)
             self%d%v(ix, iy, iz) = 1.0/ self%d%v(ix, iy, iz)
          end do
       end do
    end do
  end subroutine SetPreconditioner_DC
  
  !**
  ! apply pre-conditioner, solving lower and upper triangular systems using
  ! coefficients in db1, db2, and d.
  !*
  subroutine DivCgradILU(self, inPhi, outPhi)
    class(Preconditioner_MF_t), intent(in)    :: self
    type(cScalar3D_SG_t)      , intent(in)    :: inPhi
    type(cScalar3D_SG_t)      , intent(inout) :: outPhi
    ! Local variables
    integer :: ix, iy, iz
    
    outPhi%v = 0.0  !   should we have a routine to zero these objects?
    ! forward substitution (Solve lower triangular system)
    ! the coefficients are only for the interior nodes
    do iz = 2, inPhi%nz
       do iy = 2, inPhi%ny
          do ix = 2, inPhi%nx
             outPhi%v(ix, iy, iz) = inPhi%v(ix, iy, iz) &
                  - outPhi%v(ix-1,iy,iz)*self%ModOp%db1%x(ix,iy,iz)&
                  *self%d%v(ix-1,iy,iz) &
                  - outPhi%v(ix,iy-1,iz)*self%ModOp%db1%y(ix,iy,iz)&
                  *self%d%v(ix,iy-1,iz) &
                  - outPhi%v(ix,iy,iz-1)*self%ModOp%db1%z(ix,iy,iz)&
                  *self%d%v(ix,iy,iz-1)
          end do
       end do
    end do

    ! backward substitution (Solve upper triangular system)
    ! the coefficients are only for the interior nodes
    do iz = inPhi%nz,2,-1
       do iy = inPhi%ny,2,-1
          do ix = inPhi%nx,2,-1
             outPhi%v(ix, iy, iz) = (outPhi%v(ix, iy, iz)  &
                  - outPhi%v(ix+1, iy, iz)*self%ModOp%db2%x(ix, iy, iz)  &
                  - outPhi%v(ix, iy+1, iz)*self%ModOp%db2%y(ix, iy, iz)  &
                  - outPhi%v(ix, iy, iz+1)*self%ModOp%db2%z(ix, iy, iz)) &
                  *self%d%v(ix, iy, iz)
          end do
       end do
    end do
  end subroutine DivCgradILU
end module Preconditioner_MF

