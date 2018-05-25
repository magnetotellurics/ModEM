! *****************************************************************************
! model_data is used for sharing the data for the joint forward modeling-
! inversion scheme. Data module where the current model definition (grid,
! conductivity, frequency) is stored. This module is sued by SetUp routines
! (for initialization, modification of any parameter values), and PDE coefficient
! initialization routines.
module modelOperator3D
  !  This module merges old modules model_data, model_data_update, multA,
  !      preconditioner, divcorr
  !   ... into a single module that does all operations involving model
  !     operators including the full EM operator and operators needed for
  !     divergence correction and preconditioning.  Everything that
  !     initializes or uses the equation coefficients are now in this module
  !     allowing arrays of equation coefficients, etc. to be private
  !       and essentially global to this module
  !
  !  A goal is to update this module to fully use the Curl, Div & Grad operators
  !  in the sg_diff_oper.f90 module, and the general grid element matrices.
  !  At present we are not doing this, since we found that doing so increases
  !  run time by 20%. Once optimized, will incorporate this logic in modelOperator3D.
  !  For now, we are including the operators module but keeping the logic intact.
  !  This version is just as efficient as the original stable version.
  !  Anna Kelbert, 14 May 2018.

  use math_constants
  use utilities
  use gridcalc             ! staggered grid definitions
  use sg_diff_oper         ! will be used for operators in the future but not in this version
  use sg_vector            ! generic routines for vector operations on the
  use sg_boundary
  use ModelSpace
  use boundary_ws          ! sets the boundary conditions
  use nestedEM
  implicit none

  ! * These variables are used by model equation
  ! * and preconditioner modules;
  ! * All variables are saved until changed/deallocated

  save

  !!!!!!!>>>>>>>>> conductivities in cells and edges (private)
  type(rvector), private    :: sigma_E ! edges
  type(rscalar), private    :: sigma_C ! cells - needed for boundary conditions

  !!!!!!!>>>>>>>>> FROM model_data
  type(grid_t), private, target 	::	mGrid   ! THE model grid
  !type(rvector), public			::	volE    ! THE volume elements
  !type(rvector), private		::	condE   ! THE edge conductivities
  real(kind=prec),private	::      omega   ! THE (active) frequency

  ! NOTE: THIS VARIABLE IS TEMPORARILY REQUIRED TO SET THE BOUNDARY CONDITIONS
  !type(rscalar), private        :: Cond3D

  !!!!!!>>>>>>>>> FROM multA
  ! aBC is for the Ea equation using the b component in the c direction.
  ! The two array elements correspond to coefficients required to form
  !   the derivative using adjacent faces
  ! E.g., xXY --> coefficient for the Ex equation using Ex variables
  ! in the y direction.
  ! xY is the product of grid spacings in X and Y-directions, etc.
  ! The two array elements again correspond to adjacent faces
  ! xXO is the sum of the all the products for the spacing in all the
  ! directions for adjacent faces for Ex term at ix, ij, ik point
  ! (collection of left/ right horizontal and top/ bottom vertical faces),
  !    etc.

  real (kind=prec), pointer, dimension(:,:,:), private    :: xXYm, xXZm   ! Lana
  real (kind=prec), pointer, dimension(:,:,:), private    :: xXYp, xXZp   ! Lana
!
! Extra 8 arrays for xY*, xZ*
  real (kind=prec), pointer, dimension(:,:,:), private    :: xYmm, xZmm
  real (kind=prec), pointer, dimension(:,:,:), private    :: xYmp, xZmp
  real (kind=prec), pointer, dimension(:,:,:), private    :: xYpm, xZpm
  real (kind=prec), pointer, dimension(:,:,:), private    :: xYpp, xZpp
!
  real (kind=prec), pointer, dimension(:,:,:), private    :: xXO 
!

  real (kind=prec), pointer, dimension(:,:,:), private    :: yYO
  real (kind=prec), pointer, dimension(:,:,:), private    :: yYXm, yYZm   ! Lana
  real (kind=prec), pointer, dimension(:,:,:), private    :: yYXp, yYZp   ! Lana
!
! Extra 8 arrays for yX*, yZ*
  real (kind=prec), pointer, dimension(:,:,:), private    :: yXmm, yZmm
  real (kind=prec), pointer, dimension(:,:,:), private    :: yXmp, yZmp
  real (kind=prec), pointer, dimension(:,:,:), private    :: yXpm, yZpm
  real (kind=prec), pointer, dimension(:,:,:), private    :: yXpp, yZpp
!
  real (kind=prec), pointer, dimension(:,:,:), private    :: zZO
  real (kind=prec), pointer, dimension(:,:,:), private    :: zZXm, zZYm   ! Lana
  real (kind=prec), pointer, dimension(:,:,:), private    :: zZXp, zZYp   ! Lana

! Extra 8 arrays for zX*, zY*
  real (kind=prec), pointer, dimension(:,:,:), private    :: zXmm, zYmm
  real (kind=prec), pointer, dimension(:,:,:), private    :: zXmp, zYmp
  real (kind=prec), pointer, dimension(:,:,:), private    :: zXpm, zYpm
  real (kind=prec), pointer, dimension(:,:,:), private    :: zXpp, zYpp

  ! coefficients of diagonal of (unweighted) A operator
  type (cvector), private                                :: Adiag
  ! information about the heritage ... probably this is not needed!
  real (kind=prec), private			:: whichOmega, whichCondE

  !!!!!!>>>>>>>>> FROM preconditioner
  ! coefficients of diagonal of preconditoner for A operator
  type (cvector), private		               :: Dilu

  !!!!!!>>>>>>>>> FROM divcorr
  ! coefficients for operators div sigma grad
  !  (which acts on scalars used for corner nodes),
  !  and the diagonal of the ilu preconditoner

  type (rvector) , private	:: db1, db2
  !   db1%x contains coefficients of the stencil for shift -1 of ix index
  !    (%y,%z give corresponding coefficients for iy, iz)
  !   db2  contains coefficients for corresponding shift of +1

  type (rscalar) , private	:: c, d
  ! c contains the coefficients for div sigma grad operator diagonal
  ! d contains the inverse of diagonal elements of D-ILU used for
  !  preconditioner


  ! *****************************************************************************
  !  routines from model_data_update:
  public                             	:: UpdateFreq, UpdateCond
  public                                :: UpdateFreqCond
  public                                :: ModelDataInit
  !   These are used to initialize the modules grid, and to set/modify
  !     conductivity and/or frequency

  !  routine to set the boundary conditions (a wrapper for BC_x0_WS for now)
  public                                :: ComputeBC

  !  routines from multA
  public                             	:: CurlcurleSetUp, CurlcurlE, CurlcurleCleanUp
  public                                :: AdiagInit, AdiagSetUp, deall_Adiag, Maxwell
  public                                :: MultA_O, MultA_N, AdjtBC
  !   These are used to initialize differential equation coefficients, and
  !    then to apply the differential operator

  ! routines from precondtioner
  public                      :: DiluInit, DiluSetUp, DeallocateDilu
  public                      :: M1Solve, M2Solve

  ! routines from divcorr
  public                :: DivCorrInit, DivCorrSetUp, Deallocate_DivCorr
  public                :: DivCgradILU, DivCgrad, DivC

  ! interface for data_update  ... is this needed ?
  INTERFACE updateModelData
     module procedure UpdateFreq
     module procedure UpdateCond
     module procedure UpdateFreqCond
  END INTERFACE

Contains

!**********************************************************************
  subroutine ModelDataInit(inGrid)
  !**********************************************************************
  ! *   Copies grid to mGrid
  !      and/or compute variables stored in model_data module:
  !         create: allocate for edge conductivity
  !              volume weights;
  !         EdgeVolume:  compute volume weights for edge-centered prisms
  !
  !**********************************************************************

    implicit none
    !  INPUTS:
    type (grid_t), intent(in)		  :: inGrid

    !   copy inGrid to mGrid
    call copy_grid(mGrid,inGrid)

    ! Allocate data structure for volume elements, and compute these
    Call create(mGrid, V_E, EDGE)

    ! Use the grid (which, potentially, maybe have been updated!) to set up
    !   all the grid length, surface and volume elements stored in GridCalc.
    ! Want to initialize them here in case the grid gets updated along the way.
    ! The reason for storing them in GridCalc is that they are also used
    !   by ModelMap, EMfieldInterp, nestedEM
    Call EdgeLength(mGrid, l_E)
    Call EdgeArea(mGrid, S_E)
    Call EdgeVolume(mGrid, V_E, l_E, S_E)
    Call FaceLength(mGrid, l_F)
    Call FaceArea(mGrid, S_F)
    Call FaceVolume(mGrid, V_F, l_F, S_F)
    Call CellVolume(mGrid, V_C)
    Call NodeVolume(mGrid, V_N) ! used for divergence correction

    !  Allocate sigma_E, conductivity defined on computational cell edges
    !   sigma_E is also needed in EMfieldInterp
    !!! Call create(mGrid, sigma_E, EDGE)
    Call create_rvector(mGrid, sigma_E, EDGE)

    ! set a default omega
    omega = 0.0

  end subroutine ModelDataInit


  subroutine ModelDataCleanUp

    ! Deallocated the grid
    call deall_grid(mGrid)

    ! and the grid elements stored in GridCalc
    call deall_rvector(l_E)
    call deall_rvector(S_E)
    call deall_rvector(V_E)
    call deall_rvector(l_F)
    call deall_rvector(S_F)
    call deall_rvector(V_F)
    call deall_rscalar(V_C)
    call deall_rscalar(V_N)

    ! and the edge conductivities
    call deall_rvector(sigma_E)

    ! and the cell conductivities
    ! note that sigma_C is only needed to set up boundary conditions
    !   until we set up a better way to do this
    call deall_rscalar(sigma_C)

  end subroutine ModelDataCleanUp

  ! **************************************************************************
  ! * UpdateFreq updates the frequency that is currently being use
  subroutine UpdateFreq(inOmega)

    implicit none
    real (kind=prec), intent (in)             :: inOmega

    omega = inOmega
    Call AdiagSetUp()
    Call DiluSetUp()

  end subroutine UpdateFreq  ! UpdateFreq

  ! ***************************************************************************
  ! * UpdateCond _updates the conductivity values on the edges
  subroutine UpdateCond(CondParam)

    implicit none
    type(modelParam_t), intent(in)      :: CondParam      ! input conductivity

    ! structure on the center of the grid

    !  ModelParamToEdge is to be interpreted as an abstract routine
    !    that maps from the external conductivity parameter to the
    !    internal edge representation  ... the type of CondParam
    !    is now fixed as rscalar;  if a different representation is
    !    to be used changes to the declarations in this routine will
    !    be required, along with changes in the module interface
    Call ModelParamToEdge(CondParam, sigma_E)

    Call DivCorrSetUp()

    ! TEMPORARY; REQUIRED FOR BOUNDARY CONDITIONS
    !  set static array for cell conductivities
    call ModelParamToCell(CondParam, sigma_C)

  end subroutine UpdateCond  ! UpdateCond

  ! ***************************************************************************
  ! * UpdateFreqCond updates the frequency that is currently being use and
  !*  conductivity values on the edges
  subroutine UpdateFreqCond(inOmega, CondParam)

    implicit none
    real(kind=prec)                 :: inOmega
    type(modelParam_t), intent(in)            :: CondParam      ! input conductivity
    ! structure on the center of the grid

    omega = inOmega

    !  ModelParamToEdge is to be interpreted as an abstract routine
    !    that maps from the external conductivity parameter to the
    !    internal edge representation  ...
    Call ModelParamToEdge(CondParam, sigma_E)

    Call AdiagSetUp()
    Call DiluSetUp()
    Call DivCorrSetUp()

    ! TEMPORARY; REQUIRED FOR BOUNDARY CONDITIONS
    !  set static array for cell conductivities
    call ModelParamToCell(CondParam,sigma_C)

  end subroutine UpdateFreqCond  ! UpdateFreqCond

!**********************************************************************
! Sets boundary conditions. Currently a wrapper for BC_x0_WS.
! Uses input 3D conductivity in cells sigma_C, that has to be initialized
! by updateCond before calling this routine. Also uses mGrid set by
! ModelDataInit. Uses omega, which is set by updateFreq.
! We always run this after setting the private variable omega, anyway.
  Subroutine ComputeBC(iTx,imode,E0,BC)

    !  Input mode, period
    integer, intent(in)     :: imode
    integer, intent(in)     :: iTx

    ! local variable
    real(kind=prec) :: period

    ! Output electric field first guess (for iterative solver)
    type(cvector), intent(inout)    :: E0
    ! Output boundary conditions
    type(cboundary), intent(inout)  :: BC

    period = (2*PI)/omega ! period is seconds

    ! Compute E0 using Weerachai 2D approach; can get BC from that
    call BC_x0_WS(imode,period,mGrid,sigma_C,E0,BC)
    
    !call getBC(E0,BC)
   
    ! Cell conductivity array is no longer needed
    ! NOT TRUE: needed for imode=2
    ! call deall_rscalar(sigma_C)

  end subroutine ComputeBC

!**********************************************************************
! Sets boundary conditions. Currently a wrapper for BC_x0_WS.
! Uses input 3D conductivity in cells sigma_C, that has to be initialized
! by updateCond before calling this routine. Also uses mGrid set by
! ModelDataInit. Uses omega, which is set by updateFreq.
! We always run this after setting the private variable omega, anyway.
!  Subroutine SetBound(imode,E0,BC,iTx)
!
!    !  Input mode, period
!    integer, intent(in)		:: imode
!    integer, intent(in)		:: iTx
!
!    ! local variable
!    real(kind=prec)	:: period
!
!    ! Output electric field first guess (for iterative solver)
!    type(cvector), intent(inout)	:: E0
!    ! Output boundary conditions
!    type(cboundary), intent(inout)	:: BC
!
!    period = (2*PI)/omega ! period is seconds
!
!    if (BC_AND_E0_FROM_FILE) then
!       ! we are going to make a huge assumption here: nPol == 1 always for this case
!       !  and of course transmitters are in same order always
!       E0 = E0_from_file(iTx)
!       call getBC(E0,BC)
!       !   do we now need to set boundary edges of E0 == 0?
!    else
!       if (BC%read_E_from_file) then
!
!          call BC_x0_WS(imode,period,mGrid,sigma_C,E0,BC)
!
!          ! The BC are already computed from a larger grid for all transmitters and modes and stored in BC_from_file.
!          ! Overwrite BC with BC_from_file.
!          ! Note: Right now we are using the same period layout for both grid.
!          ! This why, it is enough to know the period and mode index to pick up the BC from BC_from_file vector.
!          BC = BC_from_file((iTx*2)-(2-imode))
!       else
!         ! Compute the BC using Weerachai 2D approach
!          call BC_x0_WS(imode,period,mGrid,sigma_C,E0,BC)
!       end if
!   end if
!
!
!    ! Cell conductivity array is no longer needed
!    ! NOT TRUE: needed for imode=2
!    ! call deall_rscalar(sigma_C)
!
!  end subroutine SetBound

! ****************************************************************************
! Routines from multA; set up finite difference operator for quasi-static
!   frequency domain 3D EM induction equations, apply in various ways for forward,
!   adjoint Krylov solvers

  ! ***************************************************************************
  ! * CurlcurleSetUp sets up all the coefficients for finite difference
  ! * approximations for del X del X E. In SetUp routines, one may do memory
  ! * allocation inside. SetUp routines does basic calculations
  ! * (maybe one time deal or sometimes more than once)

  subroutine CurlcurleSetUp()

    implicit none
    ! Output coefficients for curlcurlE (del X del X E)
    integer                     :: status     ! for dynamic memory allocation
    integer                     :: ix, iy, iz ! dummy variables
    logical icheckarr

    ! Allocate memory for del x del operator coefficient arrays
    ! Coefficients for difference equation only uses interior
    ! nodes. however, we need boundary nodes for the adjoint problem
!
    allocate(xXYm(mGrid%nx,mGrid%ny+1,mGrid%nz), STAT=status)  
    allocate(xXZm(mGrid%nx,mGrid%ny,mGrid%nz+1), STAT=status) 
    allocate(xXYp(mGrid%nx,mGrid%ny+1,mGrid%nz), STAT=status) 
    allocate(xXZp(mGrid%nx,mGrid%ny,mGrid%nz+1), STAT=status) 
    allocate(xXO(mGrid%nx,mGrid%ny,mGrid%nz), STAT=status)  

! 8 extra arrays for xY*, xZ*
    allocate(xYmm(mGrid%nx, mGrid%ny+1,mGrid%nz), STAT=status)
    allocate(xZmm(mGrid%nx, mGrid%ny,mGrid%nz+1), STAT=status)
    allocate(xYmp(mGrid%nx, mGrid%ny+1,mGrid%nz), STAT=status)
    allocate(xZmp(mGrid%nx, mGrid%ny,mGrid%nz+1), STAT=status)
    allocate(xYpm(mGrid%nx, mGrid%ny+1,mGrid%nz), STAT=status)
    allocate(xZpm(mGrid%nx, mGrid%ny,mGrid%nz+1), STAT=status)
    allocate(xYpp(mGrid%nx, mGrid%ny+1,mGrid%nz), STAT=status)
    allocate(xZpp(mGrid%nx, mGrid%ny,mGrid%nz+1), STAT=status)
!
    allocate(yYZm(mGrid%nx,mGrid%ny,mGrid%nz+1), STAT=status)
    allocate(yYZp(mGrid%nx,mGrid%ny,mGrid%nz+1), STAT=status)  
    allocate(yYXm(mGrid%nx+1,mGrid%ny,mGrid%nz), STAT=status)   
    allocate(yYXp(mGrid%nx+1,mGrid%ny,mGrid%nz), STAT=status)  
!
! Extra 8 arrays for yX*, yZ*
    allocate(yZmm(mGrid%nx,mGrid%ny, mGrid%nz+1), STAT=status) 
    allocate(yXmm(mGrid%nx+1, mGrid%ny,mGrid%nz), STAT=status) 
    allocate(yZmp(mGrid%nx,mGrid%ny, mGrid%nz+1), STAT=status) 
    allocate(yXmp(mGrid%nx+1, mGrid%ny,mGrid%nz), STAT=status) 
    allocate(yZpm(mGrid%nx,mGrid%ny, mGrid%nz+1), STAT=status) 
    allocate(yXpm(mGrid%nx+1, mGrid%ny,mGrid%nz), STAT=status) 
    allocate(yZpp(mGrid%nx,mGrid%ny, mGrid%nz+1), STAT=status) 
    allocate(yXpp(mGrid%nx+1, mGrid%ny,mGrid%nz), STAT=status) 
    allocate(yYO(mGrid%nx, mGrid%ny, mGrid%nz), STAT=status) 
!
    allocate(zZXm(mGrid%nx+1,mGrid%ny,mGrid%nz), STAT=status)   
    allocate(zZYm(mGrid%nx,mGrid%ny+1,mGrid%nz), STAT=status)  
    allocate(zZXp(mGrid%nx+1,mGrid%ny,mGrid%nz), STAT=status) 
    allocate(zZYp(mGrid%nx,mGrid%ny+1,mGrid%nz), STAT=status) 
    allocate(zZO(mGrid%nx, mGrid%ny,mGrid%nz), STAT=status)     
! extra 8 arrays for zX*, zY*
    allocate(zXmm(mGrid%nx+1,mGrid%ny, mGrid%nz), STAT=status)     
    allocate(zYmm(mGrid%nx,mGrid%ny+1, mGrid%nz), STAT=status)    
    allocate(zXmp(mGrid%nx+1,mGrid%ny, mGrid%nz), STAT=status)     
    allocate(zYmp(mGrid%nx,mGrid%ny+1, mGrid%nz), STAT=status)    
    allocate(zXpm(mGrid%nx+1,mGrid%ny, mGrid%nz), STAT=status)     
    allocate(zYpm(mGrid%nx,mGrid%ny+1, mGrid%nz), STAT=status)    
    allocate(zXpp(mGrid%nx+1,mGrid%ny, mGrid%nz), STAT=status)     
    allocate(zYpp(mGrid%nx,mGrid%ny+1, mGrid%nz), STAT=status)    
!
    ! initalize all coefficients to zero (some remain zero)
    xXYm=0;xXZm=0;xXYp=0;xXZp=0; 
    xYmm = 0.0;xYmp = 0.0;xYpm = 0.0;xYpp = 0.0; 
    xZmm = 0.0;xZmp = 0.0;xZpm = 0.0;xZpp = 0.0; 
    xXO = 0.0
    yYXm=0;yYZm=0;yYXp=0;yYZp=0; 
!
    yXmm = 0.0;yXmp = 0.0;yXpm = 0.0;yXpp = 0.0; 
    yZmm = 0.0;yZmp = 0.0;yZpm = 0.0;yZpp = 0.0; 
    yYO = 0.0

    zZXm=0;zZYm=0;zZXp=0;zZYp=0; 
    zYmm = 0.0;zYmp = 0.0;zYpm = 0.0;zYpp = 0.0; 
    zXmm = 0.0;zXmp = 0.0;zXpm = 0.0;zXpp = 0.0; 
    zZO = 0.0
!

    ! coefficents for calculating Ex ; only loop over internal edges
    do iy = 2, mGrid%ny
     do iz = 1, mGrid%nz                                    
        xXYm(:,iy,iz) = - (l_E%x(:,iy-1,iz)*l_F%z(:,iy-1,iz)) &
                           / (S_F%z(:,iy-1,iz)*S_E%x(:,iy,iz))
        xXYp(:,iy,iz) = - (l_E%x(:,iy+1,iz)*l_F%z(:,iy,iz)) &
                         / (S_F%z(:,iy,iz)*S_E%x(:,iy,iz))
     enddo
    enddo
!
    do iz = 2, mGrid%nz
      do iy = 1, mGrid%ny             
       xXZm(:,iy,iz) = - (l_E%x(:,iy,iz-1)*l_F%y(:,iy,iz-1)) &
                  / (S_F%y(:,iy,iz-1)*S_E%x(:,iy,iz))
       xXZp(:,iy,iz) = - (l_E%x(:,iy,iz+1)*l_F%y(:,iy,iz)) &
                  / (S_F%y(:,iy,iz)*S_E%x(:,iy,iz))
      enddo
    enddo

    do iy = 2, mGrid%ny                                    
       do iz = 2, mGrid%nz                                 
         xXO(:,iy,iz) = -(xXYm(:,iy,iz)/l_E%x(:,iy-1,iz) + xXYp(:,iy,iz)/l_E%x(:,iy+1,iz) + & 
               xXZm(:,iy,iz)/l_E%x(:,iy,iz-1) + xXZp(:,iy,iz)/l_E%x(:,iy,iz+1))* &            
               l_E%x(:,iy,iz)                                                                 
       enddo                                            
    enddo

    do ix = 1, mGrid%nx
       do iy = 2, mGrid%ny
         do iz =1, mGrid%nz
          xYmm(ix,iy,iz) = + (l_E%y(ix,iy-1,iz)*l_F%z(ix,iy-1,iz)) &
                  / (S_F%z(ix,iy-1,iz)*S_E%x(ix,iy,iz))
  
          xYmp(ix,iy,iz) = + (l_E%y(ix,iy,iz)*l_F%z(ix,iy,iz)) &  
                  / (S_F%z(ix,iy,iz)*S_E%x(ix,iy,iz))

          xYpm(ix,iy,iz) = + (l_E%y(ix+1,iy-1,iz)*l_F%z(ix,iy-1,iz)) & 
                  / (S_F%z(ix,iy-1,iz)*S_E%x(ix,iy,iz))

          xYpp(ix,iy,iz) = + (l_E%y(ix+1,iy,iz)*l_F%z(ix,iy,iz)) &
                  / (S_F%z(ix,iy,iz)*S_E%x(ix,iy,iz))
         enddo
       enddo
    enddo
!
    do ix = 1, mGrid%nx
       do iz = 2, mGrid%nz
         do iy = 1, mGrid%ny
          xZmm(ix,iy,iz) = + (l_E%z(ix,iy,iz-1)*l_F%y(ix,iy,iz-1)) &  
                  / (S_F%y(ix,iy,iz-1)*S_E%x(ix,iy,iz))           

          xZmp(ix,iy,iz) = + (l_E%z(ix,iy,iz)*l_F%y(ix,iy,iz)) &
                  / (S_F%y(ix,iy,iz)*S_E%x(ix,iy,iz))

          xZpm(ix,iy,iz) = + (l_E%z(ix+1,iy,iz-1)*l_F%y(ix,iy,iz-1)) &
                  / (S_F%y(ix,iy,iz-1)*S_E%x(ix,iy,iz))

          xZpp(ix,iy,iz) = + (l_E%z(ix+1,iy,iz)*l_F%y(ix,iy,iz)) &    
                  / (S_F%y(ix,iy,iz)*S_E%x(ix,iy,iz))
         enddo
       enddo
    enddo

    ! End of Ex coefficients
!=================================================================================
    ! coefficents for calculating Ey; only loop over internal edges
    do iz = 2, mGrid%nz
     do ix = 1, mGrid%nx  
       yYZm(ix,:,iz) = - (l_E%y(ix,:,iz-1)*l_F%x(ix,:,iz-1)) & 
                  / (S_F%x(ix,:,iz-1)*S_E%y(ix,:,iz))
       yYZp(ix,:,iz) = - (l_E%y(ix,:,iz+1)*l_F%x(ix,:,iz)) &   
                  / (S_F%x(ix,:,iz)*S_E%y(ix,:,iz))
     enddo        
    enddo
!
    do ix = 2, mGrid%nx
     do iz = 1, mGrid%nz
       yYXm(ix,:,iz) = - (l_E%y(ix-1,:,iz)*l_F%z(ix-1,:,iz)) & 
                  / (S_F%z(ix-1,:,iz)*S_E%y(ix,:,iz))    
       yYXp(ix,:,iz) = - (l_E%y(ix+1,:,iz)*l_F%z(ix,:,iz)) &   
                  / (S_F%z(ix,:,iz)*S_E%y(ix,:,iz))
     enddo
    enddo

    do ix = 2, mGrid%nx
       do iz = 2, mGrid%nz
          yYO(ix,:,iz) = -(yYXm(ix,:,iz)/l_E%y(ix-1,:,iz) + yYXp(ix,:,iz)/l_E%y(ix+1,:,iz) + &  
               yYZm(ix,:,iz)/l_E%y(ix,:,iz-1) + yYZp(ix,:,iz)/l_E%y(ix,:,iz+1))*l_E%y(ix,:,iz)  
       enddo
    enddo


    do iy = 1, mGrid%ny
       do iz = 2, mGrid%nz
        do ix = 1, mGrid%nx
          yZmm(ix,iy,iz) = + (l_E%z(ix,iy,iz-1)*l_F%x(ix,iy,iz-1)) &
                  / (S_F%x(ix,iy,iz-1)*S_E%y(ix,iy,iz))
          yZmp(ix,iy,iz) = + (l_E%z(ix,iy,iz)*l_F%x(ix,iy,iz)) &
                  / (S_F%x(ix,iy,iz)*S_E%y(ix,iy,iz))
          yZpm(ix,iy,iz) = + (l_E%z(ix,iy+1,iz-1)*l_F%x(ix,iy,iz-1)) &
                  / (S_F%x(ix,iy,iz-1)*S_E%y(ix,iy,iz))
          yZpp(ix,iy,iz) = + (l_E%z(ix,iy+1,iz)*l_F%x(ix,iy,iz)) &
                  / (S_F%x(ix,iy,iz)*S_E%y(ix,iy,iz))
        enddo
       enddo
    enddo

    do ix = 2, mGrid%nx
       do iy = 1, mGrid%ny
        do iz = 1, mGrid%nz
          yXmm(ix,iy,iz) = + (l_E%x(ix-1,iy,iz)*l_F%z(ix-1,iy,iz)) &
                  / (S_F%z(ix-1,iy,iz)*S_E%y(ix,iy,iz))
          yXmp(ix,iy,iz) = + (l_E%x(ix,iy,iz)*l_F%z(ix,iy,iz)) &
                  / (S_F%z(ix,iy,iz)*S_E%y(ix,iy,iz))
          yXpm(ix,iy,iz) = + (l_E%x(ix-1,iy+1,iz)*l_F%z(ix-1,iy,iz)) &
                  / (S_F%z(ix-1,iy,iz)*S_E%y(ix,iy,iz))
          yXpp(ix,iy,iz) = + (l_E%x(ix,iy+1,iz)*l_F%z(ix,iy,iz)) &
                  / (S_F%z(ix,iy,iz)*S_E%y(ix,iy,iz))
        enddo
       enddo
    enddo
!

    ! End of Ey coefficients
!=========================================================================
    ! coefficents for calculating Ez; only loop over internal edges
    do ix = 2, mGrid%nx
       do iy = 1, mGrid%ny
       zZXm(ix,iy,:) = - (l_E%z(ix-1,iy,:)*l_F%y(ix-1,iy,:)) & 
                  / (S_F%y(ix-1,iy,:)*S_E%z(ix,iy,:)) 
       zZXp(ix,iy,:) = - (l_E%z(ix+1,iy,:)*l_F%y(ix,iy,:)) &   
                  / (S_F%y(ix,iy,:)*S_E%z(ix,iy,:))
       enddo
    enddo


    do iy = 2, mGrid%ny
     do ix = 1, mGrid%nx
       zZYm(ix,iy,:) = - (l_E%z(ix,iy-1,:)*l_F%x(ix,iy-1,:)) & 
                  / (S_F%x(ix,iy-1,:)*S_E%z(ix,iy,:))
       zZYp(ix,iy,:) = - (l_E%z(ix,iy+1,:)*l_F%x(ix,iy,:)) &   
                  / (S_F%x(ix,iy,:)*S_E%z(ix,iy,:))
     enddo
    enddo       

    do ix = 2, mGrid%nx
       do iy = 2, mGrid%ny
          zZO(ix, iy,:) = -(zZXm(ix,iy,:)/l_E%z(ix-1,iy,:) + zZXp(ix,iy,:)/l_E%z(ix+1,iy,:) + &  
               zZYm(ix,iy,:)/l_E%z(ix,iy-1,:) + zZYp(ix,iy,:)/l_E%z(ix+1,iy,:))*l_E%z(ix,iy,:)   
       enddo
    enddo

    do ix = 2, mGrid%nx
       do iz = 1, mGrid%nz
        do iy = 1, mGrid%ny
          zXmm(ix,iy,iz) = + (l_E%x(ix-1,iy,iz)*l_F%y(ix-1,iy,iz)) &
                  / (S_F%y(ix-1,iy,iz)*S_E%z(ix,iy,iz))
          zXmp(ix,iy,iz) = + (l_E%x(ix,iy,iz)*l_F%y(ix,iy,iz)) &
                  / (S_F%y(ix,iy,iz)*S_E%z(ix,iy,iz))
          zXpm(ix,iy,iz) = + (l_E%x(ix-1,iy,iz+1)*l_F%y(ix-1,iy,iz)) &
                  / (S_F%y(ix-1,iy,iz)*S_E%z(ix,iy,iz))
          zXpp(ix,iy,iz) = + (l_E%x(ix,iy,iz+1)*l_F%y(ix,iy,iz)) &
                  / (S_F%y(ix,iy,iz)*S_E%z(ix,iy,iz))
        enddo
       enddo
    enddo


    do iy = 2, mGrid%ny
       do iz = 1, mGrid%nz
         do ix = 1, mGrid%nx
          zYmm(ix,iy,iz) = + (l_E%y(ix,iy-1,iz)*l_F%x(ix,iy-1,iz)) &
                  / (S_F%x(ix,iy-1,iz)*S_E%z(ix,iy,iz))
          zYmp(ix,iy,iz) = + (l_E%y(ix,iy,iz)*l_F%x(ix,iy,iz)) &
                  / (S_F%x(ix,iy,iz)*S_E%z(ix,iy,iz))
          zYpm(ix,iy,iz) = + (l_E%y(ix,iy-1,iz+1)*l_F%x(ix,iy-1,iz)) &
                  / (S_F%x(ix,iy-1,iz)*S_E%z(ix,iy,iz))
          zYpp(ix,iy,iz) = + (l_E%y(ix,iy,iz+1)*l_F%x(ix,iy,iz)) &
                  / (S_F%x(ix,iy,iz)*S_E%z(ix,iy,iz))
         enddo
       enddo
    enddo
    ! End of Ez coefficients

  end subroutine CurlcurleSetUp   ! CurlcurleSetUp


  ! ***************************************************************************
  ! * CurlcurlE computes the finite difference equation in complex vectors
  ! * for del X del X E. Note that the difference equation is only for interior
  ! * edges. However, it does use the contribution from the boundary edges. The
  ! * coefficients are calculated in CurlcurleSetUp. Remember, in the operators
  ! * that are used iterative fashion, output is always initialized outside
  subroutine CurlcurlE(inE, outE)

    implicit none
    type (cvector), intent(in)             :: inE
    ! input electrical field as complex vector
    type (cvector), target, intent(inout)  :: outE
    ! output electrical field as complex vector
    integer                                :: ix, iy, iz
! inserted by Lana 
    type (cvector)                           :: workF
    ! workE is the complex vector that is used as a work space
    call create_cvector(mGrid,workF,FACE)
! end insert
    ! dummy variables

    if (.not.inE%allocated) then
      write(0,*) 'inE in CurlcurlE not allocated yet'
      stop
    end if

    if (.not.outE%allocated) then
      write(0,*) 'outE in CurlcurlE not allocated yet'
      stop
    end if

    ! Check whether the bounds are the same
    if ((inE%nx == outE%nx).and.&
         (inE%ny == outE%ny).and.&
         (inE%nz == outE%nz)) then

       if ((inE%gridType == outE%gridType)) then
          call Curl(inE,workF)
          call Curl(workF,outE)        
       else
          write (0, *) 'CurlcurlE: not compatible usage for existing data types'
       end if
    else
       write(0, *) 'Error-complex vectors for CurlcurlE are not of same size'
    end if
    call deall(workF) 
  end subroutine CurlcurlE        ! CurlcurlE

  ! ***************************************************************************
  ! * CurlcurlE computes the finite difference equation in complex vectors
  ! * for del X del X E. Deallocate these vectors when they are no longer needed.
  subroutine CurlcurleCleanUp()

    implicit none

    integer                     :: status     ! for dynamic memory deallocation

    ! Deallocate memory for del x del operator coefficient arrays
    ! Coefficients for difference equation only uses interior
    ! nodes. however, we need boundary nodes for the adjoint problem
!
    deallocate(xXYm, STAT=status) 
    deallocate(xXZm, STAT=status) 
    deallocate(xXYp, STAT=status) 
    deallocate(xXZp, STAT=status) 
!
    deallocate(xYmm, STAT=status)
    deallocate(xZmm, STAT=status)
    deallocate(xYmp, STAT=status)
    deallocate(xZmp, STAT=status)
    deallocate(xYpm, STAT=status)
    deallocate(xZpm, STAT=status)
    deallocate(xYpp, STAT=status)
    deallocate(xZpp, STAT=status)
!
    deallocate(xXO)
!
    deallocate(yYZm, STAT=status) 
    deallocate(yYXm, STAT=status) 
    deallocate(yYZp, STAT=status) 
    deallocate(yYXp, STAT=status) 
    deallocate(yYO, STAT=status)
!
    deallocate(yZmm, STAT=status)
    deallocate(yXmm, STAT=status)
    deallocate(yZmp, STAT=status)
    deallocate(yXmp, STAT=status)
    deallocate(yZpm, STAT=status)
    deallocate(yXpm, STAT=status)
    deallocate(yZpp, STAT=status)
    deallocate(yXpp, STAT=status)
!
    deallocate(zZXm, STAT=status) 
    deallocate(zZYm, STAT=status) 
    deallocate(zZXp, STAT=status) 
    deallocate(zZYp, STAT=status) 
    deallocate(zZO, STAT=status)
!
    deallocate(zXmm, STAT=status)
    deallocate(zYmm, STAT=status)
    deallocate(zXmp, STAT=status)
    deallocate(zYmp, STAT=status)
    deallocate(zXpm, STAT=status)
    deallocate(zYpm, STAT=status)
    deallocate(zXpp, STAT=status)
    deallocate(zYpp, STAT=status)

  end subroutine CurlcurleCleanUp

  ! ***************************************************************************
  ! * AdiagInit initializes the memory for the diagonal nodes being added
  ! * with the imaginary part. Init routines mostly do memory allocation,
  ! * reading and setting up the data
  subroutine AdiagInit()

    implicit none

    Call create_cvector(mGrid, Adiag, EDGE)

  end subroutine AdiagInit


  ! ***************************************************************************
  ! * Adiag sets up the diagonal nodes with the imaginary part added to it
  ! * SetUp routines do basic calculations (maybe one time deal or more than one)
  subroutine AdiagSetUp()

    implicit  none
    integer                   :: ix, iy, iz       ! dummy variables

    if (.not.Adiag%allocated) then
      write(0,*) 'Adiag in AdiagSetUp not allocated yet'
      stop
    end if
!
! We inserted the loops below to leave boundaries set to zero
!
    do ix = 1, mGrid%nx
     do iy = 2, mGrid%ny
      do iz = 2, mGrid%nz
       Adiag%x(ix,iy,iz) = CMPLX(0.0, 1.0, 8)*omega*MU_0*sigma_E%x(ix,iy,iz)
      enddo
     enddo
    enddo

    do iy = 1, mGrid%ny
     do ix = 2, mGrid%nx
      do iz = 2, mGrid%nz
       Adiag%y(ix,iy,iz) = CMPLX(0.0, 1.0, 8)*omega*MU_0*sigma_E%y(ix,iy,iz)
      enddo
     enddo
    enddo

    do iz = 1, mGrid%nz
     do ix = 2, mGrid%nx
      do iy = 2, mGrid%ny
       Adiag%z(ix,iy,iz) = CMPLX(0.0, 1.0, 8)*omega*MU_0*sigma_E%z(ix,iy,iz)
      enddo
     enddo
    enddo

  end subroutine AdiagSetUp

  ! ***************************************************************************
  ! * deall_Adiag deallocates the memory for the diagonal nodes being added
  ! * with the imaginary part.
  subroutine deall_Adiag()

    implicit none

    Call deall_cvector(Adiag)

  end subroutine deall_Adiag

  ! ***************************************************************************
  ! * Maxwell computes the finite difference equation in complex vectors
  ! * for del X del X E +/- i*omega*mu*conductivity*E in unsymmetrical form.
  ! * Note that the difference equation is only for interior edges. However,
  ! * it does use the contribution from the boundary edges. The coefficients
  ! * are  calculated in CurlcurleSetUp. Remember, in the operators that are
  ! * used in iterative fashion, output is always initialized outside
  subroutine Maxwell(inE, adjt, outE)

    implicit none
    type (cvector), intent(in)               :: inE
    ! input electrical field as complex vector
    logical, intent (in)                     :: adjt
    type (cvector), target, intent(inout)    :: outE
    ! output electrical field as complex vector
    complex (kind=prec)                      :: diag_sign ! changed by Lana, was integer
    integer                                  :: ix, iy, iz
! inserted by Lana 
    type (cvector)                           :: workF,workE
    ! workE is the complex vector that is used as a work space
    logical old ! switch to old code
!
    old=.false.
    call create_cvector(mGrid,workF,FACE)
    call create_cvector(mGrid,workE,EDGE)
! end insert

    if (.not.inE%allocated) then
      write(0,*) 'inE in Maxwell not allocated yet'
      stop
    end if

    if (.not.outE%allocated) then
      write(0,*) 'outE in Maxwell not allocated yet'
      stop
    end if

    ! Check whether the bounds are the same
    if ((inE%nx == outE%nx).and.&
         (inE%ny == outE%ny).and.&
         (inE%nz == outE%nz)) then

       if ((inE%gridType == outE%gridType)) then

          if (adjt) then
             diag_sign = -1*ISIGN
          else
             diag_sign = ISIGN
          end if

          if(old)then  

          !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ix,iy,iz)

          ! Apply difference equation to compute Ex (only on interior nodes)
          ! the diagonal nodes have the imaginary component added
          !$OMP DO SCHEDULE(STATIC)
	      do iz = 2, inE%nz
             do iy = 2, inE%ny
                do ix = 1, inE%nx
                   outE%x(ix,iy,iz) =  & 
                        (xYpp(ix,iy,iz)*inE%y(ix+1,iy,iz)-&
                  	xYmp(ix,iy,iz)*inE%y(ix,iy,iz)-xYpm(ix,iy,iz)*inE%y(ix+1,iy-1,iz)&
                        +xYmm(ix,iy,iz)*inE%y(ix,iy-1,iz))+&
                  	(xZpp(ix,iy,iz)*inE%z(ix+1,iy,iz)-xZmp(ix,iy,iz)*inE%z(ix,iy,iz)&
                  	-xZpm(ix,iy,iz)*inE%z(ix+1,iy,iz-1)+xZmm(ix,iy,iz)*inE%z(ix,iy,iz-1))+&
                  	xXYp(ix,iy,iz)*inE%x(ix,iy+1,iz)+& 
                  	xXYm(ix,iy,iz)*inE%x(ix,iy-1,iz)+& 
                  	xXZp(ix,iy,iz)*inE%x(ix,iy,iz+1)+& 
                  	xXZm(ix,iy,iz)*inE%x(ix,iy,iz-1)+& 
                  	(xXO(ix,iy,iz)+diag_sign*Adiag%x(ix,iy,iz))*inE%x(ix,iy,iz)
                enddo
             enddo
          enddo
          !$OMP END DO NOWAIT

          ! Apply difference equation to compute Ey (only on interior nodes)
	      ! the diagonal nodes have the imaginary component added
          !$OMP DO SCHEDULE(STATIC)
          do iz = 2, inE%nz
             do iy = 1, inE%ny
                do ix = 2, inE%nx
                   outE%y(ix,iy,iz) = &
                        (yZpp(ix,iy,iz)*inE%z(ix,iy+1,iz)-&
                  	yZmp(ix,iy,iz)*inE%z(ix,iy,iz)-yZpm(ix,iy,iz)*inE%z(ix,iy+1,iz-1) & 
                        +yZmm(ix,iy,iz)*inE%z(ix,iy,iz-1))&
                  	+(yXpp(ix,iy,iz)*inE%x(ix,iy+1,iz)-yXpm(ix,iy,iz)*inE%x(ix,iy,iz) &
                  	-yXmp(ix,iy,iz)*inE%x(ix-1,iy+1,iz)+yXmm(ix,iy,iz)*inE%x(ix-1,iy,iz))+&
                  	yYZp(ix,iy,iz)*inE%y(ix,iy,iz+1)+&
                  	yYZm(ix,iy,iz)*inE%y(ix,iy,iz-1)+&
                  	yYXp(ix,iy,iz)*inE%y(ix+1,iy,iz)+&
                  	yYXm(ix,iy,iz)*inE%y(ix-1,iy,iz)+&
                  	(yYO(ix,iy,iz)+diag_sign*Adiag%y(ix,iy,iz))*inE%y(ix,iy,iz)
                enddo
             enddo
          enddo
          !$OMP END DO NOWAIT

          ! Apply difference equation to compute Ey (only on interior nodes)
	      ! the diagonal nodes have the imaginary component added
          !$OMP DO SCHEDULE(STATIC)
          do iz = 1, inE%nz
             do iy = 2, inE%ny
                do ix = 2, inE%nx
                   outE%z(ix,iy,iz) = &
                        (zXpp(ix,iy,iz)*inE%x(ix,iy,iz+1)-&
                  	zXpm(ix,iy,iz)*inE%x(ix,iy,iz)-zXmp(ix,iy,iz)*inE%x(ix-1,iy,iz+1) &
                        +zXmm(ix,iy,iz)*inE%x(ix-1,iy,iz))&
                  	+(zYpp(ix,iy,iz)*inE%y(ix,iy,iz+1)-zYpm(ix,iy,iz)*inE%y(ix,iy,iz)&
                  	-zYmp(ix,iy,iz)*inE%y(ix,iy-1,iz+1)+zYmm(ix,iy,iz)*inE%y(ix,iy-1,iz))+&
                  	zZXp(ix,iy,iz)*inE%z(ix+1,iy,iz)+&
                  	zZXm(ix,iy,iz)*inE%z(ix-1,iy,iz)+&
                  	zZYp(ix,iy,iz)*inE%z(ix,iy+1,iz)+&
                  	zZYm(ix,iy,iz)*inE%z(ix,iy-1,iz)+&
                  	(zZO(ix,iy,iz)+diag_sign*Adiag%z(ix,iy,iz))*inE%z(ix,iy,iz)
                enddo
             enddo
          enddo

          !$OMP END DO NOWAIT

          !$OMP END PARALLEL
         else ! New programming
          call Curl(inE,workF)
          call Curl(workF,outE)
          call diagMult_cvector(Adiag,inE,workE)
          call scMultAdd_cvector(diag_sign,workE,outE)  
         endif
       else
          write (0, *) ' Maxwell: not compatible usage for existing data types'
       end if
    else
       write(0, *) 'Error-complex vectors for Maxwell are not of same size'
    end if
    call deall(workF) 
    call deall(workE) 
  end subroutine Maxwell        ! Maxwell


  ! **************************************************************************
  ! * Gets the Maxwell's equation in the complete symmetrical form,
  ! * del X del X E +/- i*omega*mu*conductivity*E. E is the complex vector
  ! * defining the electrical field _O is to denote that this is the original
  ! * subroutine. Diagonally multiplied by weights for symmetry.
  subroutine MultA_O(inE, adjt, outE)

    implicit none
    type (cvector), intent (in)              :: inE
    logical, intent (in)                     :: adjt
    type (cvector), intent (inout)           :: outE
    type (cvector)                           :: workE
    ! workE is the complex vector that is used as a work space
    integer                                  :: diag_sign
    complex(kind=prec)               :: c2
    ! a complex multiplier

    if (.not.inE%allocated) then
      write(0,*) 'inE in MultA_O not allocated yet'
      stop
    end if

    if (.not.outE%allocated) then
      write(0,*) 'outE in MultA_O not allocated yet'
      stop
    end if

    Call create_cvector(mGrid, workE, EDGE)

    ! Check whether the bounds are the same
    if ((inE%nx == outE%nx).and.&
         (inE%ny == outE%ny).and.&
         (inE%nz == outE%nz).and.&
         (inE%nx == workE%nx).and.&
         (inE%ny == workE%ny).and.&
         (inE%nz == workE%nz)) then

       if ((inE%gridType == outE%gridType).and.(inE%gridType == workE%gridType)) &
            then

          Call CurlcurlE(inE, workE)
          ! done with preparing del X del X E

          if (adjt) then
             diag_sign = -1*ISIGN
          else
             diag_sign = ISIGN
          end if

          ! now preparing +/-i*omega*mu*conductivity*E
          Call diagMult_crvector(inE, sigma_E, outE)
          c2 = diag_sign*C_ONE*omega*MU_0
          Call linComb_cvector(C_ONE, workE, c2, outE, outE)

          ! diagonally multiply the final results with weights (edge volume)
          Call diagMult_crvector(outE, V_E, outE)

       else
          write (0, *) 'MultA_O: not compatible usage for existing data types'
       end if

    else
       write(0, *) 'Error-complex vectors for MultA_O are not of same size'
    end if

    Call deall(workE)

  end subroutine MultA_O


  ! ***************************************************************************
  ! * Gets the Maxwell's equation in the complete symmetrical form,
  ! * del X del X E +/- i*omega*mu*conductivity*E. E is the complex vector
  ! * defining the electrical field _N is to denote that this is the new
  ! * subroutine where the imaginary part at the at the diagonal is inbuilt
  ! *  Diagonally multiplied by weights for symmetry.
  subroutine MultA_N(inE, adjt, outE)

    implicit none
    type (cvector), intent (in)              :: inE
    logical, intent (in)                     :: adjt
    type (cvector), intent (inout)           :: outE


    if (.not.inE%allocated) then
      write(0,*) 'inE in MultA_N not allocated yet'
      stop
    end if

    if (.not.outE%allocated) then
      write(0,*) 'outE in MultA_N not allocated yet'
      stop
    end if

    ! Check whether the bounds are the same
    if ((inE%nx == outE%nx).and.&
         (inE%ny == outE%ny).and.&
         (inE%nz == outE%nz)) then

       if ((inE%gridType == outE%gridType)) then
          !!!write(0,*)maxval(abs(xXYp-xXY2)),maxval(abs(xXYm-xXY1))
          Call Maxwell(inE, adjt, outE)

          ! done with preparing del X del X E +/- i*omega*mu*conductivity*E

          ! diagonally multiply the final results with weights (edge volume)
          Call diagMult(outE, V_E, outE)

       else
          write (0, *) 'MultA_N: not compatible usage for existing data types'
       end if

    else
       write(0, *) 'Error-complex vectors for MultA_N are not of same size'
    end if

  end subroutine MultA_N

  subroutine AdjtBC(eIn, BC)
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

    implicit none

    ! INPUT: electrical fields stored as cvector
    type (cvector), intent(in)             		:: eIn
    ! OUTPUT: boundary condition structure: should be allocated
    !   and initialized before call to this routine
    type (cboundary),intent(inout)  			:: BC

    ! local variables
    integer                   :: ix,iy,iz,nx,ny,nz

    !  Multiply FD electric field vector defined on interior nodes (eIn) by
    !  adjoint of A_IB, the interior/boundary sub-block of the differential
    !  operator.

    nx = eIn%nx
    ny = eIn%ny
    nz = eIn%nz
    !write(0,*)'AdjtBC' ! It is never used for FWD?
!  Ex components in x/z plane (iy=1; iy = ny+1)
!  NOTE: C_ZERO = (0,0) (double complex) is defined in SG_Basics/math_constants.f90
    BC%xYMax(:,1) = C_ZERO
    BC%xYmin(:,1) = C_ZERO
    BC%xYMax(:,nz+1) = C_ZERO
    BC%xyMin(:,nz+1) = C_ZERO
! CHANGES BELOW WERE NOT DEBUGGED since AdjtBC is never called by FWD
    do ix = 1, nx
       do iz = 2, nz
          BC%xYmin(ix,iz) = - yXmp(ix,1,iz)*Ein%y(ix,1,iz)       & 
                            + yXmm(ix+1,1,iz)*Ein%y(ix+1,1,iz)   & 
                            + xXYm(ix,2,iz)*Ein%x(ix,2,iz)          
          BC%xYmax(ix,iz) = + yXpp(ix,ny,iz)*Ein%y(ix,ny,iz)     &  
                            - yXpm(ix+1,ny,iz)*Ein%y(ix+1,ny,iz) &  
                            + xXYp(ix,ny,iz)*Ein%x(ix,ny,iz)           
 
        enddo
     enddo

!  Ez components in x/z plane (iy=1; iy = ny+1)
    BC%zYMin(1,:) = C_ZERO
    BC%zYmax(1,:) = C_ZERO
    BC%zYmin(nx+1,:) = C_ZERO
    BC%zYmax(nx+1,:) = C_ZERO
    do iz = 1, nz
       do ix = 2, nx
          BC%zYmin(ix,iz) = - yZmp(ix,1,iz)*Ein%y(ix,1,iz)        &  
                            + yZmm(ix,1,iz+1)*Ein%y(ix,1,iz+1)    &  
                            + zZYm(ix,2,iz)*Ein%z(ix,2,iz)            
          BC%zYmax(ix,iz) = + yZpp(ix,ny,iz)*Ein%y(ix,ny,iz)      &  
                            - yZpm(ix,ny,iz+1)*Ein%y(ix,ny,iz+1)  &  
                            + zZYp(ix,ny,iz)*Ein%z(ix,ny,iz)                
        enddo
     enddo

!  Ey components in y/z plane (ix=1; ix = nx+1)
    BC%yXmin(:,1) = C_ZERO
    BC%yXmax(:,1) = C_ZERO
    BC%yXmin(:,nz+1) = C_ZERO
    BC%yXmax(:,nz+1) = C_ZERO
    do iy = 1, ny
       do iz = 2, nz
          BC%yXmin(iy,iz) = - xYmp(1,iy,iz)*Ein%x(1,iy,iz)        & 
                            + xYmm(1,iy+1,iz)*Ein%x(1,iy+1,iz)    & 
                            + yYXm(2,iy,iz)*Ein%y(2,iy,iz)           
          BC%yXmax(iy,iz) = + xYpp(nx,iy,iz)*Ein%x(nx,iy,iz)      & 
                            - xYpm(nx,iy+1,iz)*Ein%x(nx,iy+1,iz)  & 
                            + yYXp(nx,iy,iz)*Ein%y(nx,iy,iz)         
        enddo
     enddo

!  Ez components in y/z plane (ix=1; ix = nx+1)
    BC%zXmin(1,:) = C_ZERO
    BC%zXmax(1,:) = C_ZERO
    BC%zXmin(ny+1,:) = C_ZERO
    BC%zXmax(ny+1,:) = C_ZERO
    do iz = 1, nz
       do iy = 2, ny
          BC%zXmin(iy,iz) = - xZmp(1,iy,iz)*Ein%x(1,iy,iz)       &
                            + xZmm(1,iy,iz+1)*Ein%x(1,iy,iz+1)   &  
                            + zZXm(2,iy,iz)*Ein%z(2,iy,iz)          
          BC%zXmax(iy,iz) = + xZpp(nx,iy,iz)*Ein%x(nx,iy,iz)     & 
                            - xZpm(nx,iy,iz+1)*Ein%x(nx,iy,iz+1) &  
                            + zZXp(nx,iy,iz)*Ein%z(nx,iy,iz)        

        enddo
     enddo

!  Ex components in x/y plane (iz=1; iz = nz+1)
    BC%xZmin(:,1) = C_ZERO
    BC%xZmax(:,1) = C_ZERO
    BC%xZmin(:,ny+1) = C_ZERO
    BC%xZmax(:,ny+1) = C_ZERO
    do ix = 1, nx
       do iy = 2, ny
          BC%xZmin(ix,iy) = - zXmp(ix,iy,1)*Ein%z(ix,iy,1)       & 
                            + zXmm(ix+1,iy,1)*Ein%z(ix+1,iy,1)   & 
                            + xXZm(ix,iy,2)*Ein%x(ix,iy,2)         
          BC%xZmax(ix,iy) = + zXpp(ix,iy,nz)*Ein%z(ix,iy,nz)     & 
                            - zXpm(ix+1,iy,nz)*Ein%z(ix+1,iy,nz) &  
                            + xXZp(ix,iy,nz)*Ein%x(ix,iy,nz)        
        enddo
     enddo

!  Ey components in x/y plane (iz=1; iz = nz+1)
    BC%yZmin(:,1) = C_ZERO
    BC%yZmax(:,1) = C_ZERO
    BC%yZmin(nx+1,:) = C_ZERO
    BC%yZmin(nx+1,:) = C_ZERO
    do iy = 1, ny
       do ix = 2, nx
          BC%yZmin(ix,iy) = - zYmp(ix,iy,1)*Ein%z(ix,iy,1)        & 
                            + zYmm(ix,iy+1,1)*Ein%z(ix,iy+1,1)    & 
                            + yYZm(ix,iy,2)*Ein%y(ix,iy,2)           
          BC%yZmax(ix,iy) = + zYpp(ix,iy,nz)*Ein%z(ix,iy,nz)      &
                            - zYpm(ix,iy+1,nz)*Ein%z(ix,iy+1,nz)  & 
                            + yYZp(ix,iy,nz)*Ein%y(ix,iy,nz)        
        enddo
     enddo

  end subroutine AdjtBC

! ****************************************************************************
! PRECONDITIONER ROUTINES: set up ILU-Level I preconditioner for
!     Maxwell's equation, solve lower and upper triangular systems to
!      apply preconditioner
  !****************************************************************************
  ! initializes a diagonal of preconditioner for A operator
  subroutine DiluInit()

    implicit none
    integer                                 :: status
    integer                                 :: ix, iy, iz

    if (.not.Dilu%allocated) then

       Call create(mGrid, Dilu, EDGE)

    else

       if ((Dilu%nx /= mGrid%nx).or.(Dilu%ny /= mGrid%ny).or.&
            (Dilu%nz /= mGrid%nz)) then

          deallocate(Dilu%x, Dilu%y, Dilu%z, STAT = status)
          Call create(mGrid, Dilu, EDGE)

       end if
    end if

  end subroutine DiluInit ! DiluInit

  !****************************************************************************
  ! sets up a diagonal of preconditioner for A operator
  subroutine DiluSetUp()

    implicit none
    integer                                 :: status
    integer                                 :: ix, iy, iz

    if (.not.Dilu%allocated) then
       write (0, *) 'Dilu not allocated yet'
    else

       if ((Dilu%nx /= mGrid%nx).or.(Dilu%ny /= mGrid%ny).or.&
            (Dilu%nz /= mGrid%nz)) then
         write (0, *) 'Dilu that is right now existing has the wrong size'
       end if
    end if

    ! initializing the non-interior values
    ! only the interior edge values are really used
    Dilu%x(:,1,:) = cmplx(1.0, 0.0, 8)
    Dilu%x(:,:,1) = cmplx(1.0, 0.0, 8)

    Dilu%y(1,:,:) = cmplx(1.0, 0.0, 8)
    Dilu%y(:,:,1) = cmplx(1.0, 0.0, 8)

    Dilu%z(1,:,:) = cmplx(1.0, 0.0, 8)
    Dilu%z(:,1,:) = cmplx(1.0, 0.0, 8)

    do ix = 1, mGrid%nx
       do iy = 2, mGrid%ny
          do iz = 2, mGrid%nz

             Dilu%x(ix, iy, iz) = xXO(ix,iy,iz) - & 
                  CMPLX(0.0, 1.0, 8)*omega*MU_0*sigma_E%x(ix, iy, iz)  &
                  - xXYm(ix,iy,iz)*xXYp(ix,iy-1,iz)*Dilu%x(ix,iy-1,iz) &  
                  - xXZm(ix,iy,iz)*xXZp(ix,iy,iz-1)*Dilu%x(ix,iy,iz-1)

             Dilu%x(ix, iy, iz) = 1.0/ Dilu%x(ix, iy, iz)

          enddo
       enddo
    enddo

    ! the coefficients for y are only for the interior nodes
    !  but need to initialize edges for recursive algorithm
    do iy = 1, mGrid%ny
       do iz = 2, mGrid%nz
          do ix = 2, mGrid%nx

             Dilu%y(ix, iy, iz) = yYO(ix,iy,iz) - &
                  CMPLX(0.0, 1.0, 8)*omega*MU_0*sigma_E%y(ix, iy, iz) &    
                  - yYZm(ix,iy,iz)*yYZp(ix,iy,iz-1)*Dilu%y(ix, iy, iz-1) & 
                  - yYXm(ix,iy,iz)*yYXp(ix-1,iy,iz)*Dilu%y(ix-1, iy, iz)   

             Dilu%y(ix, iy, iz) = 1.0/ Dilu%y(ix, iy, iz)

          enddo
       enddo
    enddo

    ! the coefficients for z are only for the interior nodes
    !  but need to initialize edges for recursive algorithm
    do iz = 1, mGrid%nz
       do ix = 2, mGrid%nx
          do iy = 2, mGrid%ny

             Dilu%z(ix, iy, iz) = zZO(ix,iy,iz) - &
                  CMPLX(0.0, 1.0, 8)*omega*MU_0*sigma_E%z(ix, iy, iz) &
                  - zZXm(ix,iy,iz)*zZXp(ix-1,iy,iz)*Dilu%z(ix-1, iy, iz) &
                  - zZYm(ix,iy,iz)*zZYp(ix,iy-1,iz)*Dilu%z(ix, iy-1, iz)

             Dilu%z(ix, iy, iz) = 1.0/ Dilu%z(ix, iy, iz)

          enddo
       enddo
    enddo

  end subroutine DiluSetUp  ! DiluSetUp


  !****************************************************************************
  !  To Deallocate arrays in structure Dilu
  subroutine DeallocateDilu()
    implicit none
    integer                                 :: status

	call deall_cvector(Dilu)

  end subroutine DeallocateDilu  ! DeallocateDilu


  !****************************************************************************
  ! Purpose: to solve the lower triangular system (or it's adjoint);
  ! for the d-ilu pre-condtioner.
  subroutine M1solve(inE, adjt, outE)

    implicit none
    type (cvector), intent(in)	        :: inE
    logical, intent(in)		        :: adjt
    type (cvector), intent(inout) 	:: outE
    integer                             :: ix, iy, iz

    if (.not.inE%allocated) then
      write(0,*) 'inE in M1solve not allocated yet'
      stop
    end if

    if (.not.outE%allocated) then
      write(0,*) 'outE in M1solve not allocated yet'
      stop
    end if

    ! Check whether all the vector nodes are of the same size
    if((inE%nx == outE%nx).and.(inE%ny == outE%ny).and.&
         (inE%nz == outE%nz)) then

       if (inE%gridType == outE%gridType) then

          if (.not.adjt) then
	     ! adjoint = .false.
             Call diagDiv(inE, V_E, outE)

             ! ... note that we only parallelize the outer loops
             !$OMP PARALLEL DEFAULT(SHARED)
             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
             do ix = 1, inE%nx
                !$OMP ORDERED
                do iz = 2, inE%nz
                   do iy = 2, inE%ny

                      outE%x(ix, iy, iz) = (outE%x(ix, iy, iz) - &
                           outE%x(ix, iy-1, iz)*xXYm(ix,iy,iz) - & 
                           outE%x(ix, iy, iz-1)*xXZm(ix,iy,iz))* & 
                           Dilu%x(ix, iy, iz)                      

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO


             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
             do iy = 1, inE%ny
                !$OMP ORDERED
                do iz = 2, inE%nz
                   do ix = 2, inE%nx

                      outE%y(ix, iy, iz) = (outE%y(ix, iy, iz) - & 
                           outE%y(ix, iy, iz-1)*yYZm(ix,iy,iz) - & !
                           outE%y(ix-1, iy, iz)*yYXm(ix,iy,iz))* & !
                           Dilu%y(ix, iy, iz)

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO

             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
             do iz = 1, inE%nz
                !$OMP ORDERED
                do iy = 2, inE%ny
                   do ix = 2, inE%nx

                      outE%z(ix, iy, iz) = (outE%z(ix, iy, iz) - & 
                           outE%z(ix-1, iy, iz)*zZXm(ix,iy,iz) - & 
                           outE%z(ix, iy-1, iz)*zZYm(ix,iy,iz))* & 
                           Dilu%z(ix, iy, iz)

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO
             !$OMP END PARALLEL
             ! adjoint = .true.
          else

             ! ... note that we only parallelize the outer loops
             !$OMP PARALLEL DEFAULT(SHARED)
             ! the coefficients for x are only for the interior nodes
             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
             do ix = 1, inE%nx
                !$OMP ORDERED
                do iy = inE%ny, 2, -1
                   do iz = inE%nz, 2, -1

                      outE%x(ix, iy, iz) = (inE%x(ix, iy, iz) - &
                           outE%x(ix, iy+1, iz)*xXYm(ix,iy+1,iz) - & 
                           outE%x(ix, iy, iz+1)*xXZm(ix,iy,iz+1))* & 
                           conjg(Dilu%x(ix, iy, iz))                 

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO

             ! the coefficients for y are only for the interior nodes
             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
             do iy = 1, inE%ny
                !$OMP ORDERED
                do ix = inE%nx, 2, -1
                   do iz = inE%nz, 2, -1

                      outE%y(ix, iy, iz) = (inE%y(ix, iy, iz) - &      
                           outE%y(ix, iy, iz+1)*yYZm(ix,iy,iz+1) - &   
                           outE%y(ix+1, iy, iz)*yYXm(ix+1,iy,iz))* &   
                           conjg(Dilu%y(ix, iy, iz))                   

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO

             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
             do iz = 1, inE%nz
                !$OMP ORDERED
                do ix = inE%nx, 2, -1
                   do iy = inE%ny, 2, -1

                      outE%z(ix, iy, iz) = (inE%z(ix, iy, iz) - &     
                           outE%z(ix+1, iy, iz)*zZXm(ix+1,iy,iz) - &  
                           outE%z(ix, iy+1, iz)*zZYm(ix,iy+1,iz))* &  
                           conjg(Dilu%z(ix, iy, iz))                  

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO
             !$OMP END PARALLEL

             Call diagDiv(outE, V_E, outE)

          end if

       else
          write (0, *) 'not compatible usage for M1solve'
       end if

    else

       write(0, *) 'Error:lower triangular: vectors not same size'

    end if

  end subroutine M1solve ! M1solve


  !****************************************************************************
  ! Purpose: to solve the upper triangular system (or it's adjoint);
  ! for the d-ilu pre-condtioner
  subroutine M2solve(inE, adjt, outE)

    implicit none
    type (cvector), intent(in)	:: inE
    logical, intent(in)		:: adjt
    type (cvector), intent(inout) 	:: outE
    integer                         :: ix, iy, iz

    if (.not.inE%allocated) then
      write(0,*) 'inE in M2solve not allocated yet'
      stop
    end if

    if (.not.outE%allocated) then
      write(0,*) 'outE in M2solve not allocated yet'
      stop
    end if

    ! Check whether all the vector nodes are of the same size
    if((inE%nx == outE%nx).and.(inE%ny == outE%ny).and.&
         (inE%nz == outE%nz)) then

       if (inE%gridType == outE%gridType) then

          ! adjoint = .false.
          if (.not.adjt) then
             ! for standard upper triangular solution

             ! ... note that we only parallelize the outer loops
             !$OMP PARALLEL DEFAULT(SHARED)
             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
             do ix = 1, inE%nx
                !$OMP ORDERED
                do iz = inE%nz, 2, -1
                   do iy = inE%ny, 2, -1

                      outE%x(ix, iy, iz) = inE%x(ix, iy, iz) - &
                           ( outE%x(ix, iy+1, iz)*xXYp(ix,iy,iz) &   
                           + outE%x(ix, iy, iz+1)*xXZp(ix,iy,iz))* & 
                           Dilu%x(ix, iy, iz)                        

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO

             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
             do iy = 1, inE%ny
                !$OMP ORDERED
                do iz = inE%nz, 2, -1
                   do ix = inE%nx, 2, -1

                      outE%y(ix, iy, iz) = inE%y(ix, iy, iz) - &     
                           ( outE%y(ix, iy, iz+1)*yYZp(ix,iy,iz) &   
                           + outE%y(ix+1, iy, iz)*yYXp(ix,iy,iz))* & 
                           Dilu%y(ix, iy, iz)

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO

             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
             do iz = 1, inE%nz
                !$OMP ORDERED
                do iy = inE%ny, 2, -1
                   do ix = inE%nx, 2, -1

                      outE%z(ix, iy, iz) = inE%z(ix, iy, iz) - &     
                           ( outE%z(ix+1, iy, iz)*zZXp(ix,iy,iz) &   
                           + outE%z(ix, iy+1, iz)*zZYp(ix,iy,iz))* & 
                           Dilu%z(ix, iy, iz)

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO
             !$OMP END PARALLEL
          ! adjoint = .true.
          else

             ! ... note that we only parallelize the outer loops
             !$OMP PARALLEL DEFAULT(SHARED)
             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
             do ix = 1, inE%nx
                !$OMP ORDERED
                do iz = 2, inE%nz
                   do iy = 2, inE%ny

                      outE%x(ix, iy, iz) = inE%x(ix, iy, iz) &
                           - outE%x(ix, iy-1, iz)*xXYp(ix,iy-1,iz) & 
                           * conjg(Dilu%x(ix,iy-1,iz))   &           
                           - outE%x(ix, iy, iz-1)*xXZp(ix,iy,iz-1) & 
                           * conjg(Dilu%x(ix, iy, iz-1))             

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO

             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
             do iy = 1, inE%ny
                !$OMP ORDERED
                do iz = 2, inE%nz
                   do ix = 2, inE%nx

                      outE%y(ix, iy, iz) = inE%y(ix, iy, iz) &       
                           - outE%y(ix, iy, iz-1)*yYZp(ix,iy,iz-1) & 
                           * conjg(Dilu%y(ix,iy,iz-1)) &             
                           - outE%y(ix-1, iy, iz)*yYXp(ix-1,iy,iz) & 
                           * conjg(Dilu%y(ix-1, iy, iz))

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO

             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
             do iz = 1, inE%nz
                !$OMP ORDERED
                do iy = 2, inE%ny
                   do ix = 2, inE%nx

                      outE%z(ix, iy, iz) = inE%z(ix, iy, iz) &       
                           - outE%z(ix-1, iy, iz)*zZXp(ix-1,iy,iz) & 
                           * conjg(Dilu%z(ix-1,iy,iz)) &             
                           - outE%z(ix, iy-1, iz)*zZYp(ix,iy-1,iz) & 
                           * conjg(Dilu%z(ix, iy-1, iz))             

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO
             !$OMP END PARALLEL
          end if

       else
          write (0, *) 'not compatible usage for M2solve'
       end if

    else

       write(0, *) 'Error:lower triangular: vectors not same size'

    end if

  end subroutine M2solve ! M2solve

! *****************************************************************************
! Routines used by divergence correction for 3D finite difference
! EM modeling code. Initialization and application of operator used
! for divergence correction. These routines are used by the divergence
! correction driver routine to (1) compute divergence of currents
! rho =  div sigma E ; (2) set up coefficient matrix for the PDE
! div sigma grad phi = rho ; (3) apply the operator div sigma grad
! The PDE is solved using conjugate gradients with a D-ILU preconditoner.
! The inverse of the pre-conditioner diagonal is set up at the
! same time as the coefficients.  Note that the potential phi that is
! solved for should be constant (phi=0) on the boundary, so that
! tangential components of the correction E-field are zero on the bounary


  !**********************************************************************
  ! to initialize the operator and preconditioner coefficients. Init routines
  ! do memory allocation, reading and setting up arrays
  subroutine DivCorrInit()

    implicit none

    Call create_rvector(mGrid, db1, EDGE)
    Call create_rvector(mGrid, db2, EDGE)
    Call create_rscalar(mGrid, c, CORNER)
    ! d contains the inEerse of diagonal elements of ILU of divCgrad
    Call create_rscalar(mGrid, d, CORNER)
    ! set top nodes to 1.0
    d%v(1,:,:) = 1.0
    d%v(:,1,:) = 1.0
    d%v(:,:,1) = 1.0

    ! initialize volume weights centered at corners
    ! commented out - already initialized in ModelDataInit
    !Call create_rscalar(mGrid, V_N, CORNER)
    !Call NodeVolume(mGrid, V_N)

   end subroutine DivCorrInit  ! DivCorrInit


   !**********************************************************************
   ! SetUp routines do calculations (maybe once; possibly more than once)
   ! DivCorrSetup must be called once for each conductivity distribuition
   !  (i.e., before first forward run; after any change to conductivity)
  subroutine DivCorrSetUp()

    implicit none

    integer                               :: ix, iy, iz
    character*20 ModeName

    !type (cvector):: wE
    !call create_cvector(mGrid,wE,EDGE)

    IF(.not.sigma_E%allocated) THEN
 	WRITE(0,*) 'sigma_E not allocated yet: DivCorrSetUp'
 	STOP
    ENDIF

    IF(.not.db1%allocated) THEN
 	WRITE(0,*) 'db1 not allocated yet: DivCorrSetUp'
 	STOP
    ENDIF

    IF(.not.db2%allocated) THEN
 	WRITE(0,*) 'db2 not allocated yet: DivCorrSetUp'
 	STOP
    ENDIF

        IF(.not.c%allocated) THEN
 	WRITE(0,*) 'c not allocated yet: DivCorrSetUp'
 	STOP
    ENDIF

    IF(.not.d%allocated) THEN
 	WRITE(0,*) 'd not allocated yet: DivCorrSetUp'
 	STOP
    ENDIF

    ! conductivity of air is modified for computing divergence correction
    ! operator coefficients ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do iz = 1, mGrid%nzAir
       do iy = 1, mGrid%ny
          do ix = 1, mGrid%nx
             sigma_E%x(ix, iy, iz) = SIGMA_AIR
             sigma_E%y(ix, iy, iz) = SIGMA_AIR
             sigma_E%z(ix, iy, iz) = SIGMA_AIR
          enddo
       enddo
    enddo    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    ! the coefficients are only for the interior nodes
    ! these coefficients have not been multiplied by volume elements
    ! yet
    do iz = 2, mGrid%nz
       do iy = 2, mGrid%ny
          do ix = 2, mGrid%nx

             db1%x(ix, iy, iz) = sigma_E%x(ix-1,iy,iz)*S_E%x(ix-1,iy,iz)/l_E%x(ix-1,iy,iz)

             db2%x(ix, iy, iz) = sigma_E%x(ix,iy,iz)*S_E%x(ix,iy,iz)/l_E%x(ix,iy,iz)

             db1%y(ix, iy, iz) = sigma_E%y(ix, iy-1,iz)*S_E%y(ix,iy-1,iz)/l_E%y(ix,iy-1,iz)
 
             db2%y(ix, iy, iz) = sigma_E%y(ix, iy, iz)*S_E%y(ix,iy,iz)/l_E%y(ix,iy,iz)
 
             db1%z(ix, iy, iz) = sigma_E%z(ix, iy, iz-1)*S_E%z(ix,iy,iz-1)/l_E%z(ix,iy,iz-1)

             db2%z(ix, iy, iz) = sigma_E%z(ix, iy, iz)*S_E%z(ix,iy,iz)/l_E%z(ix,iy,iz)

             c%v(ix, iy, iz) = - (db1%x(ix, iy, iz) + &
                  db2%x(ix, iy, iz) + &
                  db1%y(ix, iy, iz) + &
                  db2%y(ix, iy, iz) + &
                  db1%z(ix, iy, iz) + &
                  db2%z(ix, iy, iz)   &
                  )

          enddo
       enddo
    enddo
!!!!!!!
    ! change conductivity of air back to zero
    do iz = 1, mGrid%nzAir
       do iy = 1, mGrid%ny
          do ix = 1, mGrid%nx
             sigma_E%x(ix, iy, iz) = R_ZERO
             sigma_E%y(ix, iy, iz) = R_ZERO
             sigma_E%z(ix, iy, iz) = R_ZERO
          enddo
       enddo
    enddo
!!!!!!!!!!!!!!!
    !  To be explicit about forcing coefficients that multiply boundary
    !    nodes to be zero (this gaurantees that the BC on the potential
    !    is phi = 0):
    db1%x(2,:,:) = R_ZERO
    db1%y(:,2,:) = R_ZERO
    db1%z(:,:,2) = R_ZERO
    db2%x(mGrid%nx,:,:) = R_ZERO
    db2%y(:,mGrid%ny,:) = R_ZERO
    db2%z(:,:,mGrid%nz) = R_ZERO

    ! Compute inverse diagonal elements for D-ILU (interior nodes only)
    ! set top nodes to 1.0
    d%v(1,:,:) = 1.0
    d%v(:,1,:) = 1.0
    d%v(:,:,1) = 1.0
    do iz = 2, mGrid%nz
       do iy = 2, mGrid%ny
          do ix = 2, mGrid%nx

             d%v(ix, iy, iz) = c%v(ix, iy, iz) - &
                  db1%x(ix,iy,iz)*db2%x(ix-1,iy,iz)*d%v(ix-1,iy,iz)-&
                  db1%y(ix,iy,iz)*db2%y(ix,iy-1,iz)*d%v(ix,iy-1,iz)-&
                  db1%z(ix,iy,iz)*db2%z(ix,iy,iz-1)*d%v(ix,iy,iz-1)

	     d%v(ix, iy, iz) = 1.0/ d%v(ix, iy, iz)

          enddo
       enddo
    enddo
!
  end subroutine DivCorrSetUp	! DivCorrSetUp


  !**********************************************************************
  ! to deallocate the coefficients used for divergence correction
  subroutine Deallocate_DivCorr()

    implicit none

    Call deall_rvector(db1)
    Call deall_rvector(db2)
    Call deall_rscalar(c)
    Call deall_rscalar(d)
    ! corner volumes
    !Call deall_rscalar(V_N)

  end subroutine Deallocate_DivCorr	! Deallocate_DivCorr


  !**********************************************************************
  ! apply pre-conditioner, solving lower and upper triangular systems using
  ! coefficients in db1, db2, and d.
  subroutine DivCgradILU(inPhi, outPhi)

    implicit none
    type (cscalar), intent(in)                :: inPhi
    type (cscalar), intent(inout)             :: outPhi
    integer                                   :: ix, iy, iz

    IF(.not.inPhi%allocated) THEN
 	WRITE(0,*) 'inPhi not allocated in DivCgradILU'
 	STOP
    ENDIF

    IF(.not.outPhi%allocated) THEN
 	WRITE(0,*) 'outPhi not allocated in DivCgradILU'
 	STOP
    ENDIF

    if (outPhi%allocated) then

       ! Check whether all the inputs/ outputs involved are even of the same
       ! size
       if ((inPhi%nx == outPhi%nx).and.&
            (inPhi%ny == outPhi%ny).and.&
            (inPhi%nz == outPhi%nz)) then

          if (inPhi%gridType == outPhi%gridType) then

             outPhi%v = 0.0
             ! forward substitution (Solve lower triangular system)
             ! the coefficients are only for the interior nodes
             do iz = 2, inPhi%nz
                do iy = 2, inPhi%ny
                   do ix = 2, inPhi%nx

                      outPhi%v(ix, iy, iz) = inPhi%v(ix, iy, iz) &
                           - outPhi%v(ix-1,iy,iz)*db1%x(ix,iy,iz)&
                           *d%v(ix-1,iy,iz) &
                           - outPhi%v(ix,iy-1,iz)*db1%y(ix,iy,iz)&
                           *d%v(ix,iy-1,iz) &
                           - outPhi%v(ix,iy,iz-1)*db1%z(ix,iy,iz)&
                           *d%v(ix,iy,iz-1)

                   enddo
                enddo
             enddo

             ! backward substitution (Solve upper triangular system)
             ! the coefficients are only for the interior nodes
             do iz = inPhi%nz,2,-1
                do iy = inPhi%ny,2,-1
                   do ix = inPhi%nx,2,-1

                      outPhi%v(ix, iy, iz) = (outPhi%v(ix, iy, iz)  &
                           - outPhi%v(ix+1, iy, iz)*db2%x(ix, iy, iz)  &
                           - outPhi%v(ix, iy+1, iz)*db2%y(ix, iy, iz)  &
                           - outPhi%v(ix, iy, iz+1)*db2%z(ix, iy, iz)) &
                           *d%v(ix, iy, iz)

                   enddo
                enddo
             enddo

          else
             write (0, *) 'DivCgradILU: not compatible existing data types'
          end if

       else
          write(0, *) 'Error: DivCgradILU: scalars not same size'
       end if

    else
       write(0, *) 'Error: DivCgradILU: output scalar not even allocated yet'
    end if

  end subroutine DivCgradILU  		! DivCgradILU


  !**********************************************************************
  ! Apply operator div sigma grad to a scalar field (used for corners)
  !  called by PCG for iterative solution of divergence correction equation
  subroutine DivCgrad(inPhi, outPhi)

    implicit none
    type (cscalar), intent(in)                :: inPhi
    type (cscalar), intent(inout)             :: outPhi
    integer                                   :: ix, iy, iz

   IF(.not.inPhi%allocated) THEN
 	WRITE(0,*) 'inPhi not allocated in DivCgrad'
 	STOP
    ENDIF

    IF(.not.outPhi%allocated) THEN
 	WRITE(0,*) 'outPhi not allocated in DivCgrad'
 	STOP
    ENDIF

    if (outPhi%allocated) then

       ! Check whether all the inputs/ outputs involved are even of the same
       ! size
       if ((inPhi%nx == outPhi%nx).and.&
            (inPhi%ny == outPhi%ny).and.&
            (inPhi%nz == outPhi%nz)) then

          if (inPhi%gridType == outPhi%gridType) then

             ! the coefficients are only for interior nodes
             !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz)
             do iz = 2, inPhi%nz
                do iy = 2, inPhi%ny
                   do ix = 2, inPhi%nx

                      outPhi%v(ix,iy,iz) = inPhi%v(ix+1,iy,iz)&
                           *db2%x(ix,iy,iz)+&
                           inPhi%v(ix-1, iy, iz)*db1%x(ix, iy, iz) + &
                           inPhi%v(ix, iy+1, iz)*db2%y(ix, iy, iz) + &
                           inPhi%v(ix, iy-1, iz)*db1%y(ix, iy, iz) + &
                           inPhi%v(ix, iy, iz+1)*db2%z(ix, iy, iz) + &
                           inPhi%v(ix, iy, iz-1)*db1%z(ix, iy, iz) + &
                           inPhi%v(ix, iy, iz)*c%v(ix, iy, iz)

                   enddo
                enddo
             enddo
             !$OMP END PARALLEL DO

          else
             write (0, *) 'DivCgrad: not compatible existing data types'
          end if

       else
          write(0, *) 'Error: DivCgrad: scalars not same size'
       end if

    else
       write(0, *) 'Error: DivCgrad: output scalar not even allocated yet'
    end if

  end subroutine DivCgrad	! DivCgrad


  !**********************************************************************
  ! Purpose is to compute div sigma inE (input electrical field)
  ! where sigma is the edge conductivity. Thus, in practice, this computes
  ! divergence of currents.
  ! NOTE that conductivity in air is modified to SIGMA_AIR for this
  ! subroutine.
  ! This is done as a separate specialized routine to avoid carrying
  ! around multiple edge conductivities
  subroutine DivC(inE, outSc)

    implicit none
    type (cvector), intent(in)		         :: inE
    type (cscalar), intent(inout)		 :: outSc
    integer                                      :: ix, iy, iz

    IF(.not.inE%allocated) THEN
 	WRITE(0,*) 'inE not allocated in DivC'
 	STOP
    ENDIF

    IF(.not.outSc%allocated) THEN
 	WRITE(0,*) 'outSc not allocated in DivC'
 	STOP
    ENDIF

    if (outSc%gridType == CORNER) then

       ! Check whether all the inputs/ outputs involved are even of the same
       ! size
       if ((inE%nx == outSc%nx).and.&
            (inE%ny == outSc%ny).and.&
            (inE%nz == outSc%nz)) then

          ! computation done only for internal nodes
          do ix = 2, outSc%nx
             do iy = 2, outSc%ny

	        ! FOR NODES IN THE AIR ONLY
                do iz = 2,outSc%grid%nzAir
                   outSc%v(ix, iy, iz) = &
                        SIGMA_AIR*(inE%x(ix,iy,iz)*S_E%x(ix,iy,iz)-inE%x(ix-1,iy,iz)*S_E%x(ix-1,iy,iz))   & 
                        + SIGMA_AIR*(inE%y(ix,iy,iz)*S_E%y(ix,iy,iz)-inE%y(ix,iy-1,iz)*S_E%y(ix,iy-1,iz)) & 
                        + SIGMA_AIR*(inE%z(ix,iy,iz)*S_E%z(ix,iy,iz)-inE%z(ix,iy,iz-1)*S_E%z(ix,iy,iz-1))  
                enddo   ! iz

	        ! FOR NODES AT THE AIR-EARTH INTERFACE
                iz = outSc%grid%nzAir+1
                   outSc%v(ix, iy, iz) = &
                        (sigma_E%x(ix,iy,iz)*inE%x(ix, iy, iz)*S_E%x(ix,iy,iz) -         & 
                        sigma_E%x(ix - 1,iy,iz)*inE%x(ix - 1, iy, iz)*S_E%x(ix-1,iy,iz)) & 
                        +  (sigma_E%y(ix,iy,iz)*inE%y(ix, iy, iz)*S_E%y(ix,iy,iz) -      & 
                        sigma_E%y(ix,iy - 1,iz)*inE%y(ix, iy - 1, iz)*S_E%y(ix,iy-1,iz)) & 
                        +  (sigma_E%z(ix,iy,iz)*inE%z(ix, iy, iz)*S_E%z(ix,iy,iz) -      & 
                        SIGMA_AIR*inE%z(ix, iy, iz - 1)*S_E%z(ix,iy,iz-1))                 


                ! FOR NODES INSIDE THE EARTH ONLY
		! THE TOP MOST EARTH NODE HAS AN INTERFACE WITH
		! AIR, THEREFORE THAT ONE IS SKIPPED HERE

                do iz = outSc%grid%nzAir+2, outSc%nz
                   outSc%v(ix, iy, iz) = &
                        (sigma_E%x(ix,iy,iz)*inE%x(ix, iy, iz)*S_E%x(ix,iy,iz) -       & 
                        sigma_E%x(ix-1,iy,iz)*inE%x(ix-1, iy, iz)*S_E%x(ix-1 ,iy,iz))  & 
                        +  (sigma_E%y(ix,iy,iz)*inE%y(ix, iy, iz)*S_E%y(ix,iy,iz) -    & 
                        sigma_E%y(ix,iy-1,iz)*inE%y(ix, iy-1,iz)*S_E%y(ix,iy-1,iz))    & 
                        +  (sigma_E%z(ix,iy,iz)*inE%z(ix, iy, iz)*S_E%z(ix, iy, iz) -  & 
                        sigma_E%z(ix,iy,iz-1)*inE%z(ix,iy,iz-1)*S_E%z(ix,iy,iz-1))       

                enddo   ! iz

             enddo      ! iy
          enddo         ! ix

       else
          write(0, *) 'Error: DivC: scalars not same size'
       end if
       Call diagDiv_crscalar(outSc, V_N, outSc)
    else
       write(0, *) 'Error: DivC: output scalar not compatible use'
    end if

  end subroutine DivC	! DivC

end  module modelOperator3D
