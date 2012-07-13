
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

  use math_constants
  use utilities
  use gridcalc             ! staggered grid definitions
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
  !!!!!!!>>>>>>>>> FROM model_data
  type(grid_t), private, target 	::	mGrid   ! THE model grid
  type(rvector), public			::	volE    ! THE volume elements
  type(rvector), private		::	condE   ! THE edge conductivities
  real(kind=prec),private	::      omega   ! THE (active) frequency

  ! NOTE: THIS VARIABLE IS TEMPORARILY REQUIRED TO SET THE BOUNDARY CONDITIONS
  type(rscalar), private        :: Cond3D

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
  real (kind=prec), pointer, dimension(:,:), private    :: xXY, xXZ
  real (kind=prec), pointer, dimension(:,:), private    :: xY, xZ
  real (kind=prec), pointer, dimension(:,:), private    :: xXO
  real (kind=prec), pointer, dimension(:,:), private    :: yYX, yYZ
  real (kind=prec), pointer, dimension(:,:), private    :: yX, yZ
  real (kind=prec), pointer, dimension(:,:), private    :: yYO
  real (kind=prec), pointer, dimension(:,:), private    :: zZX, zZY
  real (kind=prec), pointer, dimension(:,:), private    :: zX, zY
  real (kind=prec), pointer, dimension(:,:), private    :: zZO

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

  type (rscalar) , public	:: volC
  ! volume contains weights for corners

  ! *****************************************************************************
  !  routines from model_data_update:
  public                             	:: UpdateFreq, UpdateCond
  public                                :: UpdateFreqCond
  public                                :: ModelDataInit
  !   These are used to initialize the modules grid, and to set/modify
  !     conductivity and/or frequency

  !  routine to set the boundary conditions (a wrapper for BC_x0_WS for now)
  public                                :: SetBound

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
    Call create(mGrid, volE, EDGE)
    Call EdgeVolume(mGrid, volE)

    !  Allocate condE, conductivity defined on computational cell edges
    !   condE is also kept in module model_data
    Call create(mGrid, condE, EDGE)

    ! set a default omega
    omega = 0.0

  end subroutine ModelDataInit


  subroutine ModelDataCleanUp

    call deall_grid(mGrid)
    call deall_rvector(volE)
    call deall_rvector(condE)

    ! Cond3D is temporary
    call deall_rscalar(Cond3D)

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
    Call ModelParamToEdge(CondParam, condE)

    Call DivCorrSetUp()

    ! TEMPORARY; REQUIRED FOR BOUNDARY CONDITIONS
    !  set static array for cell conductivities
    call ModelParamToCell(CondParam,Cond3D)

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
    Call ModelParamToEdge(CondParam, condE)

    Call AdiagSetUp()
    Call DiluSetUp()
    Call DivCorrSetUp()

    ! TEMPORARY; REQUIRED FOR BOUNDARY CONDITIONS
    !  set static array for cell conductivities
    call ModelParamToCell(CondParam,Cond3D)

  end subroutine UpdateFreqCond  ! UpdateFreqCond

!**********************************************************************
! Sets boundary conditions. Currently a wrapper for BC_x0_WS.
! Uses input 3D conductivity in cells Cond3D, that has to be initialized
! by updateCond before calling this routine. Also uses mGrid set by
! ModelDataInit. Could use omega, which is set by updateFreq.
  Subroutine SetBound(imode,period,E0,BC,iTx)

    !  Input mode, period
    integer, intent(in)		:: imode
    integer, intent(in)		:: iTx
    real(kind=prec)	:: period

    ! Output electric field first guess (for iterative solver)
    type(cvector), intent(inout)	:: E0
    ! Output boundary conditions
    type(cboundary), intent(inout)	:: BC

    if (BC%read_E_from_file) then
          call BC_x0_WS(imode,period,mGrid,Cond3D,E0,BC)  
          ! The BC are already computed from a larger grid for all transmitters and modes and stored in BC_from_file.
          ! Overwrite BC with BC_from_file.
          ! Note: Right now we are using the same period layout for both grid. 
          ! This why, it is enough to know the period and mode index to pick up the BC from BC_from_file vector.
          BC = BC_from_file((iTx*2)-(2-imode))  
    else
         ! Compute the BC using Weerachai 2D approach 
          call BC_x0_WS(imode,period,mGrid,Cond3D,E0,BC)                          
    end if
    
    
    ! Cell conductivity array is no longer needed
    ! NOT TRUE: needed for imode=2
    ! call deall_rscalar(Cond3D)

  end subroutine SetBound

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


    ! Allocate memory for del x del operator coefficient arrays
    ! Coefficients for difference equation only uses interior
    ! nodes. however, we need boundary nodes for the adjoint problem
    allocate(xXY(mGrid%ny+1, 2), STAT=status)   ! Allocate memory
    allocate(xXZ(mGrid%nz+1, 2), STAT=status)   ! Allocate memory
    allocate(xY(mGrid%nx, mGrid%ny+1), STAT=status)
    allocate(xZ(mGrid%nx, mGrid%nz+1), STAT=status)
    allocate(xXO(mGrid%ny, mGrid%nz), STAT=status)

    allocate(yYZ(mGrid%nz+1, 2), STAT=status)   ! Allocate memory
    allocate(yYX(mGrid%nx+1, 2), STAT=status)   ! Allocate memory
    allocate(yZ(mGrid%ny, mGrid%nz+1), STAT=status)
    allocate(yX(mGrid%nx+1, mGrid%ny), STAT=status)
    allocate(yYO(mGrid%nx, mGrid%nz), STAT=status)

    allocate(zZX(mGrid%nx+1, 2), STAT=status)   ! Allocate memory
    allocate(zZY(mGrid%ny+1, 2), STAT=status)   ! Allocate memory
    allocate(zX(mGrid%nx+1, mGrid%nz), STAT=status)
    allocate(zY(mGrid%ny+1, mGrid%nz), STAT=status)
    allocate(zZO(mGrid%nx, mGrid%ny), STAT=status)


    ! initalize all coefficients to zero (some remain zero)
    xXY = 0.0
    xXZ = 0.0
    xY = 0.0
    xZ = 0.0
    xXO = 0.0
    yYX = 0.0
    yYZ = 0.0
    yX = 0.0
    yZ = 0.0
    zZX = 0.0
    zZY = 0.0
    zX = 0.0
    zY = 0.0
    zZO = 0.0

    ! coefficents for calculating Ex ; only loop over internal edges
    do iy = 2, mGrid%ny
       xXY(iy, 2) = -1.0/ (mGrid%delY(iy) * mGrid%dy(iy))
       xXY(iy, 1) = -1.0/ (mGrid%delY(iy) * mGrid%dy(iy-1))
    enddo

    do iz = 2, mGrid%nz
       xXZ(iz, 2) = -1.0/ (mGrid%delZ(iz) * mGrid%dz(iz))
       xXZ(iz, 1) = -1.0/ (mGrid%delZ(iz) * mGrid%dz(iz-1))
    enddo

    do iy = 2, mGrid%ny
       do iz = 2, mGrid%nz
          xXO(iy, iz) = -(xXY(iy,1) + xXY(iy,2) + &
               xXZ(iz,1) + xXZ(iz,2))

       enddo
    enddo

    do ix = 1, mGrid%nx
       do iy = 2, mGrid%ny
          xY(ix, iy) = 1.0/ (mGrid%delY(iy)*mGrid%dx(ix))
       enddo
    enddo

    do ix = 1, mGrid%nx
       do iz = 2, mGrid%nz
          xZ(ix, iz) = 1.0/ (mGrid%delZ(iz)*mGrid%dx(ix))
       enddo
    enddo
    ! End of Ex coefficients

    ! coefficents for calculating Ey; only loop over internal edges
    do iz = 2, mGrid%nz
       yYZ(iz, 2) = -1.0/ (mGrid%delZ(iz)*mGrid%dz(iz))
       yYZ(iz, 1) = -1.0/ (mGrid%delZ(iz)*mGrid%dz(iz-1))
    enddo

    do ix = 2, mGrid%nx
       yYX(ix, 2) = -1.0/ (mGrid%delX(ix)*mGrid%dx(ix))
       yYX(ix, 1) = -1.0/ (mGrid%delX(ix)*mGrid%dx(ix-1))
    enddo

    do ix = 2, mGrid%nx
       do iz = 2, mGrid%nz
          yYO(ix, iz) = -(yYX(ix,1) + yYX(ix,2) + &
               yYZ(iz,1) + yYZ(iz,2))
       enddo
    enddo

    do iy = 1, mGrid%ny
       do iz = 2, mGrid%nz
          yZ(iy, iz) = 1.0/ (mGrid%delZ(iz)*mGrid%dy(iy))
       enddo
    enddo

    do ix = 2, mGrid%nx
       do iy = 1, mGrid%ny
          yX(ix, iy) = 1.0/ (mGrid%delX(ix)*mGrid%dy(iy))
       enddo
    enddo
    ! End of Ey coefficients

    ! coefficents for calculating Ez; only loop over internal edges
    do ix = 2, mGrid%nx
       zZX(ix, 2) = -1.0/ (mGrid%delX(ix)*mGrid%dx(ix))
       zZX(ix, 1) = -1.0/ (mGrid%delX(ix)*mGrid%dx(ix-1))
    enddo

    do iy = 2, mGrid%ny
       zZY(iy, 2) = -1.0/ (mGrid%delY(iy)*mGrid%dy(iy))
       zZY(iy, 1) = -1.0/ (mGrid%delY(iy)*mGrid%dy(iy-1))
    enddo

    do ix = 2, mGrid%nx
       do iy = 2, mGrid%ny
          zZO(ix, iy) = -(zZX(ix,1) + zZX(ix,2) + &
               zZY(iy,1) + zZY(iy,2))
       enddo
    enddo

    do ix = 2, mGrid%nx
       do iz = 1, mGrid%nz
          zX(ix, iz) = 1.0/ (mGrid%delX(ix)*mGrid%dz(iz))
       enddo
    enddo

    do iy = 2, mGrid%ny
       do iz = 1, mGrid%nz
          zY(iy, iz) = 1.0/ (mGrid%delY(iy)*mGrid%dz(iz))
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

          ! Apply difference equation to compute Ex (only on interior nodes)
          do iz = 2, inE%nz
             do iy = 2, inE%ny
                do ix = 1, inE%nx

                   outE%x(ix,iy,iz) = xY(ix,iy)*(inE%y(ix+1,iy,iz)-&
                  	inE%y(ix,iy,iz)-inE%y(ix+1,iy-1,iz)&
                        +inE%y(ix,iy-1,iz))+&
                  	xZ(ix,iz)*(inE%z(ix+1,iy,iz)-inE%z(ix,iy,iz)&
                  	-inE%z(ix+1,iy,iz-1)+inE%z(ix,iy,iz-1))+&
                  	xXY(iy,2)*inE%x(ix,iy+1,iz)+&
                  	xXY(iy,1)*inE%x(ix,iy-1,iz)+&
                  	xXZ(iz,2)*inE%x(ix,iy,iz+1)+&
                  	xXZ(iz,1)*inE%x(ix,iy,iz-1)+&
                  	xXO(iy,iz)*inE%x(ix,iy,iz)

                enddo
             enddo
          enddo

          ! Apply difference equation to compute Ey (only on interior nodes)
          do iz = 2, inE%nz
             do iy = 1, inE%ny
                do ix = 2, inE%nx

                   outE%y(ix,iy,iz) = yZ(iy,iz)*(inE%z(ix,iy+1,iz)-&
                  	inE%z(ix,iy,iz)-inE%z(ix,iy+1,iz-1)+inE%z(ix,iy,iz-1))&
                  	+yX(ix,iy)*(inE%x(ix,iy+1,iz)-inE%x(ix,iy,iz)&
                  	-inE%x(ix-1,iy+1,iz)+inE%x(ix-1,iy,iz))+&
                  	yYZ(iz,2)*inE%y(ix,iy,iz+1)+&
                  	yYZ(iz,1)*inE%y(ix,iy,iz-1)+&
                  	yYX(ix,2)*inE%y(ix+1,iy,iz)+&
                  	yYX(ix,1)*inE%y(ix-1,iy,iz)+&
                  	yYO(ix,iz)*inE%y(ix,iy,iz)

                enddo
             enddo
          enddo

          ! Apply difference equation to compute Ez (only on interior nodes)
          do iz = 1, inE%nz
             do iy = 2, inE%ny
                do ix = 2, inE%nx

                   outE%z(ix,iy,iz) = zX(ix,iz)*(inE%x(ix,iy,iz+1)-&
                  	inE%x(ix,iy,iz)-inE%x(ix-1,iy,iz+1)+inE%x(ix-1,iy,iz))&
                  	+zY(iy,iz)*(inE%y(ix,iy,iz+1)-inE%y(ix,iy,iz)&
                  	-inE%y(ix,iy-1,iz+1)+inE%y(ix,iy-1,iz))+&
                  	zZX(ix,2)*inE%z(ix+1,iy,iz)+&
                  	zZX(ix,1)*inE%z(ix-1,iy,iz)+&
                  	zZY(iy,2)*inE%z(ix,iy+1,iz)+&
                  	zZY(iy,1)*inE%z(ix,iy-1,iz)+&
                  	zZO(ix,iy)*inE%z(ix,iy,iz)

                enddo
             enddo
          enddo

       else
          write (0, *) 'CurlcurlE: not compatible usage for existing data types'
       end if

    else
       write(0, *) 'Error-complex vectors for CurlcurlE are not of same size'
    end if

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
    deallocate(xXY, STAT=status)
    deallocate(xXZ, STAT=status)
    deallocate(xY, STAT=status)
    deallocate(xZ, STAT=status)
    deallocate(xXO, STAT=status)

    deallocate(yYZ, STAT=status)
    deallocate(yYX, STAT=status)
    deallocate(yZ, STAT=status)
    deallocate(yX, STAT=status)
    deallocate(yYO, STAT=status)

    deallocate(zZX, STAT=status)
    deallocate(zZY, STAT=status)
    deallocate(zX, STAT=status)
    deallocate(zY, STAT=status)
    deallocate(zZO, STAT=status)

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

    do ix = 1, mGrid%nx
       Adiag%x(ix,:,:) = CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%x(ix,:,:)
    enddo

    do iy = 1, mGrid%ny
       Adiag%y(:,iy,:) = CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%y(:,iy,:)
    enddo

    do iz = 1, mGrid%nz
       Adiag%z(:,:,iz) = CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%z(:,:,iz)
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
    integer                                  :: diag_sign
    integer                                  :: ix, iy, iz
    ! dummy variables

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

          !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ix,iy,iz)

          ! Apply difference equation to compute Ex (only on interior nodes)
          ! the diagonal nodes have the imaginary component added
          !$OMP DO SCHEDULE(STATIC)
	      do iz = 2, inE%nz
             do iy = 2, inE%ny
                do ix = 1, inE%nx
                   outE%x(ix,iy,iz) = xY(ix,iy)*(inE%y(ix+1,iy,iz)-&
                  	inE%y(ix,iy,iz)-inE%y(ix+1,iy-1,iz)&
                        +inE%y(ix,iy-1,iz))+&
                  	xZ(ix,iz)*(inE%z(ix+1,iy,iz)-inE%z(ix,iy,iz)&
                  	-inE%z(ix+1,iy,iz-1)+inE%z(ix,iy,iz-1))+&
                  	xXY(iy,2)*inE%x(ix,iy+1,iz)+&
                  	xXY(iy,1)*inE%x(ix,iy-1,iz)+&
                  	xXZ(iz,2)*inE%x(ix,iy,iz+1)+&
                  	xXZ(iz,1)*inE%x(ix,iy,iz-1)+&
                  	(xXO(iy,iz)+diag_sign*Adiag%x(ix,iy,iz))*inE%x(ix,iy,iz)
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
                   outE%y(ix,iy,iz) = yZ(iy,iz)*(inE%z(ix,iy+1,iz)-&
                  	inE%z(ix,iy,iz)-inE%z(ix,iy+1,iz-1)+inE%z(ix,iy,iz-1))&
                  	+yX(ix,iy)*(inE%x(ix,iy+1,iz)-inE%x(ix,iy,iz)&
                  	-inE%x(ix-1,iy+1,iz)+inE%x(ix-1,iy,iz))+&
                  	yYZ(iz,2)*inE%y(ix,iy,iz+1)+&
                  	yYZ(iz,1)*inE%y(ix,iy,iz-1)+&
                  	yYX(ix,2)*inE%y(ix+1,iy,iz)+&
                  	yYX(ix,1)*inE%y(ix-1,iy,iz)+&
                  	(yYO(ix,iz)+diag_sign*Adiag%y(ix,iy,iz))*inE%y(ix,iy,iz)
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
                   outE%z(ix,iy,iz) = zX(ix,iz)*(inE%x(ix,iy,iz+1)-&
                  	inE%x(ix,iy,iz)-inE%x(ix-1,iy,iz+1)+inE%x(ix-1,iy,iz))&
                  	+zY(iy,iz)*(inE%y(ix,iy,iz+1)-inE%y(ix,iy,iz)&
                  	-inE%y(ix,iy-1,iz+1)+inE%y(ix,iy-1,iz))+&
                  	zZX(ix,2)*inE%z(ix+1,iy,iz)+&
                  	zZX(ix,1)*inE%z(ix-1,iy,iz)+&
                  	zZY(iy,2)*inE%z(ix,iy+1,iz)+&
                  	zZY(iy,1)*inE%z(ix,iy-1,iz)+&
                  	(zZO(ix,iy)+diag_sign*Adiag%z(ix,iy,iz))*inE%z(ix,iy,iz)
                enddo
             enddo
          enddo
          !$OMP END DO NOWAIT

          !$OMP END PARALLEL

       else
          write (0, *) ' Maxwell: not compatible usage for existing data types'
       end if

    else
       write(0, *) 'Error-complex vectors for Maxwell are not of same size'
    end if

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
          Call diagMult_crvector(inE, condE, outE)
          c2 = diag_sign*C_ONE*omega*MU_0
          Call linComb_cvector(C_ONE, workE, c2, outE, outE)

          ! diagonally multiply the final results with weights (edge volume)
          Call diagMult_crvector(outE, volE, outE)

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

write(0,*) 'DEBUG 2-1'
          Call Maxwell(inE, adjt, outE)
          ! done with preparing del X del X E +/- i*omega*mu*conductivity*E
write(0,*) 'DEBUG 2-2'

          ! diagonally multiply the final results with weights (edge volume)
          Call diagMult(outE, volE, outE)
write(0,*) 'DEBUG 2-3'

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

!  Ex components in x/z plane (iy=1; iy = ny+1)
!  NOTE: C_ZERO = (0,0) (double complex) is defined in SG_Basics/math_constants.f90
    BC%xYMax(:,1) = C_ZERO
    BC%xYmin(:,1) = C_ZERO
    BC%xYMax(:,nz+1) = C_ZERO
    BC%xyMin(:,nz+1) = C_ZERO
    do ix = 1, nx
       do iz = 2, nz
          BC%xYmin(ix,iz) = - yX(ix,1)*Ein%y(ix,1,iz)       &
                            + yX(ix+1,1)*Ein%y(ix+1,1,iz)   &
                            + xXY(2,1)*Ein%x(ix,2,iz)
          BC%xYmax(ix,iz) = + yX(ix,ny)*Ein%y(ix,ny,iz)     &
                            - yX(ix+1,ny)*Ein%y(ix+1,ny,iz) &
                            + xXY(ny,2)*Ein%x(ix,ny,iz)
        enddo
     enddo

!  Ez components in x/z plane (iy=1; iy = ny+1)
    BC%zYMin(1,:) = C_ZERO
    BC%zYmax(1,:) = C_ZERO
    BC%zYmin(nx+1,:) = C_ZERO
    BC%zYmax(nx+1,:) = C_ZERO
    do iz = 1, nz
       do ix = 2, nx
          BC%zYmin(ix,iz) = - yZ(1,iz)*Ein%y(ix,1,iz)        &
                            + yZ(1,iz+1)*Ein%y(ix,1,iz+1)    &
                            + zZY(2,1)*Ein%z(ix,2,iz)
          BC%zYmax(ix,iz) = + yZ(ny,iz)*Ein%y(ix,ny,iz)      &
                            - yZ(ny,iz+1)*Ein%y(ix,ny,iz+1)  &
                            + zZY(ny,2)*Ein%z(ix,ny,iz)
        enddo
     enddo

!  Ey components in y/z plane (ix=1; ix = nx+1)
    BC%yXmin(:,1) = C_ZERO
    BC%yXmax(:,1) = C_ZERO
    BC%yXmin(:,nz+1) = C_ZERO
    BC%yXmax(:,nz+1) = C_ZERO
    do iy = 1, ny
       do iz = 2, nz
          BC%yXmin(iy,iz) = - xY(1,iy)*Ein%x(1,iy,iz)        &
                            + xY(1,iy+1)*Ein%x(1,iy+1,iz)    &
                            + yYX(2,1)*Ein%y(2,iy,iz)
          BC%yXmax(iy,iz) = + xY(nx,iy)*Ein%x(nx,iy,iz)      &
                            - xY(nx,iy+1)*Ein%x(nx,iy+1,iz)  &
                            + yYX(nx,2)*Ein%y(nx,iy,iz)
        enddo
     enddo

!  Ez components in y/z plane (ix=1; ix = nx+1)
    BC%zXmin(1,:) = C_ZERO
    BC%zXmax(1,:) = C_ZERO
    BC%zXmin(ny+1,:) = C_ZERO
    BC%zXmax(ny+1,:) = C_ZERO
    do iz = 1, nz
       do iy = 2, ny
          BC%zXmin(iy,iz) = - xZ(1,iz)*Ein%x(1,iy,iz)       &
                            + xZ(1,iz+1)*Ein%x(1,iy,iz+1)   &
                            + zZX(2,1)*Ein%z(2,iy,iz)
          BC%zXmax(iy,iz) = + xZ(nx,iz)*Ein%x(nx,iy,iz)     &
                            - xZ(nx,iz+1)*Ein%x(nx,iy,iz+1) &
                            + zZX(nx,2)*Ein%z(nx,iy,iz)
        enddo
     enddo

!  Ex components in x/y plane (iz=1; iz = nz+1)
    BC%xZmin(:,1) = C_ZERO
    BC%xZmax(:,1) = C_ZERO
    BC%xZmin(:,ny+1) = C_ZERO
    BC%xZmax(:,ny+1) = C_ZERO
    do ix = 1, nx
       do iy = 2, ny
          BC%xZmin(ix,iy) = - zX(ix,1)*Ein%z(ix,iy,1)       &
                            + zX(ix+1,1)*Ein%z(ix+1,iy,1)   &
                            + xXZ(2,1)*Ein%x(ix,iy,2)
          BC%xZmax(ix,iy) = + zX(ix,nz)*Ein%z(ix,iy,nz)     &
                            - zX(ix+1,nz)*Ein%z(ix+1,iy,nz) &
                            + xXZ(nz,2)*Ein%x(ix,iy,nz)
        enddo
     enddo

!  Ey components in x/y plane (iz=1; iz = nz+1)
    BC%yZmin(:,1) = C_ZERO
    BC%yZmax(:,1) = C_ZERO
    BC%yZmin(nx+1,:) = C_ZERO
    BC%yZmin(nx+1,:) = C_ZERO
    do iy = 1, ny
       do ix = 2, nx
          BC%yZmin(ix,iy) = - zY(iy,1)*Ein%z(ix,iy,1)        &
                            + zY(iy+1,1)*Ein%z(ix,iy+1,1)    &
                            + yYZ(2,1)*Ein%y(ix,iy,2)
          BC%yZmax(ix,iy) = + zY(iy,nz)*Ein%z(ix,iy,nz)      &
                            - zY(iy+1,nz)*Ein%z(ix,iy+1,nz)  &
                            + yYZ(nz,2)*Ein%y(ix,iy,nz)
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

             Dilu%x(ix, iy, iz) = xXO(iy,iz) - &
                  CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%x(ix, iy, iz)  &
                  - xXY(iy, 1)*xXY(iy-1, 2)*Dilu%x(ix,iy-1,iz) &
                  - xXZ(iz, 1)*xXZ(iz-1, 2)*Dilu%x(ix,iy,iz-1)
             Dilu%x(ix, iy, iz) = 1.0/ Dilu%x(ix, iy, iz)

          enddo
       enddo
    enddo

    ! the coefficients for y are only for the interior nodes
    !  but need to initialize edges for recursive algorithm
    do iy = 1, mGrid%ny
       do iz = 2, mGrid%nz
          do ix = 2, mGrid%nx

             Dilu%y(ix, iy, iz) = yYO(ix,iz) - &
                  CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%y(ix, iy, iz) &
                  - yYZ(iz, 1)*yYZ(iz-1, 2)*Dilu%y(ix, iy, iz-1) &
                  - yYX(ix, 1)*yYX(ix-1, 2)*Dilu%y(ix-1, iy, iz)
             Dilu%y(ix, iy, iz) = 1.0/ Dilu%y(ix, iy, iz)

          enddo
       enddo
    enddo

    ! the coefficients for z are only for the interior nodes
    !  but need to initialize edges for recursive algorithm
    do iz = 1, mGrid%nz
       do ix = 2, mGrid%nx
          do iy = 2, mGrid%ny

             Dilu%z(ix, iy, iz) = zZO(ix,iy) - &
                  CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%z(ix, iy, iz) &
                  - zZX(ix, 1)*zZX(ix-1, 2)*Dilu%z(ix-1, iy, iz) &
                  - zZY(iy, 1)*zZY(iy-1, 2)*Dilu%z(ix, iy-1, iz)
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
             Call diagDiv(inE, volE, outE)

             ! ... note that we only parallelize the outer loops
             !$OMP PARALLEL DEFAULT(SHARED)
             !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
             do ix = 1, inE%nx
                !$OMP ORDERED
                do iz = 2, inE%nz
                   do iy = 2, inE%ny

                      outE%x(ix, iy, iz) = (outE%x(ix, iy, iz) - &
                           outE%x(ix, iy-1, iz)*xXY(iy, 1) - &
                           outE%x(ix, iy, iz-1)*xXZ(iz, 1))* &
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
                           outE%y(ix, iy, iz-1)*yYZ(iz, 1) - &
                           outE%y(ix-1, iy, iz)*yYX(ix, 1))* &
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
                           outE%z(ix-1, iy, iz)*zZX(ix, 1) - &
                           outE%z(ix, iy-1, iz)*zZY(iy, 1))* &
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
                           outE%x(ix, iy+1, iz)*xXY(iy+1, 1) - &
                           outE%x(ix, iy, iz+1)*xXZ(iz+1, 1))* &
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
                           outE%y(ix, iy, iz+1)*yYZ(iz+1, 1) - &
                           outE%y(ix+1, iy, iz)*yYX(ix+1, 1))* &
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
                           outE%z(ix+1, iy, iz)*zZX(ix+1, 1) - &
                           outE%z(ix, iy+1, iz)*zZY(iy+1, 1))* &
                           conjg(Dilu%z(ix, iy, iz))

                   enddo
                enddo
                !$OMP END ORDERED
             enddo
             !$OMP END DO
             !$OMP END PARALLEL

             Call diagDiv(outE, volE, outE)

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
                           ( outE%x(ix, iy+1, iz)*xXY(iy, 2) &
                           + outE%x(ix, iy, iz+1)*xXZ(iz, 2))* &
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
                           ( outE%y(ix, iy, iz+1)*yYZ(iz, 2) &
                           + outE%y(ix+1, iy, iz)*yYX(ix, 2))* &
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
                           ( outE%z(ix+1, iy, iz)*zZX(ix, 2) &
                           + outE%z(ix, iy+1, iz)*zZY(iy, 2))* &
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
                           - outE%x(ix, iy-1, iz)*xXY(iy-1, 2) &
                           * conjg(Dilu%x(ix,iy-1,iz))   &
                           - outE%x(ix, iy, iz-1)*xXZ(iz-1, 2) &
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
                           - outE%y(ix, iy, iz-1)*yYZ(iz-1, 2) &
                           * conjg(Dilu%y(ix,iy,iz-1)) &
                           - outE%y(ix-1, iy, iz)*yYX(ix-1, 2) &
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
                           - outE%z(ix-1, iy, iz)*zZX(ix-1, 2) &
                           * conjg(Dilu%z(ix-1,iy,iz)) &
                           - outE%z(ix, iy-1, iz)*zZY(iy-1, 2) &
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
    Call create_rscalar(mGrid, volC, CORNER)
    Call CornerVolume(mGrid, volC)

   end subroutine DivCorrInit  ! DivCorrInit


   !**********************************************************************
   ! SetUp routines do calculations (maybe once; possibly more than once)
   ! DivCorrSetup must be called once for each conductivity distribuition
   !  (i.e., before first forward run; after any change to conductivity)
  subroutine DivCorrSetUp()

    implicit none
    integer                               :: ix, iy, iz

    IF(.not.condE%allocated) THEN
 	WRITE(0,*) 'condE not allocated yet: DivCorrSetUp'
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
    do iz = 1, mGrid%nzAir
       do iy = 1, mGrid%ny
          do ix = 1, mGrid%nx
             condE%x(ix, iy, iz) = SIGMA_AIR
             condE%y(ix, iy, iz) = SIGMA_AIR
             condE%z(ix, iy, iz) = SIGMA_AIR
          enddo
       enddo
    enddo

    ! the coefficients are only for the interior nodes
    ! these coefficients have not been multiplied by volume elements
    ! yet
    do iz = 2, mGrid%nz
       do iy = 2, mGrid%ny
          do ix = 2, mGrid%nx

             db1%x(ix, iy, iz) = condE%x(ix-1, iy, iz)/ &
                  (mGrid%dx(ix-1)*mGrid%delX(ix))
             db2%x(ix, iy, iz) = condE%x(ix, iy, iz)/ &
                  (mGrid%dx(ix)*mGrid%delX(ix))
             db1%y(ix, iy, iz) = condE%y(ix, iy-1, iz)/ &
                  (mGrid%dy(iy-1)*mGrid%delY(iy))
             db2%y(ix, iy, iz) = condE%y(ix, iy, iz)/ &
                  (mGrid%dy(iy)*mGrid%delY(iy))
             db1%z(ix, iy, iz) = condE%z(ix, iy, iz-1)/ &
                  (mGrid%dz(iz-1)*mGrid%delZ(iz))
             db2%z(ix, iy, iz) = condE%z(ix, iy, iz)/ &
                  (mGrid%dz(iz)*mGrid%delZ(iz))
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

    ! change conductivity of air back to zero
    do iz = 1, mGrid%nzAir
       do iy = 1, mGrid%ny
          do ix = 1, mGrid%nx
             condE%x(ix, iy, iz) = R_ZERO
             condE%y(ix, iy, iz) = R_ZERO
             condE%z(ix, iy, iz) = R_ZERO
          enddo
       enddo
    enddo

    ! Multiply by corner volume elements to make operator symmetric
    do iz = 2, mGrid%nz
       do iy = 2, mGrid%ny
          do ix = 2, mGrid%nx

             db1%x(ix, iy, iz) = db1%x(ix, iy, iz)*volC%v(ix, iy, iz)
	     db1%y(ix, iy, iz) = db1%y(ix, iy, iz)*volC%v(ix, iy, iz)
	     db1%z(ix, iy, iz) = db1%z(ix, iy, iz)*volC%v(ix, iy, iz)
	     db2%x(ix, iy, iz) = db2%x(ix, iy, iz)*volC%v(ix, iy, iz)
	     db2%y(ix, iy, iz) = db2%y(ix, iy, iz)*volC%v(ix, iy, iz)
	     db2%z(ix, iy, iz) = db2%z(ix, iy, iz)*volC%v(ix, iy, iz)

          enddo
       enddo
    enddo
    Call diagMult_rscalar(c, volC, c)

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
    Call deall_rscalar(volC)

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
                        SIGMA_AIR*(inE%x(ix,iy,iz)-inE%x(ix - 1,iy,iz)) * &
                        inE%grid%delXinv(ix)    &
                        + SIGMA_AIR*(inE%y(ix,iy,iz)-inE%y(ix,iy - 1,iz)) * &
                        inE%grid%delYinv(iy)    &
                        + SIGMA_AIR*(inE%z(ix,iy,iz)-inE%z(ix,iy,iz - 1)) * &
                        inE%grid%delZinv(iz)
                enddo   ! iz

	        ! FOR NODES AT THE AIR-EARTH INTERFACE
                iz = outSc%grid%nzAir+1
                   outSc%v(ix, iy, iz) = &
                        (condE%x(ix,iy,iz)*inE%x(ix, iy, iz) -         &
                        condE%x(ix - 1,iy,iz)*inE%x(ix - 1, iy, iz)) * &
                        inE%grid%delXinv(ix)      &
                        +  (condE%y(ix,iy,iz)*inE%y(ix, iy, iz) -      &
                        condE%y(ix,iy - 1,iz)*inE%y(ix, iy - 1, iz)) * &
                        inE%grid%delYinv(iy)      &
                        +  (condE%z(ix,iy,iz)*inE%z(ix, iy, iz) -      &
                        SIGMA_AIR*inE%z(ix, iy, iz - 1)) * &
                        inE%grid%delZinv(iz)

                ! FOR NODES INSIDE THE EARTH ONLY
		! THE TOP MOST EARTH NODE HAS AN INTERFACE WITH
		! AIR, THEREFORE THAT ONE IS SKIPPED HERE
                do iz = outSc%grid%nzAir+2, outSc%nz
                   outSc%v(ix, iy, iz) = &
                        (condE%x(ix,iy,iz)*inE%x(ix, iy, iz) -         &
                        condE%x(ix - 1,iy,iz)*inE%x(ix - 1, iy, iz)) * &
                        inE%grid%delXinv(ix)      &
                        +  (condE%y(ix,iy,iz)*inE%y(ix, iy, iz) -      &
                        condE%y(ix,iy - 1,iz)*inE%y(ix, iy - 1, iz)) * &
                        inE%grid%delYinv(iy)      &
                        +  (condE%z(ix,iy,iz)*inE%z(ix, iy, iz) -      &
                        condE%z(ix,iy,iz - 1)*inE%z(ix, iy, iz - 1)) * &
                        inE%grid%delZinv(iz)
                enddo   ! iz

             enddo      ! iy
          enddo         ! ix

       else
          write(0, *) 'Error: DivC: scalars not same size'
       end if

    else
       write(0, *) 'Error: DivC: output scalar not compatible use'
    end if

  end subroutine DivC	! DivC

end  module modelOperator3D
