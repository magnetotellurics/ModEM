
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
      ! modified by Cherevatova (2012) for the multi-grid

      use math_constants
      use utilities
      use gridcalc             ! staggered grid definitions
      use sg_vector            ! generic routines for vector operations on the
      use sg_vector_mg
      use sg_boundary
      use ModelSpace
      use boundary_ws          ! sets the boundary conditions
      implicit none

      ! * These variables are used by model equation
      ! * and preconditioner modules;
      ! * All variables are saved until changed/deallocated

      save
      !!!!!!!>>>>>>>>> FROM model_data
      type(grid_t), private, target 	::	mGrid   ! THE model grid
      type(rvector_mg), public  ::  volE    ! THE volume elements !!! MULTIGRID!!!!
      type(rvector_mg), private      ::  condE   ! THE edge conductivities on MULTIGRID
      real(kind=prec),private	::      omega   ! THE (active) frequency
      !Vector to hold the BC interpolated from a file E solution file
      type(cboundary), pointer, dimension(:)	:: BC_from_file
      type(cvector), pointer, dimension(:)	    :: E0_from_file

      ! NOTE: THIS VARIABLE IS TEMPORARILY REQUIRED TO SET THE BOUNDARY CONDITIONS
      type(rscalar), private  :: Cond3D
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
      real (kind=prec), pointer, dimension(:,:,:), private    :: xXY, xXZ
      real (kind=prec), pointer, dimension(:,:,:), private    :: xY, xZ
      real (kind=prec), pointer, dimension(:,:,:), private    :: xXO
      real (kind=prec), pointer, dimension(:,:,:), private    :: yYX, yYZ
      real (kind=prec), pointer, dimension(:,:,:), private    :: yX, yZ
      real (kind=prec), pointer, dimension(:,:,:), private    :: yYO
      real (kind=prec), pointer, dimension(:,:,:), private    :: zZX, zZY
      real (kind=prec), pointer, dimension(:,:,:), private    :: zX, zY
      real (kind=prec), pointer, dimension(:,:,:), private    :: zZO

      ! coefficients of diagonal of (unweighted) A operator
      type (cvector_mg), private   :: Adiag
      ! information about the heritage ... probably this is not needed!
      real (kind=prec), private			:: whichOmega, whichCondE

      !!!!!!>>>>>>>>> FROM preconditioner
      ! coefficients of diagonal of preconditoner for A operator
      type (cvector_mg), private  :: Dilu

      !!!!!!>>>>>>>>> FROM divcorr
      ! coefficients for operators div sigma grad
      !  (which acts on scalars used for corner nodes),
      !  and the diagonal of the ilu preconditoner

      type (rvector_mg) , private	:: db1, db2
      !   db1%x contains coefficients of the stencil for shift -1 of ix index
      !    (%y,%z give corresponding coefficients for iy, iz)
      !   db2  contains coefficients for corresponding shift of +1

      type (rscalar_mg) , private	:: c, d
      ! c contains the coefficients for div sigma grad operator diagonal
      ! d contains the inverse of diagonal elements of D-ILU used for
      !  preconditioner
      type (rscalar_mg) , public   :: volC
      ! volume contains weights for corners

      ! which layer update
      character(len=10),parameter  :: first = 'first'
      character(len=10),parameter  :: last = 'last'
      character(len=10),parameter  :: lastDirect = 'lastDirect'

      type (timer_t), private        :: timer

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
      public                             	:: CurlCurlInit, CurlcurleCleanUp, Div, UpdateDivC, DivC,Grad
      public                                :: AdiagInit, AdiagSetUp, deall_Adiag, Maxwell,UpdateCurlCurl
      public                                :: MultA_N, AdjtBC
      !   These are used to initialize differential equation coefficients, and
      !    then to apply the differential operator

      ! routines from precondtioner
      public                      :: DiluInit, DiluSetUp, DeallocateDilu
      public                      :: M1Solve, M2Solve

      ! routines from divcorr
      public                :: DivCorrInit, DivCorrSetUp, UpdateDivCorr, Deallocate_DivCorr
      public                :: DivCgradILU, UpdateDivCgradILU, DivCgrad, UpdateDivCgrad

      ! interface for data_update  ... is this needed ?
      INTERFACE updateModelData
         module procedure UpdateFreq
         module procedure UpdateCond
         module procedure UpdateFreqCond
      END INTERFACE

      INTERFACE UpdateOperators
        module procedure UpdateCurlCurl
        module procedure UpdateDivC
        module procedure UpdateDivCorr
        module procedure UpdateDivCgradILU
      END INTERFACE

    Contains

      subroutine ModelDataInit(inGrid)
      !**********************************************************************
      ! * Copies multigrid to mgrid
      !   and/or compute variables stored in model_data module:
      !   create: allocate for edge conductivity volume weights;
      !   EdgeVolumeMG:  compute volume weights for edge-centered prisms
      !
      !**********************************************************************

      implicit none
        !  INPUTS:
        type (grid_t), intent(in)  :: inGrid  ! this is multigrid

        !   copy inGrid to mGrid
        call copy_mgrid(mGrid,inGrid)

        ! Allocate data structure for volume elements, and compute these
         call create(mGrid, volE, EDGE)      ! creates VolE as multigrid cvector

        call EdgeVolume(mGrid, volE)      ! in some cases we need volE defined on multigrid

        !  Allocate condE, conductivity defined on computational multigrid cell edges
        call create(mGrid, condE, EDGE)     ! creates condE defined on multigrid cell edges
        ! set a default omega
        omega = 0.0

      end subroutine ModelDataInit

    ! ********************************************************************************************************************8
      subroutine ModelDataCleanUp

        call deall_grid(mGrid)
        call deall(volE) !deall_rvector_mg
        call deall(condE) !deall_rvector_mg

        ! Cond3D is temporary
        call deall(Cond3D) !deall_rscalar

      end subroutine ModelDataCleanUp

      ! **************************************************************************
      ! * UpdateFreq updates the frequency that is currently being use
      subroutine UpdateFreq(inOmega)

        implicit none
        real (kind=prec), intent (in)             :: inOmega

        omega = inOmega

        Call AdiagSetUp()
#IFDEF Timer
      write (*,'(a12,a40,f12.6)')    node_info, 'AdiagSetUp: elapsed time (sec)', elapsed_time(Globaltimer)
#ENDIF
        Call DiluSetUp()
#IFDEF Timer
      write (*,'(a12,a40,f12.6)')    node_info, 'DiluSetUp: elapsed time (sec)', elapsed_time(Globaltimer)
#ENDIF

      end subroutine UpdateFreq  ! UpdateFreq

      ! ***************************************************************************
      ! * UpdateCond _updates the conductivity values on the edges
      subroutine UpdateCond(CondParam)

        implicit none
        type(modelParam_t), intent(in)      :: CondParam      ! input conductivity
        ! structure on the center of the grid

        !  ModelParamToEdge is to be interpreted as an abstract routine
        !    that maps from the external conductivity parameter to the
        !    internal edge representation (multi-grid)... the type of CondParam
        !    is now fixed as rscalar;  if a different representation is
        !    to be used changes to the declarations in this routine will
        !    be required, along with changes in the module interface
        Call ModelParamToEdge(CondParam, condE)
#IFDEF Timer
      write (*,'(a12,a40,f12.6)')    node_info, 'ModelParamtoEdge: elapsed time (sec)', elapsed_time(Globaltimer)
#ENDIF
        Call DivCorrSetUp()
#IFDEF Timer
      write (*,'(a12,a40,f12.6)')    node_info, 'DivCorrSetUp: elapsed time (sec)', elapsed_time(Globaltimer)
#ENDIF
        ! TEMPORARY; REQUIRED FOR BOUNDARY CONDITIONS
        !  set static array for cell conductivities
        call ModelParamToCell(CondParam,Cond3D)
#IFDEF Timer
      write (*,'(a12,a40,f12.6)')    node_info, 'ModelParamToCell: elapsed time (sec)', elapsed_time(Globaltimer)
#ENDIF
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
        !    internal edge representation (multi-grid) ...
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
      Subroutine SetBound(imode,period,E0mg,BC,iTx)

        !  Input mode, period
        integer, intent(in)		:: imode
        integer, intent(in)		:: iTx
        real(kind=prec)	:: period

        ! Output electric field first guess (for iterative solver)
        type(cvector_mg), intent(inout)	:: E0mg
        ! Output boundary conditions
        type(cboundary), intent(inout)	:: BC

        ! temporal cvector
        type(cvector)  :: TempE0

        call create(mgrid, TempE0, EDGE)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! FOR FWD PROBLEM (WITHOUT SENS) DECIDED TO CREATE E0 AS CVECTOR
       ! NEED TO THINK HOW TO WORK IN CASE of SENSITIVITIES !!!!!!!!!!

        if (BC%read_E_from_file) then

              call BC_x0_WS(imode,period,mGrid,Cond3D,TempE0,BC)
              ! The BC are already computed from a larger grid for all transmitters and modes and stored in BC_from_file.
              ! Overwrite BC with BC_from_file.
              ! Note: Right now we are using the same period layout for both grid.
              ! This why, it is enough to know the period and mode index to pick up the BC from BC_from_file vector.
              BC = BC_from_file((iTx*2)-(2-imode))
        else
             ! Compute the BC using Weerachai 2D approach
              call BC_x0_WS(imode,period,mGrid,Cond3D,TempE0,BC)
        end if

        ! need to average and copy fields from E0 cvector to cvector_mg
        e0mg = TempE0 ! c2mg

        ! Cell conductivity array is no longer needed
        ! NOT TRUE: needed for imode=2
        ! call deall_rscalar(Cond3D)

        call deall(TempE0)

      end subroutine SetBound
      ! *********************************************************************************************88

      ! * CurlcurlSetUp sets up all the coefficients for finite difference
      ! * approximations for del X del X E. In SetUp routines, one may do memory
      ! * allocation inside. SetUp routines does basic calculations
      ! * (maybe one time deal or sometimes more than once)

      subroutine CurlCurlInit()

      implicit none
        ! Output coefficients for curlcurlE (del X del X E)
        integer  :: ix, iy, iz, izv ! dummy variables
        integer  :: imgrid
        integer  :: nzCum, nx, ny,nz

        ! Output coefficients for curlcurlE (del X del X E)
        integer :: status     ! for dynamic memory allocation

        ! Allocate memory for del x del operator coefficient arrays
        ! Coefficients for difference equation only uses interior
        ! nodes. however, we need boundary nodes for the adjoint problem
        allocate(xXY(mgrid%mgridSize,mGrid%ny+1, 2), STAT=status)   ! Allocate memory
        allocate(xXZ(mgrid%mgridSize,mGrid%nz+1, 2), STAT=status)   ! Allocate memory
        allocate(xY(mgrid%mgridSize,mGrid%nx, mGrid%ny+1), STAT=status)
        allocate(xZ(mgrid%mgridSize,mGrid%nx, mGrid%nz+1), STAT=status)
        allocate(xXO(mgrid%mgridSize,mGrid%ny, mGrid%nz+1), STAT=status)

        allocate(yYZ(mgrid%mgridSize,mGrid%nz+1, 2), STAT=status)   ! Allocate memory
        allocate(yYX(mgrid%mgridSize,mGrid%nx+1, 2), STAT=status)   ! Allocate memory
        allocate(yZ(mgrid%mgridSize,mGrid%ny, mGrid%nz+1), STAT=status)
        allocate(yX(mgrid%mgridSize,mGrid%nx+1, mGrid%ny), STAT=status)
        allocate(yYO(mgrid%mgridSize,mGrid%nx, mGrid%nz+1), STAT=status)

        allocate(zZX(mgrid%mgridSize,mGrid%nx+1, 2), STAT=status)   ! Allocate memory
        allocate(zZY(mgrid%mgridSize,mGrid%ny+1, 2), STAT=status)   ! Allocate memory
        allocate(zX(mgrid%mgridSize,mGrid%nx+1, mGrid%nz), STAT=status)
        allocate(zY(mgrid%mgridSize,mGrid%ny+1, mGrid%nz), STAT=status)
        allocate(zZO(mgrid%mgridSize,mGrid%nx, mGrid%ny), STAT=status)

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

        nzCum = 0
        do imgrid = 1,mgrid%mgridSize  ! Global loop on multigrid
          nx = mgrid%gridArray(imgrid)%nx
          ny = mgrid%gridArray(imgrid)%ny
          nz = mgrid%gridArray(imgrid)%nz
          ! coefficents for calculating Ex ; only loop over internal edges
          do iy = 2, ny
             xXY(imgrid,iy,2) = -1.0/ (mgrid%gridArray(imgrid)%delY(iy) * mgrid%gridArray(imgrid)%dy(iy))
             xXY(imgrid,iy,1) = -1.0/ (mgrid%gridArray(imgrid)%delY(iy) * mgrid%gridArray(imgrid)%dy(iy-1))
          enddo

          if(imgrid == 1) then
            do iz = 2, nz+1
               izv = iz + nzCum
               xXZ(imgrid, iz, 2) = -1.0/ (mgrid%delZ(izv) * mgrid%dz(izv))
               xXZ(imgrid,iz, 1) = -1.0/ (mgrid%delZ(izv) * mgrid%dz(izv-1))
            enddo
          else
             do iz = 1, nz+1
               izv = iz + nzCum
               xXZ(imgrid, iz, 2) = -1.0/ (mgrid%delZ(izv) * mgrid%dz(izv))
               xXZ(imgrid,iz, 1) = -1.0/ (mgrid%delZ(izv) * mgrid%dz(izv-1))
            enddo
          endif

          do iy = 2, ny
            do iz = 1, nz+1
              xXO(imgrid, iy, iz) = -(xXY(imgrid, iy,1) + xXY(imgrid, iy,2) + &
                   xXZ(imgrid, iz,1) + xXZ(imgrid, iz,2))
           enddo
        enddo

        do ix = 1, nx
           do iy = 2, ny
              xY(imgrid, ix, iy) = 1.0/ (mgrid%gridArray(imgrid)%delY(iy)*mgrid%gridArray(imgrid)%dx(ix))
           enddo
        enddo

        do ix = 1, nx
           do iz = 1, nz+1
              xZ(imgrid, ix, iz) = 1.0/ (mgrid%gridArray(imgrid)%delZ(iz)*mgrid%gridArray(imgrid)%dx(ix))
           enddo
        enddo
        ! End of Ex coefficients

        ! coefficents for calculating Ey; only loop over internal edges
        if(imgrid == 1) then
          do iz = 2, nz+1
             izv = iz + nzCum
             yYZ(imgrid, iz, 2) = -1.0/ (mgrid%delZ(izv)*mgrid%dz(izv))
             yYZ(imgrid, iz, 1) = -1.0/ (mgrid%delZ(izv)*mgrid%dz(izv-1))
          enddo
        else
          do iz = 1, nz+1
             izv = iz + nzCum
             yYZ(imgrid, iz, 2) = -1.0/ (mgrid%delZ(izv)*mgrid%dz(izv))
             yYZ(imgrid, iz, 1) = -1.0/ (mgrid%delZ(izv)*mgrid%dz(izv-1))
          enddo
        endif

        do ix = 2, nx
           yYX(imgrid, ix, 2) = -1.0/ (mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%dx(ix))
           yYX(imgrid, ix, 1) = -1.0/ (mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%dx(ix-1))
        enddo

        do ix = 2, nx
           do iz = 1, nz+1
              yYO(imgrid, ix, iz) = -(yYX(imgrid,ix,1) + yYX(imgrid, ix,2) + &
                   yYZ(imgrid, iz,1) + yYZ(imgrid, iz,2))

           enddo
        enddo

        do iy = 1, ny
           do iz = 1, nz+1
              yZ(imgrid, iy, iz) = 1.0/ (mgrid%gridArray(imgrid)%delZ(iz)*mgrid%gridArray(imgrid)%dy(iy))
           enddo
        enddo


        do ix = 2, nx
           do iy = 1, ny
              yX(imgrid, ix, iy) = 1.0/ (mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%dy(iy))
           enddo
        enddo
        ! End of Ey coefficients

        ! coefficents for calculating Ez; only loop over internal edges
        do ix = 2, nx
           zZX(imgrid, ix, 2) = -1.0/ (mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%dx(ix))
           zZX(imgrid, ix, 1) = -1.0/ (mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%dx(ix-1))
        enddo

        do iy = 2, ny
           zZY(imgrid, iy, 2) = -1.0/ (mgrid%gridArray(imgrid)%delY(iy)*mgrid%gridArray(imgrid)%dy(iy))
           zZY(imgrid, iy, 1) = -1.0/ (mgrid%gridArray(imgrid)%delY(iy)*mgrid%gridArray(imgrid)%dy(iy-1))
        enddo

        do ix = 2, nx
           do iy = 2, ny
              zZO(imgrid, ix, iy) = -(zZX(imgrid,ix,1) + zZX(imgrid,ix,2) + &
                   zZY(imgrid,iy,1) + zZY(imgrid,iy,2))
           enddo
        enddo

        do ix = 2, nx
           do iz = 1, nz
              zX(imgrid,ix, iz) = 1.0/ (mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%dz(iz))
           enddo
        enddo

        do iy = 2, ny
           do iz = 1, nz
              zY(imgrid,iy, iz) = 1.0/ (mgrid%gridArray(imgrid)%delY(iy)*mgrid%gridArray(imgrid)%dz(iz))
           enddo
        enddo

        ! End of Ez coefficients
        nzCum = nz + nzCum
        enddo ! loop on subgrids

      end subroutine CurlCurlInit
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
      ! * Div computes the divergence for a complex vector define on multi-grid
      subroutine Div(inV, outSc)

        implicit none
        type (cvector_mg), intent(in)  :: inV
        type (cscalar_mg), intent(inout)  :: outSc
        ! local variables
        integer  :: ix, iy, iz, imgrid

        if(.not.inV%allocated) then
          print*, 'inV not allocated in Div'
          stop
        endif

        if(.not.outSc%allocated) then
          print*, 'outSc not allocated in Div'
          stop
        endif

        ! Check size
        if(inV%mgridSize /= outSc%mgridSize) then
          print*,  'Error-mgridSize input/ output in Div are not same size'
          stop
        endif

        do imgrid = 1, inV%mgridSize  ! Global loop over sub-grids
          ! Check whether all the inputs/ outputs involved
          ! are even of the same size
          if ((inV%cvArray(imgrid)%nx == outSc%csArray(imgrid)%nx).and.&
             (inV%cvArray(imgrid)%ny == outSc%csArray(imgrid)%ny).and.&
             (inV%cvArray(imgrid)%nz == outSc%csArray(imgrid)%nz)) then

           if ((inV%gridType == EDGE).and.(outSc%gridType == CORNER)) then
              call UpdateDiv(inV,outSc,imgrid)
              ! computation done only for internal nodes
              do ix = 2, outSc%csArray(imgrid)%nx
                 do iy = 2, outSc%csArray(imgrid)%ny
                    do iz = 2, outSc%csArray(imgrid)%nz

                       outSc%csArray(imgrid)%v(ix, iy, iz) = (inV%cvArray(imgrid)%x(ix, iy, iz) - &
                            inV%cvArray(imgrid)%x(ix - 1, iy, iz))/ mgrid%gridArray(imgrid)%delX(ix) + &
                            (inV%cvArray(imgrid)%y(ix, iy, iz) - inV%cvArray(imgrid)%y(ix, iy - 1, iz))/&
                            mgrid%gridArray(imgrid)%delY(iy) + &
                            (inV%cvArray(imgrid)%z(ix, iy, iz) - inV%cvArray(imgrid)%z(ix, iy, iz - 1))/&
                            mgrid%gridArray(imgrid)%delZ(iz)

                    enddo   ! iz
                 enddo      ! iy
              enddo         ! ix

           else if ((inV%gridType == FACE).and.(outSc%gridType == CENTER)) then
              call UpdateDiv(inV,outSc,imgrid)
              ! computation done only for internal nodes
              ! there us no problems with update z in this case,
              ! because cvector_mg allocated nz+1 in gridType == FACE
              do ix = 1, outSc%csArray(imgrid)%nx
                 do iy = 1, outSc%csArray(imgrid)%ny
                    do iz = 1, outSc%csArray(imgrid)%nz

                       outSc%csArray(imgrid)%v(ix, iy, iz) = (inV%cvArray(imgrid)%x(ix+1, iy, iz) - &
                            inV%cvArray(imgrid)%x(ix, iy, iz))/ mgrid%gridArray(imgrid)%dx(ix) + &
                            (inV%cvArray(imgrid)%y(ix, iy+1, iz) - inV%cvArray(imgrid)%y(ix, iy, iz))/&
                            mgrid%gridArray(imgrid)%dy(iy) + &
                            (inV%cvArray(imgrid)%z(ix, iy, iz+1) - inV%cvArray(imgrid)%z(ix, iy, iz))/&
                            mgrid%gridArray(imgrid)%dz(iz)

                    enddo   ! iz
                 enddo      ! iy
              enddo         ! ix

           else
              print*, 'Div: not compatible usage for existing data types'
           end if

         else
           print*, 'Error-all input/ output in Div are not same size'
         endif
       enddo ! Global loop over sub-grids

      end subroutine Div  ! Div
      ! *********************************************************************
      !Update first layer in Div
      subroutine UpdateDiv(inV,outSc,imgrid)

        implicit none

        type (cscalar_mg), intent(inout)   :: outSc
        type (cvector_mg), intent(in)  :: inV
        integer, intent(in)  :: imgrid
        ! local variables
        integer  :: ix, iy, iz

          ! Since computations are done for the internal nodes
          ! we do not need to update the first sub-grid
          if(imgrid == 1) then
            return
          endif

          ! the basic strategy is to fill in the interface layer with values from the finer grid!
          if (mGrid%coarseness(imgrid).gt.mGrid%coarseness(imgrid-1))then
            ! interface layer: finer to coarser
            ! current sub-grid is coarser
            ! take values from nz of the previous sub-grid
           do iy = 2, mGrid%gridArray(imgrid)%ny
              do ix = 2, mGrid%gridArray(imgrid)%nx
                outSc%csArray(imgrid)%v(ix, iy, 1) = (inV%cvArray(imgrid)%x(ix, iy, 1) - &
                            inV%cvArray(imgrid)%x(ix-1, iy, 1))/ mgrid%gridArray(imgrid)%delX(ix) + &
                            (inV%cvArray(imgrid)%y(ix, iy, 1) - inV%cvArray(imgrid)%y(ix, iy-1, 1))/&
                            mgrid%gridArray(imgrid)%delY(iy) + &
                            (inV%cvArray(imgrid)%z(ix, iy, 1) - inV%cvArray(imgrid-1)%z(2*ix-1, 2*iy-1, mGrid%gridArray(imgrid-1)%nz))/&
                            mgrid%gridArray(imgrid)%delZ(1)
               enddo
           enddo
          else if (mGrid%coarseness(imgrid).lt.mGrid%coarseness(imgrid-1)) then
           ! interface : coarser to finer
           ! current grid is finer
           ! vice versa
           do iy = 2, mGrid%gridArray(imgrid-1)%ny
             do ix = 2, mGrid%gridArray(imgrid-1)%nx
               outSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, 1) = (inV%cvArray(imgrid)%x(2*ix-1, 2*iy-1, 1) - &
                            inV%cvArray(imgrid)%x(2*(ix-1), 2*iy-1, 1))/ mgrid%gridArray(imgrid)%delX(2*ix-1) + &
                            (inV%cvArray(imgrid)%y(2*ix-1, 2*iy-1, 1) - inV%cvArray(imgrid)%y(2*ix-1, 2*(iy-1), 1))/&
                            mgrid%gridArray(imgrid)%delY(2*iy-1) + &
                            (inV%cvArray(imgrid)%z(2*ix-1, 2*iy-1, 1) - inV%cvArray(imgrid-1)%z(ix, iy, mGrid%gridArray(imgrid-1)%nz))/&
                            mgrid%gridArray(imgrid)%delZ(1)
              enddo
            enddo

          else if (mGrid%coarseness(imgrid).eq.mGrid%coarseness(imgrid-1)) then
           do iy = 2, mGrid%gridArray(imgrid)%ny
             do ix = 2, mGrid%gridArray(imgrid)%nx
               outSc%csArray(imgrid)%v(ix, iy, 1) = (inV%cvArray(imgrid)%x(ix, iy, 1) - &
                            inV%cvArray(imgrid)%x(ix - 1, iy, 1))/ mgrid%gridArray(imgrid)%delX(ix) + &
                            (inV%cvArray(imgrid)%y(ix, iy, 1) - inV%cvArray(imgrid)%y(ix, iy - 1, 1))/&
                            mgrid%gridArray(imgrid)%delY(iy) + &
                            (inV%cvArray(imgrid)%z(ix, iy, 1) - inV%cvArray(imgrid-1)%z(ix, iy, mGrid%gridArray(imgrid-1)%nz))/&
                            mgrid%gridArray(imgrid)%delZ(1)
             enddo
           enddo
         endif

      end subroutine UpdateDiv
      ! ***************************************************************************
      ! * Grad computes the gradient for a complex scalar

      subroutine Grad(inSc, outV)

      implicit none
        type (cscalar_mg), intent(in)  :: inSc
        type (cvector_mg), intent(inout)  :: outV
        integer  :: ix, iy, iz, imgrid
        integer  :: nx,ny,nz

        if(.not.inSc%allocated) then
          print *, 'inSc not allocated in Grad'
          stop
        endif

        if(.not.outV%allocated) then
          print *,'outV not allocated in Grad'
          stop
        endif

          do imgrid = 1, inSc%mgridSize  ! Global loop over subgrids
            nx = inSc%csArray(imgrid)%nx
            ny = inSc%csArray(imgrid)%ny
            nz = inSc%csArray(imgrid)%nz

            if ((nx == outV%cvArray(imgrid)%nx).and.&
                (ny == outV%cvArray(imgrid)%ny).and.&
                (nz == outV%cvArray(imgrid)%nz)) then

              if ((outV%gridType == EDGE).and.(inSc%gridType == CORNER).and.&
                                       (outV%mgridSize == inSc%mgridSize)) then

                ! the conversion in Grad is only done for interior nodes
                do ix = 1, nx
                 do iy = 2, ny
                    do iz = 1, nz+1

                       outV%cvArray(imgrid)%x(ix, iy, iz) = (inSc%csArray(imgrid)%v(ix+1, iy, iz) - &
                            inSc%csArray(imgrid)%v(ix, iy, iz))/ inSc%csArray(imgrid)%grid%dx(ix)

                    enddo
                 enddo
                enddo

                do ix = 2, nx
                  do iy = 1, ny
                   do iz = 1, nz+1

                    outV%cvArray(imgrid)%y(ix, iy, iz) = (inSc%csArray(imgrid)%v(ix, iy+1, iz) - &
                         inSc%csArray(imgrid)%v(ix, iy, iz))/ inSc%csArray(imgrid)%grid%dy(iy)

                   enddo
                  enddo
                enddo

                do ix = 2, nx
                  do iy = 2, ny
                   do iz = 1, nz

                     outV%cvArray(imgrid)%z(ix, iy, iz) = (inSc%csArray(imgrid)%v(ix, iy, iz+1) - &
                            inSc%csArray(imgrid)%v(ix, iy, iz))/ inSc%csArray(imgrid)%grid%dz(iz)

                    enddo
                   enddo
                enddo

              else if ((outV%gridType == FACE).and.(inSc%gridType == CENTER)) then

                ! the conversion in Grad is only done for interior nodes
                do ix = 2, nx
                 do iy = 1, ny
                    do iz = 1, nz+1

                       outV%cvArray(imgrid)%x(ix, iy, iz) = (inSc%csArray(imgrid)%v(ix, iy, iz) - &
                            inSc%csArray(imgrid)%v(ix-1, iy, iz))/ inSc%csArray(imgrid)%grid%delX(ix)

                    enddo
                 enddo
                enddo

                do ix = 1, nx
                 do iy = 2, ny
                    do iz = 1, nz+1

                       outV%cvArray(imgrid)%y(ix, iy, iz) = (inSc%csArray(imgrid)%v(ix, iy, iz) - &
                            inSc%csArray(imgrid)%v(ix, iy-1, iz))/ inSc%csArray(imgrid)%grid%delY(iy)

                    enddo
                 enddo
                enddo

                do ix = 1, nx
                 do iy = 1, ny
                    do iz = 2, nz

                       outV%cvArray(imgrid)%z(ix, iy, iz) = (inSc%csArray(imgrid)%v(ix, iy, iz) - &
                            inSc%csArray(imgrid)%v(ix, iy, iz-1))/ inSc%csArray(imgrid)%grid%delZ(iz)

                    enddo
                 enddo
                enddo

            else
              print *, 'Error-all input/ output in Grad are not same size'
            endif


          else
            print *, 'Grad: not compatible usage for existing data types'
          endif
        enddo ! Global loop over subgrids
      end subroutine Grad
      ! ***************************************************************************
      ! * AdiagInit initializes the memory for the diagonal nodes being added
      ! * with the imaginary part. Init routines mostly do memory allocation,
      ! * reading and setting up the data
      subroutine AdiagInit()

        implicit none

       ! Call create_cvector(mGrid, Adiag, EDGE)
        Call create(mGrid, Adiag, EDGE)

      end subroutine AdiagInit

      ! ***************************************************************************
      ! * Adiag sets up the diagonal nodes with the imaginary part added to it
      ! * SetUp routines do basic calculations (maybe one time deal or more than one)
      subroutine AdiagSetUp()

        implicit  none
        integer                   :: ix, iy, iz , imgrid      ! dummy variables
        integer  :: nx,ny,nz

        if (.not.Adiag%allocated) then
          print *, 'Adiag in AdiagSetUp not allocated yet'
          stop
        end if

          ! check size
          if(mgrid%mgridSize /= Adiag%mgridSize) then
            print *, 'AdiagSetUp error: mgridSize'
            else
              do imgrid = 1, mgrid%mgridSize
                nx = mgrid%gridArray(imgrid)%nx
                ny = mgrid%gridArray(imgrid)%ny
                nz = mgrid%gridArray(imgrid)%nz
                if ((nx /= Adiag%cvArray(imgrid)%nx.or.ny /= Adiag%cvArray(imgrid)%ny.or.&
                                ny /= Adiag%cvArray(imgrid)%ny).and.(nx /= condE%rvArray(imgrid)%nx.or.ny /= condE%rvArray(imgrid)%ny.or.&
                                ny /= condE%rvArray(imgrid)%ny)) then
                   print *, 'AdiagSetUp error: grids are not the same size'
                else

                  do ix = 1, nx
                    Adiag%cvArray(imgrid)%x(ix,:,:) = CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%rvArray(imgrid)%x(ix,:,:)
                  enddo
                  do iy = 1, ny
                    Adiag%cvArray(imgrid)%y(:,iy,:) = CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%rvArray(imgrid)%y(:,iy,:)
                  enddo
                  do iz = 1, nz
                    Adiag%cvArray(imgrid)%z(:,:,iz) = CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%rvArray(imgrid)%z(:,:,iz)
                  enddo

                endif
               enddo
            endif

      end subroutine AdiagSetUp

    ! *************************************************************************************
      ! * deall_Adiag deallocates the memory for the diagonal nodes being added
      ! * with the imaginary part.
      subroutine deall_Adiag()

        implicit none
        Call deall(Adiag)
      end subroutine deall_Adiag

        ! ***************************************************************************
      ! * Maxwell computes the finite difference equation in complex vectors
      ! * for del X del X E +/- i*omega*mu*conductivity*E in unsymmetrical form.
      ! * Note that the difference equation is only for interior edges. However,
      ! * it does use the contribution from the boundary edges. The coefficients
      ! * are  calculated in CurlcurleSetUp. Remember, in the operators that are
      ! * used in iterative fashion, output is always initialized outside
       subroutine Maxwell(inE, adjt, outE)   ! Multigrid

        implicit none
        type (cvector_mg), intent(in)               :: inE
        ! input electrical field as complex vector
        logical, intent (in)                     :: adjt
        type (cvector_mg), target, intent(inout)    :: outE
        ! output electrical field as complex vector
        integer  :: diag_sign
        integer  :: ix, iy, iz, imgrid
        integer  :: nx, ny,nz
        ! dummy variables

        if (.not.inE%allocated) then
          write(0,*) 'inE in Maxwell not allocated yet'
          stop
        end if

        if (.not.outE%allocated) then
          write(0,*) 'outE in Maxwell not allocated yet'
          stop
        end if

        if (inE%mgridSize == outE%mgridSize) then
          ! Check grid type also
          if ((inE%gridType == outE%gridType)) then
            do imgrid = 1, inE%mgridSize    ! Global loop on subgrids
              ! Check whether the bounds are the same
              if ((inE%cvArray(imgrid)%nx == outE%cvArray(imgrid)%nx).and.&
               (inE%cvArray(imgrid)%ny == outE%cvArray(imgrid)%ny).and.&
                (inE%cvArray(imgrid)%nz == outE%cvArray(imgrid)%nz)) then
               nx = inE%cvArray(imgrid)%nx
               ny = inE%cvArray(imgrid)%ny
               nz = inE%cvArray(imgrid)%nz
               if (adjt) then
                  diag_sign = -1*ISIGN
               else
                  diag_sign = ISIGN
               end if

              !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ix,iy,iz)
              ! Apply difference equation to compute Ex (only on interior nodes)
              ! the diagonal nodes have the imaginary component added
              !$OMP DO SCHEDULE(STATIC)
              do iz = 2, nz
                 do iy = 2, ny
                    do ix = 1, nx
                       outE%cvArray(imgrid)%x(ix,iy,iz) = xY(imgrid,ix,iy)*(inE%cvArray(imgrid)%y(ix+1,iy,iz)-&
                        inE%cvArray(imgrid)%y(ix,iy,iz)-inE%cvArray(imgrid)%y(ix+1,iy-1,iz)&
                            +inE%cvArray(imgrid)%y(ix,iy-1,iz))+&
                        xZ(imgrid,ix,iz)*(inE%cvArray(imgrid)%z(ix+1,iy,iz)-inE%cvArray(imgrid)%z(ix,iy,iz)&
                        -inE%cvArray(imgrid)%z(ix+1,iy,iz-1)+inE%cvArray(imgrid)%z(ix,iy,iz-1))+&
                        xXY(imgrid,iy,2)*inE%cvArray(imgrid)%x(ix,iy+1,iz)+&
                        xXY(imgrid,iy,1)*inE%cvArray(imgrid)%x(ix,iy-1,iz)+&
                        xXZ(imgrid,iz ,2)*inE%cvArray(imgrid)%x(ix,iy,iz+1)+&
                        xXZ(imgrid,iz,1)*inE%cvArray(imgrid)%x(ix,iy,iz-1)+&
                        (xXO(imgrid,iy,iz)+diag_sign*Adiag%cvArray(imgrid)%x(ix,iy,iz))*inE%cvArray(imgrid)%x(ix,iy,iz)
                    enddo
                 enddo
              enddo
              !$OMP END DO NOWAIT

              ! Apply difference equation to compute Ey (only on interior nodes)
              ! the diagonal nodes have the imaginary component added
              !$OMP DO SCHEDULE(STATIC)
              do iz = 2, nz
                 do iy = 1, ny
                    do ix = 2, nx
                       outE%cvArray(imgrid)%y(ix,iy,iz) = yZ(imgrid,iy,iz)*(inE%cvArray(imgrid)%z(ix,iy+1,iz)-&
                        inE%cvArray(imgrid)%z(ix,iy,iz)-inE%cvArray(imgrid)%z(ix,iy+1,iz-1)+inE%cvArray(imgrid)%z(ix,iy,iz-1))&
                        +yX(imgrid,ix,iy)*(inE%cvArray(imgrid)%x(ix,iy+1,iz)-inE%cvArray(imgrid)%x(ix,iy,iz)&
                        -inE%cvArray(imgrid)%x(ix-1,iy+1,iz)+inE%cvArray(imgrid)%x(ix-1,iy,iz))+&
                        yYZ(imgrid,iz,2)*inE%cvArray(imgrid)%y(ix,iy,iz+1)+&
                        yYZ(imgrid,iz,1)*inE%cvArray(imgrid)%y(ix,iy,iz-1)+&
                        yYX(imgrid,ix,2)*inE%cvArray(imgrid)%y(ix+1,iy,iz)+&
                        yYX(imgrid,ix,1)*inE%cvArray(imgrid)%y(ix-1,iy,iz)+&
                        (yYO(imgrid,ix,iz)+diag_sign*Adiag%cvArray(imgrid)%y(ix,iy,iz))*inE%cvArray(imgrid)%y(ix,iy,iz)
                    enddo
                 enddo
              enddo
              !$OMP END DO NOWAIT

              ! Apply difference equation to compute Ey (only on interior nodes)
              ! the diagonal nodes have the imaginary component added
              !$OMP DO SCHEDULE(STATIC)
              do iz = 1, nz
                 do iy = 2, ny
                    do ix = 2, nx
                       outE%cvArray(imgrid)%z(ix,iy,iz) = zX(imgrid,ix,iz)*(inE%cvArray(imgrid)%x(ix,iy,iz+1)-&
                          inE%cvArray(imgrid)%x(ix,iy,iz)-inE%cvArray(imgrid)%x(ix-1,iy,iz+1)+inE%cvArray(imgrid)%x(ix-1,iy,iz))&
                        +zY(imgrid,iy,iz)*(inE%cvArray(imgrid)%y(ix,iy,iz+1)-inE%cvArray(imgrid)%y(ix,iy,iz)&
                        -inE%cvArray(imgrid)%y(ix,iy-1,iz+1)+inE%cvArray(imgrid)%y(ix,iy-1,iz))+&
                        zZX(imgrid,ix,2)*inE%cvArray(imgrid)%z(ix+1,iy,iz)+&
                        zZX(imgrid,ix,1)*inE%cvArray(imgrid)%z(ix-1,iy,iz)+&
                        zZY(imgrid,iy,2)*inE%cvArray(imgrid)%z(ix,iy+1,iz)+&
                        zZY(imgrid,iy,1)*inE%cvArray(imgrid)%z(ix,iy-1,iz)+&
                        (zZO(imgrid,ix,iy)+diag_sign*Adiag%cvArray(imgrid)%z(ix,iy,iz))*inE%cvArray(imgrid)%z(ix,iy,iz)
                    enddo
                 enddo
              enddo
              !$OMP END DO NOWAIT
              !$OMP END PARALLEL

              ! update first and last layers
              call UpdateCurlCurl(inE,outE,diag_sign,imgrid)

              else
                print *, 'Error-complex vectors for Maxwell are not of same size'
              end if
            enddo ! Global loop on subgrids

          else
            print *, 'Maxwell: not compatible usage for existing data types'
           end if

        else
          print *, 'Error complex vectors for Maxwell'
        endif

      end subroutine Maxwell        ! Maxwell  Multigrid

      ! ******************************************************************************************************
      subroutine UpdateCurlCurl(inE, outE,diag_sign,imgrid)

      ! Updates first z layer in curl curl
               !last
      implicit none

      type(cvector_mg), intent(inout)  :: outE
      type(cvector_mg), intent(in)  :: inE
      integer, intent(in)  :: imgrid
      integer  :: diag_sign

      ! local variables
      integer  :: ix, iy, iz
      integer  :: nx, ny, nz

        nx = inE%cvArray(imgrid)%nx
        ny = inE%cvArray(imgrid)%ny
        nz = inE%cvArray(imgrid)%nz

       if(imgrid == inE%mgridSize) then
       ! the last multigrid layer
       ! z=1 must be filled in
       ! nz+1 should be zero
       ! nothing to do
         return
       endif

       ! the basic strategy is to fill in the interface layer with values from the finer grid!
       if (inE%coarseness(imgrid).lt.inE%coarseness(imgrid+1))then
          ! interface layer: finer -> coarser
          ! Update nz+1 curl curl operator

          do iy = 2, inE%cvarray(imgrid+1)%ny
            do ix = 1, inE%cvarray(imgrid+1)%nx
              ! Update Ex
              ! 'odd values'
              outE%cvArray(imgrid)%x(2*ix-1,2*iy-1,nz+1) = xY(imgrid,2*ix-1,2*iy-1)*(inE%cvArray(imgrid)%y(2*ix,2*iy-1,nz+1)-&
                        inE%cvArray(imgrid)%y(2*ix-1,2*iy-1,nz+1)-inE%cvArray(imgrid)%y(2*ix,2*(iy-1),nz+1)&
                            +inE%cvArray(imgrid)%y(2*ix-1,2*(iy-1),nz+1))+&
                        xZ(imgrid,2*ix-1,nz)*(inE%cvArray(imgrid)%z(2*ix,2*iy-1,nz)-inE%cvArray(imgrid)%z(2*ix-1,2*iy-1,nz))&
                        -xZ(imgrid+1,ix,1)*(inE%cvArray(imgrid+1)%z(ix+1,iy,1)-inE%cvArray(imgrid+1)%z(ix,iy,1))+&
                        xXY(imgrid,2*iy-1,2)*inE%cvArray(imgrid)%x(2*ix-1,2*iy,nz+1)+&
                        xXY(imgrid,2*iy-1,1)*inE%cvArray(imgrid)%x(2*ix-1,2*(iy-1),nz+1)+&
                        xXZ(imgrid,nz+1,1)*inE%cvArray(imgrid)%x(2*ix-1,2*iy-1,nz)+&
                        xXZ(imgrid,nz+1,2)*(inE%cvArray(imgrid+1)%x(ix,iy,2)-inE%cvArray(imgrid+1)%x(ix,iy,1))+&
                        (-(xXY(imgrid,2*iy-1,1)+xXY(imgrid,2*iy-1,2)+xXZ(imgrid,nz+1,1))+diag_sign*Adiag%cvArray(imgrid)%x(2*ix-1,2*iy-1,nz+1))*inE%cvArray(imgrid)%x(2*ix-1,2*iy-1,nz+1)
              ! 'even values'
              outE%cvArray(imgrid)%x(2*ix,2*iy-1,nz+1) = xY(imgrid,2*ix,2*iy-1)*(inE%cvArray(imgrid)%y(2*ix+1,2*iy-1,nz+1)-&
                        inE%cvArray(imgrid)%y(2*ix,2*iy-1,nz+1)-inE%cvArray(imgrid)%y(2*ix+1,2*(iy-1),nz+1)&
                            +inE%cvArray(imgrid)%y(2*ix,2*(iy-1),nz+1))+&
                        xZ(imgrid,2*ix,nz)*(inE%cvArray(imgrid)%z(2*ix+1,2*iy-1,nz)-inE%cvArray(imgrid)%z(2*ix,2*iy-1,nz))&
                        -xZ(imgrid+1,ix,1)*(inE%cvArray(imgrid+1)%z(ix+1,iy,1)-inE%cvArray(imgrid+1)%z(ix,iy,1))+&
                        xXY(imgrid,2*iy-1,2)*inE%cvArray(imgrid)%x(2*ix,2*iy,nz+1)+&
                        xXY(imgrid,2*iy-1,1)*inE%cvArray(imgrid)%x(2*ix,2*(iy-1),nz+1)+&
                        xXZ(imgrid,nz+1,1)*inE%cvArray(imgrid)%x(2*ix,2*iy-1,nz)+&
                        xXZ(imgrid,nz+1,2)*(inE%cvArray(imgrid+1)%x(ix,iy,2)-inE%cvArray(imgrid+1)%x(ix,iy,1))+&
                        (-(xXY(imgrid,2*iy-1,1)+xXY(imgrid,2*iy-1,2)+xXZ(imgrid,nz+1,1))+diag_sign*Adiag%cvArray(imgrid)%x(2*ix,2*iy-1,nz+1))*inE%cvArray(imgrid)%x(2*ix,2*iy-1,nz+1)

              ! fill in first layer of the next multigrid layer
              outE%cvArray(imgrid+1)%x(ix,iy,1) =  (outE%cvArray(imgrid)%x(2*ix-1,2*iy-1,nz+1)*mGrid%gridArray(imgrid)%dx(2*ix-1)&
                                                  +outE%cvArray(imgrid)%x(2*ix,2*iy-1,nz+1)*mGrid%gridArray(imgrid)%dx(2*ix))/&
                                                  (mGrid%gridArray(imgrid)%dx(2*ix-1)+mGrid%gridArray(imgrid)%dx(2*ix))
            enddo
          enddo

         ! Update Ey
         do iy = 1, inE%cvarray(imgrid+1)%ny
           do ix = 2, inE%cvarray(imgrid+1)%nx

             ! 'odd values'
                       outE%cvArray(imgrid)%y(2*ix-1,2*iy-1,nz+1) = yZ(imgrid,2*iy-1,nz)*(inE%cvArray(imgrid)%z(2*ix-1,2*iy,nz)-&
                        inE%cvArray(imgrid)%z(2*ix-1,2*iy-1,nz))&
                        -yZ(imgrid+1,iy,1)*(inE%cvArray(imgrid+1)%z(ix,iy+1,1)-inE%cvArray(imgrid+1)%z(ix,iy,1))&
                        +yX(imgrid,2*ix-1,2*iy-1)*(inE%cvArray(imgrid)%x(2*ix-1,2*iy,nz+1)-inE%cvArray(imgrid)%x(2*ix-1,2*iy-1,nz+1)&
                        -inE%cvArray(imgrid)%x(2*(ix-1),2*iy,nz+1)+inE%cvArray(imgrid)%x(2*(ix-1),2*iy-1,nz+1))+&
                        yYZ(imgrid,nz+1,1)*inE%cvArray(imgrid)%y(2*ix-1,2*iy-1,nz)+&
                        yYZ(imgrid,nz+1,2)*(inE%cvArray(imgrid+1)%y(ix,iy,2)-inE%cvArray(imgrid+1)%y(ix,iy,1))+&
                        yYX(imgrid,2*ix-1,2)*inE%cvArray(imgrid)%y(2*ix,2*iy-1,nz+1)+&
                        yYX(imgrid,2*ix-1,1)*inE%cvArray(imgrid)%y(2*(ix-1),2*iy-1,nz+1)+&
                        (-(yYX(imgrid,2*ix-1,1)+yYZ(imgrid,nz+1,1)+yYX(imgrid,2*ix-1,2))+diag_sign*Adiag%cvArray(imgrid)%y(2*ix-1,2*iy-1,nz+1))*inE%cvArray(imgrid)%y(2*ix-1,2*iy-1,nz+1)

             ! 'even values'
                       outE%cvArray(imgrid)%y(2*ix-1,2*iy,nz+1) = yZ(imgrid,2*iy,nz)*(inE%cvArray(imgrid)%z(2*ix-1,2*iy+1,nz)-&
                        inE%cvArray(imgrid)%z(2*ix-1,2*iy,nz))&
                        -yZ(imgrid+1,iy,1)*(inE%cvArray(imgrid+1)%z(ix,iy+1,1)-inE%cvArray(imgrid+1)%z(ix,iy,1))&
                        +yX(imgrid,2*ix-1,2*iy)*(inE%cvArray(imgrid)%x(2*ix-1,2*iy+1,nz+1)-inE%cvArray(imgrid)%x(2*ix-1,2*iy,nz+1)&
                        -inE%cvArray(imgrid)%x(2*(ix-1),2*iy+1,nz+1)+inE%cvArray(imgrid)%x(2*(ix-1),2*iy,nz+1))+&
                        yYZ(imgrid,nz+1,1)*inE%cvArray(imgrid)%y(2*ix-1,2*iy,nz)+&
                        yYZ(imgrid,nz+1,2)*(inE%cvArray(imgrid+1)%y(ix,iy,2)-inE%cvArray(imgrid+1)%y(ix,iy,1))+&
                        yYX(imgrid,2*ix-1,2)*inE%cvArray(imgrid)%y(2*ix,2*iy,nz+1)+&
                        yYX(imgrid,2*ix-1,1)*inE%cvArray(imgrid)%y(2*(ix-1),2*iy,nz+1)+&
                        (-(yYX(imgrid,2*ix-1,1)+yYZ(imgrid,nz+1,1)+yYX(imgrid,2*ix-1,2))+diag_sign*Adiag%cvArray(imgrid)%y(2*ix-1,2*iy,nz+1))*inE%cvArray(imgrid)%y(2*ix-1,2*iy,nz+1)

              ! fill in zero layer of the next multigrid layer
              outE%cvArray(imgrid+1)%y(ix,iy,1) =  (outE%cvArray(imgrid)%y(2*ix-1,2*iy-1,nz+1)*mGrid%gridArray(imgrid)%dy(2*iy-1)&
                                                  +outE%cvArray(imgrid)%y(2*ix-1,2*iy,nz+1)*mGrid%gridArray(imgrid)%dy(2*iy))/&
                                                  (mGrid%gridArray(imgrid)%dy(2*iy-1)+mGrid%gridArray(imgrid)%dy(2*iy))
            enddo
          enddo
       else if (inE%coarseness(imgrid).gt.inE%coarseness(imgrid+1))then

         ! interface layer: coarse to fine
          do iy = 2, ny
            do ix = 1, nx
              ! Update Ex
              ! 'odd values'
              outE%cvArray(imgrid+1)%x(2*ix-1,2*iy-1,1) = xY(imgrid+1,2*ix-1,2*iy-1)*(inE%cvArray(imgrid+1)%y(2*ix,2*iy-1,1)-&
                        inE%cvArray(imgrid+1)%y(2*ix-1,2*iy-1,1)-inE%cvArray(imgrid+1)%y(2*ix,2*(iy-1),1)&
                            +inE%cvArray(imgrid+1)%y(2*ix-1,2*(iy-1),1))+&
                        xZ(imgrid+1,2*ix-1,1)*(inE%cvArray(imgrid+1)%z(2*ix,2*iy-1,1)-inE%cvArray(imgrid+1)%z(2*ix-1,2*iy-1,1))&
                        -xZ(imgrid,ix,nz)*(inE%cvArray(imgrid)%z(ix+1,iy,nz)-inE%cvArray(imgrid)%z(ix,iy,nz))+&
                        xXY(imgrid+1,2*iy-1,2)*inE%cvArray(imgrid+1)%x(2*ix-1,2*iy,1)+&
                        xXY(imgrid+1,2*iy-1,1)*inE%cvArray(imgrid+1)%x(2*ix-1,2*(iy-1),1)+&
                        xXZ(imgrid+1,1,2)*inE%cvArray(imgrid+1)%x(2*ix-1,2*iy-1,2)+&
                        xXZ(imgrid+1,1,1)*(inE%cvArray(imgrid)%x(ix,iy,nz)-inE%cvArray(imgrid)%x(ix,iy,nz+1))+&
                        (-(xXY(imgrid+1,2*iy-1,1)+xXY(imgrid+1,2*iy-1,2)+xXZ(imgrid+1,1,2))+diag_sign*Adiag%cvArray(imgrid+1)%x(2*ix-1,2*iy-1,1))*inE%cvArray(imgrid+1)%x(2*ix-1,2*iy-1,1)
              ! 'even values'
              outE%cvArray(imgrid+1)%x(2*ix,2*iy-1,1) = xY(imgrid+1,2*ix,2*iy-1)*(inE%cvArray(imgrid+1)%y(2*ix+1,2*iy-1,1)-&
                        inE%cvArray(imgrid+1)%y(2*ix,2*iy-1,1)-inE%cvArray(imgrid+1)%y(2*ix+1,2*(iy-1),1)&
                            +inE%cvArray(imgrid+1)%y(2*ix,2*(iy-1),1))+&
                        xZ(imgrid+1,2*ix,1)*(inE%cvArray(imgrid+1)%z(2*ix+1,2*iy-1,1)-inE%cvArray(imgrid+1)%z(2*ix,2*iy-1,1))&
                        -xZ(imgrid,ix,nz)*(inE%cvArray(imgrid)%z(ix+1,iy,nz)-inE%cvArray(imgrid)%z(ix,iy,nz))+&
                        xXY(imgrid+1,2*iy-1,2)*inE%cvArray(imgrid+1)%x(2*ix,2*iy,1)+&
                        xXY(imgrid+1,2*iy-1,1)*inE%cvArray(imgrid+1)%x(2*ix,2*(iy-1),1)+&
                        xXZ(imgrid+1,1,2)*inE%cvArray(imgrid+1)%x(2*ix,2*iy-1,2)+&
                        xXZ(imgrid+1,1,1)*(inE%cvArray(imgrid)%x(ix,iy,nz)-inE%cvArray(imgrid)%x(ix,iy,nz+1))+&
                        (-(xXY(imgrid+1,2*iy-1,1)+xXY(imgrid+1,2*iy-1,2)+xXZ(imgrid+1,1,2))+diag_sign*Adiag%cvArray(imgrid+1)%x(2*ix,2*iy-1,1))*inE%cvArray(imgrid+1)%x(2*ix,2*iy-1,1)

              ! fill in zero layer of the next multigrid layer
              outE%cvArray(imgrid)%x(ix,iy,nz+1) =  (outE%cvArray(imgrid+1)%x(2*ix-1,2*iy-1,1)*mGrid%gridArray(imgrid+1)%dx(2*ix-1)&
                                                   +outE%cvArray(imgrid+1)%x(2*ix,2*iy-1,1)*mGrid%gridArray(imgrid+1)%dx(2*ix))/&
                                                   (mGrid%gridArray(imgrid+1)%dx(2*ix-1)+mGrid%gridArray(imgrid+1)%dx(2*ix))
            enddo
          enddo

          ! Update Ey
          do iy = 1, ny
            do ix = 2, nx
             ! 'odd values'
                       outE%cvArray(imgrid+1)%y(2*ix-1,2*iy-1,1) = yZ(imgrid+1,2*iy-1,1)*(inE%cvArray(imgrid+1)%z(2*ix-1,2*iy,1)-&
                        inE%cvArray(imgrid+1)%z(2*ix-1,2*iy-1,1))&
                        -yZ(imgrid,iy,nz)*(inE%cvArray(imgrid)%z(ix,iy+1,nz)-inE%cvArray(imgrid)%z(ix,iy,nz))&
                        +yX(imgrid+1,2*ix-1,2*iy-1)*(inE%cvArray(imgrid+1)%x(2*ix-1,2*iy,1)-inE%cvArray(imgrid+1)%x(2*ix-1,2*iy-1,1)&
                        -inE%cvArray(imgrid+1)%x(2*(ix-1),2*iy,1)+inE%cvArray(imgrid+1)%x(2*(ix-1),2*iy-1,1))+&
                        yYZ(imgrid+1,1,2)*inE%cvArray(imgrid+1)%y(2*ix-1,2*iy-1,2)+&
                        yYZ(imgrid+1,1,1)*(inE%cvArray(imgrid)%y(ix,iy,nz)-inE%cvArray(imgrid)%y(ix,iy,nz+1))+&
                        yYX(imgrid+1,2*ix-1,2)*inE%cvArray(imgrid+1)%y(2*ix,2*iy-1,1)+&
                        yYX(imgrid+1,2*ix-1,1)*inE%cvArray(imgrid+1)%y(2*(ix-1),2*iy-1,1)+&
                        (-(yYX(imgrid+1,2*ix-1,1)+yYX(imgrid+1,2*ix-1,2)+yYZ(imgrid+1,1,2))+diag_sign*Adiag%cvArray(imgrid+1)%y(2*ix-1,2*iy-1,1))*inE%cvArray(imgrid+1)%y(2*ix-1,2*iy-1,1)
             ! 'even values'
                       outE%cvArray(imgrid+1)%y(2*ix-1,2*iy,1) = yZ(imgrid+1,2*iy,1)*(inE%cvArray(imgrid+1)%z(2*ix-1,2*iy+1,1)-&
                        inE%cvArray(imgrid+1)%z(2*ix-1,2*iy,1))&
                        -yZ(imgrid,iy,nz)*(inE%cvArray(imgrid)%z(ix,iy+1,nz)-inE%cvArray(imgrid)%z(ix,iy,nz))&
                        +yX(imgrid+1,2*ix-1,2*iy)*(inE%cvArray(imgrid+1)%x(2*ix-1,2*iy+1,1)-inE%cvArray(imgrid+1)%x(2*ix-1,2*iy,1)&
                        -inE%cvArray(imgrid+1)%x(2*(ix-1),2*iy+1,1)+inE%cvArray(imgrid+1)%x(2*(ix-1),2*iy,1))+&
                        yYZ(imgrid+1,1,2)*inE%cvArray(imgrid+1)%y(2*ix-1,2*iy,2)+&
                        yYZ(imgrid+1,1,1)*(inE%cvArray(imgrid)%y(ix,iy,nz)-inE%cvArray(imgrid)%y(ix,iy,nz+1))+&
                        yYX(imgrid+1,2*ix-1,2)*inE%cvArray(imgrid+1)%y(2*ix,2*iy,1)+&
                        yYX(imgrid+1,2*ix-1,1)*inE%cvArray(imgrid+1)%y(2*(ix-1),2*iy,1)+&
                        (-(yYX(imgrid+1,2*ix-1,1)+yYX(imgrid+1,2*ix-1,2)+yYZ(imgrid+1,1,2))+diag_sign*Adiag%cvArray(imgrid+1)%y(2*ix-1,2*iy,1))*inE%cvArray(imgrid+1)%y(2*ix-1,2*iy,1)

              ! fill in zero layer of the next multigrid layer
              outE%cvArray(imgrid)%y(ix,iy,nz+1) =  (outE%cvArray(imgrid+1)%y(2*ix-1,2*iy-1,1)*mGrid%gridArray(imgrid+1)%dy(2*iy-1)&
                                                   +outE%cvArray(imgrid+1)%y(2*ix-1,2*iy,1)*mGrid%gridArray(imgrid+1)%dy(2*iy))/&
                                                   (mGrid%gridArray(imgrid+1)%dy(2*iy-1)+mGrid%gridArray(imgrid+1)%dy(2*iy))
            enddo
          enddo

         ! if coarseness == 0,0,0,0,0 ...
         ! we get the regular grid
         ! it was needed for debugging
         ! may also be needed ?!

         else if (inE%coarseness(imgrid).eq.inE%coarseness(imgrid+1))then
          do iy = 2, ny
            do ix= 1, nx
              ! standard curl curl for interface
              ! curl curl Ex
              outE%cvArray(imgrid)%x(ix,iy,nz+1) = &
                  xY(imgrid,ix,iy)*(inE%cvArray(imgrid)%y(ix+1,iy,nz+1)-&
                                    inE%cvArray(imgrid)%y(ix,iy,nz+1)-inE%cvArray(imgrid)%y(ix+1,iy-1,nz+1)&
                                    +inE%cvArray(imgrid)%y(ix,iy-1,nz+1))+&
                  xZ(imgrid,ix,nz+1)*(inE%cvArray(imgrid+1)%z(ix+1,iy,1)-inE%cvArray(imgrid+1)%z(ix,iy,1)&
                        -inE%cvArray(imgrid)%z(ix+1,iy,nz)+inE%cvArray(imgrid)%z(ix,iy,nz))+&
                  xXY(imgrid,iy,2)*inE%cvArray(imgrid)%x(ix,iy+1,nz+1)+&
                  xXY(imgrid,iy,1)*inE%cvArray(imgrid)%x(ix,iy-1,nz+1)+&
                  xXZ(imgrid,nz+1,2)*inE%cvArray(imgrid+1)%x(ix,iy,2)+&
                  xXZ(imgrid,nz+1,1)*inE%cvArray(imgrid)%x(ix,iy,nz)+&
                  (xXO(imgrid,iy,nz+1)+diag_sign*Adiag%cvArray(imgrid)%x(ix,iy,nz+1))*inE%cvArray(imgrid)%x(ix,iy,nz+1)

              outE%cvArray(imgrid+1)%x(ix,iy,1) = outE%cvArray(imgrid)%x(ix,iy,nz+1)
            enddo
          enddo
          do iy = 1, ny
            do ix = 2, nx
              ! standard curl curl for interface
              ! curl curl Ey
                outE%cvArray(imgrid)%y(ix,iy,nz+1) =  &
                    yZ(imgrid,iy,nz+1)*(inE%cvArray(imgrid+1)%z(ix,iy+1,1)-inE%cvArray(imgrid+1)%z(ix,iy,1)-&
                                        inE%cvArray(imgrid)%z(ix,iy+1,nz)+inE%cvArray(imgrid)%z(ix,iy,nz))+&
                    yX(imgrid,ix,iy)*(inE%cvArray(imgrid)%x(ix,iy+1,nz+1)-inE%cvArray(imgrid)%x(ix,iy,nz+1)&
                                     -inE%cvArray(imgrid)%x(ix-1,iy+1,nz+1)+inE%cvArray(imgrid)%x(ix-1,iy,nz+1))+&
                    yYZ(imgrid,nz+1,2)*inE%cvArray(imgrid+1)%y(ix,iy,2)+&
                    yYZ(imgrid,nz+1,1)*inE%cvArray(imgrid)%y(ix,iy,nz)+&
                    yYX(imgrid,ix,2)*inE%cvArray(imgrid)%y(ix+1,iy,nz+1)+&
                    yYX(imgrid,ix,1)*inE%cvArray(imgrid)%y(ix-1,iy,nz+1)+&
                    (yYO(imgrid,ix,nz+1)+diag_sign*Adiag%cvArray(imgrid)%y(ix,iy,nz+1))*inE%cvArray(imgrid)%y(ix,iy,nz+1)

              outE%cvArray(imgrid+1)%y(ix,iy,1) = outE%cvArray(imgrid)%y(ix,iy,nz+1)
            enddo
          enddo
        endif

      end subroutine UpdateCurlCurl
      ! ***************************************************************************
      ! * Gets the Maxwell's equation in the complete symmetrical form,
      ! * del X del X E +/- i*omega*mu*conductivity*E. E is the complex vector
      ! * defining the electrical field _N is to denote that this is the new
      ! * subroutine where the imaginary part at the at the diagonal is inbuilt
      ! *  Diagonally multiplied by weights for symmetry.
      subroutine MultA_N(inE, adjt, outE)

        implicit none
        type (cvector_mg), intent (in)  :: inE    !  inE(:) array of cvectors on multigrid
        logical, intent (in)  :: adjt
        type (cvector_mg), intent (inout) :: outE !  outE(:) array of cvectors on multigird

        !local variables
        integer  :: imgrid
        integer  :: nx,ny,nz

        if (.not.inE%allocated) then
          print *, 'inE in MultA_N not allocated yet'
          stop
        end if

        if (.not.outE%allocated) then
          print *, 'out in MultA_N not allocated yet'
          stop
        end if

         Call Maxwell(inE, adjt, outE)
         ! done with preparing del X del X E +/- i*omega*mu*conductivity*E
         ! diagonally multiply the final results with weights (edge volume)
         Call diagMult(outE, volE, outE)


      end subroutine MultA_N
    ! *************************************************************************************
      subroutine AdjtBC(eIn, BC)

     ! modified by Cherevatova (Aug,2012)
     ! for the multi-grid
     ! BC here are defined on the multi-grid
     ! and belong to extended type cboundary_mg


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

        ! INPUT: electrical fields stored as cvector_mg
        type (cvector_mg), intent(in)             		:: eIn
        ! OUTPUT: boundary condition structure: should be allocated
        !   and initialized before call to this routine
        type (cboundary_mg),intent(inout)  		    	:: BC
        !local
        integer                   :: imgrid,ix,iy,iz,nx,ny,nz

        !  Multiply FD electric field vector defined on interior nodes (eIn) by
        !  adjoint of A_IB, the interior/boundary sub-block of the differential
        !  operator.

        !Ex components in x/z plane (iy=1; iy = ny+1)
        !NOTE: C_ZERO = (0,0) (double complex) is defined in SG_Basics/math_constants.f90

        do imgrid = 1,BC%mgridSize ! Main loop over sub-grid

          nx = BC%bcArray(imgrid)%nx
          ny = BC%bcArray(imgrid)%ny
          nz = BC%bcArray(imgrid)%nz

            do ix = 1, nx
               do iz = 1, nz+1
                  BC%bcArray(imgrid)%xYmin(ix,iz) = - yX(imgrid,ix,1)*eIn%cvarray(imgrid)%y(ix,1,iz)       &
                                    + yX(imgrid,ix+1,1)*eIn%cvArray(imgrid)%y(ix+1,1,iz)   &
                                    + xXY(imgrid,2,1)*eIn%cvArray(imgrid)%x(ix,2,iz)
                  BC%bcArray(imgrid)%xYmax(ix,iz) = + yX(imgrid,ix,ny)*eIn%cvArray(imgrid)%y(ix,ny,iz)     &
                                    - yX(imgrid,ix+1,ny)*eIn%cvArray(imgrid)%y(ix+1,ny,iz) &
                                    + xXY(imgrid,ny,2)*eIn%cvArray(imgrid)%x(ix,ny,iz)
                enddo ! iz
             enddo ! ix

            !Ez components in x/z plane (iy=1; iy = ny+1)
            BC%bcArray(imgrid)%zYMin(1,:) = C_ZERO
            BC%bcArray(imgrid)%zYmax(1,:) = C_ZERO
            BC%bcArray(imgrid)%zYmin(nx+1,:) = C_ZERO
            BC%bcArray(imgrid)%zYmax(nx+1,:) = C_ZERO
            do iz = 1, nz
               do ix = 2, nx
                  BC%bcArray(imgrid)%zYmin(ix,iz) = - yZ(imgrid,1,iz)*eIn%cvArray(imgrid)%y(ix,1,iz)        &
                                    + yZ(imgrid,1,iz+1)*eIn%cvArray(imgrid)%y(ix,1,iz+1)    &
                                    + zZY(imgrid,2,1)*eIn%cvArray(imgrid)%z(ix,2,iz)
                  BC%bcArray(imgrid)%zYmax(ix,iz) = + yZ(imgrid,ny,iz)*eIn%cvArray(imgrid)%y(ix,ny,iz)      &
                                    - yZ(imgrid,ny,iz+1)*eIn%cvArray(imgrid)%y(ix,ny,iz+1)  &
                                    + zZY(imgrid,ny,2)*eIn%cvArray(imgrid)%z(ix,ny,iz)
                enddo
             enddo

            !Ey components in y/z plane (ix=1; ix = nx+1)
            do iy = 1, ny
               do iz = 1, nz+1
                  BC%bcArray(imgrid)%yXmin(iy,iz) = - xY(imgrid,1,iy)*eIn%cvArray(imgrid)%x(1,iy,iz)        &
                                    + xY(imgrid,1,iy+1)*eIn%cvArray(imgrid)%x(1,iy+1,iz)    &
                                    + yYX(imgrid,2,1)*eIn%cvArray(imgrid)%y(2,iy,iz)
                  BC%bcArray(imgrid)%yXmax(iy,iz) = + xY(imgrid,nx,iy)*eIn%cvArray(imgrid)%x(nx,iy,iz)      &
                                    - xY(imgrid,nx,iy+1)*eIn%cvArray(imgrid)%x(nx,iy+1,iz)  &
                                    + yYX(imgrid,nx,2)*eIn%cvArray(imgrid)%y(nx,iy,iz)
                enddo
             enddo

            !Ez components in y/z plane (ix=1; ix = nx+1)
            BC%bcArray(imgrid)%zXmin(1,:) = C_ZERO
            BC%bcArray(imgrid)%zXmax(1,:) = C_ZERO
            BC%bcArray(imgrid)%zXmin(ny+1,:) = C_ZERO
            BC%bcArray(imgrid)%zXmax(ny+1,:) = C_ZERO
            do iz = 1, nz
               do iy = 2, ny
                  BC%bcArray(imgrid)%zXmin(iy,iz) = - xZ(imgrid,1,iz)*eIn%cvArray(imgrid)%x(1,iy,iz)       &
                                    + xZ(imgrid,1,iz+1)*eIn%cvArray(imgrid)%x(1,iy,iz+1)   &
                                    + zZX(imgrid,2,1)*eIn%cvArray(imgrid)%z(2,iy,iz)
                  BC%bcArray(imgrid)%zXmax(iy,iz) = + xZ(imgrid,nx,iz)*eIn%cvArray(imgrid)%x(nx,iy,iz)     &
                                    - xZ(imgrid,nx,iz+1)*eIn%cvArray(imgrid)%x(nx,iy,iz+1) &
                                    + zZX(imgrid,nx,2)*eIn%cvArray(imgrid)%z(nx,iy,iz)
                enddo
             enddo

            !Ex components in x/y plane (iz=1; iz = nz+1)
            BC%bcArray(imgrid)%xZmin(:,1) = C_ZERO
            BC%bcArray(imgrid)%xZmax(:,1) = C_ZERO
            BC%bcArray(imgrid)%xZmin(:,ny+1) = C_ZERO
            BC%bcArray(imgrid)%xZmax(:,ny+1) = C_ZERO
            do ix = 1, nx
               do iy = 2, ny
                  BC%bcArray(imgrid)%xZmin(ix,iy) = - zX(imgrid,ix,1)*eIn%cvArray(imgrid)%z(ix,iy,1)       &
                                    + zX(imgrid,ix+1,1)*eIn%cvArray(imgrid)%z(ix+1,iy,1)   &
                                    + xXZ(imgrid,2,1)*eIn%cvArray(imgrid)%x(ix,iy,2)
                  BC%bcArray(imgrid)%xZmax(ix,iy) = + zX(imgrid,ix,nz)*eIn%cvArray(imgrid)%z(ix,iy,nz)     &
                                    - zX(imgrid,ix+1,nz)*eIn%cvArray(imgrid)%z(ix+1,iy,nz) &
                                    + xXZ(imgrid,nz,2)*eIn%cvArray(imgrid)%x(ix,iy,nz)
                enddo
             enddo

            !Ey components in x/y plane (iz=1; iz = nz+1)
            BC%yZmin(1,:) = C_ZERO
            BC%yZmax(1,:) = C_ZERO
            BC%yZmin(nx+1,:) = C_ZERO
            BC%yZmin(nx+1,:) = C_ZERO
            do iy = 1, ny
               do ix = 2, nx
                  BC%bcArray(imgrid)%yZmin(ix,iy) = - zY(imgrid,iy,1)*eIn%cvArray(imgrid)%z(ix,iy,1)        &
                                    + zY(imgrid,iy+1,1)*eIn%cvArray(imgrid)%z(ix,iy+1,1)    &
                                    + yYZ(imgrid,2,1)*eIn%cvArray(imgrid)%y(ix,iy,2)
                  BC%bcArray(imgrid)%yZmax(ix,iy) = + zY(imgrid,iy,nz)*eIn%cvArray(imgrid)%z(ix,iy,nz)      &
                                    - zY(imgrid,iy+1,nz)*eIn%cvArray(imgrid)%z(ix,iy+1,nz)  &
                                    + yYZ(imgrid,nz,2)*eIn%cvArray(imgrid)%y(ix,iy,nz)
                enddo
             enddo
      enddo ! end Main loop over sub-grids
        BC%bcArray(1)%xYMax(:,1) = C_ZERO
        BC%bcArray(1)%xYmin(:,1) = C_ZERO
        BC%bcArray(BC%mgridSize)%xYMax(:,BC%bcArray(BC%mgridSize)%nz+1) = C_ZERO
        BC%bcArray(BC%mgridSize)%xyMin(:,BC%bcArray(BC%mgridSize)%nz+1) = C_ZERO

        BC%bcArray(BC%mgridSize)%yXmin(:,1) = C_ZERO
        BC%bcArray(BC%mgridSize)%yXmax(:,1) = C_ZERO
        BC%bcArray(BC%mgridSize)%yXmin(:,BC%bcArray(BC%mgridSize)%nz+1) = C_ZERO
        BC%bcArray(BC%mgridSize)%yXmax(:,BC%bcArray(BC%mgridSize)%nz+1) = C_ZERO

      end subroutine AdjtBC

      ! ****************************************************************************
      ! PRECONDITIONER ROUTINES: set up ILU-Level I preconditioner for
      ! Maxwell's equation, solve lower and upper triangular systems to
      ! apply preconditioner
      !****************************************************************************
      ! initializes a diagonal of preconditioner for A operator
      subroutine DiluInit()

        implicit none
        integer  :: status
        integer  :: ix, iy, iz, imgrid


        if (.not.Dilu%allocated) then

           Call create(mGrid, Dilu, EDGE) ! allocate cvector_mg

        else

           if (Dilu%mgridSize /= mGrid%mgridSize) then
             deallocate(Dilu%cvArray, STAT = status)
             Call create(mGrid, Dilu, EDGE)
           else
             do imgrid = 1, mgrid%mgridSize
                if ((Dilu%cvArray(imgrid)%nx /= mGrid%gridArray(imgrid)%nx).or. &
                (Dilu%cvArray(imgrid)%ny /= mGrid%gridArray(imgrid)%ny).or.&
                (Dilu%cvArray(imgrid)%nz /= mGrid%gridArray(imgrid)%nz)) exit
                deallocate(Dilu%cvArray, STAT = status)
                Call create(mGrid, Dilu, EDGE)
             enddo
           end if
        end if

      end subroutine DiluInit ! DiluInit

      !****************************************************************************
      ! sets up a diagonal of preconditioner for A operator
      subroutine DiluSetUp()

        implicit none
        integer   :: status
        integer   :: ix, iy, iz, imgrid
        integer  :: nx, ny, nz

        if (.not.Dilu%allocated) then
          print *, 'Dilu not allocated yet'
        else

        ! initializing the non-interior values
        ! only the interior edge values are really used
        Dilu%cvArray(1)%x(:,:,1) = cmplx(1.0, 0.0, 8)
        Dilu%cvArray(1)%y(:,:,1) = cmplx(1.0, 0.0, 8)

        do imgrid = 1, mgrid%mgridSize  ! Global loop on subgrids

          Dilu%cvArray(imgrid)%x(:,1,:) = cmplx(1.0, 0.0, 8)
          Dilu%cvArray(imgrid)%y(1,:,:) = cmplx(1.0, 0.0, 8)
          Dilu%cvArray(imgrid)%z(1,:,:) = cmplx(1.0, 0.0, 8)
          Dilu%cvArray(imgrid)%z(:,1,:) = cmplx(1.0, 0.0, 8)

          if ((Dilu%cvArray(imgrid)%nx /= mGrid%gridArray(imgrid)%nx).or. &
             (Dilu%cvArray(imgrid)%ny /= mGrid%gridArray(imgrid)%ny).or.&
             (Dilu%cvArray(imgrid)%nz /= mGrid%gridArray(imgrid)%nz)) then
             print *,'Dilu that is right now existing has the wrong size'
          endif

          nx = mGrid%gridArray(imgrid)%nx
          ny = mGrid%gridArray(imgrid)%ny
          nz = mGrid%gridArray(imgrid)%nz

          ! update first layer in Dilu
          call UpdateZ(Dilu,first,imgrid)

          do ix = 1, nx
            do iz = 2, nz+1
              do iy = 2, ny

               Dilu%cvArray(imgrid)%x(ix, iy, iz) = xXO(imgrid,iy,iz) - &
                      CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%rvArray(imgrid)%x(ix, iy, iz)  &
                      - xXY(imgrid,iy, 1)*xXY(imgrid,iy-1, 2)*Dilu%cvArray(imgrid)%x(ix,iy-1,iz) &
                      - xXZ(imgrid,iz, 1)*xXZ(imgrid,iz-1, 2)*Dilu%cvArray(imgrid)%x(ix,iy,iz-1)
               Dilu%cvArray(imgrid)%x(ix, iy, iz) = 1.0/ Dilu%cvArray(imgrid)%x(ix, iy, iz)
              enddo
            enddo
          enddo

          ! the coefficients for y are only for the interior nodes
          !  but need to initialize edges for recursive algorithm
          do iy = 1, ny
            do iz = 2, nz+1
              do ix = 2, nx

                 Dilu%cvArray(imgrid)%y(ix, iy, iz) = yYO(imgrid,ix,iz) - &
                     CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%rvArray(imgrid)%y(ix, iy, iz) &
                     - yYZ(imgrid,iz, 1)*yYZ(imgrid,iz-1, 2)*Dilu%cvArray(imgrid)%y(ix, iy, iz-1) &
                     - yYX(imgrid,ix, 1)*yYX(imgrid,ix-1, 2)*Dilu%cvArray(imgrid)%y(ix-1, iy, iz)
                 Dilu%cvArray(imgrid)%y(ix, iy, iz) = 1.0/ Dilu%cvArray(imgrid)%y(ix, iy, iz)

              enddo
            enddo
          enddo

          ! the coefficients for z are only for the interior nodes
          !  but need to initialize edges for recursive algorithm
          do iz = 1, nz
            do ix = 2, ny
              do iy = 2, nx

               Dilu%cvArray(imgrid)%z(ix, iy, iz) = zZO(imgrid,ix,iy) - &
                    CMPLX(0.0, 1.0, 8)*omega*MU_0*condE%rvArray(imgrid)%z(ix, iy, iz) &
                    - zZX(imgrid,ix, 1)*zZX(imgrid,ix-1, 2)*Dilu%cvArray(imgrid)%z(ix-1, iy, iz) &
                    - zZY(imgrid,iy, 1)*zZY(imgrid,iy-1, 2)*Dilu%cvArray(imgrid)%z(ix, iy-1, iz)
               Dilu%cvArray(imgrid)%z(ix, iy, iz) = 1.0/ Dilu%cvArray(imgrid)%z(ix, iy, iz)

              enddo
            enddo
          enddo

          enddo ! Global loop on subgrid
          ! last nz+1 must be 1 and 0
          Dilu%cvArray(mGrid%mgridSize)%x(:,:,mGrid%gridArray(mGrid%mgridSize)%nz+1) = C_ZERO
          Dilu%cvArray(mGrid%mgridSize)%y(:,:,mGrid%gridArray(mGrid%mgridSize)%nz+1) = C_ZERO
          Dilu%cvArray(mGrid%mgridSize)%x(:,1,mGrid%gridArray(mGrid%mgridSize)%nz+1) = cmplx(1.0, 0.0, 8)
          Dilu%cvArray(mGrid%mgridSize)%y(1,:,mGrid%gridArray(mGrid%mgridSize)%nz+1) = cmplx(1.0, 0.0, 8)

       endif
      end subroutine DiluSetUp  ! DiluSetUp

      !****************************************************************************
      !  To Deallocate arrays in structure Dilu
      subroutine DeallocateDilu()
        implicit none

        call deall(Dilu) !deall_cvector_mg
      end subroutine DeallocateDilu  ! DeallocateDilu


      !****************************************************************************
      ! Purpose: to solve the lower triangular system (or it's adjoint);
      ! for the d-ilu pre-condtioner.

      subroutine M1solve(inE, adjt, outE)

        implicit none
        type (cvector_mg), intent(in)  :: inE
        logical, intent(in)  :: adjt
        type (cvector_mg), intent(inout)  :: outE
        ! local variables
        integer  :: ix, iy, iz, imgrid
        integer  :: nx, ny, nz

        if (.not.inE%allocated) then
          write(0,*) 'inE in M1solve not allocated yet'
          stop
        end if

        if (.not.outE%allocated) then
          write(0,*) 'outE in M1solve not allocated yet'
          stop
        end if

        ! Check whether the multigrid sizes are the same
        if (inE%mgridSize == outE%mgridSize) then
          ! Check grid type also
          if ((inE%gridType == outE%gridType)) then

          if (.not.adjt) then
          !adjoint = .false.

            call diagDiv(inE, volE, outE)
              !$OMP PARALLEL DEFAULT(SHARED)
              do imgrid = 1, inE%mgridSize    ! Direct loop on subgrids
                ! Check whether the bounds are the same
                if ((inE%cvArray(imgrid)%nx == outE%cvArray(imgrid)%nx).and.&
                   (inE%cvArray(imgrid)%ny == outE%cvArray(imgrid)%ny).and.&
                   (inE%cvArray(imgrid)%nz == outE%cvArray(imgrid)%nz)) then
                  nx = inE%cvArray(imgrid)%nx
                  ny = inE%cvArray(imgrid)%ny
                  nz = inE%cvArray(imgrid)%nz

                  ! update first layer
                  call UpdateZ(outE,first,imgrid)
                 ! ... note that we only parallelize the outer loops
                 !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
                  do ix = 1, nx
                    !$OMP ORDERED
                    do iz = 2, nz+1
                       do iy = 2, ny

                          outE%cvArray(imgrid)%x(ix, iy, iz) = (outE%cvArray(imgrid)%x(ix, iy, iz) - &
                               outE%cvArray(imgrid)%x(ix, iy-1, iz)*xXY(imgrid,iy, 1) - &
                               outE%cvArray(imgrid)%x(ix, iy, iz-1)*xXZ(imgrid,iz, 1))* &
                               Dilu%cvArray(imgrid)%x(ix, iy, iz)

                        enddo
                     enddo
                    !$OMP END ORDERED
                   enddo
                  !$OMP END DO

                 !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
                 do iy = 1, ny
                    !$OMP ORDERED
                    do iz = 2, nz+1
                       do ix = 2, nx

                          outE%cvArray(imgrid)%y(ix, iy, iz) = (outE%cvArray(imgrid)%y(ix, iy, iz) - &
                               outE%cvArray(imgrid)%y(ix, iy, iz-1)*yYZ(imgrid,iz, 1) - &
                               outE%cvArray(imgrid)%y(ix-1, iy, iz)*yYX(imgrid,ix, 1))* &
                               Dilu%cvArray(imgrid)%y(ix, iy, iz)

                       enddo
                    enddo
                    !$OMP END ORDERED
                 enddo
                 !$OMP END DO

                 !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
                 do iz = 1, nz
                    !$OMP ORDERED
                    do iy = 2, ny
                       do ix = 2, nx
                          outE%cvArray(imgrid)%z(ix, iy, iz) = (outE%cvArray(imgrid)%z(ix, iy, iz) - &
                               outE%cvArray(imgrid)%z(ix-1, iy, iz)*zZX(imgrid,ix, 1) - &
                               outE%cvArray(imgrid)%z(ix, iy-1, iz)*zZY(imgrid,iy, 1))* &
                               Dilu%cvArray(imgrid)%z(ix, iy, iz)
                       enddo
                    enddo
                    !$OMP END ORDERED
                 enddo
                 !$OMP END DO
                else
                  print *, 'Error-complex vectors for M1Solve are not of same size'
                endif
              enddo ! Direct loop on subgrids
              !$OMP END PARALLEL
           ! to be sure that the last layer is equal zero
           outE%cvArray(inE%mgridSize)%x(:,:,inE%cvArray(inE%mgridSize)%nz+1) = 0.0
           outE%cvArray(inE%mgridSize)%y(:,:,inE%cvArray(inE%mgridSize)%nz+1) = 0.0

           else
             ! adjoint = .true.
              !$OMP PARALLEL DEFAULT(SHARED)
              do imgrid = inE%mgridSize, 1, -1 ! Reverse loop on subgrids
                nx = inE%cvArray(imgrid)%nx
                ny = inE%cvArray(imgrid)%ny
                nz = inE%cvArray(imgrid)%nz

                call UpdateZ(outE, last, imgrid)

                ! ... note that we only parallelize the outer loops
                ! the coefficients for x are only for the interior nodes
                !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
                do ix = 1, nx
                  !$OMP ORDERED
                  do iy = ny, 2, -1
                     do iz = nz, 1, -1

                          outE%cvArray(imgrid)%x(ix, iy, iz) = (inE%cvArray(imgrid)%x(ix, iy, iz) - &
                               outE%cvArray(imgrid)%x(ix, iy+1, iz)*xXY(imgrid,iy+1, 1) - &
                               outE%cvArray(imgrid)%x(ix, iy, iz+1)*xXZ(imgrid,iz+1, 1))* &
                               conjg(Dilu%cvArray(imgrid)%x(ix, iy, iz))

                     enddo
                  enddo
                  !$OMP END ORDERED
                enddo
                !$OMP END DO

                ! the coefficients for y are only for the interior nodes
                !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
                do iy = 1, ny
                  !$OMP ORDERED
                  do ix = nx, 2, -1
                     do iz = nz, 1, -1

                          outE%cvArray(imgrid)%y(ix, iy, iz) = (inE%cvArray(imgrid)%y(ix, iy, iz) - &
                               outE%cvArray(imgrid)%y(ix, iy, iz+1)*yYZ(imgrid,iz+1, 1) - &
                               outE%cvArray(imgrid)%y(ix+1, iy, iz)*yYX(imgrid,ix+1, 1))* &
                               conjg(Dilu%cvArray(imgrid)%y(ix, iy, iz))

                     enddo
                  enddo
                  !$OMP END ORDERED
                enddo
                !$OMP END DO
             enddo !  Reverse loop on subgrids
                !$OMP END PARALLEL

             do imgrid = 1, inE%mgridSize
                nx = inE%cvArray(imgrid)%nx
                ny = inE%cvArray(imgrid)%ny
                nz = inE%cvArray(imgrid)%nz
                !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
                do iz = 1, nz
                  !$OMP ORDERED
                  do ix = nx, 2, -1
                    do iy = ny, 2, -1

                          outE%cvArray(imgrid)%z(ix, iy, iz) = (inE%cvArray(imgrid)%z(ix, iy, iz) - &
                               outE%cvArray(imgrid)%z(ix+1, iy, iz)*zZX(imgrid,ix+1, 1) - &
                               outE%cvArray(imgrid)%z(ix, iy+1, iz)*zZY(imgrid,iy+1, 1))* &
                               conjg(Dilu%cvArray(imgrid)%z(ix, iy, iz))

                     enddo
                  enddo
                  !$OMP END ORDERED
                enddo
                !$OMP END DO

    enddo

             call diagDiv(outE, volE, outE)

           endif

          else
            print *, 'M1 solver: not compatible usage for existing data types'
          end if

        else
          print *, 'Error complex vectors for M1 solver; multigridSize are not of same size'
        endif

      end subroutine M1solve ! M1solve

      !****************************************************************************
      ! Purpose: to solve the upper triangular system (or it's adjoint);
      ! for the d-ilu pre-condtioner
      subroutine M2solve(inE, adjt, outE)

        implicit none
        type (cvector_mg), intent(in)	:: inE
        logical, intent(in)		:: adjt
        type (cvector_mg), intent(inout) 	:: outE
        integer  :: ix, iy, iz, imgrid
        integer  :: nx, ny, nz

        if (.not.inE%allocated) then
          write(0,*) 'inE in M2solve not allocated yet'
          stop
        end if

        if (.not.outE%allocated) then
          write(0,*) 'outE in M2solve not allocated yet'
          stop
        end if

        ! Check whether the multigrid sizes are the same
        if (inE%mgridSize == outE%mgridSize) then
          ! Check grid type also
          if ((inE%gridType == outE%gridType)) then
            if (.not.adjt) then
            !adjoint = .false.
              ! ... note that we only parallelize the outer loops
              !$OMP PARALLEL DEFAULT(SHARED)
              do imgrid = inE%mgridSize, 1, -1    ! Reverse loop on subgrids
                ! Check whether the bounds are the same
                if ((inE%cvArray(imgrid)%nx == outE%cvArray(imgrid)%nx).and.&
                (inE%cvArray(imgrid)%ny == outE%cvArray(imgrid)%ny).and.&
                (inE%cvArray(imgrid)%nz == outE%cvArray(imgrid)%nz)) then

                ! update nz+1 layer
                call UpdateZ(outE, last, imgrid)

                 !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
                   do ix = 1, inE%cvArray(imgrid)%nx
                      !$OMP ORDERED
                      do iz = inE%cvArray(imgrid)%nz, 1, -1
                         do iy = inE%cvArray(imgrid)%ny, 2, -1

                            outE%cvArray(imgrid)%x(ix, iy, iz) = inE%cvArray(imgrid)%x(ix, iy, iz) - &
                               ( outE%cvArray(imgrid)%x(ix, iy+1, iz)*xXY(imgrid,iy, 2) &
                               + outE%cvArray(imgrid)%x(ix, iy, iz+1)*xXZ(imgrid,iz, 2))* &
                               Dilu%cvArray(imgrid)%x(ix, iy, iz)
                         enddo
                      enddo
                    !$OMP END ORDERED
                   enddo
                  !$OMP END DO

                   !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
                   do iy = 1, inE%cvArray(imgrid)%ny
                      !$OMP ORDERED
                      do iz = inE%cvArray(imgrid)%nz, 1, -1
                         do ix = inE%cvArray(imgrid)%nx, 2, -1

                            outE%cvArray(imgrid)%y(ix, iy, iz) = inE%cvArray(imgrid)%y(ix, iy, iz) - &
                                 ( outE%cvArray(imgrid)%y(ix, iy, iz+1)*yYZ(imgrid,iz, 2) &
                                 + outE%cvArray(imgrid)%y(ix+1, iy, iz)*yYX(imgrid,ix, 2))* &
                                 Dilu%cvArray(imgrid)%y(ix, iy, iz)
                         enddo
                      enddo
                     !$OMP END ORDERED
                   enddo
                   !$OMP END DO
                else
                  print *, 'Error-complex vectors for M2 solver are not of same size'
                end if
              enddo ! Reverse loop on subgrids
              !$OMP END PARALLEL
              !$OMP PARALLEL DEFAULT(SHARED)
              do imgrid = 1, inE%mgridSize ! Direct loop on subgrids
                   !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
                   do iz = 1, inE%cvArray(imgrid)%nz
                      !$OMP ORDERED
                      do iy = inE%cvArray(imgrid)%ny, 2, -1
                         do ix = inE%cvArray(imgrid)%nx, 2, -1

                            outE%cvArray(imgrid)%z(ix, iy, iz) = inE%cvArray(imgrid)%z(ix, iy, iz) - &
                                 ( outE%cvArray(imgrid)%z(ix+1, iy, iz)*zZX(imgrid,ix, 2) &
                                 + outE%cvArray(imgrid)%z(ix, iy+1, iz)*zZY(imgrid,iy, 2))* &
                                 Dilu%cvArray(imgrid)%z(ix, iy, iz)

                         enddo
                      enddo
                      !$OMP END ORDERED
                   enddo
                   !$OMP END DO
              enddo ! Direct loop subgrids
              !$OMP END PARALLEL
            ! adjoint = .true.
            else

              !$OMP PARALLEL DEFAULT(SHARED)
              do imgrid = 1, inE%mgridSize ! Direct loop on sub-grid

                ! ... note that we only parallelize the outer loops

              ! update 1 layer
              call UpdateZ(outE, first, imgrid)
                !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(ix)
                do ix = 1, inE%cvArray(imgrid)%nx
                  !$OMP ORDERED
                  do iz = 2, inE%cvArray(imgrid)%nz+1
                     do iy = 2, inE%cvArray(imgrid)%ny

                            outE%cvArray(imgrid)%x(ix, iy, iz) = inE%cvArray(imgrid)%x(ix, iy, iz) &
                                 - outE%cvArray(imgrid)%x(ix, iy-1, iz)*xXY(imgrid,iy-1, 2) &
                                 * conjg(Dilu%cvArray(imgrid)%x(ix,iy-1,iz))   &
                                 - outE%cvArray(imgrid)%x(ix, iy, iz-1)*xXZ(imgrid,iz-1, 2) &
                                 * conjg(Dilu%cvArray(imgrid)%x(ix, iy, iz-1))
                     enddo
                  enddo
                 !$OMP END ORDERED
               enddo
               !$OMP END DO

               !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iy)
               do iy = 1, inE%cvArray(imgrid)%ny
                !$OMP ORDERED
                  do iz = 2, inE%cvArray(imgrid)%nz+1
                    do ix = 2, inE%cvArray(imgrid)%nx

                          outE%cvArray(imgrid)%y(ix, iy, iz) = inE%cvArray(imgrid)%y(ix, iy, iz) &
                               - outE%cvArray(imgrid)%y(ix, iy, iz-1)*yYZ(imgrid,iz-1, 2) &
                               * conjg(Dilu%cvArray(imgrid)%y(ix,iy,iz-1)) &
                               - outE%cvArray(imgrid)%y(ix-1, iy, iz)*yYX(imgrid,ix-1, 2) &
                               * conjg(Dilu%cvArray(imgrid)%y(ix-1, iy, iz))
                     enddo
                  enddo
                  !$OMP END ORDERED
               enddo
               !$OMP END DO

               !$OMP DO ORDERED SCHEDULE(STATIC) PRIVATE(iz)
               do iz = 1, inE%cvArray(imgrid)%nz
                 !$OMP ORDERED
                 do iy = 2, inE%cvArray(imgrid)%ny
                   do ix = 2, inE%cvArray(imgrid)%nx

                          outE%cvArray(imgrid)%z(ix, iy, iz) = inE%cvArray(imgrid)%z(ix, iy, iz) &
                               - outE%cvArray(imgrid)%z(ix-1, iy, iz)*zZX(imgrid,ix-1, 2) &
                               * conjg(Dilu%cvArray(imgrid)%z(ix-1,iy,iz)) &
                               - outE%cvArray(imgrid)%z(ix, iy-1, iz)*zZY(imgrid,iy-1, 2) &
                               * conjg(Dilu%cvArray(imgrid)%z(ix, iy-1, iz))

                   enddo
                 enddo
                !$OMP END ORDERED
               enddo
               !$OMP END DO
              enddo !  Direct loop subgrids
               !$OMP END PARALLEL

              ! values in the last layer must be zero
              outE%cvArray(outE%mgridSize)%x(:, :, outE%cvarray(outE%mgridSize)%nz+1) = C_ZERO
              outE%cvArray(outE%mgridSize)%y(:, :, outE%cvarray(outE%mgridSize)%nz+1) = C_ZERO

           end if

          else
            print *, 'M2 solver: not compatible usage for existing data types'
          end if

        else
           print *, 'Error complex vectors for M2 solver; multigridSize are not of same size'
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

        integer  :: imgrid

        call create(mGrid, db1, EDGE)
        call create(mGrid, db2, EDGE)
        call create(mGrid, c, CORNER)
        ! d contains the inEerse of diagonal elements of ILU of divCgrad
        call create(mGrid, d, CORNER)
        ! set top nodes to 1.0
        do imgrid=1,mgrid%mgridSize

          d%rscArray(imgrid)%v(1,:,:) = 1.0
          d%rscArray(imgrid)%v(:,1,:) = 1.0
        enddo
        d%rscArray(1)%v(:,:,1) = 1.0

        ! initialize volume weights centered at corners
        Call create(mGrid, volC, CORNER)
        Call CornerVolume(mGrid, volC)

       end subroutine DivCorrInit  ! DivCorrInit Multigrid case

      !**********************************************************************
      ! SetUp routines do calculations (maybe once; possibly more than once)
      ! DivCorrSetup must be called once for each conductivity distribuition
      !  (i.e., before first forward run; after any change to conductivity)

       subroutine DivCorrSetUp()

      implicit none
        integer  :: ix, iy, iz, imgrid
        integer  :: nx, ny, nz, nzAir, nzCum

        if(.not.condE%allocated) then
          print *, 'condE not allocated yet: DivCorrSetUp'
          stop
        endif

        if(.not.db1%allocated) then
          print *,'db1 not allocated yet: DivCorrSetUp'
          stop
        endif

      if(.not.db2%allocated) then
        print *, 'db2 not allocated yet: DivCorrSetUp'
        stop
      endif

      if(.not.c%allocated) then
        print *, 'c not allocated yet: DivCorrSetUp'
        stop
      endif

      if(.not.d%allocated) then
        print *, 'd not allocated yet: DivCorrSetUp'
        stop
      endif

      d%rscArray(1)%v(:,:,1) = 1.0

      do imgrid = 1, mgrid%mgridSize  ! Global loop on subgrids
        nx = mgrid%gridArray(imgrid)%nx
        ny = mgrid%gridArray(imgrid)%ny
        nz = mgrid%gridArray(imgrid)%nz
        nzAir = mgrid%gridArray(imgrid)%nzAir

        ! conductivity of air is modified for computing divergence correction
        ! operator coefficients ...
        if(mgrid%gridArray(imgrid)%flag == 1) then

           condE%rvArray(imgrid)%x(:, :, :) = SIGMA_AIR
           condE%rvArray(imgrid)%y(:, :, :) = SIGMA_AIR
           condE%rvArray(imgrid)%z(:, :, :) = SIGMA_AIR

        else

                 condE%rvArray(imgrid)%x(:, :, 1:nzAir) = SIGMA_AIR
                 condE%rvArray(imgrid)%y(:, :, 1:nzAir) = SIGMA_AIR
                 condE%rvArray(imgrid)%z(:, :, 1:nzAir) = SIGMA_AIR
        endif

        ! update first Z
        call UpdateDivCorr(imgrid)
        ! the coefficients are only for the interior nodes
        ! these coefficients have not been multiplied by volume elements
        ! yet
        do iz = 2, nz
           do iy = 2, ny
              do ix = 2, nx

                 db1%rvArray(imgrid)%x(ix, iy, iz) = condE%rvArray(imgrid)%x(ix-1, iy, iz)/ &
                      (mGrid%gridArray(imgrid)%dx(ix-1)*mGrid%gridArray(imgrid)%delX(ix))
                 db2%rvArray(imgrid)%x(ix, iy, iz) = condE%rvArray(imgrid)%x(ix, iy, iz)/ &
                      (mGrid%gridArray(imgrid)%dx(ix)*mGrid%gridArray(imgrid)%delX(ix))
                 db1%rvArray(imgrid)%y(ix, iy, iz) = condE%rvArray(imgrid)%y(ix, iy-1, iz)/ &
                      (mGrid%gridArray(imgrid)%dy(iy-1)*mGrid%gridArray(imgrid)%delY(iy))
                 db2%rvArray(imgrid)%y(ix, iy, iz) = condE%rvArray(imgrid)%y(ix, iy, iz)/ &
                      (mGrid%gridArray(imgrid)%dy(iy)*mGrid%gridArray(imgrid)%delY(iy))
                 db1%rvArray(imgrid)%z(ix, iy, iz) = condE%rvArray(imgrid)%z(ix, iy, iz-1)/ &
                      (mGrid%gridArray(imgrid)%dz(iz-1)*mGrid%gridarray(imgrid)%delZ(iz))
                 db2%rvArray(imgrid)%z(ix, iy, iz) = condE%rvArray(imgrid)%z(ix, iy, iz)/ &
                      (mGrid%gridArray(imgrid)%dz(iz)*mGrid%gridArray(imgrid)%delZ(iz))
                 c%rscArray(imgrid)%v(ix, iy, iz) = - (db1%rvArray(imgrid)%x(ix, iy, iz) + &
                      db2%rvArray(imgrid)%x(ix, iy, iz) + &
                      db1%rvArray(imgrid)%y(ix, iy, iz) + &
                      db2%rvArray(imgrid)%y(ix, iy, iz) + &
                      db1%rvArray(imgrid)%z(ix, iy, iz) + &
                      db2%rvArray(imgrid)%z(ix, iy, iz))
              enddo
           enddo
        enddo

        ! Multiply by corner volume elements to make operator symmetric
        do iz = 1, nz
           do iy = 2, ny
              do ix = 2, nx

             db1%rvArray(imgrid)%x(ix, iy, iz) = db1%rvArray(imgrid)%x(ix, iy, iz)*volC%rscArray(imgrid)%v(ix, iy, iz)
             db1%rvArray(imgrid)%y(ix, iy, iz) = db1%rvArray(imgrid)%y(ix, iy, iz)*volC%rscArray(imgrid)%v(ix, iy, iz)
             db1%rvArray(imgrid)%z(ix, iy, iz) = db1%rvArray(imgrid)%z(ix, iy, iz)*volC%rscArray(imgrid)%v(ix, iy, iz)
             db2%rvArray(imgrid)%x(ix, iy, iz) = db2%rvArray(imgrid)%x(ix, iy, iz)*volC%rscArray(imgrid)%v(ix, iy, iz)
             db2%rvArray(imgrid)%y(ix, iy, iz) = db2%rvArray(imgrid)%y(ix, iy, iz)*volC%rscArray(imgrid)%v(ix, iy, iz)
             db2%rvArray(imgrid)%z(ix, iy, iz) = db2%rvArray(imgrid)%z(ix, iy, iz)*volC%rscArray(imgrid)%v(ix, iy, iz)

              enddo
           enddo
        enddo

        Call diagMult(c%rscArray(imgrid), volC%rscArray(imgrid), c%rscArray(imgrid)) !diagMult_rscalar

        !  To be explicit about forcing coefficients that multiply boundary
        !    nodes to be zero (this gaurantees that the BC on the potential
        !    is phi = 0):

        db1%rvArray(imgrid)%x(2,:,:) = R_ZERO
        db1%rvArray(imgrid)%y(:,2,:) = R_ZERO
        db1%rvArray(1)%z(:,:,2) = R_ZERO
        db2%rvArray(imgrid)%x(mGrid%gridArray(imgrid)%nx,:,:) = R_ZERO
        db2%rvArray(imgrid)%y(:,mGrid%gridArray(imgrid)%ny,:) = R_ZERO
        db2%rvArray(mgrid%mgridSize)%z(:,:,mGrid%gridArray(mgrid%mgridSize)%nz) = R_ZERO

        ! Compute inverse diagonal elements for D-ILU (interior nodes only)
        ! set top nodes to 1.0
        d%rscArray(imgrid)%v(1,:,:) = 1.0
        d%rscArray(imgrid)%v(:,1,:) = 1.0

        ! update first Z in d array
        call UpdateDivCorr(imgrid,d)

        do iz = 2, nz
          do iy = 2, ny
            do ix = 2, nx

                 d%rscArray(imgrid)%v(ix, iy, iz) = c%rscArray(imgrid)%v(ix, iy, iz) - &
                      db1%rvArray(imgrid)%x(ix,iy,iz)*db2%rvArray(imgrid)%x(ix-1,iy,iz)*d%rscArray(imgrid)%v(ix-1,iy,iz)-&
                      db1%rvArray(imgrid)%y(ix,iy,iz)*db2%rvArray(imgrid)%y(ix,iy-1,iz)*d%rscArray(imgrid)%v(ix,iy-1,iz)-&
                      db1%rvArray(imgrid)%z(ix,iy,iz)*db2%rvArray(imgrid)%z(ix,iy,iz-1)*d%rscArray(imgrid)%v(ix,iy,iz-1)
                 d%rscArray(imgrid)%v(ix, iy, iz) = 1.0/ d%rscArray(imgrid)%v(ix, iy, iz)

            enddo
          enddo
        enddo

        enddo ! Global loop over subgrids

        do imgrid = 1, mgrid%mgridSize
          ! change conductivity of air back to zero
          if(mgrid%gridArray(imgrid)%flag == 1) then
                 condE%rvArray(imgrid)%x(:, :, :) = R_ZERO
                 condE%rvArray(imgrid)%y(:, :, :) = R_ZERO
                 condE%rvArray(imgrid)%z(:, :, :) = R_ZERO
          else
               condE%rvArray(imgrid)%x(:, :, 1:mgrid%gridArray(imgrid)%nzAir) = R_ZERO
               condE%rvArray(imgrid)%y(:, :, 1:mgrid%gridArray(imgrid)%nzAir) = R_ZERO
               condE%rvArray(imgrid)%z(:, :, 1:mgrid%gridArray(imgrid)%nzAir) = R_ZERO
          endif
        enddo

      end subroutine DivCorrSetUp   ! DivCorrSetUp
    ! ********************************************************************************************************************
      subroutine UpdateDivCorr(imgrid,d)
      ! created by Cherevatova (May, 2012)
      ! Computes db1, db2, d, c on the interfaces between sub-grids
      implicit none

        integer, intent(in)  :: imgrid
        type (rscalar_mg) ,optional   :: d
        ! local variables
        integer  :: ix, iy, iz
        integer  :: nz

          if(imgrid == 1) then
          ! in the first sub-grid
          ! z=1 should not be changed
          ! nothing to do
          return
          endif
          nz = mGrid%gridArray(imgrid-1)%nz
          if (mGrid%coarseness(imgrid).gt.mGrid%coarseness(imgrid-1))then
            ! interface layer: finer -> coarser
            ! current sub-grid is coarser
            ! Update 1 taking values from nz of the previous sub-grid
             do iy = 2, mGrid%gridArray(imgrid)%ny
                do ix = 2,  mGrid%gridArray(imgrid)%nx
                  if(present(d))then
                    d%rscArray(imgrid)%v(ix, iy, 1) = c%rscArray(imgrid)%v(ix, iy, 1)- &
                      db1%rvArray(imgrid)%x(ix,iy,1)*db2%rvArray(imgrid)%x(ix-1,iy,1)*d%rscArray(imgrid)%v(ix-1,iy,1)-&
                      db1%rvArray(imgrid)%y(ix,iy,1)*db2%rvArray(imgrid)%y(ix,iy-1,1)*d%rscArray(imgrid)%v(ix,iy-1,1)-&
                      db1%rvArray(imgrid)%z(ix,iy,1)*db2%rvArray(imgrid-1)%z(2*ix-1,2*iy-1,nz)*d%rscArray(imgrid-1)%v(2*ix-1,2*iy-1,nz)

                    d%rscArray(imgrid)%v(ix, iy, 1) = 1.0/ d%rscArray(imgrid)%v(ix, iy, 1)
                  else
                    db1%rvArray(imgrid)%x(ix, iy, 1) = condE%rvArray(imgrid)%x(ix-1,iy,1)/ &
                      (mGrid%gridArray(imgrid)%dx(ix-1)*mGrid%gridArray(imgrid)%delX(ix))
                    db2%rvArray(imgrid)%x(ix, iy, 1) = condE%rvArray(imgrid)%x(ix,iy,1)/&
                      (mGrid%gridArray(imgrid)%dx(ix)*mGrid%gridArray(imgrid)%delX(ix))
                    db1%rvArray(imgrid)%y(ix, iy, 1) = condE%rvArray(imgrid)%y(ix, iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dy(iy-1)*mGrid%gridArray(imgrid)%delY(iy))
                    db2%rvArray(imgrid)%y(ix, iy, 1) = condE%rvArray(imgrid)%y(ix, iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dy(iy)*mGrid%gridArray(imgrid)%delY(iy))
                    db1%rvArray(imgrid)%z(ix, iy, 1) = condE%rvArray(imgrid-1)%z(2*ix-1,2*iy-1,nz)/ &
                      (mGrid%gridArray(imgrid-1)%dz(nz)*mGrid%gridArray(imgrid)%delZ(1))
                    db2%rvArray(imgrid)%z(ix, iy, 1) = condE%rvArray(imgrid)%z(ix, iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dz(1)*mGrid%gridArray(imgrid)%delZ(1))
                    c%rscArray(imgrid)%v(ix, iy, 1) = - (db1%rvArray(imgrid)%x(ix, iy, 1) + &
                      db2%rvArray(imgrid)%x(ix, iy, 1) + &
                      db1%rvArray(imgrid)%y(ix, iy, 1) + &
                      db2%rvArray(imgrid)%y(ix, iy, 1) + &
                      db1%rvArray(imgrid)%z(ix, iy, 1) + &
                      db2%rvArray(imgrid)%z(ix, iy, 1))
                  endif
                enddo
             enddo
          else if (mGrid%coarseness(imgrid).lt.mGrid%coarseness(imgrid-1)) then
           ! interface : coarser to finer
           ! current grid is finer
           do iy = 2, mGrid%gridArray(imgrid-1)%ny
             do ix = 2, mGrid%gridArray(imgrid-1)%nx
               if(present(d))then
                 d%rscArray(imgrid)%v(2*ix-1, 2*iy-1, 1) = c%rscArray(imgrid)%v(2*ix-1, 2*iy-1, 1) - &
                      db1%rvArray(imgrid)%x(2*ix-1, 2*iy-1, 1)*db2%rvArray(imgrid)%x(2*(ix-1), 2*iy-1, 1)*d%rscArray(imgrid)%v(2*(ix-1), 2*iy-1, 1)-&
                      db1%rvArray(imgrid)%y(2*ix-1, 2*iy-1, 1)*db2%rvArray(imgrid)%y(2*ix-1, 2*(iy-1), 1)*d%rscArray(imgrid)%v(2*ix-1, 2*(iy-1), 1)-&
                      db1%rvArray(imgrid)%z(2*ix-1, 2*iy-1, 1)*db2%rvArray(imgrid-1)%z(ix,iy,nz)*d%rscArray(imgrid-1)%v(ix, iy, nz)
                 d%rscArray(imgrid)%v(2*ix-1, 2*iy-1, 1) = 1.0/ d%rscArray(imgrid)%v(2*ix-1, 2*iy-1, 1)

                 d%rscArray(imgrid)%v(2*ix, 2*iy-1, 1) = c%rscArray(imgrid)%v(2*ix, 2*iy-1, 1) - &
                      db1%rvArray(imgrid)%x(2*ix,2*iy-1,1)*db2%rvArray(imgrid)%x(2*ix-1,2*iy-1,1)*d%rscArray(imgrid)%v(2*ix-1,2*iy-1,1)-&
                      db1%rvArray(imgrid)%y(2*ix,2*iy-1,1)*db2%rvArray(imgrid)%y(2*ix,2*(iy-1),1)*d%rscArray(imgrid)%v(2*ix,2*(iy-1),1)-&
                      db1%rvArray(imgrid)%z(2*ix,2*iy-1,1)*db2%rvArray(imgrid-1)%z(ix,iy,nz)*d%rscArray(imgrid-1)%v(ix,iy,nz)
                 d%rscArray(imgrid)%v(2*ix, 2*iy-1, 1) = 1.0/ d%rscArray(imgrid)%v(2*ix, 2*iy-1, 1)

                 d%rscArray(imgrid)%v(2*ix, 2*iy, 1) = c%rscArray(imgrid)%v(2*ix, 2*iy, 1) - &
                      db1%rvArray(imgrid)%x(2*ix,2*iy,1)*db2%rvArray(imgrid)%x(2*ix-1,2*iy,1)*d%rscArray(imgrid)%v(2*ix-1,2*iy,1)-&
                      db1%rvArray(imgrid)%y(2*ix,2*iy,1)*db2%rvArray(imgrid)%y(2*ix,2*iy-1,1)*d%rscArray(imgrid)%v(2*ix,2*iy-1,1)-&
                      db1%rvArray(imgrid)%z(2*ix,2*iy,1)*db2%rvArray(imgrid-1)%z(ix,iy,nz) *d%rscArray(imgrid-1)%v(ix,iy,nz)
                 d%rscArray(imgrid)%v(2*ix, 2*iy, 1) = 1.0/ d%rscArray(imgrid)%v(2*ix, 2*iy, 1)

                 d%rscArray(imgrid)%v(2*ix-1, 2*iy, 1) = c%rscArray(imgrid)%v(2*ix-1, 2*iy, 1) - &
                      db1%rvArray(imgrid)%x(2*ix-1,2*iy,1)*db2%rvArray(imgrid)%x(2*(ix-1),2*iy,1)*d%rscArray(imgrid)%v(2*(ix-1),2*iy,1)-&
                      db1%rvArray(imgrid)%y(2*ix-1,2*iy,1)*db2%rvArray(imgrid)%y(2*ix-1,2*iy-1,1)*d%rscArray(imgrid)%v(2*ix-1,2*iy-1,1)-&
                      db1%rvArray(imgrid)%z(2*ix-1,2*iy,1)*db2%rvArray(imgrid-1)%z(ix,iy,nz)*d%rscArray(imgrid-1)%v(ix,iy,nz)
                 d%rscArray(imgrid)%v(2*ix-1, 2*iy, 1) = 1.0/ d%rscArray(imgrid)%v(2*ix-1, 2*iy, 1)
               else
                 ! X
                 !odd/odd
                  db1%rvArray(imgrid)%x(2*ix-1, 2*iy-1, 1) = condE%rvArray(imgrid)%x(2*(ix-1), 2*iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dx(2*(ix-1))*mGrid%gridArray(imgrid)%delX(2*ix-1))
                 ! even/odd
                  db1%rvArray(imgrid)%x(2*ix, 2*iy-1, 1) = condE%rvArray(imgrid)%x(2*ix-1, 2*iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dx(2*ix-1)*mGrid%gridArray(imgrid)%delX(2*ix))
                 ! even/even
                  db1%rvArray(imgrid)%x(2*ix, 2*iy, 1) = condE%rvArray(imgrid)%x(2*ix-1, 2*iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dx(2*ix-1)*mGrid%gridArray(imgrid)%delX(2*ix))
                 ! odd/even
                  db1%rvArray(imgrid)%x(2*ix-1, 2*iy, 1) = condE%rvArray(imgrid)%x(2*(ix-1), 2*iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dx(2*(ix-1))*mGrid%gridArray(imgrid)%delX(2*ix))

                  db2%rvArray(imgrid)%x(2*ix-1, 2*iy-1, 1) = condE%rvArray(imgrid)%x(2*ix-1, 2*iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dx(2*ix-1)*mGrid%gridArray(imgrid)%delX(2*ix-1))
                  db2%rvArray(imgrid)%x(2*ix, 2*iy-1, 1) = condE%rvArray(imgrid)%x(2*ix, 2*iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dx(2*ix)*mGrid%gridArray(imgrid)%delX(2*ix))
                  db2%rvArray(imgrid)%x(2*ix, 2*iy, 1) = condE%rvArray(imgrid)%x(2*ix, 2*iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dx(2*ix)*mGrid%gridArray(imgrid)%delX(2*ix))
                  db2%rvArray(imgrid)%x(2*ix-1, 2*iy, 1) = condE%rvArray(imgrid)%x(2*ix-1, 2*iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dx(2*ix-1)*mGrid%gridArray(imgrid)%delX(2*ix-1))
                 ! Y
                  db1%rvArray(imgrid)%y(2*ix-1, 2*iy-1, 1)= condE%rvArray(imgrid)%y(2*ix-1, 2*(iy-1), 1)/ &
                      (mGrid%gridArray(imgrid)%dy(2*(iy-1))*mGrid%gridArray(imgrid)%delY(2*iy-1))
                  db1%rvArray(imgrid)%y(2*ix, 2*iy-1, 1) = condE%rvArray(imgrid)%y(2*ix, 2*(iy-1), 1)/ &
                      (mGrid%gridArray(imgrid)%dy(2*(iy-1))*mGrid%gridArray(imgrid)%delY(2*iy-1))
                  db1%rvArray(imgrid)%y(2*ix, 2*iy, 1) = condE%rvArray(imgrid)%y(2*ix, 2*iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dy(2*iy-1)*mGrid%gridArray(imgrid)%delY(2*iy))
                  db1%rvArray(imgrid)%y(2*ix-1, 2*iy, 1) = condE%rvArray(imgrid)%y(2*ix-1, 2*iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dy(2*iy-1)*mGrid%gridArray(imgrid)%delY(2*iy))

                  db2%rvArray(imgrid)%y(2*ix-1, 2*iy-1, 1) = condE%rvArray(imgrid)%y(2*ix-1, 2*iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dy(2*iy-1)*mGrid%gridArray(imgrid)%delY(2*iy-1))
                  db2%rvArray(imgrid)%y(2*ix, 2*iy-1, 1) = condE%rvArray(imgrid)%y(2*ix, 2*iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dy(2*iy-1)*mGrid%gridArray(imgrid)%delY(2*iy-1))
                  db2%rvArray(imgrid)%y(2*ix, 2*iy, 1) = condE%rvArray(imgrid)%y(2*ix, 2*iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dy(2*iy)*mGrid%gridArray(imgrid)%delY(2*iy))
                  db2%rvArray(imgrid)%y(2*ix-1, 2*iy, 1) = condE%rvArray(imgrid)%y(2*ix-1, 2*iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dy(2*iy)*mGrid%gridArray(imgrid)%delY(2*iy))
                 ! Z
                  db1%rvArray(imgrid)%z(2*ix-1, 2*iy-1, 1) = condE%rvArray(imgrid-1)%z(ix, iy, nz)/ &
                      (mGrid%gridArray(imgrid-1)%dz(nz)*mGrid%gridarray(imgrid)%delZ(1))
                  db1%rvArray(imgrid)%z(2*ix, 2*iy-1, 1) = condE%rvArray(imgrid-1)%z(ix, iy, nz)/ &
                      (mGrid%gridArray(imgrid-1)%dz(nz)*mGrid%gridarray(imgrid)%delZ(1))
                  db1%rvArray(imgrid)%z(2*ix, 2*iy, 1) = condE%rvArray(imgrid-1)%z(ix, iy, nz)/ &
                      (mGrid%gridArray(imgrid-1)%dz(nz)*mGrid%gridarray(imgrid)%delZ(1))
                  db1%rvArray(imgrid)%z(2*ix-1, 2*iy, 1) = condE%rvArray(imgrid-1)%z(ix, iy, nz)/ &
                      (mGrid%gridArray(imgrid-1)%dz(nz)*mGrid%gridarray(imgrid)%delZ(1))

                  db2%rvArray(imgrid)%z(2*ix-1, 2*iy-1, 1) = condE%rvArray(imgrid)%z(2*ix-1, 2*iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dz(1)*mGrid%gridArray(imgrid)%delZ(1))
                  db2%rvArray(imgrid)%z(2*ix, 2*iy-1, 1) = condE%rvArray(imgrid)%z(2*ix, 2*iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dz(1)*mGrid%gridArray(imgrid)%delZ(1))
                  db2%rvArray(imgrid)%z(2*ix, 2*iy, 1) = condE%rvArray(imgrid)%z(2*ix, 2*iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dz(1)*mGrid%gridArray(imgrid)%delZ(1))
                  db2%rvArray(imgrid)%z(2*ix-1, 2*iy, 1) = condE%rvArray(imgrid)%z(2*ix-1, 2*iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dz(1)*mGrid%gridArray(imgrid)%delZ(1))
                  ! C
                  c%rscArray(imgrid)%v(2*ix-1, 2*iy-1, 1) = - (db1%rvArray(imgrid)%x(2*ix-1, 2*iy-1, 1) + &
                      db2%rvArray(imgrid)%x(2*ix-1, 2*iy-1, 1) + &
                      db1%rvArray(imgrid)%y(2*ix-1, 2*iy-1, 1) + &
                      db2%rvArray(imgrid)%y(2*ix-1, 2*iy-1, 1) + &
                      db1%rvArray(imgrid)%z(2*ix-1, 2*iy-1, 1) + &
                      db2%rvArray(imgrid)%z(2*ix-1, 2*iy-1, 1))
                  c%rscArray(imgrid)%v(2*ix, 2*iy-1, 1) = - (db1%rvArray(imgrid)%x(2*ix, 2*iy-1, 1) + &
                      db2%rvArray(imgrid)%x(2*ix, 2*iy-1, 1) + &
                      db1%rvArray(imgrid)%y(2*ix, 2*iy-1, 1) + &
                      db2%rvArray(imgrid)%y(2*ix, 2*iy-1, 1) + &
                      db1%rvArray(imgrid)%z(2*ix, 2*iy-1, 1) + &
                      db2%rvArray(imgrid)%z(2*ix, 2*iy-1, 1))
                  c%rscArray(imgrid)%v(2*ix, 2*iy, 1) = - (db1%rvArray(imgrid)%x(2*ix, 2*iy, 1) + &
                      db2%rvArray(imgrid)%x(2*ix, 2*iy, 1) + &
                      db1%rvArray(imgrid)%y(2*ix, 2*iy, 1) + &
                      db2%rvArray(imgrid)%y(2*ix, 2*iy, 1) + &
                      db1%rvArray(imgrid)%z(2*ix, 2*iy, 1) + &
                      db2%rvArray(imgrid)%z(2*ix, 2*iy, 1))
                  c%rscArray(imgrid)%v(2*ix-1, 2*iy, 1) = - (db1%rvArray(imgrid)%x(2*ix-1, 2*iy, 1) + &
                      db2%rvArray(imgrid)%x(2*ix-1, 2*iy, 1) + &
                      db1%rvArray(imgrid)%y(2*ix-1, 2*iy, 1) + &
                      db2%rvArray(imgrid)%y(2*ix-1, 2*iy, 1) + &
                      db1%rvArray(imgrid)%z(2*ix-1, 2*iy, 1) + &
                      db2%rvArray(imgrid)%z(2*ix-1, 2*iy, 1))
               endif
             enddo !ix
           enddo !iy

          else if (mGrid%coarseness(imgrid).eq.mGrid%coarseness(imgrid-1)) then
           do iy = 2, mGrid%gridArray(imgrid)%ny
             do ix = 2, mGrid%gridArray(imgrid)%nx
               if(present(d))then
               d%rscArray(imgrid)%v(ix, iy, 1) = c%rscArray(imgrid)%v(ix, iy, 1) - &
                      db1%rvArray(imgrid)%x(ix,iy,1)*db2%rvArray(imgrid)%x(ix-1,iy,1)*d%rscArray(imgrid)%v(ix-1,iy,1)-&
                      db1%rvArray(imgrid)%y(ix,iy,1)*db2%rvArray(imgrid)%y(ix,iy-1,1)*d%rscArray(imgrid)%v(ix,iy-1,1)-&
                      db1%rvArray(imgrid)%z(ix,iy,1)*db2%rvArray(imgrid-1)%z(ix,iy,nz)*d%rscArray(imgrid-1)%v(ix,iy,nz)
               d%rscArray(imgrid)%v(ix, iy, 1) = 1.0/ d%rscArray(imgrid)%v(ix, iy, 1)

               else
                db1%rvArray(imgrid)%x(ix, iy, 1) = condE%rvArray(imgrid)%x(ix-1, iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dx(ix-1)*mGrid%gridArray(imgrid)%delX(ix))
                 db2%rvArray(imgrid)%x(ix, iy, 1) = condE%rvArray(imgrid)%x(ix, iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dx(ix)*mGrid%gridArray(imgrid)%delX(ix))
                 db1%rvArray(imgrid)%y(ix, iy, 1) = condE%rvArray(imgrid)%y(ix, iy-1, 1)/ &
                      (mGrid%gridArray(imgrid)%dy(iy-1)*mGrid%gridArray(imgrid)%delY(iy))
                 db2%rvArray(imgrid)%y(ix, iy, 1) = condE%rvArray(imgrid)%y(ix, iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dy(iy)*mGrid%gridArray(imgrid)%delY(iy))
                 db1%rvArray(imgrid)%z(ix, iy, 1) = condE%rvArray(imgrid-1)%z(ix, iy, nz)/ &
                      (mGrid%gridArray(imgrid-1)%dz(nz)*mGrid%gridArray(imgrid)%delZ(1))
                 db2%rvArray(imgrid)%z(ix, iy, 1) = condE%rvArray(imgrid)%z(ix, iy, 1)/ &
                      (mGrid%gridArray(imgrid)%dz(1)*mGrid%gridArray(imgrid)%delZ(1))
                 c%rscArray(imgrid)%v(ix, iy, 1) = - (db1%rvArray(imgrid)%x(ix, iy, 1) + &
                      db2%rvArray(imgrid)%x(ix, iy, 1) + &
                      db1%rvArray(imgrid)%y(ix, iy, 1) + &
                      db2%rvArray(imgrid)%y(ix, iy, 1) + &
                      db1%rvArray(imgrid)%z(ix, iy, 1) + &
                      db2%rvArray(imgrid)%z(ix, iy, 1))
               endif
             enddo
           enddo
         endif

      end subroutine UpdateDivCorr
      !**********************************************************************
      ! to deallocate the coefficients used for divergence correction
      subroutine Deallocate_DivCorr()

        implicit none

        Call deall(db1) !deall_rvector_mg
        Call deall(db2)
        Call deall(c) !deall_rscalar_mg
        Call deall(d)
        ! corner volumes
        Call deall(volC)

      end subroutine Deallocate_DivCorr	! Deallocate_DivCorr

      !**********************************************************************
      ! apply pre-conditioner, solving lower and upper triangular systems using
      ! coefficients in db1, db2, and d.

      subroutine DivCgradILU(inPhi, outPhi)

      implicit none

        type (cscalar_mg), intent(in)   :: inPhi
        type (cscalar_mg), intent(inout)  :: outPhi
        ! local variables
        integer  :: ix, iy, iz, imgrid
        integer  :: nx, ny ,nz
        character(len=7)  :: reverse='reverse'

        if(.not.inPhi%allocated) then
          print *, 'inPhi not allocated in DivCgradILU'
        stop
        endif

        if(.not.outPhi%allocated) then
         print *, 'outPhi not allocated in DivCgradILU'
         stop
        endif

        if (outPhi%allocated) then
          if ((inPhi%gridType == outPhi%gridType).and.(inPhi%mgridSize == outPhi%mgridSize)) then

            do imgrid= 1, mGrid%mgridSize ! Direct loop on multigrid
              nx = inPhi%csArray(imgrid)%nx
              ny = inPhi%csArray(imgrid)%ny
              nz = inPhi%csArray(imgrid)%nz

              ! Check whether all the inputs/ outputs involved are even of the same size
              if ((nx == outPhi%csArray(imgrid)%nx).and.&
                (ny == outPhi%csArray(imgrid)%ny).and.&
                (nz == outPhi%csArray(imgrid)%nz)) then

              ! update first layer
              call UpdateDivCgradILU(inPhi, outPhi, imgrid)

                do iz = 2, nz
                  do iy = 2, ny
                    do ix = 2, nx

                    outPhi%csArray(imgrid)%v(ix, iy, iz) = inPhi%csArray(imgrid)%v(ix, iy, iz) &
                      - outPhi%csArray(imgrid)%v(ix-1,iy,iz)*db1%rvArray(imgrid)%x(ix,iy,iz)&
                        *d%rscArray(imgrid)%v(ix-1,iy,iz) &
                       - outPhi%csArray(imgrid)%v(ix,iy-1,iz)*db1%rvArray(imgrid)%y(ix,iy,iz)&
                         *d%rscArray(imgrid)%v(ix,iy-1,iz) &
                        - outPhi%csArray(imgrid)%v(ix,iy,iz-1)*db1%rvArray(imgrid)%z(ix,iy,iz)&
                          *d%rscArray(imgrid)%v(ix,iy,iz-1)
                    enddo
                  enddo
               enddo
             else
               print *, 'DivCgradILU: not compatible existing data types'
             endif
          enddo !Loop on subgrids

          ! backward substitution (Solve upper triangular system)
          ! the coefficients are only for the interior nodes
          do imgrid = mGrid%mgridSize, 1, -1 ! Reverse loop on subgrids

             nx = inPhi%csArray(imgrid)%nx
             ny = inPhi%csArray(imgrid)%ny
             nz = inPhi%csArray(imgrid)%nz

             ! update nz+1 layer
             call UpdateDivCgradILU(inPhi, outPhi, imgrid, reverse)

             do iz = nz,1,-1
               do iy = ny,2,-1
                  do ix = nx,2,-1

                    outPhi%csArray(imgrid)%v(ix, iy, iz) = (outPhi%csArray(imgrid)%v(ix, iy, iz) &
                      - outPhi%csArray(imgrid)%v(ix+1, iy, iz)*db2%rvArray(imgrid)%x(ix, iy, iz)  &
                      - outPhi%csArray(imgrid)%v(ix, iy+1, iz)*db2%rvArray(imgrid)%y(ix, iy, iz)  &
                      - outPhi%csArray(imgrid)%v(ix, iy, iz+1)*db2%rvArray(imgrid)%z(ix, iy, iz)) &
                            *d%rscArray(imgrid)%v(ix, iy, iz)
                       enddo
                    enddo
                 enddo
          enddo ! GLobal loop on subgrids (reverse)

        else
          print *, 'Error: DivCgradILU: scalars not same size'
        endif

      else
        print *, 'Error: DivCgradILU: output scalar not even allocated yet'
      endif

      end subroutine DivCgradILU  		! DivCgradILU
      ! ****************************************************************************************************
      subroutine UpdateDivCgradILU (inSc, outSc, imgrid,reverse)
      ! created by Cherevatova (May, 2012)
      ! Computes DivCgradILU on the interfaces between sub-grids

        implicit none

        type (cscalar_mg), intent(in)   :: inSc
        type (cscalar_mg), intent(inout)  :: outSc
        integer, intent(in)  :: imgrid
        ! local variables
        integer  :: ix, iy, iz
        integer  :: nz
        character(len=7),optional  :: reverse

         if (present(reverse)) then
            if(imgrid==mgrid%mgridSize)then
             ! do not compute in the last sub-grid
             !nothing to do
              return
            endif
            nz = mgrid%gridArray(imgrid)%nz
            if (mGrid%coarseness(imgrid).lt.mGrid%coarseness(imgrid+1))then
            ! coarser --> finer
              do iy = 2, mGrid%gridArray(imgrid+1)%ny
                do ix = 2, mGrid%gridArray(imgrid+1)%nx
                  outSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, nz+1) = outSc%csArray(imgrid+1)%v(ix, iy, 1)
                  outSc%csArray(imgrid)%v(2*ix, 2*iy-1, nz+1) =  outSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, nz+1)
                  outSc%csArray(imgrid)%v(2*ix, 2*iy, nz+1) = outSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, nz+1)
                  outSc%csArray(imgrid)%v(2*ix-1, 2*iy, nz+1) = outSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, nz+1)
                enddo !ix
              enddo !iy
            else if  (mGrid%coarseness(imgrid).gt.mGrid%coarseness(imgrid+1)) then
            ! finer --> coarser
              do iy = 2, mGrid%gridArray(imgrid)%ny
                do ix = 2, mGrid%gridArray(imgrid)%nx
                  outSc%csArray(imgrid)%v(ix, iy, nz+1) = (outSc%csArray(imgrid+1)%v(2*ix-1, 2*iy-1, 1)+&
                                                           outSc%csArray(imgrid+1)%v(2*ix, 2*iy-1, 1))/2
                enddo
              enddo
            else if (mGrid%coarseness(imgrid).eq.mGrid%coarseness(imgrid+1)) then
              do iy = 2, mGrid%gridArray(imgrid)%ny
                do ix = 2, mGrid%gridArray(imgrid)%nx
                  outSc%csArray(imgrid)%v(ix, iy, nz+1) =  outSc%csArray(imgrid+1)%v(ix, iy, 1)
                enddo
              enddo
            endif
          else ! if direct loop

            if(imgrid == 1) then
              ! again, do not compute in the first sub-grid
              ! nothing to do
              return
            endif
            nz = mGrid%gridArray(imgrid-1)%nz
            if (mGrid%coarseness(imgrid).gt.mGrid%coarseness(imgrid-1))then
              ! interface layer: finer -> coarser
              ! current sub-grid is coarser
              ! Update 1 taking valies from nz of the previous sub-grid
              do iy = 2, mGrid%gridArray(imgrid)%ny
                do ix = 2, mGrid%gridArray(imgrid)%nx
                  outSc%csArray(imgrid)%v(ix, iy, 1) = inSc%csArray(imgrid)%v(ix, iy, 1) &
                  - outSc%csArray(imgrid)%v(ix-1,iy,1)*db1%rvArray(imgrid)%x(ix,iy,1)&
                                                       *d%rscArray(imgrid)%v(ix-1,iy,1) &
                       - outSc%csArray(imgrid)%v(ix,iy-1,1)*db1%rvArray(imgrid)%y(ix,iy,1)&
                         *d%rscArray(imgrid)%v(ix,iy-1,1) &
                        - outSc%csArray(imgrid-1)%v(2*ix-1,2*iy-1,nz)*db1%rvArray(imgrid)%z(ix,iy,1)&
                          *d%rscArray(imgrid-1)%v(2*ix-1,2*iy-1,nz)
                enddo !ix
              enddo !iy
            else if (mGrid%coarseness(imgrid).lt.mGrid%coarseness(imgrid-1)) then
              ! interface : coarser to finer
              do iy = 2, mGrid%gridArray(imgrid-1)%ny
                do ix = 2, mGrid%gridArray(imgrid-1)%nx
                  ! odd/odd nodes
                  outSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, 1) = inSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, 1) &
                      - outSc%csArray(imgrid)%v(2*(ix-1), 2*iy-1, 1)*db1%rvArray(imgrid)%x(2*ix-1, 2*iy-1, 1)&
                        *d%rscArray(imgrid)%v(2*(ix-1), 2*iy-1, 1) &
                       - outSc%csArray(imgrid)%v(2*ix-1, 2*(iy-1), 1)*db1%rvArray(imgrid)%y(2*ix-1, 2*iy-1, 1)&
                         *d%rscArray(imgrid)%v(2*ix-1, 2*(iy-1), 1) &
                        - outSc%csArray(imgrid-1)%v(ix,iy,nz)*db1%rvArray(imgrid)%z(2*ix-1, 2*iy-1, 1)&
                          *d%rscArray(imgrid-1)%v(ix,iy,nz)
    ! simple copy
                  outSc%csArray(imgrid)%v(2*ix, 2*iy-1, 1) = outSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, 1)
                  outSc%csArray(imgrid)%v(2*ix-1, 2*iy, 1) = outSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, 1)
                  outSc%csArray(imgrid)%v(2*ix, 2*iy, 1) = outSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, 1)

                enddo
              enddo

            else if (mGrid%coarseness(imgrid).eq.mGrid%coarseness(imgrid-1)) then
              do iy = 2, mGrid%gridArray(imgrid)%ny
                do ix = 2, mGrid%gridArray(imgrid)%nx
                 outSc%csArray(imgrid)%v(ix, iy, 1) = inSc%csArray(imgrid)%v(ix, iy, 1) &
                      - outSc%csArray(imgrid)%v(ix-1,iy,1)*db1%rvArray(imgrid)%x(ix,iy,1)&
                        *d%rscArray(imgrid)%v(ix-1,iy,1) &
                       - outSc%csArray(imgrid)%v(ix,iy-1,1)*db1%rvArray(imgrid)%y(ix,iy,1)&
                         *d%rscArray(imgrid)%v(ix,iy-1,1) &
                        - outSc%csArray(imgrid-1)%v(ix,iy,nz)*db1%rvArray(imgrid)%z(ix,iy,1)&
                          *d%rscArray(imgrid-1)%v(ix,iy,nz)
                enddo
              enddo
            endif
          endif
      end subroutine UpdateDivCgradILU
      !**********************************************************************
      ! Apply operator div sigma grad to a scalar field (used for corners)
      !  called by PCG for iterative solution of divergence correction equation
      subroutine DivCgrad(inPhi, outPhi)

      implicit none
        type (cscalar_mg), intent(in)  :: inPhi
        type (cscalar_mg), intent(inout)  :: outPhi
        integer  :: ix, iy, iz,imgrid

        if(.not.inPhi%allocated) then
          print *, 'inPhi not allocated in DivCgrad'
          stop
        endif

        if(.not.outPhi%allocated) then
         print *, 'outPhi not allocated in DivCgrad'
         stop
        endif

        if (outPhi%allocated .and. inPhi%mgridSize == outPhi%mgridSize) then

          if (inPhi%gridType == outPhi%gridType) then

           do imgrid = 1, inPhi%mgridSize  ! loop on subgrids

            ! Check whether all the inputs/ outputs involved are even of the same size
            if ((inPhi%csArray(imgrid)%nx == outPhi%csArray(imgrid)%nx).and.&
                (inPhi%csArray(imgrid)%ny == outPhi%csArray(imgrid)%ny).and.&
                (inPhi%csArray(imgrid)%nz == outPhi%csArray(imgrid)%nz)) then

           ! update first layer
           call UpdateDivCgrad(inPhi, outPhi, imgrid)

              ! the coefficients are only for interior nodes
              !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz)
              do iz = 2, inPhi%csArray(imgrid)%nz
                do iy = 2, inPhi%csArray(imgrid)%ny
                  do ix = 2, inPhi%csArray(imgrid)%nx

                    outPhi%csArray(imgrid)%v(ix,iy,iz) = inPhi%csArray(imgrid)%v(ix+1,iy,iz)&
                                                      *db2%rvArray(imgrid)%x(ix,iy,iz)+&
                               inPhi%csArray(imgrid)%v(ix-1, iy, iz)*db1%rvArray(imgrid)%x(ix, iy, iz) + &
                               inPhi%csArray(imgrid)%v(ix, iy+1, iz)*db2%rvArray(imgrid)%y(ix, iy, iz) + &
                               inPhi%csArray(imgrid)%v(ix, iy-1, iz)*db1%rvArray(imgrid)%y(ix, iy, iz) + &
                               inPhi%csArray(imgrid)%v(ix, iy, iz+1)*db2%rvArray(imgrid)%z(ix, iy, iz) + &
                               inPhi%csArray(imgrid)%v(ix, iy, iz-1)*db1%rvArray(imgrid)%z(ix, iy, iz) + &
                               inPhi%csArray(imgrid)%v(ix, iy, iz)*c%rscArray(imgrid)%v(ix, iy, iz)

                  enddo
                enddo
              enddo
            !$OMP END PARALLEL DO

            else
              print *, 'Error: DivCgrad: scalars not same size'
            end if

        enddo ! loop on subgrids

        else
          print *, 'DivCgrad: not compatible existing data types'
        endif

      else
        print *, 'Error: DivCgrad: output scalar not even allocated yet'
      endif

      end subroutine DivCgrad	! DivCgrad

      ! *********************************************************************
      !Update first layer in DivCgrad
      subroutine UpdateDivCgrad(inSc,outSc,imgrid)
        implicit none

        type (cscalar_mg), intent(in)   :: inSc
        type (cscalar_mg), intent(inout)  :: outSc
        integer, intent(in)  :: imgrid
        ! local variables
        integer  :: ix, iy, iz
        integer  :: nz


          if(imgrid == 1) then
          return
          endif

          ! the basic strategy is to fill in the interface layer with values from the finer grid!
          if (mGrid%coarseness(imgrid).gt.mGrid%coarseness(imgrid-1))then
            ! interface layer: finer to coarser
            ! current sub-grid is coarser
            ! Update 1 taking values from nz of the previous sub-grid
           do iy = 2, mGrid%gridArray(imgrid)%ny
              do ix = 2, mGrid%gridArray(imgrid)%nx
                outSc%csArray(imgrid)%v(ix,iy,1) = inSc%csArray(imgrid)%v(ix+1,iy,1)&
                                                     *db2%rvArray(imgrid)%x(ix,iy,1)+&
                               inSc%csArray(imgrid)%v(ix-1, iy, 1)*db1%rvArray(imgrid)%x(ix, iy, 1) + &
                               inSc%csArray(imgrid)%v(ix, iy+1, 1)*db2%rvArray(imgrid)%y(ix, iy, 1) + &
                               inSc%csArray(imgrid)%v(ix, iy-1, 1)*db1%rvArray(imgrid)%y(ix, iy, 1) + &
                               inSc%csArray(imgrid)%v(ix, iy, 2)*db2%rvArray(imgrid)%z(ix, iy, 1) + &
                               inSc%csArray(imgrid-1)%v(2*ix-1, 2*iy-1,mGrid%gridArray(imgrid-1)%nz )*db1%rvArray(imgrid)%z(ix, iy, 1) + &
                               inSc%csArray(imgrid)%v(ix, iy, 1)*c%rscArray(imgrid)%v(ix, iy, 1)
               enddo
           enddo
         else if (mGrid%coarseness(imgrid).lt.mGrid%coarseness(imgrid-1)) then
           ! interface : coarser to finer
           ! current grid is finer
           ! vice versa
           nz = mGrid%gridArray(imgrid-1)%nz
           do iy = 2, mGrid%gridArray(imgrid-1)%ny
             do ix = 2, mGrid%gridArray(imgrid-1)%nx
                outSc%csArray(imgrid)%v(2*ix-1,2*iy-1,1) = inSc%csArray(imgrid)%v(2*ix,2*iy-1,1)&
                                                      *db2%rvArray(imgrid)%x(2*ix-1,2*iy-1,1)+&
                               inSc%csArray(imgrid)%v(2*(ix-1), 2*iy-1, 1)*db1%rvArray(imgrid)%x(2*ix-1, 2*iy-1, 1) + &
                               inSc%csArray(imgrid)%v(2*ix-1, 2*iy, 1)*db2%rvArray(imgrid)%y(2*ix-1, 2*iy-1, 1) + &
                               inSc%csArray(imgrid)%v(2*ix-1, 2*(iy-1), 1)*db1%rvArray(imgrid)%y(2*ix-1, 2*iy-1, 1) + &
                               inSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, 2)*db2%rvArray(imgrid)%z(2*ix-1, 2*iy-1, 1) + &
                               inSc%csArray(imgrid-1)%v(ix, iy, nz)*db1%rvArray(imgrid)%z(2*ix-1, 2*iy-1, 1) + &
                               inSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, 1)*c%rscArray(imgrid)%v(2*ix-1, 2*iy-1, 1)
              enddo
            enddo

          else if (mGrid%coarseness(imgrid).eq.mGrid%coarseness(imgrid-1)) then
           do iy = 2, mGrid%gridArray(imgrid)%ny
             do ix = 2, mGrid%gridArray(imgrid)%nx
                    outSc%csArray(imgrid)%v(ix,iy,1) = inSc%csArray(imgrid)%v(ix+1,iy,1)&
                                                      *db2%rvArray(imgrid)%x(ix,iy,1)+&
                               inSc%csArray(imgrid)%v(ix-1, iy, 1)*db1%rvArray(imgrid)%x(ix, iy, 1) + &
                               inSc%csArray(imgrid)%v(ix, iy+1, 1)*db2%rvArray(imgrid)%y(ix, iy, 1) + &
                               inSc%csArray(imgrid)%v(ix, iy-1, 1)*db1%rvArray(imgrid)%y(ix, iy, 1) + &
                               inSc%csArray(imgrid)%v(ix, iy, 2)*db2%rvArray(imgrid)%z(ix, iy, 1) + &
                               inSc%csArray(imgrid-1)%v(ix, iy, mGrid%gridArray(imgrid-1)%nz)*db1%rvArray(imgrid)%z(ix, iy, 1) + &
                               inSc%csArray(imgrid)%v(ix, iy, 1)*c%rscArray(imgrid)%v(ix, iy, 1)
             enddo
           enddo
         endif
      end subroutine UpdateDivCgrad
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
          type (cvector_mg), intent(in)  :: inE
          type (cscalar_mg), intent(inout)   :: outSc
          integer   :: ix, iy, iz,imgrid
          integer  :: nz, nzAir, nx, ny

          if(.not.inE%allocated) then
            print *, 'inE not allocated in DivC'
            stop
          endif

          if(.not.outSc%allocated) then
            print *, 'outSc not allocated in DivC'
            stop
          endif

          if (outSc%gridType == CORNER) then

          do imgrid = 1, inE%mgridSize ! Loop over subgrids
            nx = inE%cvArray(imgrid)%nx
            ny = inE%cvArray(imgrid)%ny
            nz = inE%cvArray(imgrid)%nz
            nzAir = inE%cvArray(imgrid)%grid%nzAir

            ! Check whether all the inputs/ outputs involved are even of the same
            ! size
            if ((nx == outSc%csArray(imgrid)%nx).and.&
                (ny == outSc%csArray(imgrid)%ny).and.&
                (nz == outSc%csArray(imgrid)%nz)) then
              ! computation done only for internal nodes

           ! update first layer
             call UpdateDivC(inE,outSc,imgrid) ! call UpdateZ_DivC subroutine

             do ix = 2, nx
               do iy = 2, ny

                ! FOR NODES IN THE AIR ONLY
                   do iz = 2,nzAir

                   outSc%csArray(imgrid)%v(ix, iy, iz) = &
                      SIGMA_AIR*(inE%cvArray(imgrid)%x(ix,iy,iz)-inE%cvArray(imgrid)%x(ix-1,iy,iz)) * &
                        inE%cvArray(imgrid)%grid%delXinv(ix)    &
                        + SIGMA_AIR*(inE%cvArray(imgrid)%y(ix,iy,iz)-inE%cvArray(imgrid)%y(ix,iy-1,iz)) * &
                        inE%cvArray(imgrid)%grid%delYinv(iy)    &
                        + SIGMA_AIR*(inE%cvArray(imgrid)%z(ix,iy,iz)-inE%cvArray(imgrid)%z(ix,iy,iz-1)) * &
                       inE%cvArray(imgrid)%grid%delZinv(iz)
                   enddo   ! iz

                   ! FOR NODES AT THE AIR-EARTH INTERFACE

                   if(inE%cvArray(imgrid)%grid%flag == 0) then
                     iz = nzAir+1

                       outSc%csArray(imgrid)%v(ix, iy, iz) = &
                            (condE%rvArray(imgrid)%x(ix,iy,iz)*inE%cvArray(imgrid)%x(ix, iy, iz) -         &
                            condE%rvArray(imgrid)%x(ix-1,iy,iz)*inE%cvArray(imgrid)%x(ix-1, iy, iz)) * &
                            inE%cvArray(imgrid)%grid%delXinv(ix)      &
                            +  (condE%rvArray(imgrid)%y(ix,iy,iz)*inE%cvArray(imgrid)%y(ix, iy, iz) -      &
                            condE%rvArray(imgrid)%y(ix,iy - 1,iz)*inE%cvArray(imgrid)%y(ix, iy - 1, iz)) * &
                            inE%cvArray(imgrid)%grid%delYinv(iy)      &
                            +  (condE%rvArray(imgrid)%z(ix,iy,iz)*inE%cvArray(imgrid)%z(ix, iy, iz) -      &
                            SIGMA_AIR*inE%cvArray(imgrid)%z(ix, iy, iz-1)) * &
                            inE%cvArray(imgrid)%grid%delZinv(iz)
                   endif
                    ! FOR NODES INSIDE THE EARTH ONLY
                    ! THE TOP MOST EARTH NODE HAS AN INTERFACE WITH
                    ! AIR, THEREFORE THAT ONE IS SKIPPED HERE
                    do iz = nzAir+2, nz
                       outSc%csArray(imgrid)%v(ix, iy, iz) = &
                            (condE%rvArray(imgrid)%x(ix,iy,iz)*inE%cvArray(imgrid)%x(ix, iy, iz) -         &
                            condE%rvArray(imgrid)%x(ix - 1,iy,iz)*inE%cvArray(imgrid)%x(ix - 1, iy, iz)) * &
                            inE%cvArray(imgrid)%grid%delXinv(ix)      &
                            +  (condE%rvArray(imgrid)%y(ix,iy,iz)*inE%cvArray(imgrid)%y(ix, iy, iz) -      &
                            condE%rvArray(imgrid)%y(ix,iy - 1,iz)*inE%cvArray(imgrid)%y(ix, iy - 1, iz)) * &
                            inE%cvArray(imgrid)%grid%delYinv(iy)      &
                            +  (condE%rvArray(imgrid)%z(ix,iy,iz)*inE%cvArray(imgrid)%z(ix, iy, iz) -      &
                            condE%rvArray(imgrid)%z(ix,iy,iz - 1)*inE%cvArray(imgrid)%z(ix, iy, iz - 1)) * &
                            inE%cvArray(imgrid)%grid%delZinv(iz)
                    enddo   ! iz
                 enddo      ! iy
              enddo         ! ix

           else
              write(0, *) 'Error: DivC: scalars not same size'
           end if

          enddo ! loop over sub-grids
        else
           write(0, *) 'Error: DivC: output scalar not compatible use'
        end if
      end subroutine DivC	! DivC
      ! ************************************************************************************************************

      subroutine UpdateDivC(inE, outSc, imgrid)
      ! created by Cherevatova (May,2012)
      ! compute DivC on the interfaces between sub-grids
        implicit none
        type (cvector_mg), intent(in)   :: inE
        type (cscalar_mg), intent(inout)   :: outSc
        integer, intent(in)  :: imgrid
        ! local variables
        integer  :: ix, iy
        integer  :: nz

          if(imgrid == 1) then
            ! do not need to compute DivC for z=1 in 1st sub-grid
            ! nothing to do
            return
          endif

          nz = mGrid%gridArray(imgrid-1)%nz
          if (mGrid%coarseness(imgrid).gt.mGrid%coarseness(imgrid-1))then
            ! interface layer: finer to coarser
            ! current sub-grid is coarser
            ! Update iz=1 taking values from iz=nz of the previous sub-grid
           do iy = 2, mGrid%gridArray(imgrid)%ny
              do ix = 2, mGrid%gridArray(imgrid)%nx
                ! FOR NODES IN THE AIR
                if(mGrid%gridArray(imgrid)%flag == 1 .or. mGrid%gridArray(imgrid)%flag == 0) then

                   outSc%csArray(imgrid)%v(ix, iy, 1) = &
                      SIGMA_AIR*(inE%cvArray(imgrid)%x(ix,iy,1)-inE%cvArray(imgrid)%x(ix-1,iy,1)) * &
                        inE%cvArray(imgrid)%grid%delXinv(ix)    &
                        + SIGMA_AIR*(inE%cvArray(imgrid)%y(ix,iy,1)-inE%cvArray(imgrid)%y(ix,iy-1,1)) * &
                        inE%cvArray(imgrid)%grid%delYinv(iy)    &
                        + SIGMA_AIR*(inE%cvArray(imgrid)%z(ix,iy,1)-inE%cvArray(imgrid-1)%z(2*ix-1,2*iy-1,nz)) * &
                       inE%cvArray(imgrid)%grid%delZinv(1)

               ! FOR NODES INSIDE THE EARTH ONLY
               else if(mGrid%gridArray(imgrid)%flag == 2) then
                   outSc%csArray(imgrid)%v(ix, iy, 1) = &
                          (condE%rvArray(imgrid)%x(ix,iy, 1)*inE%cvArray(imgrid)%x(ix, iy, 1) -        &
                            condE%rvArray(imgrid)%x(ix-1, iy, 1)*inE%cvArray(imgrid)%x(ix-1, iy, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(ix)      &
                            +  (condE%rvArray(imgrid)%y(ix, iy, 1)*inE%cvArray(imgrid)%y(ix, iy, 1) -  &
                            condE%rvArray(imgrid)%y(ix,iy-1, 1)*inE%cvArray(imgrid)%y(ix, iy-1, 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(iy)      &
                            +  (condE%rvArray(imgrid)%z(ix,iy,1)*inE%cvArray(imgrid)%z(ix, iy, 1)-&
                                condE%rvArray(imgrid-1)%z(2*ix-1, 2*iy-1 , nz)*inE%cvArray(imgrid-1)%z(2*ix-1, 2*iy-1, nz)) * &
                            inE%cvArray(imgrid)%grid%delZinv(1)
                endif
             enddo ! ix
           enddo  ! iy

         else if (mGrid%coarseness(imgrid).lt.mGrid%coarseness(imgrid-1)) then
           ! interface : coarser to finer
           ! current grid is finer
           do iy = 2, mGrid%gridArray(imgrid-1)%ny
             do ix = 2, mGrid%gridArray(imgrid-1)%nx
                ! for nodes in the air
                  if(mGrid%gridArray(imgrid)%flag == 1 .or. mGrid%gridArray(imgrid)%flag == 0) then
                  ! odd/odd nodes
                    outSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, 1) = &
                      SIGMA_AIR*(inE%cvArray(imgrid)%x(2*ix-1, 2*iy-1, 1)-inE%cvArray(imgrid)%x(2*(ix-1),2*iy-1, 1)) * &
                        inE%cvArray(imgrid)%grid%delXinv(2*ix-1)    &
                        + SIGMA_AIR*(inE%cvArray(imgrid)%y(2*ix-1,2*iy-1,1)-inE%cvArray(imgrid)%y(2*ix-1,2*(iy-1), 1)) * &
                        inE%cvArray(imgrid)%grid%delYinv(2*iy-1)    &
                        + SIGMA_AIR*(inE%cvArray(imgrid)%z(2*ix-1,2*iy-1,1)-inE%cvArray(imgrid-1)%z(ix, iy, nz)) * &
                       inE%cvArray(imgrid)%grid%delZinv(1)
                  ! Since we do not know E(z-1) fields in many nodes of the finer grid
                  ! we suppose that these are E(z-1) = E(z)
                  ! thus Ex term in the DivC are zero
                  ! even/odd
                    outSc%csArray(imgrid)%v(2*ix, 2*iy-1, 1) = &
                      SIGMA_AIR*(inE%cvArray(imgrid)%x(2*ix, 2*iy-1, 1)-inE%cvArray(imgrid)%x(2*ix-1,2*iy-1, 1)) * &
                        inE%cvArray(imgrid)%grid%delXinv(2*ix)    &
                        + SIGMA_AIR*(inE%cvArray(imgrid)%y(2*ix,2*iy-1,1)-inE%cvArray(imgrid)%y(2*ix,2*(iy-1), 1)) * &
                        inE%cvArray(imgrid)%grid%delYinv(2*iy-1)
                  ! even/even
                   outSc%csArray(imgrid)%v(2*ix, 2*iy, 1) = &
                      SIGMA_AIR*(inE%cvArray(imgrid)%x(2*ix, 2*iy, 1)-inE%cvArray(imgrid)%x(2*ix-1,2*iy, 1)) * &
                        inE%cvArray(imgrid)%grid%delXinv(2*ix)    &
                        + SIGMA_AIR*(inE%cvArray(imgrid)%y(2*ix,2*iy,1)-inE%cvArray(imgrid)%y(2*ix,2*iy-1, 1)) * &
                        inE%cvArray(imgrid)%grid%delYinv(2*iy)
                  ! odd/even
                  outSc%csArray(imgrid)%v(2*ix-1, 2*iy, 1) = &
                      SIGMA_AIR*(inE%cvArray(imgrid)%x(2*ix-1, 2*iy, 1)-inE%cvArray(imgrid)%x(2*(ix-1),2*iy, 1)) * &
                        inE%cvArray(imgrid)%grid%delXinv(2*ix-1)    &
                        + SIGMA_AIR*(inE%cvArray(imgrid)%y(2*ix-1,2*iy,1)-inE%cvArray(imgrid)%y(2*ix-1,2*iy-1, 1)) * &
                        inE%cvArray(imgrid)%grid%delYinv(2*iy)

                  else if(mGrid%gridArray(imgrid)%flag == 2) then
                  ! FOR NODES INSIDE THE EARTH ONLY
                      ! odd/odd
                       outSc%csArray(imgrid)%v(2*ix-1, 2*iy-1, 1) = &
                            (condE%rvArray(imgrid)%x(2*ix-1, 2*iy-1, 1)*inE%cvArray(imgrid)%x(2*ix-1, 2*iy-1, 1) -     &
                            condE%rvArray(imgrid)%x(2*(ix-1), 2*iy-1, 1)*inE%cvArray(imgrid)%x(2*(ix-1), 2*iy-1, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2*ix-1)      &
                            +  (condE%rvArray(imgrid)%y(2*ix-1, 2*iy-1,1)*inE%cvArray(imgrid)%y(2*ix-1, 2*iy-1, 1) -   &
                            condE%rvArray(imgrid)%y(2*ix-1, 2*(iy-1), 1)*inE%cvArray(imgrid)%y(2*ix-1, 2*(iy-1), 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2*iy-1)      &
                            +  (condE%rvArray(imgrid)%z(2*ix-1,2*iy-1,1)*inE%cvArray(imgrid)%z(2*ix-1, 2*iy-1, 1) -    &
                            condE%rvArray(imgrid-1)%z(ix,iy,nz)*inE%cvArray(imgrid-1)%z(ix, iy, nz)) * &
                            inE%cvArray(imgrid)%grid%delZinv(1)
                      ! even/odd
                      outSc%csArray(imgrid)%v(2*ix, 2*iy-1, 1) = &
                            (condE%rvArray(imgrid)%x(2*ix, 2*iy-1, 1)*inE%cvArray(imgrid)%x(2*ix, 2*iy-1, 1) -     &
                            condE%rvArray(imgrid)%x(2*ix-1, 2*iy-1, 1)*inE%cvArray(imgrid)%x(2*ix-1, 2*iy-1, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2*ix)      &
                            +  (condE%rvArray(imgrid)%y(2*ix, 2*iy-1,1)*inE%cvArray(imgrid)%y(2*ix, 2*iy-1, 1) -   &
                            condE%rvArray(imgrid)%y(2*ix, 2*(iy-1), 1)*inE%cvArray(imgrid)%y(2*ix, 2*(iy-1), 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2*iy-1)
                      ! even/even
                      outSc%csArray(imgrid)%v(2*ix, 2*iy, 1) = &
                            (condE%rvArray(imgrid)%x(2*ix, 2*iy, 1)*inE%cvArray(imgrid)%x(2*ix, 2*iy, 1) -     &
                            condE%rvArray(imgrid)%x(2*ix-1, 2*iy, 1)*inE%cvArray(imgrid)%x(2*ix-1, 2*iy, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2*ix)      &
                            +  (condE%rvArray(imgrid)%y(2*ix, 2*iy,1)*inE%cvArray(imgrid)%y(2*ix, 2*iy, 1) -   &
                            condE%rvArray(imgrid)%y(2*ix, 2*iy-1, 1)*inE%cvArray(imgrid)%y(2*ix, 2*iy-1, 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2*iy)
                      ! odd/even
                      outSc%csArray(imgrid)%v(2*ix-1, 2*iy, 1) = &
                            (condE%rvArray(imgrid)%x(2*ix-1, 2*iy, 1)*inE%cvArray(imgrid)%x(2*ix-1, 2*iy, 1) -     &
                            condE%rvArray(imgrid)%x(2*(ix-1), 2*iy, 1)*inE%cvArray(imgrid)%x(2*(ix-1), 2*iy, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2*ix-1)      &
                            +  (condE%rvArray(imgrid)%y(2*ix-1, 2*iy,1)*inE%cvArray(imgrid)%y(2*ix-1, 2*iy, 1) -   &
                            condE%rvArray(imgrid)%y(2*ix-1, 2*iy-1, 1)*inE%cvArray(imgrid)%y(2*ix-1, 2*iy-1, 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2*iy)
                  endif
              enddo ! ix
            enddo ! iy

            ! need to fill in additional nodes
            if(mGrid%gridArray(imgrid)%flag == 1 .or. mGrid%gridArray(imgrid)%flag == 0) then
                ! iy= 1
                do ix= 2, mGrid%gridArray(imgrid-1)%nx
                   outSc%csArray(imgrid)%v(2*ix-1, 2, 1) = &
                          SIGMA_AIR*(inE%cvArray(imgrid)%x(2*ix-1, 2, 1)-inE%cvArray(imgrid)%x(2*(ix-1),2, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2*ix-1)    &
                            + SIGMA_AIR*(inE%cvArray(imgrid)%y(2*ix-1,2,1)-inE%cvArray(imgrid)%y(2*ix-1,1, 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2)

                   outSc%csArray(imgrid)%v(2*ix, 2, 1) = &
                          SIGMA_AIR*(inE%cvArray(imgrid)%x(2*ix, 2, 1)-inE%cvArray(imgrid)%x(2*ix-1,2, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2*ix)    &
                            + SIGMA_AIR*(inE%cvArray(imgrid)%y(2*ix,2,1)-inE%cvArray(imgrid)%y(2*ix,1, 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2)
                enddo ! ix
                outSc%csArray(imgrid)%v(2, 2, 1) =  outSc%csArray(imgrid)%v(3, 2, 1)
                ! ix = 1
                do iy = 2, mGrid%gridArray(imgrid-1)%ny
                  outSc%csArray(imgrid)%v(2, 2*iy-1, 1) = &
                          SIGMA_AIR*(inE%cvArray(imgrid)%x(2, 2*iy-1, 1)-inE%cvArray(imgrid)%x(1,2*iy-1, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2)    &
                            + SIGMA_AIR*(inE%cvArray(imgrid)%y(2,2*iy-1,1)-inE%cvArray(imgrid)%y(2,2*(iy-1), 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2*iy-1)
                  outSc%csArray(imgrid)%v(2, 2*iy, 1) = &
                          SIGMA_AIR*(inE%cvArray(imgrid)%x(2, 2*iy, 1)-inE%cvArray(imgrid)%x(1,2*iy, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2)    &
                            + SIGMA_AIR*(inE%cvArray(imgrid)%y(2,2*iy,1)-inE%cvArray(imgrid)%y(2,2*iy-1, 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2*iy)
                enddo !iy
            elseif(mGrid%gridArray(imgrid)%flag == 2) then
                 ! iy= 1
                do ix= 2, mGrid%gridArray(imgrid-1)%nx
                       outSc%csArray(imgrid)%v(2*ix-1, 2, 1) = &
                            (condE%rvArray(imgrid)%x(2*ix-1, 2, 1)*inE%cvArray(imgrid)%x(2*ix-1, 2, 1) -     &
                            condE%rvArray(imgrid)%x(2*(ix-1), 2, 1)*inE%cvArray(imgrid)%x(2*(ix-1), 2, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2*ix-1)      &
                            +  (condE%rvArray(imgrid)%y(2*ix-1, 2,1)*inE%cvArray(imgrid)%y(2*ix-1, 2, 1) -   &
                            condE%rvArray(imgrid)%y(2*ix-1, 1, 1)*inE%cvArray(imgrid)%y(2*ix-1, 1, 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2)

                       outSc%csArray(imgrid)%v(2*ix, 2, 1) = &
                            (condE%rvArray(imgrid)%x(2*ix, 2, 1)*inE%cvArray(imgrid)%x(2*ix, 2, 1) -     &
                            condE%rvArray(imgrid)%x(2*ix-1, 2, 1)*inE%cvArray(imgrid)%x(2*ix-1, 2, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2*ix)      &
                            +  (condE%rvArray(imgrid)%y(2*ix, 2,1)*inE%cvArray(imgrid)%y(2*ix, 2, 1) -   &
                            condE%rvArray(imgrid)%y(2*ix, 1, 1)*inE%cvArray(imgrid)%y(2*ix, 1, 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2)
                enddo ! ix
                outSc%csArray(imgrid)%v(2, 2, 1) =  outSc%csArray(imgrid)%v(3, 2, 1)
                ! ix = 1
                do iy = 2, mGrid%gridArray(imgrid-1)%ny
                       outSc%csArray(imgrid)%v(2, 2*iy-1, 1) = &
                            (condE%rvArray(imgrid)%x(2, 2*iy-1, 1)*inE%cvArray(imgrid)%x(2, 2*iy-1, 1) -     &
                            condE%rvArray(imgrid)%x(1, 2*iy-1, 1)*inE%cvArray(imgrid)%x(1, 2*iy-1, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2)      &
                            +  (condE%rvArray(imgrid)%y(2, 2*iy-1,1)*inE%cvArray(imgrid)%y(2, 2*iy-1, 1) -   &
                            condE%rvArray(imgrid)%y(2, 2*(iy-1), 1)*inE%cvArray(imgrid)%y(2, 2*(iy-1), 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2*iy-1)

                        outSc%csArray(imgrid)%v(2, 2*iy, 1) = &
                            (condE%rvArray(imgrid)%x(2, 2*iy, 1)*inE%cvArray(imgrid)%x(2, 2*iy, 1) -     &
                            condE%rvArray(imgrid)%x(1, 2*iy, 1)*inE%cvArray(imgrid)%x(1, 2*iy, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(2)      &
                            +  (condE%rvArray(imgrid)%y(2, 2*iy,1)*inE%cvArray(imgrid)%y(2, 2*iy, 1) -   &
                            condE%rvArray(imgrid)%y(2, 2*iy-1, 1)*inE%cvArray(imgrid)%y(2, 2*iy-1, 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(2*iy)
                enddo !iy
            endif

         else if (mGrid%coarseness(imgrid).eq.mGrid%coarseness(imgrid-1)) then
           do ix = 2, mGrid%gridArray(imgrid)%ny
             do iy = 2, mGrid%gridArray(imgrid)%nx

                ! FOR NODES IN THE AIR
                  if(mGrid%gridArray(imgrid)%flag == 1 .or. mGrid%gridArray(imgrid)%flag == 0) then

                   outSc%csArray(imgrid)%v(ix, iy, 1) = &
                      SIGMA_AIR*(inE%cvArray(imgrid)%x(ix,iy,1)-inE%cvArray(imgrid)%x(ix-1,iy,1)) * &
                        inE%cvArray(imgrid)%grid%delXinv(ix)    &
                        + SIGMA_AIR*(inE%cvArray(imgrid)%y(ix,iy,1)-inE%cvArray(imgrid)%y(ix,iy-1,1)) * &
                        inE%cvArray(imgrid)%grid%delYinv(iy)    &
                        + SIGMA_AIR*(inE%cvArray(imgrid)%z(ix,iy,1)-inE%cvArray(imgrid-1)%z(ix,iy,nz)) * &
                       inE%cvArray(imgrid)%grid%delZinv(1)

                   else if(mGrid%gridArray(imgrid)%flag == 2) then
                    ! FOR NODES INSIDE THE EARTH ONLY
                    ! THE TOP MOST EARTH NODE HAS AN INTERFACE WITH
                    ! AIR, THEREFORE THAT ONE IS SKIPPED HERE
                       outSc%csArray(imgrid)%v(ix, iy,1) = &
                            (condE%rvArray(imgrid)%x(ix,iy,1)*inE%cvArray(imgrid)%x(ix, iy, 1) -         &
                            condE%rvArray(imgrid)%x(ix - 1,iy,1)*inE%cvArray(imgrid)%x(ix - 1, iy, 1)) * &
                            inE%cvArray(imgrid)%grid%delXinv(ix)      &
                            +  (condE%rvArray(imgrid)%y(ix,iy,1)*inE%cvArray(imgrid)%y(ix, iy, 1) -      &
                            condE%rvArray(imgrid)%y(ix,iy - 1,1)*inE%cvArray(imgrid)%y(ix, iy - 1, 1)) * &
                            inE%cvArray(imgrid)%grid%delYinv(iy)      &
                            +  (condE%rvArray(imgrid)%z(ix,iy,1)*inE%cvArray(imgrid)%z(ix, iy, 1) -      &
                            condE%rvArray(imgrid-1)%z(ix,iy,nz)*inE%cvArray(imgrid-1)%z(ix, iy,nz)) * &
                            inE%cvArray(imgrid)%grid%delZinv(1)
                   endif
                 enddo      ! iy
              enddo         ! ix
         endif
      end  subroutine UpdateDivC

    ! *********************************************************************************************************
    end  module modelOperator3D
