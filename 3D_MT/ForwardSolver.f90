module ForwardSolver

!  High level interface/control module used by top level routines
!   for initializing and using the solver.  The key public routines
!   in this module have only generic (abstract data type) arguments
!   and can thus be called from the top-level inversion routines.
!  A similar solver interface will be required for any specific use
!   of the top-level inversion modules
!  All public routines must have the same generic names, parameter lists
!   and abstract functionality; this implementation is for 3D-MT
!
!  Main difference from 2D: no need to keep track of TE/TM modes

use math_constants
use datafunc
use dataspace
use solnspace
use emsolve3d
use transmitters

implicit none

! options for boundary conditions input or computation
logical, save, public   :: COMPUTE_BC = .true.
logical, save, public   :: NESTED_BC = .false.
logical, save, public   :: BC_FROM_RHS_FILE = .false.
logical, save, public   :: BC_FROM_E0_FILE = .false.

! option to read a primary solution from file for SFF
logical, save, public   :: PRIMARY_E_FROM_FILE = .false.

!=======================================================================
!Mar. 13, 2011============================== Use only for MT Calculation
!=======================================================================

!! Used for interpolating the BC from a larger grid.
  type(grid_t),save,public                   ::  Larg_Grid   
  type(solnVectorMTX_t),save,public          ::  eAll_larg
  integer              ,save,public          ::  nTx_nPol
  logical              ,save,public          ::  nestedEM_initialized

!=======================================================================
!May 15, 2018== AK == New general RHS for all transmitters now stored
!=======================================================================
  type(rhsVectorMTX_t),save,public           ::  bAll

!=======================================================================
!Aug 18, 2021== AK == Also store an array of primary fields (not good
!in terms of memory usage, will read each from file as needed ...)
!=======================================================================
  type(solnVectorMTX_t),save,public           ::  eAllPrimary
  type(modelParam_t),save,public              ::  sigmaPrimary

!  initialization routines (call Fwd version if no sensitivities are
!     are calculated).  Note that these routines are set up to
!    automatically manage memory and to figure out which initialization
!    (or reinitialization) steps are required (e.g., when the frequency
!    changes from the previous solver call, appropriate solver
!    coefficients are updated, matrices factored, etc.).  This
!    functionality needs to be maintained in implementations for new
!    problems!

!=======================================================================
!Aug 18, 2021== AK == Local variables used to SFF (originally CSEM only) 
!=======================================================================
   type(rvector), save, private :: edgeCond ! Full conductivity (on edge)
   type(rvector), save, private :: condAnomaly ! Anomalous conductivity (on edge)
   type(cvector), save, private :: E_p         ! Primary field
public initSolver

!  cleanup/deallocation routines
public exitSolver

! solver routines
public fwdSolve, sensSolve

logical, save, private		:: modelDataInitialized = .false.
logical, save, private		:: BC_from_file_Initialized = .false.
!  logical, save, private		:: sigmaNotCurrent = .true.

Contains

!**********************************************************************
subroutine init_BC_from_file(grid,nTx_nPol)
  type(grid_t), intent(in)        :: grid
  integer, intent(in)             :: nTx_nPol

  ! local
  type(cboundary)   :: BC
  integer           :: iTx,j,status
  

     call create_cboundary(grid,BC)
     
     
     allocate (BC_from_file(nTx_nPol), STAT=status)
     do j=1,nTx_nPol
       BC_from_file(j)=BC
     end do
     BC_from_file_Initialized=.true.
     

end subroutine init_BC_from_file

!**********************************************************************
subroutine pack_BC_from_file(grid,bAll,nTx_nPol)
  ! The boundary conditions are sent around through the MPI messaging
  ! using the linear BC_from_file array (defined in nestedEM). Here,
  ! we pack this into an array of BCs. Used when bAll are directly
  ! read or extracted from an E file. [AK 2/24/2023]

  type(grid_t), intent(in)                :: grid
  type(rhsVectorMTX_t), intent(in)        :: bAll
  integer, intent(inout)                  :: nTx_nPol

  ! local
  integer           :: iTx,iMode,counter

    if (.not. BC_from_file_Initialized) then
       call count_rhsVectorMTX(bAll,nTx_nPol)
       call init_BC_from_file(grid,nTx_nPol)
    else if (.not. bAll%allocated) then
       call warning('bAll needs to be allocated before calling pack_BC_from_file in ForwardSolver')
    end if

    counter=0
    do iTx = 1,bAll%nTx
       do iMode = 1,bAll%combs(iTx)%nPol
         counter=counter+1
         BC_from_file(counter) = bAll%combs(iTx)%b(iMode)%bc
       end do
    end do

    if (counter .ne. nTx_nPol) then
        call warning('Boundary conditions are corrupted in the MPI message transfer')
    end if

end subroutine pack_BC_from_file

!**********************************************************************
subroutine unpack_BC_from_file(bAll)
  ! The boundary conditions are sent around through the MPI messaging
  ! using the linear BC_from_file array (defined in nestedEM). Here,
  ! we unpack this back into a bAll structure to use for specific
  ! iTx and iMode combination. This is needed to enable flexibility
  ! in the number of modes between transmitters... which is otherwise
  ! eliminated by storing the BCs in an array. [AK 2/24/2023]

  type(rhsVectorMTX_t), intent(inout)        :: bAll

  ! local
  integer           :: iTx,iMode,counter

    if (.not. BC_from_file_Initialized) then
       call warning('BC_from_file not initialized: cannot create bAll in ForwardSolver')
    else if (.not. bAll%allocated) then
       call warning('bAll needs to be allocated before calling unpack_BC_from_file in ForwardSolver')
    end if

    counter=0
    do iTx = 1,bAll%nTx
       do iMode = 1,bAll%combs(iTx)%nPol
         counter=counter+1
         bAll%combs(iTx)%b(iMode)%bc = BC_from_file(counter)
       end do
    end do

    if (counter .ne. nTx_nPol) then
        call warning('Boundary conditions are corrupted in the MPI message transfer')
    end if

end subroutine unpack_BC_from_file


   !**********************************************************************
   subroutine initSolver(iTx,sigma,grid,e0,e,comb)
   !   Initializes forward solver for transmitter iTx.
   !     Idea is to call this before calling fwdSolve or sensSolve,
   !     in particular before the first solution for each transmitter
   !     (frequency).  If called for the first time (in a program run,
   !     or after a call to exitSolver), full initialization
   !     (after deallocation/cleanup if required) is performed.
   !
   !   iTx defines transmitter: for 2D MT, this provides info about
   !       frequency and TE/TM mode; for 3D frequency and number
   !       of polarizations
   !
   !   This now does all setup (including matrix factorization) for
   !     the appropriate mode/frequency
   !   NOTE: e and comb are optional calling arguments;
   !     both should be present if one is

   integer, intent(in)				            :: iTx
   type(modelParam_t),intent(in), target		:: sigma
   type(grid_t), intent(in), target             :: grid
   !  following structures are initialized
   !	solution vector for forward problem
   type(solnVector_t), intent(inout)			:: e0
   !	solution vector for sensitivity
   type(solnVector_t), intent(inout), optional	:: e
   !	forcing for sensitivity
   type(rhsVector_t), intent(inout), optional   :: comb

   !  local variables
   integer		:: IER,k
   character*80 :: gridType
   logical		:: initForSens,sigmaNotCurrent


   initForSens = present(comb)

   !  allocate for background solution
   call create_solnVector(grid,iTx,e0)

   if(initForSens) then
      !  allocate for sensitivity solution, RHS - same for all TX types
      !  assuming here that we don't use sparse storage ... we could!
      call create_solnVector(grid,iTx,e)
      comb%nonzero_source = .true.
      comb%sparse_source = .false.
      comb%nonzero_bc = .false.
      call create_rhsVector(grid,iTx,comb)
!      do k = 1,comb%nPol
!        comb%b(k)%sparse_Source = .false.
!        comb%b(k)%adj = ''
!        !  using all this information, reallocate storage for each polarization
!        call create_RHS(grid,iTx,comb%b(k))
!      enddo
   endif

   if(.NOT.modelDataInitialized) then
   !   Initialize modelData, setup model operators
      call ModelDataInit(grid)
      call ModelOperatorSetup()
      modelDataInitialized = .true.
	  
   !=======================================================================
   !Mar. 13, 2011================================= Use for CSEM Calculation
   !=======================================================================
   !  call setPrimaryCond(sigma, .true.) ! Originally created by Aihua;
   ! This subroutine can create the primary model by averaging Sigma
   ! or read it from file. If the second argument is true, this routine
   ! will read Primary model from file.
   endif
   !if (txDict(iTx)%Tx_type=='CSEM') then	
   !   xTx1D = txDict(iTx)%xyzTx(1)
   !   yTx1D = txDict(iTx)%xyzTx(2)   
   !   Call set1DModel(sigma,xTx1D,yTx1D)
   !end if
   

!    the following needs work ... want to avoid reinitializing
!     operator coefficients when conductivity does not change;
!     need to have a way to reset sigmaNotCurrent to false when
!     conductivity changes (one idea: assign a random number whenever
!     a conductivity parameter is modified (by any of the routines in
!     module ModelSpace); store this in the modelOperator module (which
!     is where updateCond sits) and have updateCond compare the random key
!     with what is stored)
!  if(sigmaNotCurrent) then
       call updateCond(sigma)
!      sigmaNotCurrent = .false.
!   endif

   if (txDict(iTx)%Tx_type=='SFF') then
      ! compute sigma-sigma1D for the source... NOT PHYSICAL!
      !Call linComb_modelParam(ONE,sigma,MinusONE,sigmaPrimary,sigmaTemp)
      ! sigmaTemp is the anomalous conductivity, map it onto edges
      !Call ModelParamToEdge(sigmaTemp,condAnomaly)
      ! sigmaTemp is the anomalous conductivity, map it onto edges
      Call ModelParamToEdge(sigma,edgeCond)
      Call ModelParamToEdge(sigmaPrimary,condAnomaly)
      ! compute sigma-sigma1D for the source... NOT PHYSICAL!
      Call linComb_rvector(ONE,edgeCond,MinusONE,condAnomaly,condAnomaly)

      write(0,'(a35,3i5,a2)') 'Conductivity anomaly dimensions: [',condAnomaly%nx,condAnomaly%ny,condAnomaly%nz,' ]'
   end if
   ! This needs to be called before solving for a different frequency
   !!!!!!!  BUT AFTER UPDATECOND !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call UpdateFreq(txDict(iTx)%omega)

   end subroutine initSolver

   !**********************************************************************
   subroutine exitSolver(e0,e,comb)
   !   deallocates b0, comb, e0, e and solver arrays
   type(solnVector_t), intent(inout), optional  :: e0
   type(solnVector_t), intent(inout), optional	::e
   type(rhsVector_t), intent(inout), optional	::comb

   ! local variables
   logical			:: initForSens

   initForSens = present(comb)

   if(present(e0)) then
      call deall_solnVector(e0)
   endif

   if(initForSens) then
      call deall_rhsVector(comb)
      call deall_solnVector(e)
   endif

   if(modelDataInitialized) then
      ! cleanup/deallocation routines for model operators
      call ModelDataCleanUp() ! FWD/modelOperator3D.f90
      call ModelOperatorCleanUp() ! FWD/EMsolve3D.f90
      modelDataInitialized = .false.
   endif

   end subroutine exitSolver


!**********************************************************************
! Sets boundary conditions in RHS b0, and the initial conditions in e0,
! for a transmitter index iTx. Typically, this is done for all polarizations.
! However in the MPI case, e0%nPol is artificially set to 1 to ensure optimal
! load distribution, so need to be careful with the indexing here.
! In all cases, this is meant to be called just before fwdSolve ONLY.
! The difference with fwdSolve is that here b0 is computed; in fwdSolve
! it is input only.
!
! Boundary conditions would have already been initialized from file
! into bAll for option BC_FROM_RHS_FILE. Otherwise, need to be computed.
! At present, we are not using this routine to set up the initial value
! of e0, but we might do so in the future.
!
! We are doing all this in a separate routine because we want this to be
! a high-level function that can be called from the upper level.
!
! A. Kelbert, 24 May 2018; last edited 18 Aug 2021
  Subroutine fwdSetup(iTx,e0,b0)

    !  Input mode, period
    integer, intent(in)     :: iTx
    ! Output electric field first guess (for iterative solver)
    type(solnVector_t), intent(inout) :: e0
    ! Output boundary conditions
    type(rhsVector_t), intent(inout)  :: b0

    ! local
    type(cboundary)     :: BC
    integer             :: iMode,j
    real(kind=prec)     :: omega
    complex(kind=prec)  :: i_omega_mu

    ! local variables for TIDE
   integer		 :: ios,istat
   character*80 :: file_name,comment
   character*2  :: tidal_component
   logical		 :: exists
   type(sparsevecc) :: jInt

    ! For RHS vector, the major differences between MT and CSEM are
    ! (1) Boundary condition for MT isn't zeros while CSEM is zeros
    ! (2) CSEM has a source term while that of MT is zeros
    ! Initialize the RHS vector; should we always clean it up on input?
    if (.not. b0%allocated) then
      select case (txDict(iTx)%Tx_type)
      case ('CSEM','SFF')
        b0%nonzero_Source = .true.
        b0%sparse_Source = .false.
        b0%nonzero_BC = .false.
      case ('MT')
        b0%nonzero_Source = .false.
        b0%sparse_Source = .false.
        b0%nonzero_BC = .true.
      case('TIDE')
        b0%nonzero_Source = .true.
        b0%sparse_Source = .false.
        b0%nonzero_BC = .true.
      case default
        write(0,*) node_info,'Unknown FWD problem type',trim(txDict(iTx)%Tx_type),'; unable to initialize RHS'
      end select
      call create_rhsVector(e0%grid,iTx,b0)
    end if


    ! careful here with imode indexing. In MPI, we trick the code into thinking
    ! that there is only one mode for each processor. So instead of using plain indexing
    ! to determine the mode, we use Pol_index integer variable.
    do j = 1,e0%nPol
        iMode = e0%Pol_index(j)

        select case (txDict(iTx)%Tx_type)

            case ('CSEM')
                !  set period, complete setup of 3D EM equation system
                omega = txDict(iTx)%omega
                i_omega_mu = cmplx(0.,ISIGN*MU_0*omega,kind=prec)
                ! Now finish up the computation of the general b0%s = - ISIGN * i\omega\mu_0 j
                call diagMult(condAnomaly,E_P,b0%b(j)%s)
                call scMult(-i_omega_mu,b0%b(j)%s,b0%b(j)%s)

            case ('SFF')
                ! this is currently implemented only for 1 mode - check for this...
                !if (iMode .ne. 1) then
                !  write(0,*) 'ERROR: SFF only implemented for one mode at present. Exiting...'
                !  stop
                !end if
                ! we've read eAllPrimary from EM soln file already
                E_P = eAllPrimary%solns(iTx)%pol(iMode)
                !  set period, complete setup of 3D EM equation system
                omega = txDict(iTx)%omega
                i_omega_mu = cmplx(0.,ISIGN*MU_0*omega,kind=prec)
                ! Now finish up the computation of the general b0%s = - ISIGN * i\omega\mu_0 j
                call diagMult(condAnomaly,E_P,b0%b(j)%s)
                call scMult(-i_omega_mu,b0%b(j)%s,b0%b(j)%s)

            case ('MT')
                if (BC_FROM_RHS_FILE .or. BC_FROM_E0_FILE .or. NESTED_BC) then
                    ! in all these cases, we've read or interpolated BCs from a file; stored them in BC_from_file
                    ! which was MPI messaged to all nodes, then unpacked back into bAll. This allows an arbitrary
                    ! number of modes per transmitter (2 modes used to be hardcoded) but the period and mode
                    ! layout needs to be consistent throughout. [AK 2/24/2023]
                    write (*,'(a12,a29,a12,i4,a15,i2)') node_info, 'Setting the BC from RHS or E-field file ', &
                        ' for period ',iTx,' & mode # ',iMode
                    BC = bAll%combs(iTx)%b(iMode)%bc

                elseif (COMPUTE_BC) then
                    ! For e0 and b0, use the same fake polarization index j for MPI modeling context
                    write (*,'(a12,a28,a12,i4,a15,i2)') node_info, 'Computing the BC internally ', &
                        ' for period ',iTx,' & mode # ',iMode
                    BC = b0%b(j)%bc
                    call ComputeBC(iTx,iMode,e0%pol(j),BC)

                end if
                ! store the BC in b0 and set up the forward problem - use fake indexing in MPI
                b0%b(j)%adj = 'FWD'
                b0%b(j)%bc = BC
                call deall_cboundary(BC)

            case ('TIDE') ! in the future, may use the BC options from MT block above

               file_name = trim(txDict(iTx)%id)//'.source'
               inquire(FILE=file_name,EXIST=exists)
               if (exists) then
                  write(*,*) node_info,'Reading source - i \omega \mu \sigma_E (v x B) from interior source file: ',trim(file_name)
                  open(ioREAD,file=file_name,status='unknown',form='formatted',iostat=ios)
                  read(ioREAD,'(a35)',iostat=istat) comment
                  read(ioREAD,'(a2)',iostat=istat) tidal_component
                  if (tidal_component .ne. trim(txDict(iTx)%id)) then
                     write(0,*) node_info,'Warning: tidal component ',tidal_component,' is read from file ',trim(file_name)
                  end if
                  call read_sparsevecc(ioREAD,jInt)
                  close(ioREAD)
               end if

               ! Assume that the source - i \omega \mu \sigma_E (v x B) in on the edges
               if(jInt%allocated) then
                  write(0,*) node_info,'Using interior forcing to compute the RHS for the FWD problem'
                  call add_scvector(ISIGN*C_ONE,jInt,b0%b(j)%s)
                  call deall_sparsevecc(jInt)
               end if

            case default
                write(0,*) node_info,'Unknown FWD problem type',trim(txDict(iTx)%Tx_type),'; unable to compute RHS'
        end select
    end do

  end subroutine fwdSetup

   !**********************************************************************
   subroutine fwdSolve(iTx,e0,b0,device_id,comm_local)

   !  driver for 3d forward solver; sets up for transmitter iTx, returns
   !   solution in e0 ; rhs vector (b0) is generated locally--i.e.
   !   boundary conditions are set internally (NOTE: could use transmitter
   !   dictionary to indireclty provide information about boundary
   !    conditions.  Presently we set BC using WS approach.
   !
   ! this *should* works with the SP/SP2/SPETSc2 versions as well

   integer, intent(in)		:: iTx
   type(solnVector_t), intent(inout)	:: e0
   type(rhsVector_t), intent(in)        :: b0
   ! NOTE: this is only needed for two layer parallelization
   ! normal (none-petsc) programmes have no need of this
   integer, intent(in), optional        :: device_id
   integer, intent(in), optional        :: comm_local

   ! local variables
   real(kind=prec)	:: omega
   integer			:: IER,iMode
   complex(kind=prec)	:: i_omega_mu

   omega = txDict(iTx)%omega
   !  set period, complete setup of 3D EM equation system
   i_omega_mu = cmplx(0.,ISIGN*MU_0*omega,kind=prec)

   !  complete operator intialization, for this frequency
   if (txDict(iTx)%Tx_type=='CSEM') then 

      ! THIS IS THE TEMPORARY SETUP FOR CSEM - NOT ACTIVE AT PRESENT [AK] 
      ! Now finish up the computation of the general b0%s = - ISIGN * i\omega\mu_0 j
      !call diagMult(condAnomaly,E_P,b0%s)
      !call scMult(-i_omega_mu,b0%s,b0%s)
      !   call forward solver, compute secondary field
      write(*,'(a12,a3,a25,i3,a14,es15.7)') 'Solving the ','FWD', &
                  ' problem for transmitter ',iTx,' at frequency ',txDict(iTx)%PERIOD
      call zero_solnVector(e0)
      call FWDsolve3D(b0%b(1),omega,e0%pol(1))

      !   add primary field to secondary field
      !e0%pol(1)=E_p
      call add(E_p,e0%pol(1),e0%pol(1))
      !term=1.0/10.0 ! txDict(iTx)%Moment  
      !call scMult(term,e0%pol(1),e0%pol(1))

   elseif (txDict(iTx)%Tx_type=='SFF') then 
 
      ! General b0%s = - ISIGN * i\omega\mu_0 (sigma-sigma1d) E1D already computed
      do iMode = 1,e0%nPol
         ! Extract primary solution again...
         E_P = eAllPrimary%solns(iTx)%pol(iMode)
		   ! call forward solver, compute secondary field
         ! set the starting solution to zero
		   ! NOTE that in the MPI parallelization, e0 may only contain a single mode;
		   ! mode number is determined by Pol_index, NOT by its order index in e0
		   ! ... but b0 uses the same fake indexing as e0
		   write (*,'(a12,a12,a3,a20,i4,a2,es13.6,a15,i2)') node_info, 'Solving the ','SFF', &
			   	' problem for period ',iTx,': ',(2*PI)/omega,' secs & mode # ',e0%Pol_index(iMode)
         call zero(e0%pol(iMode))
		   call FWDsolve3D(b0%b(iMode),omega,e0%pol(iMode))
		   write (6,*)node_info,'FINISHED solve, nPol',e0%nPol
         ! now add primary field to secondary field
         call add(E_p,e0%pol(iMode),e0%pol(iMode))
      enddo

   elseif ((txDict(iTx)%Tx_type=='MT') .or. (txDict(iTx)%Tx_type=='TIDE')) then
       !  loop over polarizations
       do iMode = 1,e0%nPol
       ! compute boundary conditions for polarization iMode
       !   uses cell conductivity already set by updateCond
       ! call setBound(iTx,e0%Pol_index(iMode),e0%pol(imode),b0%bc)
       ! NOTE that in the MPI parallelization, e0 may only contain a single mode;
       ! mode number is determined by Pol_index, NOT by its order index in e0
       ! ... but b0 uses the same fake indexing as e0
         write (*,'(a12,a12,a3,a20,i4,a2,es13.6,a15,i2)') node_info, &
                 'Solving the ','FWD', ' problem for period ',iTx,   &
                 ': ',(2*PI)/omega,' secs & mode # ',e0%Pol_index(iMode)
         if (present(device_id)) then
             if (present(comm_local)) then
                 call FWDsolve3D(b0%b(iMode),omega,e0%pol(iMode), device_id, &
    &               comm_local)
             else
                 call FWDsolve3D(b0%b(iMode),omega,e0%pol(iMode), device_id)
             end if
         else
             call FWDsolve3D(b0%b(iMode),omega,e0%pol(iMode))
         end if
         ! write (6,*) node_info,' finished solving, nPol', e0%nPol
         write (*,'(a12,a24,i4)') node_info, &
                 'finished solving, nPol' , e0%nPol
       enddo
   else
       write(0,*) node_info,'Unknown FWD problem type',trim(txDict(iTx)%Tx_type),'; unable to run fwdSolve'
   endif

   ! update pointer to the transmitter in solnVector
   e0%tx = iTx

   end subroutine fwdSolve

   !**********************************************************************
   subroutine sensSolve(iTx,FWDorADJ,e,comb,device_id,comm_local)
   !   Uses forcing input from comb, which must be set before calling
   !    solves forward or adjoint problem, depending on comb%ADJ
   !  NOTE that this routine DOES NOT call UpdateFreq routine to complete
   !   initialization of solver for a particular frequency.
   !  This final initialization step must (at present) be done by
   !    calling fwdSolve before calling this routine.
   !
   ! this *should* works with the SP/SP2/SPETSc2 versions as well

   integer, intent(in)          	:: iTx
   character*3, intent(in)		:: FWDorADJ
   type(solnVector_t), intent(inout)		:: e
   type(rhsVector_t), intent(inout)		:: comb
   ! NOTE: this is only needed for two layer parallelization
   ! normal (none-petsc) programmes have no need of this
   integer, intent(in), optional        :: device_id
   integer, intent(in), optional        :: comm_local

   ! local variables
   integer      			:: IER,iMode
   real(kind=prec) 		:: omega, period

!  zero starting solution, solve for all modes
   call zero_solnVector(e)
   
   if (txDict(iTx)%Tx_type=='MT' .or. txDict(iTx)%Tx_type=='CSEM' ) then 
   	omega = txDict(iTx)%omega
   	period = txDict(iTx)%period
      do iMode = 1,e%nPol
      	comb%b(e%Pol_index(iMode))%adj = FWDorADJ
      	write (*,'(a12,a12,a3,a20,i4,a2,es13.6,a15,i2)') node_info, &
     &          'Solving the ',FWDorADJ, ' problem for period ',iTx,&
     &          ': ',(2*PI)/omega,' secs & mode # ',e%Pol_index(iMode)
        if (present(device_id)) then
            if (present(comm_local)) then
                call FWDsolve3d(comb%b(e%Pol_index(iMode)),omega,       &
     &           e%pol(imode), device_id, comm_local)
            else
                call FWDsolve3d(comb%b(e%Pol_index(iMode)),omega,       &
     &           e%pol(imode), device_id)
            end if
        else
            call FWDsolve3d(comb%b(e%Pol_index(iMode)),omega,e%pol(imode))
        end if
      enddo
   else
       write(0,*) node_info,'Unknown FWD problem type',trim(txDict(iTx)%Tx_type),'; unable to run sensSolve'
   endif

   ! update pointer to the transmitter in solnVector
   e%tx = iTx

   end subroutine sensSolve


  !**********************************************************************
  ! uses nestedEM module to extract the boundary conditions directly from
  ! a full EMsolnMTX vector on a larger (and coarser) grid
  ! NOTE AK 24 May 2018:
  ! As coded by NM, it looks like the BC from large E-solution
  ! are the same for ALL periods and modes. This is actually not true.
  ! Here, BC were merely input used to create BC_from_file in nestedEM.f90.
  ! I now initialize them there.
  ! NOTE: at present, this sets up BC_from_file array that is stored in
  !       nestedEM module. No interaction with the bAll variable needed.
  !       Used to require b0 as input but that was merely for initialization.
  !       Anyway, we might want to rewrite this somewhat. [AK 20 May 2019]
  ! NOTE AK 24 Feb 2023: some cleaning; cannot replace BC_from_file
  ! completely with bAll because the linear array is needed for MPI messaging
  ! but we are now unpacking it before we run fwdSetup on worker nodes.

  subroutine Interpolate_BC_from_E_soln(eAll_larg,Larg_Grid,grid,nTx_nPol)

  type(grid_t)  ,intent(inout)              ::  Larg_Grid
  type(solnVectorMTX_t),intent(inout)       ::  eAll_larg
  type(grid_t),intent(in)                   ::  Grid
  integer, intent(inout)                    ::  nTx_nPol

    ! local variables needed for nesting calculations
  integer                      :: status,iMode,iTx,counter

  ! If eAll_larg solution is already allocated, use that to interpolate the BCs to smaller grid
  if (.not. eAll_larg%allocated) then
    call warning('Unable to interpolate BCs from the large grid to smaller: large grid not allocated')
  end if

  write(15,*) ' Start interpolating',grid%nx,grid%ny,grid%nz

  call count_solnVectorMTX(eAll_larg,nTx_nPol)

  call setup_BC_from_file(Grid,nTx_nPol)

  ! For now extract the BC from eAll_larg
    counter=0
    do iTx = 1, eAll_larg%nTx
       do iMode=1, eAll_larg%solns(iTx)%nPol
         counter=counter+1
         Call compute_BC_from_file(Larg_Grid,eAll_larg%solns(iTx)%pol(iMode),Grid,counter)
       end do
    end do

  ! Done with eAll_larg, deallocate it, and keep track that BC_from_file are already Initialized.
   call deall(eAll_larg)
   call deall_grid(Larg_Grid)
   BC_from_file_Initialized=.true.

   write(15,*) ' End interpolating',BC_from_file(1)%yXMin(10,11)

  end subroutine Interpolate_BC_from_E_soln


end module ForwardSolver
