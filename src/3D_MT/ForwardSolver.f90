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
use CSEM_module
implicit none

! options for boundary conditions input or computation
logical, save, public   :: COMPUTE_BC = .true.
logical, save, public   :: NESTED_BC = .false.
logical, save, public   :: BC_FROM_RHS_FILE = .false.
logical, save, public   :: BC_FROM_E0_FILE = .false.
type(modelParam_t),save,private		:: sigma_temp

!=======================================================================
!Mar. 13, 2011============================== Use only for MT Calculation
!=======================================================================

!! Used for interpolating the BC from a larger grid.
  type(grid_t),save,public               	 ::	 Larg_Grid   
  type(solnVectorMTX_t),save,public          ::  eAll_larg
  integer              ,save,public          ::   nTx_nPol
  logical              ,save,public          ::  nestedEM_initialized

!=======================================================================
!May 15, 2018== AK == New general RHS for all transmitters now stored
!=======================================================================
  type(rhsVectorMTX_t),save,public           ::  bAll

!  initialization routines (call Fwd version if no sensitivities are
!     are calculated).  Note that these routines are set up to
!    automatically manage memory and to figure out which initialization
!    (or reinitialization) steps are required (e.g., when the frequency
!    changes from the previous solver call, appropriate solver
!    coefficients are updated, matrices factored, etc.).  This
!    functionality needs to be maintained in implementations for new
!    problems!
public initSolver

!  cleanup/deallocation routines
public exitSolver

! solver routines
public fwdSolve, sensSolve
#ifdef PETSC
public fwdSolve_petsc,sensSolve_petsc
#endif

logical, save, private		:: modelDataInitialized = .false.
logical, save, private		:: BC_from_file_Initialized = .false.
!  logical, save, private		:: sigmaNotCurrent = .true.

Contains
   !**********************************************************************
   ! edited by AK 24 May 2018 to store the output in bAll. Strongly
   ! suspect that BC_from_file can be completely replaced by using bAll
   ! but for now, trying to minimize changes to the extent possible...
   ! cleaning up for nested modeling is not complete yet.
   ! Note that as coded by NM, it looks like the BC from large E-solution
   ! are the same for ALL periods and modes. This is actually not true.
   ! Here, BC were merely input used to create BC_from_file in nestedEM.f90.
   ! I now initialize them there.
subroutine Interpolate_BC(grid)
  type(grid_t), intent(in)        :: grid

  ! local
  integer       :: iTx,ix,iy,iz
  

  !In case of interpolating the BC from eAll_larg   
  ! If eAll_ larg solution is already allocated, then use that to interpolate the BC from it
     
  if (eAll_larg%allocated) then
  	write(15,*) ' Start interpolating',grid%nx,grid%ny,grid%nz
	call Interpolate_BC_from_E_soln (eAll_larg,Larg_Grid,Grid)
        !Once we are ready from eAll_larg, deallocate it, and keep track, that BC_from_file are already Initialized.
        call deall(eAll_larg)
        call deall_grid(Larg_Grid)
        BC_from_file_Initialized=.true.
  end if
  write(15,*) ' End interpolating',BC_from_file(1)%yXMin(10,11)
  
      
end subroutine Interpolate_BC
!**********************************************************************
subroutine ini_BC_from_file(grid)
  type(grid_t), intent(in)        :: grid

  ! local
  type(cboundary)   :: BC
  integer           :: iTx,j,status
  

     call create_cboundary(grid,BC)
     
     
     allocate (BC_from_file(nTx_nPol), STAT=status)
     do j=1,nTx_nPol
       BC_from_file(j)=BC
     end do
     BC_from_file_Initialized=.true.
     

end subroutine ini_BC_from_file     

!**********************************************************************
subroutine copyE0fromFile()
   !   this is just another kluge --- eAll_larg is not available to SetBound,
   !         a routine in ModelOperator3D; copy to E0_from_file (in NestedEM) which is
   !   now have moved the logic out of ModelOperator3D but keeping this temporarily
   !   for historic reasons, until we can revisit [AK]
  integer   :: counter,j,k,status
     allocate (E0_from_file(nTx_nPol), STAT=status)
     counter = 0
     do j=1,eAll_larg%nTx
        !  hard to imagine anything but 1 polarization per transmitter here!
        do k = 1,eAll_larg%solns(1)%nPol
           counter = counter+1
           E0_from_file(counter)=eAll_larg%solns(j)%pol(k)
        end do
     end do

end subroutine copyE0fromFile

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

   integer, intent(in)				:: iTx
   type(modelParam_t),intent(in), target		:: sigma
   type(grid_t), intent(in), target         :: grid
   !  following structures are initialized
   !	solution vector for forward problem
   type(solnVector_t), intent(inout)			:: e0
   !	solution vector for sensitivity
   type(solnVector_t), intent(inout), optional	:: e
   !	forcing for sensitivity
   type(rhsVector_t), intent(inout), optional		:: comb

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
   endif

   ! now calls this all-at-once
   ! call UpdateFreqCond(txDict(iTx)%omega, sigma)
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
	   call copy_modelParam(sigma_temp,sigma) 
!      sigmaNotCurrent = .false.
!   endif

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
   character(10)    :: txType

   initForSens = present(comb)

   if(present(e0)) then
      if(e0%allocated) then
        txType = txDict(e0%tx)%tx_type
      endif
   endif

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
! A. Kelbert, 24 May 2018
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


    ! Initialize the RHS vector; should we always clean it up on input?
    if (.not. b0%allocated) then
      select case (txDict(iTx)%Tx_type)
      case ('MT')
        b0%nonzero_Source = .false.
        b0%sparse_Source = .false.
        b0%nonzero_BC = .true.
	  case ('CSEM')
      	b0%nonzero_Source = .true.
		b0%sparse_Source = .false.
        b0%nonzero_BC = .false.  	  
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

            case ('MT')
                if (BC_FROM_RHS_FILE) then
                    ! in this case, we've read bAll from RHS file already
                    write (*,'(a12,a29,a12,i4,a15,i2)') node_info, 'Setting the BC from RHS file ', &
                        ' for period ',iTx,' & mode # ',iMode
                    BC = bAll%combs(iTx)%b(iMode)%bc

                elseif (BC_FROM_E0_FILE) then
                    ! THIS OPTION SHOULD BE CLEANED UP
                    ! we are going to make a huge assumption here: nPol == 1 always for this case
                    !  and of course transmitters are in same order always
                    write (*,'(a12,a28,a12,i4,a15,i2)') node_info, 'Setting the BC from E0 file ', &
                        ' for period ',iTx,' & mode # ',iMode
                    e0%pol(j) = E0_from_file(iTx)
                    call getBC(e0%pol(j),BC)
                    !   do we now need to set boundary edges of E0 == 0?

                elseif (NESTED_BC) then
                    ! The BC are already computed from a larger grid for all transmitters and modes and stored in BC_from_file.
                    ! Overwrite BC with BC_from_file.
                    ! Note [NM]: Right now we are using the same period layout for both grid.
                    ! This why, it is enough to know the period and mode index to pick up the BC from BC_from_file vector.
                    write (*,'(a12,a35,a12,i4,a15,i2)') node_info, 'Setting the BC from nested E0 file ', &
                        ' for period ',iTx,' & mode # ',iMode
                    BC = BC_from_file((iTx*2)-(2-iMode))

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
             case ('CSEM')
			     ! Naser Meqbel 04.10.2019
			     ! As  b0%nonzero_BC=true for CSEM, NO need to setup the BC for it. 
				 ! However, b0%s is needed to be filled later with Secondary field obtained:
				 ! b0%s=i_omega_mu*(sigma-sigma1d)*Ep
				 ! Ep is the primary E-field obtained from the 1D solution

                 call get_source_for_csem(sigma_temp,e0%grid,iTx,b0%b(j)%s)
				 
				 
            case default
                write(0,*) node_info,'Unknown FWD problem type',trim(txDict(iTx)%Tx_type),'; unable to compute RHS'
        end select
    end do

  end subroutine fwdSetup

   !**********************************************************************
   subroutine fwdSolve(iTx,e0,b0)

   !  driver for 3d forward solver; sets up for transmitter iTx, returns
   !   solution in e0 ; rhs vector (b0) is generated locally--i.e.
   !   boundary conditions are set internally (NOTE: could use transmitter
   !   dictionary to indireclty provide information about boundary
   !    conditions.  Presently we set BC using WS approach.
   !
   ! this *should* works with the SP versions as well

   integer, intent(in)		:: iTx
   type(solnVector_t), intent(inout)	:: e0
   type(rhsVector_t), intent(in)        :: b0

   ! local variables
   real(kind=prec)	:: omega
   integer			:: IER,iMode
   complex(kind=prec)	:: i_omega_mu

   omega = txDict(iTx)%omega
   !  set period, complete setup of 3D EM equation system
   i_omega_mu = cmplx(0.,ISIGN*MU_0*omega,kind=prec)

   !  complete operator intialization, for this frequency
   !  call UpdateFreq(txDict(iTx)%omega)
   !  loop over polarizations
   if (txDict(iTx)%Tx_type=='MT') then
       do iMode = 1,e0%nPol
		! compute boundary conditions for polarization iMode
		!   uses cell conductivity already set by updateCond
		!call setBound(iTx,e0%Pol_index(iMode),e0%pol(imode),b0%bc)
		! NOTE that in the MPI parallelization, e0 may only contain a single mode;
		! mode number is determined by Pol_index, NOT by its order index in e0
		! ... but b0 uses the same fake indexing as e0
		write (*,'(a12,a12,a3,a20,i4,a2,es13.6,a15,i2)') node_info, 'Solving the ','FWD', &
				' problem for period ',iTx,': ',(2*PI)/omega,' secs & mode # ',e0%Pol_index(iMode)
		call FWDsolve3D(b0%b(iMode),omega,e0%pol(iMode))
		write (6,*)node_info,'FINISHED solve, nPol',e0%nPol
   	   enddo
   else
       write(0,*) node_info,'Unknown FWD problem type',trim(txDict(iTx)%Tx_type),'; unable to run fwdSolve'
   endif

   ! update pointer to the transmitter in solnVector
   e0%tx = iTx

   end subroutine fwdSolve

   !**********************************************************************
   subroutine sensSolve(iTx,FWDorADJ,e,comb)
   !   Uses forcing input from comb, which must be set before calling
   !    solves forward or adjoint problem, depending on comb%ADJ
   !  NOTE that this routine DOES NOT call UpdateFreq routine to complete
   !   initialization of solver for a particular frequency.
   !  This final initialization step must (at present) be done by
   !    calling fwdSolve before calling this routine.

   integer, intent(in)          	:: iTx
   character*3, intent(in)		:: FWDorADJ
   type(solnVector_t), intent(inout)		:: e
   type(rhsVector_t), intent(inout)		:: comb

   ! local variables
   integer      			:: IER,iMode
   real(kind=prec) 		:: omega, period

!  zero starting solution, solve for all modes
   call zero_solnVector(e)
   
   if (txDict(iTx)%Tx_type=='MT') then 
   	omega = txDict(iTx)%omega
   	period = txDict(iTx)%period
   	do iMode = 1,e%nPol
      	comb%b(e%Pol_index(iMode))%adj = FWDorADJ
      	write (*,'(a12,a12,a3,a20,i4,a2,es13.6,a15,i2)') node_info,'Solving the ',FWDorADJ, &
		' problem for period ',iTx,': ',(2*PI)/omega,' secs & mode # ',e%Pol_index(iMode)
      	call FWDsolve3d(comb%b(e%Pol_index(iMode)),omega,e%pol(imode))
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
  ! NOTE: at present, this sets up BC_from_file array that is stored in
  !       nestedEM module. No interaction with the bAll variable needed.
  !       Used to require b0 as input but that was merely for initialization.
  !       Anyway, we might want to rewrite this somewhat. [AK 20 May 2019] 

  subroutine Interpolate_BC_from_E_soln(eAll_larg,Larg_Grid,grid)

  type(grid_t)  ,intent(in)                 ::  Larg_Grid
  type(solnVectorMTX_t),intent(in)          ::  eAll_larg
  type(grid_t),intent(in)                   ::  Grid

    ! local variables needed for nesting calculations
  integer                      :: status,iMode,iTx,counter

  Call setup_BC_from_file(Grid,eAll_larg%nTx,eAll_larg%solns(1)%nPol)

  ! For now extract the BC from eAll_larg
    counter=0
    do iTx = 1, eAll_larg%nTx
       do iMode=1, eAll_larg%solns(iTx)%nPol
         counter=counter+1
         Call compute_BC_from_file(Larg_Grid,eAll_larg%solns(iTx)%pol(iMode),Grid,counter)
       end do
    end do

  end subroutine Interpolate_BC_from_E_soln


end module ForwardSolver
