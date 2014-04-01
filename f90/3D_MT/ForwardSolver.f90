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

!    keep data structures used only by
!    routines in this module private
!   rhs data structures for solving forward, sensitivity probs
type(RHS_t), save, private			:: b0
! Used for interpolating the BC from a larger grid.
  type(grid_t),save,public               	 ::	 Larg_Grid   
  type(solnVectorMTX_t),save,public          ::  eAll_larg
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

logical, save, private		:: modelDataInitialized = .false.
logical, save, private		:: BC_from_file_Initialized = .false.
!  logical, save, private		:: sigmaNotCurrent = .true.

Contains
   !**********************************************************************
subroutine Interpolate_BC(grid)
  type(grid_t), intent(in)        :: grid
  integer                                  :: iTx,ix,iy,iz

  


	 b0%nonzero_Source = .false.
     b0%nonzero_bc = .true.
	 
	 
     b0%adj = 'FWD'
     b0%sparse_Source = .false.
     iTx=1
     call create_RHS(grid,iTx,b0)
	!In case of interpolating the BC from eAll_larg   
	! If eAll_ larg solution is already allocated, then use that to interpolate the BC from it
     
	  if (eAll_larg%allocated) then
           write(15,*) ' Start interploating',grid%nx,grid%ny,grid%nz,b0%bc%nx,b0%bc%ny,b0%bc%nz
		call Interpolate_BC_from_E_soln (eAll_larg,Larg_Grid,Grid,b0)
        !Once we are ready from eAll_larg, deallocate it, and keep track, that BC_from_file are already Initialized.
        call deall(eAll_larg)
        call deall_grid(Larg_Grid)
         BC_from_file_Initialized=.true.
      end if
        write(15,*) ' End interploating',BC_from_file(1)%yXMin(10,11)
  


        
      
end subroutine Interpolate_BC
!**********************************************************************
subroutine ini_BC_from_file(grid)
  type(grid_t), intent(in)        :: grid    
  integer                                  :: iTx,j,status
  
	 b0%nonzero_Source = .false.
     b0%nonzero_bc = .true.
     b0%adj = 'FWD'
     b0%sparse_Source = .false.
     iTx=1
     call create_RHS(grid,iTx,b0)
     
     
     allocate (BC_from_file(nTx_nPol), STAT=status)
     do j=1,nTx_nPol
       BC_from_file(j)=b0%bc
     end do
     BC_from_file_Initialized=.true.
     

end subroutine ini_BC_from_file     
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

   !  allocate for scratch rhsVector structure for background, sensitivity
   b0%nonzero_Source = .false.
   b0%nonzero_bc = .true.
   b0%adj = 'FWD'
   b0%sparse_Source = .false.
   call create_RHS(grid,iTx,b0)

   !In case of interpolating the BC from eAll_larg   
   ! If eAll_ larg solution is already allocated, then use that to interpolate the BC from it
   !if (eAll_larg%allocated) then
   !  call Interpolate_BC_from_E_soln (eAll_larg,Larg_Grid,grid,b0)
   !  !Once we are ready from eAll_larg, deallocate it, and keep track, that BC_from_file are already Initialized.
   !  call deall(eAll_larg)
   !  call deall_grid(Larg_Grid)
   !  BC_from_file_Initialized=.true.
   !end if
   
   if (BC_from_file_Initialized) then
     b0%bc%read_E_from_file=.true.
   end if

   !  allocate for background solution
   call create_solnVector(grid,iTx,e0)

   if(initForSens) then
      !  allocate for sensitivity solution, RHS
      call create_solnVector(grid,iTx,e)
      call create_rhsVector(grid,iTx,comb)
      do k = 1,comb%nPol
        comb%b(k)%nonzero_source = .true.
        comb%b(k)%nonzero_bc = .false.
        !  assuming here that we don't use sparse storage ... we could!
        comb%b(k)%sparse_Source = .false.
        comb%b(k)%adj = ''
        !  using all this information, reallocate storage for each polarization
        call create_RHS(grid,iTx,comb%b(k))
      enddo
   endif

   if(.NOT.modelDataInitialized) then
   !   Initialize modelData, setup model operators
      call ModelDataInit(grid)
      call ModelOperatorSetup()
      modelDataInitialized = .true.
   endif

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

   call deall_RHS(b0)
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
   subroutine fwdSolve(iTx,e0)

   !  driver for 3d forward solver; sets up for transmitter iTx, returns
   !   solution in e0 ; rhs vector (b0) is generated locally--i.e.
   !   boundary conditions are set internally (NOTE: could use transmitter
   !   dictionary to indireclty provide information about boundary
   !    conditions.  Presently we set BC using WS approach.
   !  NOTE that this routine calls UpdateFreq routine to complete
   !   initialization of solver for a particular frequency.

   integer, intent(in)		:: iTx
   type(solnVector_t), intent(inout)	:: e0

   ! local variables
   real(kind=prec)	:: period, omega
   integer			:: IER,iMode
   complex(kind=prec)	:: i_omega_mu

   omega = txDict(iTx)%omega
   period = txDict(iTx)%period
   !  set period, complete setup of 3D EM equation system
   i_omega_mu = cmplx(0.,ISIGN*MU_0*omega,kind=prec)

   !  complete operator intialization, for this frequency
   !  call UpdateFreq(txDict(iTx)%omega)
   !  loop over polarizations
   do iMode = 1,e0%nPol
      ! compute boundary conditions for polarization iMode
      !   uses cell conductivity already set by updateCond
      call SetBound(e0%Pol_index(iMode),period,e0%pol(imode),b0%bc,iTx)
      write (*,'(a12,a12,a3,a20,i4,a2,es12.6,a15,i2)') node_info, 'Solving the ','FWD', &
				' problem for period ',iTx,': ',(2*PI)/omega,' secs & mode # ',e0%Pol_index(iMode)
      call FWDsolve3D(b0,omega,e0%pol(imode))
   enddo

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

   omega = txDict(iTx)%omega
   period = txDict(iTx)%period
   !  zero starting solution, solve for all modes
   call zero_solnVector(e)
   do iMode = 1,e%nPol
      comb%b(e%Pol_index(iMode))%adj = FWDorADJ
      write (*,'(a12,a12,a3,a20,i4,a2,es12.6,a15,i2)') node_info,'Solving the ',FWDorADJ, &
				' problem for period ',iTx,': ',(2*PI)/omega,' secs & mode # ',e%Pol_index(iMode)
      call FWDsolve3d(comb%b(e%Pol_index(iMode)),omega,e%pol(imode))
   enddo

   ! update pointer to the transmitter in solnVector
   e%tx = iTx

   end subroutine sensSolve


  !**********************************************************************
  ! uses nestedEM module to extract the boundary conditions directly from
  ! a full EMsolnMTX vector on a larger (and coarser) grid

  subroutine Interpolate_BC_from_E_soln(eAll_larg,Larg_Grid,grid,b0)

  type(grid_t)  ,intent(in)                 ::  Larg_Grid
  type(solnVectorMTX_t),intent(in)          ::  eAll_larg
  type(grid_t),intent(in)                   ::  Grid
  type(RHS_t),intent(in)                    ::  b0

    ! local variables needed for nesting calculations
  integer                      :: status,iMode,iTx,counter

  Call setup_BC_from_file(Grid,b0%bc,eAll_larg%nTx,eAll_larg%solns(1)%nPol)

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
