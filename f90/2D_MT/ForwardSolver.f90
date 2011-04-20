module ForwardSolver

!  High level interface/control module used by top level routines
!   for initializing and using the solver.  The key public routines
!   in this module have only generic (abstract data type) arguments
!   and can thus be called from the top-level inversion routines.
!  A similar solver interface will be required for any specific use
!   of the top-level inversion modules
!  All public routines must have the same generic names, parameter lists
!   and abstract functionality as in this example for 2D MT.

use math_constants
use utilities
use datafunc ! inherit SolnSpace soln2d ModelSpace
use dataspace
use fwdtemod
use fwdtmmod
use SolnSpace
use transmitters

implicit none

!    keep data structures used only by
!    routines in this module private
!   rhs data structures for solving forward, sensitivity probs
type(rhsVector_t), save, private		:: b0
!    keep track of which mode was solved for most recently
!    (to minimize reinitialization ... not clear this is needed!)
character*2, save, public		:: currentMode = '  '


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


Contains

   !**********************************************************************
   subroutine initSolver(iTx,sigma,grid,e0,e,comb)
   !   Initializes forward solver for
   !    transmitter iTx: in this instance TE or TM mode solvers
   !    for the appropriate frequency depending on mode
   !   Idea is to call this before calling fwdSolve or sensSolve,
   !     in particular before the first solution for each transmitter
   !     (frequency).  If called for the first time (in a program run,
   !     or after a call to exitSolver), or if the previous call to
   !     this routine initialized for a different data type (TE vs. TM
   !     mode) full initialization (after deallocation/cleanup if required)
   !     is performed.
   !
   !   iTx defines transmitter: for 2D MT, this provides info about
   !       frequency and TE/TM mode; for 3D MT just frequency
   !
   !   This now does all setup (including matrix factorization) for
   !     the appropriate mode/frequency
   !   NOTE: e and comb are optional calling arguments;
   !     both should be present if one is

   integer, intent(in)				:: iTx
   type(modelParam_t),intent(in), target		:: sigma
   type(grid_t), intent(in), target         :: grid
   !  following structures are initialized/created in this routine
   !	solution vector for forward problem
   type(solnVector_t), intent(inout)			:: e0
   !	solution vector for sensitivity
   type(solnVector_t), intent(inout), optional	:: e
   !	forcing for sensitivity
   type(rhsVector_t), intent(inout), optional		:: comb

   !  local variables
   integer					:: IER
   real(kind=prec)			:: period
   character*2          	:: mode
   logical					:: initForSens

   initForSens = present(comb)

   mode = txDict(iTx)%mode
   period = txDict(iTx)%period

   if(currentMode .ne. mode) then
      if(currentMode .ne. '  ') then
         !  not inital solution, but mode has changed from last
         !  sensitivity calculated: will need to reinitialize
         !  solver + rhs/soln arrays ... first deallocate from
         !  previous mode
         call deall_rhsVector(b0)
         call deall_solnVector(e0)
         if(initForSens) then
            call deall_rhsVector(comb)
            call deall_solnVector(e)
         endif
         select case(currentMode)
            case('TE')
               call Fwd2DdeallTE()
            case('TM')
               call Fwd2DdeallTM()
            case default
         end select
      endif
      currentMode = mode

      ! initialize coefficient matrix (frequency indpendent part)
      select case(mode)
         case('TE')
            call FWD2DSetupTE(grid,sigma,IER)
            if(IER.lt.0) then
              call errStop('initializing for TE mode in initSolver')
            endif
         case('TM')
            call FWD2DSetupTM(grid,sigma,IER)
            if(IER.lt.0) then
              call errStop('initializing for TM mode in initSolver')
            endif
         case default
            call errStop('mode must be TE or TM in initSolver')
      end select

      !  allocate for rhs for background, scratch sensitivity solutions
      if (output_level > 3) then
         write(*,*) 'Initializing RHS and background solution for 2D forward solver... iTx=',iTx,' mode=',mode
      endif
      b0%nonzero_source = .false.
      b0%nonzero_bc = .true.
      b0%adj = 'FWD'
      call create_rhsVector(grid,iTx,b0)

      !  allocate for background solution
      call create_solnVector(grid,iTx,e0)

   endif

   !  allocate storage for the sensitivity soln and RHS if run for the first time
   !  or with a new mode; if the mode is new deallocate them first (above)
   if(initForSens) then
     !  allocate for sensitivity solution, RHS
     if (.not. e%allocated) then
        if (output_level > 3) then
            write(*,*) 'Initializing 2D EM soln for sensitivity... iTx=',iTx,' model=',mode
        endif
        call create_solnVector(grid,iTx,e)
     endif
     if (.not. comb%allocated) then
        if (output_level > 3) then
            write(*,*) 'Initializing the RHS for 2D sensitivity... iTx=',iTx,' model=',mode
        endif
        comb%nonzero_source = .true.
        comb%nonzero_bc = .false.
        comb%adj = ''
        call create_rhsVector(grid,iTx,comb)
     endif
   endif

   !   complete initialization of coefficient matrix, then factor
   !   (factored matrix is saved in TE/TM mode modeling module
   !   This needs to be called before solving for a different frequency
   if (output_level > 3) then
      write(*,*) 'Updating frequency for 2D forward solver... iTx=',iTx,' mode=',mode
   endif
   select case (mode)
      case ('TE')
         call UpdateFreqTE(period)
      case ('TM')
         call UpdateFreqTM(period)
      case default
   end select

   end subroutine initSolver

   !**********************************************************************
   subroutine exitSolver(e0,e,comb)
   !   deallocates rhs0, rhs and solver arrays
   type(solnVector_t), intent(inout), optional  :: e0
   type(solnVector_t), intent(inout), optional	::e
   type(rhsVector_t), intent(inout), optional   ::comb

   ! local variables
   logical			:: initForSens

   initForSens = present(comb)

   call deall_rhsVector(b0)
   if(present(e0)) then
      call deall_solnVector(e0)
   endif

   if(initForSens) then
      call deall_rhsVector(comb)
      call deall_solnVector(e)
   endif

   select case(currentMode)
      case('TE')
         call Fwd2DdeallTE()
      case('TM')
         call Fwd2DdeallTM()
      case default
   end select

   currentMode = '  '

   end subroutine exitSolver

   !**********************************************************************
   subroutine fwdSolve(iTx,e0)

   !  driver for 2d forward solver; sets up for transmitter iTx, returns
   !   solution in e0 ; rhs vector (b) is private to this module
   !  NOTE that this routine now does no initialization: this must be
   !   done with a call to initSolver before calling for the first
   !   time FOR EACH DATA TYPE/TRANSMITTER

   integer, intent(in)		:: iTx
   type(solnVector_t), intent(inout)	:: e0

   ! local variables
   real(kind=prec)	:: period,omega
   integer			:: IER
   complex(kind=prec)	:: i_omega_mu
   character*2          	:: mode

   mode = txDict(iTx)%mode
   period = txDict(iTx)%period
   omega = txDict(iTx)%omega
   !  set period, complete setup of TE mode equation system
   i_omega_mu = cmplx(0.,ISIGN*MU_0*omega,kind=prec)

   ! solve forward problem
   select case (mode)
      case ('TE')
         ! compute boundary conditions for 2D TE forward problem
         call SetBoundTE(period,b0%bc)
         !  solve 2D TE equations, return solution in e0
         call Fwd2DsolveTE(b0,e0%vec%v,IER)
         if(IER .LT. 0) then
            call errStop('solving TE mode equation in fwdSolve')
         endif
      case ('TM')
         ! compute boundary conditions for 2D TM forward problem
         call SetBoundTM(period,b0%bc)
         !  solve 2D TM equations, return solution in e0
         call Fwd2DsolveTM(b0,e0%vec%v,IER)
         if(IER .LT. 0) then
            call errStop('solving TM mode equation in fwdSolve')
         endif
      case default
   end select

   if (output_level > 3) then
      write(*,*) 'Completed: 2D forward solver... iTx=',iTx,' mode=',mode
   endif

   ! update pointer to the transmitter in solnVector
   e0%tx = iTx

   end subroutine fwdSolve

   !**********************************************************************
   subroutine sensSolve(iTx,FWDorADJ,e,comb)
   !   Uses forcing input from comb, which must be set before calling
   !    solves forward or adjoint problem, depending on comb%ADJ
   !  NOTE that this routine now does no initialization: this must be
   !   done with a call to initSolver before calling for the first
   !   time FOR EACH DATA TYPE/TRANSMITTER

   integer, intent(in)          	:: iTx
   character*3, intent(in)		    :: FWDorADJ
   type(solnVector_t), intent(inout)		:: e
   type(rhsVector_t), intent(inout)		:: comb

   ! local variables
   integer      :: IER

   comb%adj = FWDorADJ
   if(txDict(iTx)%mode.eq.'TE') then
       call Fwd2DsolveTE(comb,e%vec%v,IER)

       if(IER .LT. 0) then
          call errStop('solving TE mode equation in sensSolve')
       endif
    else
       call Fwd2DsolveTM(comb,e%vec%v,IER)
       if(IER .LT. 0) then
          call errStop('solving TM mode equation in sensSolve')
       endif
   endif

   if (output_level > 3) then
      write(*,*) 'Completed: 2D sensitivity computations... iTx=',iTx,' mode=',txDict(iTx)%mode
   endif

   ! update pointer to the transmitter in solnVector
   e%tx = iTx

   end subroutine sensSolve

end module ForwardSolver
