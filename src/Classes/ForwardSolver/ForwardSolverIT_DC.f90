!
! Please ad some fancy relevant description for this class
!
module ForwardSolverIT_DC
   !
   use ForwardSolverIT
   use DivergenceCorrection
   use ModelOperator
   use Solver_QMR
   !
   ! ModEM classic -- iterative solution with QMR or BiCG, divergence
   ! correction (DC) in usual way.   A slight variant/extension would be
   ! required for secondary field formulation.

   ! NEED TO THINK ABOUT MANAGING DEFAULTS/USER CONTROL
   ! default values for solver control 
   !   NOTE:   QMR and BiCG will have solver default set appropriate for 
   !    a general iterative solution WITHOUT DC -- these defaults will be
   !    changed for the DC case (to iter_per_div_cor)  
   !   NOT SURE we want to have these solver dependent parameters set explicitly
   !    here, but for now OK
   integer, parameter :: iter_per_div_corDefQMR = 40
   !
   integer, parameter :: iter_per_div_corDefBCG = 80
   ! maximum number of divergence correction calls allowed
   integer, parameter :: max_div_corDef = 20
   ! for DC object default max_iter_total = iter_per_div_cor*max_div_cor
   !  compute from other defaults, depending on solver
   ! maximum number of PCG iterations for divergence correction
   ! this will be default max_iter for Solver_PCG object
   integer, parameter :: max_iterDivCorDef = 100
   !
   ! misfit tolerance for convergence of divergence correction solver
   real( kind=prec ), parameter :: tolDivCorDef = 1E-5
   !   this default is used to set overall tolerance for ForwardSolver object
   real( kind=prec ), parameter :: tolCurlCurlDef = 1E-7
   !
   type, extends( ForwardSolverIT_t ), public :: ForwardSolverIT_DC_t
      !
      class( DivergenceCorrection_t ), allocatable :: divergence_correction ! pointer to divergence correction
      !
      !    array of relative residuals for PCG convergence, all DC steps
      real( kind=prec ), allocatable, dimension(:,:) :: DivCorRelErr
      !    current divergence before and after each divergence correction step
      real( kind=prec ), allocatable, dimension(:,:) :: divJ
      !
      integer :: nDivCor = 0
      !
      integer :: max_div_cor = 0
      integer :: max_iterDivCor = 0
      real( kind=prec ) :: tolDivCor = 0.0
      !
   contains
      !
      final :: ForwardSolverIT_DC_dtor
      !
      !   procedures with abstract interfaces
      procedure, public :: setPeriod => setPeriodForwardSolverIT_DC
      procedure, public :: setCond => setCondForwardSolverIT_DC
      procedure, public :: setIterControl => setIterControlForwardSolverIT_DC
      procedure, public :: initDiagnostics => initDiagnosticsForwardSolverIT_DC
      procedure, public :: zeroDiagnostics => zeroDiagnosticsForwardSolverIT_DC
      procedure, public :: getESolution => getESolutionForwardSolverIT_DC
      !
      !   procedure unique to DC
      procedure, public :: setIterDefaultsDC

   end type ForwardSolverIT_DC_t
   !
   interface ForwardSolverIT_DC_t
      module procedure ForwardSolverIT_DC_ctor
   end interface ForwardSolverIT_DC_t
   !
   contains
   !
    function ForwardSolverIT_DC_ctor( model_operator, solver_type ) result( self )
       implicit none
       class( ModelOperator_t ), target, intent( in ) :: model_operator
       character(*), intent(in)    :: solver_type
       type( ForwardSolverIT_DC_t ) :: self
 
       integer :: maxIter, maxItTotal
       real(kind=prec)  :: tol
       !
       !write(*,*) "Constructor ForwardSolverIT_DC_t"
       !
       ! DivergenceCorrection has only one type (might change components,
       !    but basic scheme implemented is not going to change)
       self%divergence_correction = DivergenceCorrection_t( model_operator )
	   !
       !   solver will soon have options
       select case( solver_type )
          case( QMR )
             self%solver = Solver_QMR_t( model_operator )
             maxIter = iter_per_div_corDefQMR
          case( BiCG )
             maxIter = iter_per_div_corDefBCG
             stop "Not yet coded for Bi-Conjugate Gradients"
       end select
       !
       !    set solver iteration control parameters using defaults
       call self%solver%setParameters(maxIter,tolCurlCurlDef)
 
       !    set remaining default iteration control for DC
       call self%setIterDefaultsDC()
       !
       !    initialize Fwd operator iteration control, diagonstic arrays
       !     using defaults from solver
       maxItTotal = self%max_div_cor * self%solver%max_iter
       tol = self%solver%tolerance
       call self%setIterControl( maxItTotal,tol )
       !
       call self%initDiagnostics()
       !
    end function ForwardSolverIT_DC_ctor
    !
    ! Destructor
    subroutine ForwardSolverIT_DC_dtor( self )
       implicit none
       !
       type( ForwardSolverIT_DC_t ), intent( in out ) :: self
       !
       !write(*,*) "Destructor ForwardSolverIT_DC_t"
       !
    end subroutine ForwardSolverIT_DC_dtor
    !
    !   perhaps this can be in base class?   I guess this is always the same
    subroutine setPeriodForwardSolverIT_DC( self, period )
       implicit none
       !
       class( ForwardSolverIT_DC_t ), intent( inout ) :: self
       real( kind=prec ), intent( in )                :: period
       !
       real( kind=prec ) :: rel_diff
       !
       rel_diff = ( self%period - period ) / period
 
       self%period = period
       !
       self%omega = 2.0 * PI / period
       !
       !   set frequency in solver object
       self%solver%omega = self%omega
       !     set preconditoner (depends on frequency in general)
       !   but only if there is a large enough change in period
       if( rel_diff .gt. TOL4 ) then
          call self%solver%preconditioner%SetPreconditioner( self%omega )
       endif
       !
    end subroutine setPeriodForwardSolverIT_DC
    !
    !    Sets Condctivity
    !
    subroutine setCondForwardSolverIT_DC( self, modPar )
       implicit none
       !
       class( ForwardSolverIT_DC_t ), intent( inout ) :: self
       class( ModelParameter_t ), intent( in )     :: modPar
       !
       !   set conductivity in model_operator object
       call self%solver%model_operator%setCond( modPar )
       !   set arrays in model_operator needed for divergence correction
       call self%divergence_correction%SetCond()
       !     set preconditoner for (PCG) solver (depends only on conductivity
       !       in this case)
       call self%solver%preconditioner%SetPreconditioner( self%omega )
       !
     end subroutine setCondForwardSolverIT_DC
    !
    !**********
    ! 
    !   This is a base class procedure, intended to make it easy to change
    !   overall solution tolerance, and overall maximum number of iterations
    !   For DC case we use these to adjust some other parameters; for full
    !   control we need a routine that can set all DC iteration control parameters
    subroutine setIterControlForwardSolverIT_DC( self, maxit, tol )
       implicit none
	   !
       class( ForwardSolverIT_DC_t ), intent( inout ) :: self
       integer, intent(in)                         :: maxit
       real(kind=prec), intent(in)                 :: tol
       !
       integer  :: itPerDC
       !    
       !    tolerance is property of base class -- overall tolerance for
       !    convergence (also used to set tolerance in solver--these are always same`
       self%tolerance = tol
       !   self%solver%max_iter is number of iterations per DC -- leave
       !    this as is, and adjust max_div_cor
       itPerDC = self%solver%max_iter
       self%max_div_cor = maxit/itPerDC
       self%max_iter_total = itPerDC*self%max_div_cor
       !
       !   reset solver iteration control 
       call self%solver%setParameters(itPerDC,tol)

    end subroutine setIterControlForwardSolverIT_DC
    !
    !********
    !
    subroutine setIterDefaultsDC( self )
       implicit none
       class( ForwardSolverIT_DC_t ), intent( inout ) :: self
       ! this just sets iteration control parameters specific to DC to default
       ! values -- should be good enough to get us started!
       !  Rationale:   this keeps all iteration control parameters needed
       !    for curl-curl and DC solvers.   Defaults are parameters, but 
       !    we need to allow all of these to be changed by users, so actual
       !    values are variables stored in ForwardSolver object
       !    need a way to set these 4 parameters
       !    to what the (experienced) user asks for ...
       !
       self%max_div_cor = max_div_corDef
       !   this is max_iter for PCG in DC
       self%max_IterDivCor = max_IterDivCorDef
       !   this is tolerance for PCG in DC
       self%tolDivCor = tolDivCorDef

    end subroutine setIterDefaultsDC

    !
    !***********************************************************
    !
    subroutine initDiagnosticsForwardSolverIT_DC( self )
        implicit none
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
        ! this allocates arrays for storage of solver diagnostics

        self%n_iter_actual = 0
        self%relResFinal = R_ZERO
        !
        !   vector of all relative residuals (concatenated over all divergence
        !    correction steps)
        if(allocated(self%relResVec)) deallocate(self%relResVec)
        allocate(self%relResVec(self%max_iter_total))
        !
        ! Intermediate solution divergence before/after  -- one for each DC step
        if( allocated( self%divJ ) ) deallocate( self%divJ )
        allocate( self%divJ( 2, self%max_div_cor ) )
        !
        !   convergence of PCG solver for all divergence corrections
        !  (probably not needed?  Do we ever really look at this?  Naser?)
        if( allocated( self%DivCorRelErr ) ) deallocate( self%DivCorRelErr )
        allocate( self%DivCorRelErr( self%max_IterDivCor, self%max_div_cor ) )
        !
     end subroutine initDiagnosticsForwardSolverIT_DC
     !
     !*********
     !
     subroutine zeroDiagnosticsForwardSolverIT_DC(self)
        implicit none
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self

          self%relResVec = R_ZERO
          self%divJ = R_ZERO
          self%DivCorRelErr = R_ZERO
          call self%solver%zeroDiagnostics()

     end subroutine zeroDiagnosticsForwardSolverIT_DC
     !
     !*********
     !
     subroutine getESolutionForwardSolverIT_DC( self, source, e_solution )
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
        class( Source_t ), intent( inout )                :: source
        !integer, intent( in )                            :: polarization
        class( cVector_t ), intent( inout )               :: e_solution
        ! local variables
        class( cVector_t ), allocatable :: b    ! copy of RHS--do we really need?
        class( cVector_t ), allocatable :: temp
        class( cScalar_t ), allocatable :: phi0
        integer :: iter
        !
        ! zero solver diagnostic arrays
        call self%solver%zeroDiagnostics()
        !
        ! not sure about allocation here -- solution will exist
        ! (and might be allocated) in calling routine, but b is local
        ! Note that rhs will be a TVector of same type as solution
        allocate( b, source = source%rhs )
        !
        ! set up source term for divergence correction equations
        if( source%non_zero_source ) then
          !
          ! make a copy of TScalar using model_operator template
          allocate( phi0, source = self%solver%model_operator%createScalar() )
          !
          call self%divergence_correction%rhsDivCor( self%omega, source, phi0 )
          !
        endif
        !
        allocate( temp, source = e_solution )    ! copy of solution for input to DC
        !
        loop: do while ( ( .not.self%solver%converged ).and.( .not.self%solver%failed ) )
           !
           ! TEMPORARY SELECT CASE
           !
           select type( solver => self%solver )
              class is( Solver_QMR_t )
                 call solver%solve( b, e_solution )
              class default
                 write(*, *) "ERROR:ForwardSolverIT_DC::getESolutionForwardSolverIT_DC:"
                 STOP        "         Unknow solver type."
           end select
           !
           ! I am just copying this -- while we work on this should
           !  reconsider implementation
           ! solver%converged when the relative error is less than tolerance
           self%solver%converged = self%solver%n_iter .lt. self%solver%max_iter
           !
           ! there are two ways of failing: 1) QMR did not work or
           !     2) total number of divergence corrections exceeded
           self%solver%failed = self%solver%failed .or. self%failed
           !
           ! update solver diagnostics 
           do iter = 1, self%solver%n_iter
              ! why are we using an explicit loop here?
              self%relResVec( self%n_iter_actual + iter ) = self%solver%relErr(iter)
           enddo
           self%n_iter_actual = self%n_iter_actual + self%solver%n_iter
           self%nDivCor = self%nDivCor+1
           !
           if( self%nDivCor < self%max_div_cor ) then
              !  copy current e_solution into temp (discuss if this is this needed?)
              temp = e_solution
              if( source%non_zero_source ) then
                 !
                 call self%divergence_correction%DivCorr( temp, e_solution, phi0 )
                 !
              else
                 !
                 call self%divergence_correction%DivCorr( temp, e_solution )
                 !
           endif
       	  !
          else
             ! max number of divergence corrections exceeded; convergence solver%failed
             self%solver%failed = .true.
          endif
	   !
       enddo loop
       !
       self%relResFinal = self%relResVec(self%n_iter_actual)
       !
       ! finish up solution--I am omitting boundary values for adjt case --
       ! we never used boundary of adjoint--sensitivity to boundary data,
       ! just use interior part of adjoint soln to compute sensitivity to model
       ! parameters

       ! note that here I am assuming things like mult and add are subroutines.
       ! we need to sort out conventions! In Solver_QMR I assumed functions,
       ! but I suspect we will be better off just using subroutines in terms
       ! of efficiency
       if( source%adjt ) then
          select type( modOp => self%solver%model_operator )
             class is ( ModelOperator_MF_t )
             !
             e_solution = e_solution * modOp%Metric%Vedge
             !
             class default
                write(*, *) "ERROR:ForwardSolverIT_DC_t::getESolutionForwardSolverIT_DC:"
                STOP        "model_operator type unknow"
          end select
         ! just leave bdry values set to 0
       else
          !
          e_solution = e_solution + source%bdry
          !
       endif
       !
       select type( e_solution )
          class is( cVector3D_SG_t )
       	     write( *, * ) "         ", e_solution%nx, e_solution%ny, e_solution%nz, e_solution%gridType
          class default
       	     stop "Unclassified ForwardSolverIT_DC e_solution"
       end select
       !
       ! deallocate local objects
       if( allocated( temp ) ) deallocate(temp)
       if( allocated( b ) )    deallocate(b)
       if( allocated( phi0 ) ) deallocate(phi0)

   end subroutine getESolutionForwardSolverIT_DC
        
end Module ForwardSolverIT_DC
 
