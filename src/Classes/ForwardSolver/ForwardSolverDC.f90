module ForwardSolverDC
   !
   use ForwardSolver
   use DivergenceCorrection
   use ModelOperator
   use Solver_QMR
   !
   !   ModEM classic -- iterative solution with QMR or BiCG, divergence
   !   correction (DC) in usual way.   A slight variant/extension would be
   !   required for secondary field formulation.

   !   default values for solver control 
   !    since reasonable defaults depend on the details of forward solver
   !    details, set the defaults for parameters declared in solver object
   !    here also
   !    this will be default max_iter for QMR/BiCG solver
   integer, parameter         :: iter_per_div_corDefQMR = 40
   integer, parameter         :: iter_per_div_corDefBCG = 80
   ! maximum number of divergence correction calls allowed
   integer, parameter         :: max_div_corDef = 20
   !   for DC object default max_iter_total = iter_per_div_cor*max_div_cor
   !     compute from other defaults, depending on solver
   ! maximum number of PCG iterations for divergence correction
   !    this will be default max_iter for Solver_PCG object
   integer, parameter         :: max_iterDivCorDef = 100
   !
   ! misfit tolerance for curl-curl (QMR or BCG)
   real(kind=prec), parameter ::   tolEMDef = 1E-7
   ! misfit tolerance for convergence of divergence correction solver
   real(kind=prec), parameter ::   tolDivCorDef = 1E-5
   !
   type, extends( ForwardSolver_t ), public :: ForwardSolverDC_t
      !
      class( DivergenceCorrection_t ), pointer :: divergence_correction !  pointer to divergence correction
      !
      real( kind=prec )                              :: omega = 0.0
      real( kind=prec ), allocatable, dimension(:)   :: EMrelErr
      real( kind=prec ), allocatable, dimension(:,:) :: divJ, DivCorRelErr
      !
      integer :: nDivCor = 0
      !   next two are not independent of max_iter_total
      integer :: max_iter_total = 0
      integer :: max_div_cor = 0
      integer :: iter_per_div_cor = 0
      !
   contains
      !
      final :: ForwardSolverDC_dtor
      !
      procedure, public :: setPeriod => setPeriodForwardSolverDC
      procedure, public :: setFrequency => setFrequencyForwardSolverDC
      procedure, public :: setIterDefaults
      procedure, public :: createDiagnosticArrays
      procedure, public :: initDiagnosticArrays
      procedure, public :: getESolution => getESolutionForwardSolverDC
      !   set routines for iteration control parameters (specific to DC)
      !procedure, public :: set_max_div_cor
      !procedure, public :: set_iter_per_div_cor
      !   get routines for diagonstics
      !procedure, public :: get_nDivCor
      !procedure, public :: get_divJ
      !procedure, public :: get_DivCorRelErr
      !
   end type ForwardSolverDC_t
   !
   interface ForwardSolverDC_t
      module procedure ForwardSolverDC_ctor
   end interface ForwardSolverDC_t
   !
   contains
   !
   function ForwardSolverDC_ctor( model_operator, divergence_correction ) result( self )
      !
	  class( ModelOperator_t ), target, intent( in ) :: model_operator
	  class( DivergenceCorrection_t ), target, intent( in ) :: divergence_correction
	  !
      type( ForwardSolverDC_t ) :: self
      !
      write(*,*) "Constructor ForwardSolverDC_t"
      !
      call self%init()
      !
	  self%solver = Solver_QMR_t( model_operator )
	  !
      self%divergence_correction => divergence_correction
      !
      call self%setIterDefaults()
      !
      call self%createDiagnosticArrays()
      !
   end function ForwardSolverDC_ctor
   !
   ! Destructor
   subroutine ForwardSolverDC_dtor( self )
      implicit none
      !
      type( ForwardSolverDC_t ), intent( in out ) :: self
      !
      !write(*,*) "Destructor ForwardSolverDC_t"
      !
      call self%dealloc()
      !
   end subroutine ForwardSolverDC_dtor
   !
   subroutine setPeriodForwardSolverDC( self, period )
      implicit none
      !
      class( ForwardSolverDC_t ), intent( inout ) :: self
      real( kind=prec ), intent( in )             :: period
      !
      self%period = period
      !
      self%omega = 2.0 * PI / period
      !
      call self%setFrequency( self%omega )
   !
   end subroutine setPeriodForwardSolverDC
   !
    !
    !   creator and destructor ... still need these
   subroutine setIterDefaults( self )
      class( ForwardSolverDC_t ), intent( inout ) :: self
      !   this just sets iteration control parameters to default
      !    values -- should be good enough to get us started!
      !    Note that some of the parameters are set in solver object 
      !
      self%max_div_cor = max_div_corDef
      self%max_iter_total = self%max_div_cor * self%solver%max_iter
     !
      !select type( solver => self%solver )
      !class is (Solver_QMR_t)
        self%iter_per_div_cor = iter_per_div_corDefQMR
     !class is (Solver_BiCG_t)
        !self%iter_per_div_cor = iter_per_div_corDefBCG
      !end select

      !   these are parameters in the solver objects, for curl-curl
      !      and for DC
      !self%solver%setParameters(self%IterDivCor,tolEMDef)
      call self%solver%setParameters( max_iterDivCorDef, tolDivCorDef )

    end subroutine setIterDefaults
      !***********************************************************
    subroutine createDiagnosticArrays( self )
        !
        class( ForwardSolverDC_t ), intent( inout ) :: self
        !   this allocates arrays for storage of solver diagnostics

        ! Forward object: EMrelErr, divJ, DivCorRelErr
        if( allocated( self%divJ ) ) deallocate( self%divJ )
        !
        allocate( self%divJ( 2, self%max_div_cor ) )
        !
        if( allocated( self%EMrelErr) ) deallocate( self%EMrelErr )
        !
        allocate( self%EMrelErr( self%max_iter_total ) )
        !
        if( allocated( self%DivCorRelErr ) ) deallocate( self%DivCorRelErr )
        !
        allocate( self%DivCorRelErr( self%solver%max_iter, self%max_div_cor ) )
        !
        !  Solver objects
        if( allocated( self%solver%relErr ) ) deallocate( self%solver%relErr )
        !
        allocate( self%solver%relErr( self%solver%max_iter ) )

     end subroutine createDiagnosticArrays
     !*****************************************************
     subroutine initDiagnosticArrays( self )
        !
        implicit none
        !
        class( ForwardSolverDC_t ), intent( inout ) :: self
        !   this zero"s diagnostic arrays for storage of solver diagnostics
        self%n_iter_total = 0
        self%nDivCor = 0
        self%EMrelErr = R_ZERO
        self%divJ = R_ZERO
        self%DivCorRelErr = R_ZERO
        self%solver%failed = .false.
        self%solver%converged = .false.
        !
     end subroutine initDiagnosticArrays    
     !*****************************************************
     subroutine setFrequencyForwardSolverDC( self, omega )
        !   this is not specific to DC solver -- can we implement
        !     in abstract class?
        class( ForwardSolverDC_t ), intent( inout ) :: self
        real( kind=prec ), intent( in )             :: omega
        if(abs(self%omega-omega) .lt. TOL8) then
          !   omega is close enough to input -- no need to reset
          !    freqeuncy dependent properties
          return
        end if
        ! otherwise need to update some things ...
        self%omega = omega
        self%solver%omega = omega
        call self%solver%preconditioner%SetPreconditioner( omega )
		!
     end subroutine setFrequencyForwardSolverDC
     !*****************************************************
     function getESolutionForwardSolverDC( self, source, polarization ) result( e_solution )
        implicit none
        !
        class( ForwardSolverDC_t ), intent( inout ) :: self
        class( Source_t ), intent( in )             :: source
        integer, intent( in )                       :: polarization
        !
        class( cVector_t ), allocatable :: e_solution
        !   local variables
        class( cVector_t ), allocatable :: b    !  copy of RHS--do we really need?
        class( cVector_t ), allocatable :: temp
        class( cScalar_t ), allocatable :: phi0
        integer :: iter
        !
        !   initialize diagnostics -- am assuming that setting of solver parameters
        !     is done in a set up step (once in the run) outside this object
        call self%initDiagnosticArrays()
        !
        call self%solver%zeroDiagnostics
        !
        !   initialize solution
        allocate( e_solution, source = source%e0 )
        !
        !   not sure about allocation here -- solution will exist
        !    (and might be allocated) in calling routine, but b is local
        !   Note that rhs will be a TVector of same type as solution
        allocate( b, source = source%rhs )
        !
        !   set up source term for divergence correction equations
        if( source%non_zero_source ) then
           !phi0 =  self%model_operator%p   !   make a copy of TScalar using
                                   !   model_operator template
            call self%divergence_correction%rhsDivCor( self%omega, source, phi0 )
           !allocate( phi0, source = self%divergence_correction%rhsDivCor( self%omega, source ) )
        endif
        !
        allocate( temp, source = e_solution )    !  copy of solution for input to DC
        !
        loop: do while ((.not.self%solver%converged).and.(.not.self%solver%failed))
           !
           ! TEMPORARY SELECT CASE
           !
           select type( solver => self%solver )
              class is( Solver_QMR_t )
                 call solver%solve( b, e_solution )
              class default
                 write(*, *) "ERROR:ForwardSolverDC::getESolutionForwardSolverDC:"
                 STOP "         Unknow solver type."
           end select
           !
           !   I am just copying this -- while we work on this should
           !     reconsider implementation
           ! solver%converged when the relative error is less than tolerance
           self%solver%converged = self%solver%n_iter .lt. self%solver%max_iter
           !
           ! there are two ways of failing: 1) QMR did not work or
           !        2) total number of divergence corrections exceeded
           self%solver%failed = self%solver%failed .or. self%failed
           !
           !  update solver diagnostics 
           do iter = 1, self%solver%n_iter
              ! why are we using an explicit loop here?
              self%EMrelErr( self%n_iter_total + iter ) = self%solver%relErr(iter)
           enddo
           !
           self%n_iter_total = self%n_iter_total + self%solver%n_iter
           !
           self%nDivCor = self%nDivCor+1
		   !
           if( self%nDivCor < self%max_div_cor ) then
              ! do divergence correction
              e_solution = temp !    assuming temp is already allocated,
                                ! don"t want to reallocate!
			  !
              if( source%non_zero_source ) then
                 !
			     call self%divergence_correction%DivCorr( e_solution, e_solution, phi0 )
                 !
              else
			     !
                 call self%divergence_correction%DivCorr( e_solution, e_solution )
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
       !   finish up solution--I am omitting boundary values for adjt case --
       !    we never used boundary of adjoint--sensitivity to boundary data,
       !    just use interior part of adjoint soln to compute sensitivity to model
       !    parameters

       ! note that here I am assuming things like mult and add are subroutines.
       !   we need to sort out conventions!   In Solver_QMR I assumed functions,
       !    but I suspect we will be better off just using subroutines in terms
       !    of efficiency
       if( source%adjt ) then
          select type( modOp => self%solver%model_operator )
             class is ( ModelOperator_MF_t )
             !
             e_solution = e_solution * modOp%Metric%Vedge
             !
             class default
                write(*, *) "ERROR:ForwardSolverDC_t::getESolutionForwardSolverDC:"
                STOP        "model_operator type unknow"
          end select
         !   just leave bdry values set to 0
       else
          !
          e_solution = e_solution + source%bdry
          !
       endif

       !  deallocate local objects
       if( allocated( temp ) ) deallocate(temp)
       if( allocated( b ) )    deallocate(b)
       if( allocated( phi0 ) ) deallocate(phi0)

   end function getESolutionForwardSolverDC
        
end Module ForwardSolverDC
 
