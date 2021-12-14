!
! Please ad some fancy relevant description for this class
!
module ForwardSolverIT
   !
   use ForwardSolver
   use DivergenceCorrection
   use ModelOperator
   use Solver_QMR
   !
   ! misfit tolerance for curl-curl (QMR or BCG)
   real( kind=prec ), parameter :: tolEMDef = 1E-7
   !
   type, extends( ForwardSolver_t ), public :: ForwardSolverIT_t
      !
      real( kind=prec )                              :: omega = 0.0
      real( kind=prec ), allocatable, dimension(:)   :: EMrelErr
      real( kind=prec ), allocatable, dimension(:,:) :: divJ
      !
      integer :: max_iter_total = 0
      !
   contains
      !
      final :: ForwardSolverIT_dtor
      !
      procedure, public :: setPeriod              => setPeriodForwardSolverIT
      procedure, public :: setFrequency           => setFrequencyForwardSolverIT
      procedure, public :: setIterDefaults        => setIterDefaultsForwardSolverIT
      procedure, public :: createDiagnosticArrays => createDiagnosticArraysForwardSolverIT
      procedure, public :: initDiagnosticArrays   => initDiagnosticArraysForwardSolverIT
      procedure, public :: getESolution           => getESolutionForwardSolverIT
      !
   end type ForwardSolverIT_t
   !
   interface ForwardSolverIT_t
      module procedure ForwardSolverIT_ctor
   end interface ForwardSolverIT_t
   !
contains
   !
   ! Derived class Constructor:
   !    Calls base init()
   !    Sends ModelOperator to solver
   !    
   !    setIterDefaults ???? is really necessary ?
   !
   function ForwardSolverIT_ctor( model_operator ) result( self )
      implicit none
      !
      class( ModelOperator_t ), intent( in ) :: model_operator
      !
      type( ForwardSolverIT_t )              :: self
      !
      !write(*,*) "Constructor ForwardSolverIT_t"
      !
      call self%init()
      !
      self%solver = Solver_QMR_t( model_operator )
      !
      call self%setIterDefaults()
      !
      call self%createDiagnosticArrays()
      !
   end function ForwardSolverIT_ctor
   !
   ! Derived class Destructor:
   !    Triggered when this object is dereferenced in memory.
   !    Deallocates memory from properties of this class.
   !    Calls base init()
   subroutine ForwardSolverIT_dtor( self )
      implicit none
      !
      type( ForwardSolverIT_t ), intent( in out ) :: self
      !
      !write(*,*) "Destructor ForwardSolverIT_t"
      !
      call self%dealloc()
      !
   end subroutine ForwardSolverIT_dtor
   !
   ! Base Interface setPeriod:
   !    Sets Period and Omega
   !    Calls base setFrequency()
   subroutine setPeriodForwardSolverIT( self, period )
      implicit none
      !
      class( ForwardSolverIT_t ), intent( inout ) :: self
      real( kind=prec ), intent( in )             :: period
      !
      self%period = period
      !
      self%omega = 2.0 * PI / period
      !
      call self%setFrequency( self%omega )
      !
   end subroutine setPeriodForwardSolverIT
   !
   ! Base Interface setPeriod:
   !    Sets Period and Omega
   !    Calls base setFrequency()
   subroutine setIterDefaultsForwardSolverIT( self )
      implicit none
      !
      class( ForwardSolverIT_t ), intent( inout ) :: self
      !
      self%max_iter_total = self%solver%max_iter
      !
      ! Need to set Solver again without DC ????
      !call self%solver%setParameters( max_iterDivCorDef, tolDivCorDef )

    end subroutine setIterDefaultsForwardSolverIT
    !
    ! ForwardSolverIT createDiagnosticArrays:
    !    Allocates the arrays and their dimensions used for diagnostic analysis.
    subroutine createDiagnosticArraysForwardSolverIT( self )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        !
        ! Forward object: EMrelErr, divJ
        if( allocated( self%divJ ) ) deallocate( self%divJ )
        allocate( self%divJ( 2, self%max_iter_total ) )
        !
        if( allocated( self%EMrelErr) ) deallocate( self%EMrelErr )
        allocate( self%EMrelErr( self%max_iter_total ) )
        !
        if( allocated( self%solver%relErr ) ) deallocate( self%solver%relErr )
        allocate( self%solver%relErr( self%solver%max_iter ) )

     end subroutine createDiagnosticArraysForwardSolverIT
     !
     ! ForwardSolverIT initDiagnostic:
     !    Init the arrays used for diagnostic analysis.
     subroutine initDiagnosticArraysForwardSolverIT( self )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        !
        self%n_iter_total = 0
        self%EMrelErr = R_ZERO
        self%divJ = R_ZERO
        self%solver%failed = .false.
        self%solver%converged = .false.
        !
     end subroutine initDiagnosticArraysForwardSolverIT 
     !
     ! ForwardSolverIT setFrequency:
     !    Set omega: For itself and its solver
     !    Calls base setPreconditioner()
     subroutine setFrequencyForwardSolverIT( self, omega )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        real( kind=prec ), intent( in )             :: omega
        !
        if( abs( self%omega-omega ) .lt. TOL8 ) then
          !
          ! Professional courtesy? Why not simply do everything if .gt. TOL8????
          return
        end if
        !
        self%omega = omega
        self%solver%omega = omega
        call self%solver%preconditioner%SetPreconditioner( omega )
       !
     end subroutine setFrequencyForwardSolverIT
     !
     ! ForwardSolverIT setFrequency:
     !    Init the arrays used for diagnostic analysis.
     !    Calls base setPreconditioner()
     function getESolutionForwardSolverIT( self, source, polarization ) result( e_solution )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        class( Source_t ), intent( inout )          :: source
        integer, intent( in )                       :: polarization
        !
        class( cVector_t ), allocatable :: e_solution, temp
		class( cScalar_t ), allocatable :: phi0
        integer :: iter
        !
        write(*,*) "getESolution ForwardSolverIT for pol:", polarization
        !
        ! initialize diagnostics -- am assuming that setting of solver parameters
        !  is done in a set up step (once in the run) outside this object
        call self%initDiagnosticArrays()
        !
        call self%solver%zeroDiagnostics()
        !
        ! initialize solution
        allocate( e_solution, source = source%e0 )
        !
        loop: do while ( ( .not. self%solver%converged ) &
		                   .and. &
                         ( .not. self%solver%failed ) )
           !
           ! Need to be discussed
           select type( solver => self%solver )
              class is( Solver_QMR_t )
                 call solver%solve( source%rhs, e_solution )
              class default
                 write(*, *) "ERROR:ForwardSolverIT::getESolutionForwardSolverIT:"
                 stop        "         Unknow solver type."
           end select
           !
           ! Gary comments:
           ! I am just copying this -- while we work on this should
           ! reconsider implementation
           ! solver%converged when the relative error is less than tolerance
           ! 
           ! Werdt comments:
           ! GOT IT, WORKING ON IT !!!
           self%solver%converged = self%solver%n_iter .lt. self%solver%max_iter
           ! 
           ! Gary comments:
           ! there are two ways of failing: 1) QMR did not work or
           !     2) total number of divergence corrections exceeded
           ! 
           ! Werdt comments:
           ! then better to implement in Solver?
           self%solver%failed = self%solver%failed .or. self%failed
           !
           ! update solver diagnostics 
           do iter = 1, self%solver%n_iter
              ! why are we using an explicit loop here?
              self%EMrelErr( self%n_iter_total + iter ) = self%solver%relErr(iter)
           enddo
           !
           self%n_iter_total = self%n_iter_total + self%solver%n_iter
           !
       enddo loop
       !
       if( source%adjt ) then
          select type( modOp => self%solver%model_operator )
             class is ( ModelOperator_MF_t )
             !
             e_solution = e_solution * modOp%Metric%Vedge
             !
             class default
                write(*, *) "ERROR:ForwardSolverIT_t::getESolutionForwardSolverIT:"
                STOP        "model_operator type unknow"
          end select
          ! just leave bdry values set to 0 ???
       else
          !
          e_solution = e_solution + source%bdry
          !
       endif
       !
       ! JUST TO SEE THE E_SOLUTION RESULT
       select type( e_solution )
          class is( cVector3D_SG_t )
                write( *, * ) "         ", e_solution%nx, e_solution%ny, e_solution%nz, e_solution%gridType
          class default
                stop "Unclassified ForwardSolverIT e_solution"
       end select
       !
       if( allocated( phi0 ) ) deallocate(phi0)
       !
   end function getESolutionForwardSolverIT
   !
end Module ForwardSolverIT
 
