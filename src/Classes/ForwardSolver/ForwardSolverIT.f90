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
   ! real( kind=prec ), parameter :: tolEMDef = 1E-7
   !    default solver control parameters are set in the Solver object
   !
   type, extends( ForwardSolver_t ), public :: ForwardSolverIT_t
      !
      real( kind=prec )                              :: omega = 0.0
      real( kind=prec ), allocatable, dimension(:)   :: EMrelErr
      integer :: max_iter_total = 0
      !
   contains
      !
      final :: ForwardSolverIT_dtor
      !
      procedure, public :: setPeriod              => setPeriodForwardSolverIT
      procedure, public :: setFrequency           => setFrequencyForwardSolverIT
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
      call self%solver%setDefaults()
      !
      !    just dircectly create the total relative error array...`
      self%max_iter_total = self%solver%max_iter    !   really don't need this?

      allocate( self%EMrelErr( self%max_iter_total ) )
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
        if( abs( self%omega-omega ) .gt. TOL8 ) then
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
     ! polarization should probably not be a property of ForwardSolver -- I am eliminating
     !  for "File" version we need to pass this via the source (makes more sense to have
     !      any information that defines the file in "source" object -- need an source type
     !      tailored to file input ...
     function getESolutionForwardSolverIT( self, source ) result( e_solution )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        class( Source_t ), intent( inout )          :: source
        !
        class( cVector_t ), allocatable :: e_solution, temp
        integer :: iter
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
        select type( solver => self%solver )
            class is( Solver_QMR_t )
               call solver%solve( source%rhs, e_solution )
            class default
               write(*, *) "ERROR:ForwardSolverIT::getESolutionForwardSolverIT:"
               stop        "         Unknow solver type."
        end select
        !
        ! update solver diagnostic array EMrelErr -- in this case just a copy of the
        !   relErr array in solver ...
        self%EMrelErr(1:self%solver%n_iter) = self%solver%relErr(1:self%solver%n_iter)
        !
        self%n_iter_total = self%solver%n_iter
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
    end function getESolutionForwardSolverIT
   !
end Module ForwardSolverIT
 
