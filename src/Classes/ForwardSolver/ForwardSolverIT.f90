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
      !
   contains
      !
      final :: ForwardSolverIT_dtor
      !
      procedure, public :: setPeriod              => setPeriodForwardSolverIT
      procedure, public :: initDiagnostics        => initDiagnosticsForwardSolverIT
      procedure, public :: zeroDiagnostics        => zeroDiagnosticsForwardSolverIT
      procedure, public :: setIterControl         => SetIterControlSolverIT
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
      !    
      function ForwardSolverIT_ctor( model_operator, solver_type ) result( self )
         implicit none
         class( ModelOperator_t ), intent( in ) :: model_operator
         character(*), intent(in)   :: solver_type
         type( ForwardSolverIT_t )              :: self
         !
         !write(*,*) "Constructor ForwardSolverIT_t"
         !
         call self%init()  !   is this really needed????
         !
         select case(solver_type)
            case(QMR)
               self%solver = Solver_QMR_t( model_operator )
            case(BiCG)
               write(*,*) 'Not yet coded for Bi-Conjugate Gradients'
               stop
         end select
         !   this sets defaults iteration control for particular solver
         call self%solver%setDefaults()
         !
         !    initialize Fwd operator iteration control, diagonstic arrays
         !     using defaults from solver
         call self%setIterControl(self%solver%max_iter,self%solver%tolerance)
         !
         call self%initDiagnostics()
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
      !    Sets Period and Omega -- separate setFrequency procedure elminated!
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
         !   set frequency in solver object
         self%solver%omega = omega
         !     set preconditoner (depends on frequency in general)
         call self%solver%preconditioner%SetPreconditioner( omega )
         !
       end subroutine setPeriodForwardSolverIT
      !
      ! ForwardSolverIT initDiagnostic:
      !    Init the arrays used for diagnostic analysis.
      subroutine setIterControlForwardSolverIT( self, maxit, tol )
         implicit none
         class( ForwardSolverIT_t ), intent( inout ) :: self
         integer, intent(in)  :: maxit
         real(kind=prec)      :: tol
         !
         self%max_iter_total = maxit
         self%tolerance = tol
         !
         !   if this is not called from ctor, input tol and maxit may
         !    not match what is set in solver -- set explicitly
         !     to make sure this is correct
         call self%solver%setParameters(maxit,tol)

      end subroutine setIterControlForwardSolverIT 
      !
      !**********
      !
      ! ForwardSolverIT initDiagnostic:
      !    Init the arrays used for diagnostic analysis.
      !   NOTE: this should be called AFTER any reset of iteration
      !    control parameters
      subroutine initDiagnosticsForwardSolverIT( self )
         implicit none
         class( ForwardSolverIT_t ), intent( inout ) :: self
         !
         self%n_iter_actual = 0
         self%relResFinal = R_ZERO
         !
         if allocated(self%relResVec) deallocate(self%relResVec)
         allocate(self%relResVec(self%max_iter_total))
         !
      end subroutine initDiagnosticsForwardSolverIT 
      !
      !*********
      !
      subroutine zeroDiagnosticsForwardSolverIT(self)
         implicit none
         class( ForwardSolverIT_DC_t ), intent( inout ) :: self

           self%relErr = R_ZERO
           self%solver%zeroDiagnostics

      end subroutine zeroDiagnosticsForwardSolverIT
      !
      !**********
      !
      subroutine getESolutionForwardSolverIT( self, source, e_solution )
         implicit none
         !
         class( ForwardSolverIT_t ), intent( inout ) :: self
         class( Source_t ), intent( inout )          :: source
         class( cVector_t ), intent(inout) :: e_solution
          
         integer :: iter
         !
         ! zero diagnostics, including for solver
         !  is done in a set up step (once in the run) outside this object
         call self%zeroDiagnostics()
         !
         select type( solver => self%solver )
            class is( Solver_QMR_t )
               call solver%solve( source%rhs, e_solution )
               write(*,*) 'n_iter = ',self%solver%n_iter
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
     end subroutine getESolutionForwardSolverIT
   !
end Module ForwardSolverIT
 
