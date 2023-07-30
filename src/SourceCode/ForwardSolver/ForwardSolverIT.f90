!
!> Derived class to define a iterative ForwardSolver using divergence correction
!
module ForwardSolverIT
    !
    use ForwardSolver
    use Solver_QMR
    !
    type, extends( ForwardSolver_t ) :: ForwardSolverIT_t
        !
        integer :: max_solver_calls
        !
        contains
            !
            procedure, public :: setFrequency => setFrequency_ForwardSolverIT
            !
            procedure, public :: setIterControl => setIterControl_ForwardSolverIT
            !
            procedure, public :: initDiagnostics => initDiagnostics_ForwardSolverIT
            !
            procedure, public :: zeroDiagnostics => zeroDiagnostics_ForwardSolverIT
            !
            procedure, public :: createESolution => createESolution_ForwardSolverIT
            !
            procedure, public :: setIterDefaults => setIterDefaults_ForwardSolverIT
            !
            procedure, public :: copyFrom => copyFrom_ForwardSolverIT
            !
    end type ForwardSolverIT_t
    !
    interface ForwardSolverIT_t
        module procedure ForwardSolverIT_ctor
    end interface ForwardSolverIT_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ForwardSolverIT_ctor( model_operator, solver_type ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        character(*), intent( in ) :: solver_type
        !
        type( ForwardSolverIT_t ) :: self
        !
        write( *, * ) "Constructor ForwardSolverIT_t"
        !
        call self%baseInit
        !
        select case( solver_type )
            !
            case( QMR )
                !
                if( allocated( self%solver )  ) deallocate( self%solver )
                allocate( self%solver, source = Solver_QMR_t( model_operator ) )
                !
            case( BiCG )
                call errStop( "ForwardSolverIT_ctor > Not yet coded for Bi-Conjugate Gradients" )
            case default
                call errStop( "ForwardSolverIT_ctor > Unknown solver" )
            !
        end select
        !
        !> Set default values for this ForwardSolver
        call self%setIterDefaults
        !
        !> Set max number of all forward solver iterations
        self%max_iter_total = self%max_solver_calls * self%solver%max_iters
        !
        call self%setIterControl
        !
        call self%initDiagnostics
        !
    end function ForwardSolverIT_ctor
    !
    !> Procedure setFrequency_ForwardSolverIT
    !> Set omega for this ForwardSolver (Called on the main transmitter loop at main program)
    !
    subroutine setFrequency_ForwardSolverIT( self, sigma, period )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( inout ) :: sigma
        real( kind=prec ), intent( in ) :: period
        !
        !> Set omega for this ForwardSolver solver
        self%solver%omega = ( 2.0 * PI / period )
        !
        !> Set conductivity for the model operator (again ????)
        call self%solver%preconditioner%model_operator%setCond( sigma, self%solver%omega )
        !
        !> Set preconditioner for this solver's preconditioner
        call self%solver%preconditioner%setPreconditioner( self%solver%omega )
        !
        call self%initDiagnostics
        !
    end subroutine setFrequency_ForwardSolverIT
    !
    !> No subroutine briefing
    subroutine setIterControl_ForwardSolverIT( self )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        !
        self%tolerance = self%solver%tolerance
        !
        self%max_solver_calls = self%max_iter_total / self%solver%max_iters
        !
        self%max_iter_total = self%solver%max_iters * self%max_solver_calls
        !
    end subroutine setIterControl_ForwardSolverIT
    !
    !> No subroutine briefing
    subroutine setIterDefaults_ForwardSolverIT( self )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        !
        self%max_solver_calls = max_solver_calls
        !
    end subroutine setIterDefaults_ForwardSolverIT
    !
    !> No subroutine briefing
    subroutine initDiagnostics_ForwardSolverIT( self )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        !
        self%n_iter_actual = 0
        !
        self%relResFinal = R_ZERO
        !
        if( .NOT. allocated( self%relResVec ) ) then
            allocate( self%relResVec( self%max_iter_total ) )
        endif
        !
    end subroutine initDiagnostics_ForwardSolverIT
    !
    !> No subroutine briefing
    subroutine zeroDiagnostics_ForwardSolverIT( self )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        !
        self%relResVec = R_ZERO
        !
        call self%solver%zeroDiagnostics
        !
    end subroutine zeroDiagnostics_ForwardSolverIT
    !
    !> No subroutine briefing
    !
    subroutine createESolution_ForwardSolverIT( self, pol, source, e_solution )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        integer, intent( in ) :: pol
        class( Source_t ), intent( in ) :: source
        class( Vector_t ), allocatable, intent( out ) :: e_solution
        !
        class( Vector_t ), allocatable :: temp_vec
        !
        integer :: iter
        !
        call self%solver%zeroDiagnostics
        !
        self%solver%converged = .FALSE.
        self%solver%failed = .FALSE.
        self%n_iter_actual = 0
        !
        !> Create e_solution Vector
        call self%solver%preconditioner%model_operator%metric%createVector( complex_t, EDGE, e_solution )
        !
        call e_solution%zeros
        !
        loop: do while ( ( .NOT. self%solver%converged ) .AND. ( .NOT. self%solver%failed ) )
            !
            select type( solver => self%solver )
                !
                class is( Solver_QMR_t )
                    call solver%solve( source%rhs( pol ), e_solution )
                class default
                    call errStop( "createESolution_ForwardSolverIT > Unknown solver type." )
                !
            end select
            !
            self%solver%converged = self%solver%n_iter .LT. self%solver%max_iters
            !
            self%solver%failed = self%solver%failed .OR. self%failed
            !
            do iter = 1, self%solver%n_iter
                !
                self%relResVec( self%n_iter_actual + iter ) = self%solver%relErr( iter )
                !
            enddo
            !
            self%n_iter_actual = self%n_iter_actual + self%solver%n_iter
            !
        enddo loop
        !
        self%relResFinal = self%relResVec( self%n_iter_actual )
        !
        !> Just for the serialJMult_T SourceInteriorForce case
        if( source%for_transpose ) then
            !
            call e_solution%mult( self%solver%preconditioner%model_operator%metric%v_edge )
            !
        endif
        !
        if( source%non_zero_bc ) then
            !
            call source%rhs( pol )%boundary( temp_vec )
            !
        else
            !
            call source%E( pol )%boundary( temp_vec )
            !
        endif
        !
        call e_solution%add( temp_vec )
        !
        deallocate( temp_vec )
        !
    end subroutine createESolution_ForwardSolverIT
    !
    !> No subroutine briefing
    subroutine copyFrom_ForwardSolverIT( self, rhs )
        implicit none
        !
        class( ForwardSolverIT_t ), intent( inout ) :: self
        class( ForwardSolver_t ), intent( in ) :: rhs
        !
        self%solver = rhs%solver
        !
        self%max_iter_total = rhs%max_iter_total
        !
        self%n_iter_actual = rhs%n_iter_actual
        !
        self%tolerance = rhs%tolerance
        !
        self%relResFinal = rhs%relResFinal
        !
        self%relResVec = rhs%relResVec
        !
        self%failed = rhs%failed
        !
        select type( rhs )
            !
            class is( ForwardSolverIT_t )
                !
                self%max_solver_calls = rhs%max_solver_calls
                !
            class default
               call errStop( "copyFrom_ForwardSolverIT > Incompatible input." )
            !
        end select
        !
    end subroutine copyFrom_ForwardSolverIT
    !
end Module ForwardSolverIT
