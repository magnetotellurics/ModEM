!
!> Derived class to define a iterative ForwardSolver using divergence correction
!
module ForwardSolverIT_DC
    !
    use ForwardSolverIT
    use DivergenceCorrection
    !
    type, extends( ForwardSolverIT_t ) :: ForwardSolverIT_DC_t
        !
        integer :: n_divcor, max_divcor_iters
        !
        real( kind=prec ) :: tol_div_cor
        !
        type( DivergenceCorrection_t ) :: divergence_correction 
        !
        contains
            !
            procedure, public :: setFrequency => setFrequency_ForwardSolverIT_DC
            !
            procedure, public :: setIterControl => setIterControl_ForwardSolverIT_DC
            !
            procedure, public :: initDiagnostics => initDiagnostics_ForwardSolverIT_DC
            !
            procedure, public :: zeroDiagnostics => zeroDiagnostics_ForwardSolverIT_DC
            !
            procedure, public :: createESolution => createESolution_ForwardSolverIT_DC
            !
            procedure, public :: copyFrom => copyFrom_ForwardSolverIT_DC
            !
    end type ForwardSolverIT_DC_t
    !
    interface ForwardSolverIT_DC_t
        module procedure ForwardSolverIT_DC_ctor
    end interface ForwardSolverIT_DC_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ForwardSolverIT_DC_ctor( model_operator, solver_type ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        character(*), intent( in ) :: solver_type
        !
        type( ForwardSolverIT_DC_t ) :: self
        !
        !write( *, * ) "Constructor ForwardSolverIT_DC_t"
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
                call errStop( "ForwardSolverIT_DC_ctor > Not yet coded for Bi-Conjugate Gradients" )
            case default
                call errStop( "ForwardSolverIT_DC_ctor > Unknown solver" )
            !
        end select
        !
        call self%setIterControl
        !
        call self%initDiagnostics
        !
        self%divergence_correction = DivergenceCorrection_t( model_operator )
        !
    end function ForwardSolverIT_DC_ctor
    !
    !> Procedure setFrequency_ForwardSolverIT_DC
    !> Set omega for this ForwardSolver (Called on the main transmitter loop at main program)
    !
    subroutine setFrequency_ForwardSolverIT_DC( self, sigma, period )
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
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
        !> Set conductivity for the model operator (again ????)
        call self%solver%preconditioner%model_operator%divCorSetUp
        !
        !> Set conductivity for the divergence_correction
        call self%divergence_correction%setCond( self%solver%omega )
        !
        call self%initDiagnostics
        !
    end subroutine setFrequency_ForwardSolverIT_DC
    !
    !> No subroutine briefing
    subroutine setIterControl_ForwardSolverIT_DC( self )
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
        !
        self%n_divcor = 0
        !
        self%max_solver_calls = max_solver_calls
        !
        self%max_divcor_iters = max_divcor_iters
        !
        self%tol_div_cor = tolerance_divcor
        !
        self%tolerance = self%solver%tolerance
        !
        self%max_iter_total = self%solver%max_iters * self%max_solver_calls
        !
    end subroutine setIterControl_ForwardSolverIT_DC
    !
    !> No subroutine briefing
    subroutine initDiagnostics_ForwardSolverIT_DC( self )
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
        !
        self%n_iter_actual = 0
        !
        self%relResFinal = R_ZERO
        !
        if( .NOT. allocated( self%relResVec ) ) then
            allocate( self%relResVec( self%max_iter_total ) )
        endif
        !
    end subroutine initDiagnostics_ForwardSolverIT_DC
    !
    !> No subroutine briefing
    subroutine zeroDiagnostics_ForwardSolverIT_DC( self )
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
        !
        self%relResVec = R_ZERO
        !
        call self%solver%zeroDiagnostics
        !
    end subroutine zeroDiagnostics_ForwardSolverIT_DC
    !
    !> No subroutine briefing
    !
    subroutine createESolution_ForwardSolverIT_DC( self, pol, source, e_solution )
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
        integer, intent( in ) :: pol
        class( Source_t ), intent( in ) :: source
        class( Vector_t ), allocatable, intent( out ) :: e_solution
        !
        class( Vector_t ), allocatable :: temp_vec
        class( Scalar_t ), allocatable :: phi0
        integer :: iter
        !
        call self%solver%zeroDiagnostics
        !
        self%solver%converged = .FALSE.
        self%solver%failed = .FALSE.
        self%n_divcor = 0
        self%n_iter_actual = 0
        !
        !> Create e_solution Vector
        call self%solver%preconditioner%model_operator%metric%createVector( complex_t, EDGE, e_solution )
        !
        call e_solution%zeros
        !
        if( source%non_zero_source ) then
            !
            !> Create phi0
            call self%divergence_correction%rhsDivCor( self%solver%omega, source%E( pol )%v, phi0 )
            !
        endif
        !
        loop: do while ( ( .NOT. self%solver%converged ) .AND. ( .NOT. self%solver%failed ) )
            !
            select type( solver => self%solver )
                !
                class is( Solver_QMR_t )
                    call solver%solve( source%rhs( pol )%v, e_solution )
                class default
                    stop "Error: getESolutionForwardSolverIT_DC > Unknown solver type."
                !
            end select
            !
            self%solver%converged = self%solver%n_iter .LT. self%solver%max_iters
            !
            self%solver%failed = self%solver%failed .OR. self%failed
            !
            !write( *, * ) "n_iter_actual+iter,     iter,     solver%relErr(iter)"
            !
            do iter = 1, self%solver%n_iter
                !
                self%relResVec( self%n_iter_actual + iter ) = self%solver%relErr( iter )
                !
                !write( *, * ) self%n_iter_actual + iter, iter, self%solver%relErr( iter )
                !
            enddo
            !
            self%n_iter_actual = self%n_iter_actual + self%solver%n_iter
            !
            self%n_divcor = self%n_divcor + 1
            !
            if( .NOT. self%solver%converged )  then
                !
                if( self%n_divcor < self%max_solver_calls ) then
                    !
                    !> USING THIS TEMPORARY VARIABLE IMPROVES THE EXECUTION TIME CONSIDERABLY...
                    allocate( temp_vec, source = e_solution )
                    !
                    if( source%non_zero_source ) then
                        !
                        call self%divergence_correction%divCorr( temp_vec, e_solution, phi0 )
                        !
                    else
                        !
                        call self%divergence_correction%divCorr( temp_vec, e_solution )
                        !
                    endif
                    !
                    deallocate( temp_vec )
                    !
                else
                    !
                    self%solver%failed = .TRUE.
                    !
                endif
                !
            endif
        !
        enddo loop
        !
        if( source%non_zero_source ) deallocate( phi0 )
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
            call source%rhs( pol )%v%boundary( temp_vec )
            !
        else
            !
            call source%E( pol )%v%boundary( temp_vec )
            !
        endif
        !
        call e_solution%add( temp_vec )
        !
        deallocate( temp_vec )
        !
    end subroutine createESolution_ForwardSolverIT_DC
    !
    !> No subroutine briefing
    subroutine copyFrom_ForwardSolverIT_DC( self, rhs )
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
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
            class is( ForwardSolverIT_DC_t )
                !
                self%n_divcor = rhs%n_divcor
                !
                self%max_solver_calls = rhs%max_solver_calls
                !
                self%max_divcor_iters = rhs%max_divcor_iters
                !
                self%tol_div_cor = rhs%tol_div_cor
                !
                self%divergence_correction = rhs%divergence_correction
                !
            class default
               call errStop( "copyFrom_ForwardSolverIT_DC > Incompatible input." )
            !
        end select
        !
    end subroutine copyFrom_ForwardSolverIT_DC
    !
end Module ForwardSolverIT_DC
