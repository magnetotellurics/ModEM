!
!> Derived class to define a iterative ForwardSolver using divergence correction
!
module ForwardSolverIT_DC
    !
    use ForwardSolver
    use DivergenceCorrection
    use Solver_QMR
    !
    type, extends( ForwardSolver_t ) :: ForwardSolverIT_DC_t
        !
        integer :: n_divcor, max_div_cor, max_divcor_iters
        !
        real( kind=prec ) :: tol_div_cor
        !
        type( DivergenceCorrection_t ) :: divergence_correction 
        !
        contains
            !
            procedure, public :: setFrequency => setFrequencyForwardSolverIT_DC
            !
            procedure, public :: setIterControl => setIterControlForwardSolverIT_DC
            !
            procedure, public :: initDiagnostics => initDiagnosticsForwardSolverIT_DC
            !
            procedure, public :: zeroDiagnostics => zeroDiagnosticsForwardSolverIT_DC
            !
            procedure, public :: createESolution => createESolutionForwardSolverIT_DC
            !
            procedure, public :: setIterDefaults => setIterDefaultsForwardSolverIT_DC
            !
    end type ForwardSolverIT_DC_t
    !
    interface ForwardSolverIT_DC_t
        module procedure ForwardSolverIT_DC_ctor
    end interface ForwardSolverIT_DC_t
    !
contains
    !
    !> No function briefing
    function ForwardSolverIT_DC_ctor( model_operator, solver_type ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        character(*), intent( in ) :: solver_type
        type( ForwardSolverIT_DC_t ) :: self
        !
        integer :: max_iter
        !
        !write( *, * ) "Constructor ForwardSolverIT_DC_t"
        !
        call self%init()
        !
        select case( solver_type )
            !
            case( QMR )
                !
                if( allocated( self%solver )  ) deallocate( self%solver )
                allocate( self%solver, source = Solver_QMR_t( model_operator ) )
                !
            case( BiCG )
                stop "ForwardSolverIT_DC_ctor > Not yet coded for Bi-Conjugate Gradients"
            case default
                stop "ForwardSolverIT_DC_ctor > Unknown solver"
            !
        end select
        !
        self%n_divcor = 0
        !
        self%max_div_cor = 0
        !
        self%max_divcor_iters = 0
        !
        self%tol_div_cor = R_ZERO
        !
        !> Set default values for this ForwardSolver
        call self%setIterDefaults()
        !
        !> Set max number of all forward solver iterations
        self%max_iter_total = self%max_div_cor * self%solver%max_iter
        !
        call self%setIterControl
        !
        call self%initDiagnostics()
        !
        self%divergence_correction = DivergenceCorrection_t( model_operator )
        !
    end function ForwardSolverIT_DC_ctor
    !
    !> Procedure setFrequencyForwardSolverIT_DC
    !> Set omega for this ForwardSolver (Called on the main transmitter loop at main program)
    subroutine setFrequencyForwardSolverIT_DC( self, sigma, period )
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: sigma
        real( kind=prec ), intent( in ) :: period
        !
        !> Set omega for this ForwardSolver solver
        self%solver%omega = ( 2.0 * PI / period )
        !
        !> Set conductivity for the model operator (again ????)
        call self%solver%preconditioner%model_operator%setCond( sigma )
        !
        !> Set omega for the divergence_correction´s solver
        self%divergence_correction%solver%omega = self%solver%omega
        !
        !> Set conductivity for the divergence_correction
        call self%divergence_correction%SetCond()
        !
        !> Set preconditioner for this solver´s preconditioner
        call self%solver%preconditioner%SetPreconditioner( self%solver%omega )
        !
        call self%initDiagnostics()
        !
    end subroutine setFrequencyForwardSolverIT_DC
    !
    !> No subroutine briefing
    subroutine setIterControlForwardSolverIT_DC( self )
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
        !
        self%tolerance = self%solver%tolerance
        !
        self%max_div_cor = self%max_iter_total / self%solver%max_iter
        !
        self%max_iter_total = self%solver%max_iter * self%max_div_cor
        !
    end subroutine setIterControlForwardSolverIT_DC
    !
    !> No subroutine briefing
    subroutine setIterDefaultsForwardSolverIT_DC( self )
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
        !
        self%max_div_cor      = max_divcor_calls
        self%max_divcor_iters = max_divcor_iters
        self%tol_div_cor      = tolerance_divcor
        !
    end subroutine setIterDefaultsForwardSolverIT_DC
    !
    !> No subroutine briefing
    subroutine initDiagnosticsForwardSolverIT_DC( self )
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
    end subroutine initDiagnosticsForwardSolverIT_DC
    !
    !> No subroutine briefing
    subroutine zeroDiagnosticsForwardSolverIT_DC(self)
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
        !
        self%relResVec = R_ZERO
        !
        call self%solver%zeroDiagnostics()
        !
    end subroutine zeroDiagnosticsForwardSolverIT_DC
    !
    !> No function briefing
    subroutine createESolutionForwardSolverIT_DC( self, pol, source, e_solution )
        implicit none
        !
        class( ForwardSolverIT_DC_t ), intent( inout ) :: self
        integer, intent( in ) :: pol
        class( Source_t ), intent( in ) :: source
        class( Vector_t ), intent( inout ) :: e_solution
        !
        class( Vector_t ), allocatable :: temp_vec
        class( Scalar_t ), allocatable :: phi0
        integer :: iter
        !
        call self%solver%zeroDiagnostics()
        !
        self%solver%converged = .FALSE.
        self%solver%failed    = .FALSE.
        self%n_divcor = 0
        self%n_iter_actual = 0
        !
        e_solution = cVector3D_SG_t( self%solver%preconditioner%model_operator%metric%grid, EDGE )
        !
        call e_solution%zeros()
        !
        if( source%non_zero_source ) then
            !
            allocate( phi0, source = cScalar3D_SG_t( self%solver%preconditioner%model_operator%metric%grid, NODE ) )
            !
            call self%divergence_correction%rhsDivCor( self%solver%omega, source%E( pol ), phi0 )
            !
            !> USING THIS TEMPORARY VARIABLE IMPROVES THE EXECUTION TIME CONSIDERABLY...
            allocate( temp_vec, source = e_solution )
            !
            call self%divergence_correction%divCorr( temp_vec, e_solution, phi0 )
            !
            deallocate( temp_vec )
            !
            self%n_divcor = 1
            !
        endif
        !
        loop: do while ( ( .NOT. self%solver%converged ) .AND. ( .NOT. self%solver%failed ) )
            !
            select type( solver => self%solver )
                !
                class is( Solver_QMR_t )
                    call solver%solve( source%rhs( pol ), e_solution )
                class default
                    stop "Error: getESolutionForwardSolverIT_DC > Unknown solver type."
                !
            end select
            !
            self%solver%converged = self%solver%n_iter .LT. self%solver%max_iter
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
                if( self%n_divcor < self%max_div_cor ) then
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
        !> Just for JMult_T SourceInteriorForce case
        if( source%trans ) then
            !
            call e_solution%mult( self%solver%preconditioner%model_operator%metric%VEdge )
            !
        else
            !
            call source%E( pol )%Boundary( temp_vec )
            !
            call e_solution%add( temp_vec )
            !
            deallocate( temp_vec )
            !
        endif
        !
    end subroutine createESolutionForwardSolverIT_DC
    !
end Module ForwardSolverIT_DC
