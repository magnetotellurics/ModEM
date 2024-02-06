!
!> Derived class to define a iterative ForwardSolver using Divergence Correction
!
module ForwardSolver_IT_DC
    !
    use ForwardSolver_IT
    use DivergenceCorrection
    !
    type, extends( ForwardSolver_IT_t ) :: ForwardSolver_IT_DC_t
        !
        type( DivergenceCorrection_t ) :: divergence_correction 
        !
        contains
            !
            procedure, public :: setFrequency => setFrequency_ForwardSolver_IT_DC
            !
            procedure, public :: createESolution => createESolution_ForwardSolver_IT_DC
            !
            procedure, public :: copyFrom => copyFrom_ForwardSolver_IT_DC
            !
    end type ForwardSolver_IT_DC_t
    !
    interface ForwardSolver_IT_DC_t
        module procedure ForwardSolver_IT_DC_ctor
    end interface ForwardSolver_IT_DC_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ForwardSolver_IT_DC_ctor( model_operator, solver_type ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        character(*), intent( in ) :: solver_type
        !
        type( ForwardSolver_IT_DC_t ) :: self
        !
        !write( *, * ) "Constructor ForwardSolver_IT_DC_t"
        !
        call self%baseInit
        !
        !> NEED THIS LINE ????
        !if( allocated( self%solver )  ) deallocate( self%solver )
        !
        select case( solver_type )
            !
            case( SLV_QMR )
                !
                allocate( self%solver, source = Solver_QMR_t( model_operator ) )
                !
            case( SLV_BICG )
                !
                allocate( self%solver, source = Solver_BICG_t( model_operator ) )
                !
            case( "" )
                !
                call warning( "ForwardSolver_IT_DC_ctor > solver_type not provided, using BICG." )
                !
                allocate( self%solver, source = Solver_BICG_t( model_operator ) )
                !
            case default
                call errStop( "ForwardSolver_IT_DC_ctor > Unknown solver ["//solver_type//"]" )
            !
        end select
        !
        call self%setIterControl
        !
        call self%initDiagnostics
        !
        self%divergence_correction = DivergenceCorrection_t( model_operator )
        !
    end function ForwardSolver_IT_DC_ctor
    !
    !> Procedure setFrequency_ForwardSolver_IT_DC
    !> Set omega for this ForwardSolver (Called on the main transmitter loop at main program)
    !
    subroutine setFrequency_ForwardSolver_IT_DC( self, sigma, period )
        implicit none
        !
        class( ForwardSolver_IT_DC_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: sigma
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
        !> Set conductivity for the model operator
        call self%solver%preconditioner%model_operator%divCorSetUp
        !
        !> Set conductivity for the divergence_correction
        call self%divergence_correction%setCond( self%solver%omega )
        !
        call self%initDiagnostics
        !
    end subroutine setFrequency_ForwardSolver_IT_DC
    !
    !> No subroutine briefing
    !
    subroutine createESolution_ForwardSolver_IT_DC( self, pol, source, e_solution )
        implicit none
        !
        class( ForwardSolver_IT_DC_t ), intent( inout ) :: self
        integer, intent( in ) :: pol
        class( Source_t ), intent( in ) :: source
        class( Vector_t ), allocatable, intent( out ) :: e_solution
        !
        class( Vector_t ), allocatable :: temp_e, temp_vec
        type( cVector3D_MR_t ) :: temp_e_mr
        class( Scalar_t ), allocatable :: phi0
        integer :: i
        !
        !> Create proper SG or MR temp source vectors
        select type( grid => source%E( pol )%grid )
            !
            class is( Grid3D_SG_t )
                !
                allocate( temp_e, source = source%E( pol ) )
                !
            class is( Grid3D_MR_t )
                !
                temp_e_mr = cVector3D_MR_t( grid, source%E( pol )%grid_type )
                call temp_e_mr%fromSG( source%E( pol ) )
                allocate( temp_e, source = temp_e_mr )
                !
            class default
               call errStop( "createESolution_ForwardSolver_IT > Unclassified Source grid." )
            !
        end select
        !
        !> Create e_solution Vector
        call self%solver%preconditioner%model_operator%metric%createVector( complex_t, EDGE, e_solution )
        !
        !> Inicialize DivergenceCorrection
        if( source%non_zero_source ) then
            !
            !> Create phi0
            call self%divergence_correction%rhsDivCor( self%solver%omega, temp_e, phi0 )
            !
        endif
        !
        !> Inicialize FWD Solver
        self%iter = 1
        !
        self%n_iter_actual = 0
        !
        call self%solver%zeroDiagnostics
        !
        fwd_solver_loop: do
            !
            !> 
            call self%solver%solve( source%rhs( pol )%v, e_solution )
            !
            do i = 1, self%solver%n_iter
                !
                self%relResVec( self%n_iter_actual + i ) = self%solver%relErr(i)
                !
                self%relResFinal = self%solver%relErr(i)
                !
            enddo
            !
            !> Apply Divergence Correction if solver not converged
            if( .NOT. self%solver%converged )  then
                !
                if( source%non_zero_source ) then
                    !
                    call self%divergence_correction%divCorr( e_solution, phi0 )
                    !
                else
                    !
                    call self%divergence_correction%divCorr( e_solution )
                    !
                endif
                !
            endif
            !
            !> Check Stop Conditions
            if( self%solver%converged .OR. ( self%iter .GE. self%max_solver_calls ) ) then
                exit
            endif
            !
            self%iter = self%iter + 1
            !
        enddo fwd_solver_loop
        !
        if( self%solver%converged ) then
            write( *, "( a49, i6, a7, es12.3 )" ) "IT_DC converged within ", self%iter, ": err= ", self%relResFinal
        else
            call warning( "createESolution_ForwardSolver_IT_DC failed to converge!" )
        endif
        !
        if( source%non_zero_source ) deallocate( phi0 )
        !
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
            call temp_e%boundary( temp_vec )
            !
        endif
        !
        call e_solution%add( temp_vec )
        !
        deallocate( temp_vec )
        !
    end subroutine createESolution_ForwardSolver_IT_DC
    !
    !> No subroutine briefing
    subroutine copyFrom_ForwardSolver_IT_DC( self, rhs )
        implicit none
        !
        class( ForwardSolver_IT_DC_t ), intent( inout ) :: self
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
        select type( rhs )
            !
            class is( ForwardSolver_IT_DC_t )
                !
                self%iter = rhs%iter
                !
                self%max_solver_calls = rhs%max_solver_calls
                !
                self%divergence_correction = rhs%divergence_correction
                !
            class default
               call errStop( "copyFrom_ForwardSolver_IT_DC > Incompatible input." )
            !
        end select
        !
    end subroutine copyFrom_ForwardSolver_IT_DC
    !
end Module ForwardSolver_IT_DC
!
