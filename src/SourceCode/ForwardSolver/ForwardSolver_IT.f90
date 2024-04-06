!
!> Derived class to define a iterative ForwardSolver
!
module ForwardSolver_IT
    !
    use ForwardSolver
    use Solver_QMR
    use Solver_BICG
    !use Solver_BICG_OMP
    !
    type, extends( ForwardSolver_t ) :: ForwardSolver_IT_t
        !
        integer :: iter, max_solver_calls
        !
        contains
            !
            procedure, public :: setFrequency => setFrequency_ForwardSolver_IT
            !
            procedure, public :: setIterControl => setIterControl_ForwardSolver_IT
            !
            procedure, public :: initDiagnostics => initDiagnostics_ForwardSolver_IT
            !
            procedure, public :: zeroDiagnostics => zeroDiagnostics_ForwardSolver_IT
            !
            procedure, public :: createESolution => createESolution_ForwardSolver_IT
            !
            procedure, public :: copyFrom => copyFrom_ForwardSolver_IT
            !
    end type ForwardSolver_IT_t
    !
    interface ForwardSolver_IT_t
        module procedure ForwardSolver_IT_ctor
    end interface ForwardSolver_IT_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ForwardSolver_IT_ctor( model_operator, solver_type ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        character(*), intent( in ) :: solver_type
        !
        type( ForwardSolver_IT_t ) :: self
        !
        !write( *, * ) "Constructor ForwardSolver_IT_t"
        !
        call self%baseInit
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
                call warning( "ForwardSolver_IT_ctor > solver_type not provided, using BICG." )
                !
                allocate( self%solver, source = Solver_BICG_t( model_operator ) )
                !
            case default
                call errStop( "ForwardSolver_IT_ctor > Unknown solver ["//solver_type//"]" )
            !
        end select
        !
        call self%setIterControl
        !
        call self%initDiagnostics
        !
    end function ForwardSolver_IT_ctor
    !
    !> Procedure setFrequency_ForwardSolver_IT
    !> Set omega for this ForwardSolver (Called on the main transmitter loop at main program)
    !
    subroutine setFrequency_ForwardSolver_IT( self, sigma, period )
        implicit none
        !
        class( ForwardSolver_IT_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: sigma
        real( kind=prec ), intent( in ) :: period
        !
        !> Set omega for this ForwardSolver solver
        self%solver%omega = ( 2.0 * PI / period )
        !
        !> Set conductivity for the model operator
        call self%solver%preconditioner%model_operator%setCond( sigma, self%solver%omega )
        !
        !> Set preconditioner for this solver's preconditioner
        call self%solver%preconditioner%setPreconditioner( self%solver%omega )
        !
        call self%initDiagnostics
        !
    end subroutine setFrequency_ForwardSolver_IT
    !
    !> No subroutine briefing
    !
    subroutine setIterControl_ForwardSolver_IT( self )
        implicit none
        !
        class( ForwardSolver_IT_t ), intent( inout ) :: self
        !
        self%max_solver_calls = max_solver_calls
        !
        self%tolerance = self%solver%tolerance
        !
        self%max_iter_total = self%solver%max_iters * self%max_solver_calls
        !
    end subroutine setIterControl_ForwardSolver_IT
    !
    !> No subroutine briefing
    !
    subroutine initDiagnostics_ForwardSolver_IT( self )
        implicit none
        !
        class( ForwardSolver_IT_t ), intent( inout ) :: self
        !
        self%n_iter_actual = 0
        !
        self%relResFinal = R_ZERO
        !
        if( .NOT. allocated( self%relResVec ) ) then
            allocate( self%relResVec( self%max_iter_total ) )
        endif
        !
    end subroutine initDiagnostics_ForwardSolver_IT
    !
    !> No subroutine briefing
    !
    subroutine zeroDiagnostics_ForwardSolver_IT( self )
        implicit none
        !
        class( ForwardSolver_IT_t ), intent( inout ) :: self
        !
        self%relResVec = R_ZERO
        !
        call self%solver%zeroDiagnostics
        !
    end subroutine zeroDiagnostics_ForwardSolver_IT
    !
    !> No subroutine briefing
    !
    subroutine createESolution_ForwardSolver_IT( self, pol, source, e_solution )
        implicit none
        !
        class( ForwardSolver_IT_t ), intent( inout ) :: self
        integer, intent( in ) :: pol
        class( Source_t ), intent( in ) :: source
        class( Vector_t ), allocatable, intent( out ) :: e_solution
        !
        class( Vector_t ), allocatable :: source_e_vec, boundary_vec
        type( cVector3D_MR_t ) :: source_e_vec_mr
        integer :: i
        !
        !> Create proper SG or MR source_e_vec
        select type( grid => source%E( pol )%grid )
            !
            class is( Grid3D_SG_t )
                !
                allocate( source_e_vec, source = source%E( pol ) )
                !
            class is( Grid3D_MR_t )
                !
                source_e_vec_mr = cVector3D_MR_t( grid, source%E( pol )%grid_type )
                call source_e_vec_mr%fromSG( source%E( pol ) )
                allocate( source_e_vec, source = source_e_vec_mr )
                !
            class default
               call errStop( "createESolution_ForwardSolver_IT > Unclassified Source grid." )
            !
        end select
        !
        !> Create e_solution Vector
        call self%solver%preconditioner%model_operator%metric%createVector( complex_t, EDGE, e_solution )
        !
        !> Initialize FWD Solver
        self%iter = 1
        !
        self%n_iter_actual = 0
        !
        call self%solver%zeroDiagnostics
        !
        !> 
        fwd_solver_loop: do
            !
            !> 
            call self%solver%solve( source%rhs( pol )%v, e_solution )
            !
            do i = 1, self%solver%n_iter
                !
                self%relResVec( self%n_iter_actual + i ) = self%solver%relErr(i)
                !
            enddo
            !
            self%n_iter_actual = self%n_iter_actual + self%solver%n_iter
            !
            if( self%solver%converged .OR. ( self%iter .GE. self%max_solver_calls ) ) then
                exit
            endif
            !
            self%iter = self%iter + 1
            !
        enddo fwd_solver_loop
        !
        self%relResFinal = self%relResVec( self%n_iter_actual )
        !
        if( self%solver%converged ) then
            write( *, "( a46, i6, a7, es12.3 )" ) "IT converged within ", self%n_iter_actual, ": err= ", self%relResFinal
        else
            call warning( "createESolution_ForwardSolver_IT failed to converge!" )
        endif
        !
        !> Just for the serialJMult_T SourceAdjoint case
        if( source%for_transpose ) then
            !
            call e_solution%mult( self%solver%preconditioner%model_operator%metric%v_edge )
            !
        endif
        !
        if( source%non_zero_bc ) then
            !
            call source%rhs( pol )%v%boundary( boundary_vec )
            !
        else
            !
            call source_e_vec%boundary( boundary_vec )
            !
        endif
        !
        call e_solution%add( boundary_vec )
        !
        deallocate( boundary_vec )
        !
    end subroutine createESolution_ForwardSolver_IT
    !
    !> No subroutine briefing
    !
    subroutine copyFrom_ForwardSolver_IT( self, rhs )
        implicit none
        !
        class( ForwardSolver_IT_t ), intent( inout ) :: self
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
            class is( ForwardSolver_IT_t )
                !
                self%iter = rhs%iter
                !
                self%max_solver_calls = rhs%max_solver_calls
                !
            class default
               call errStop( "copyFrom_ForwardSolver_IT > Incompatible input." )
            !
        end select
        !
    end subroutine copyFrom_ForwardSolver_IT
    !
end Module ForwardSolver_IT
!