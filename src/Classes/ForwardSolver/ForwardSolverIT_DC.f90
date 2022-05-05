!
! Please ad some fancy relevant description for this class
!
module ForwardSolverIT_DC
    !
    use ForwardSolver
    use DivergenceCorrection
    use ModelOperator_MF
    use Solver_QMR
    !
    !
    type, extends( ForwardSolver_t ) :: ForwardSolverIT_DC_t
        !
        type( DivergenceCorrection_t ) :: divergence_correction 
        !
        integer :: nDivCor, max_div_cor, max_iterDivCor
        !
        real( kind=prec ) :: tolDivCor
        !
        contains
            !
            final :: ForwardSolverIT_DC_dtor
            !
            procedure, public :: setFrequency    => setFrequencyForwardSolverIT_DC
            procedure, public :: setIterControl  => setIterControlForwardSolverIT_DC
            procedure, public :: initDiagnostics => initDiagnosticsForwardSolverIT_DC
            procedure, public :: zeroDiagnostics => zeroDiagnosticsForwardSolverIT_DC
            procedure, public :: getESolution    => getESolutionForwardSolverIT_DC
            !
            procedure, public :: setIterDefaultsDC
            !
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
            !
            class( ModelOperator_t ), intent( in ) :: model_operator
            character(*), intent(in)               :: solver_type
            type( ForwardSolverIT_DC_t )           :: self
            !
            integer :: max_iter
            !
            !write(*,*) "Constructor ForwardSolverIT_DC_t"
            !
            call self%init()
            !
            self%nDivCor = 0
            self%max_div_cor = 0
            self%max_iterDivCor = 0
            self%tolDivCor = 0.0
            !
            select case( solver_type )
                case( QMR )
                    !
                    self%solver = Solver_QMR_t( model_operator )
                    !
                case( BiCG )
                    stop "ForwardSolverIT_DC_ctor: Not yet coded for Bi-Conjugate Gradients"
                case default
                    stop "ForwardSolverIT_DC_ctor: Unknow solver"
            end select
            !
            ! Set default values for this ForwardSolver
            call self%setIterDefaultsDC()
            !
            ! Set max number of all forward solver iterations
            self%max_iter_total = self%max_div_cor * self%solver%max_iter
            !
            call self%setIterControl
            !
            !
            self%divergence_correction = DivergenceCorrection_t( model_operator )
            !
            !
            call self%initDiagnostics()
            !
        end function ForwardSolverIT_DC_ctor
        !
        !
        subroutine ForwardSolverIT_DC_dtor( self )
            implicit none
            !
            type( ForwardSolverIT_DC_t ), intent( inout ) :: self
            !
            !write(*,*) "Destructor ForwardSolverIT_DC_t"
            !
            call self%dealloc()
            !
        end subroutine ForwardSolverIT_DC_dtor
        !
        !
        subroutine setFrequencyForwardSolverIT_DC( self, model_parameter, period )
            implicit none
            !
            class( ForwardSolverIT_DC_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in )        :: model_parameter
            real( kind=prec ), intent( in )                :: period
            !
            ! Set omega for this solver
            self%solver%omega = ( 2.0 * PI / period )
            !
            !
            call self%solver%preconditioner%model_operator%setCond( model_parameter )
            !
            call self%divergence_correction%SetCond()
            !
            ! Set this preconditioner
            call self%solver%preconditioner%SetPreconditioner( self%solver%omega )
            !
            ! Set omega for divergence_correction´s solver
            self%divergence_correction%solver%omega = self%solver%omega
            !
        end subroutine setFrequencyForwardSolverIT_DC
        !
        subroutine setIterControlForwardSolverIT_DC( self )
            implicit none
            !
            class( ForwardSolverIT_DC_t ), intent( inout ) :: self
            !
            !
            self%tolerance = self%solver%tolerance
            !
            self%max_div_cor = self%max_iter_total / self%solver%max_iter
            !
            self%max_iter_total = self%solver%max_iter * self%max_div_cor
            !
        end subroutine setIterControlForwardSolverIT_DC
        !
        !
        subroutine setIterDefaultsDC( self )
            implicit none
            !
            class( ForwardSolverIT_DC_t ), intent( inout ) :: self
            !
            self%max_div_cor    = max_div_corDef
            self%max_IterDivCor = max_IterDivCorDef
            self%tolDivCor      = tolDivCorDef
            !
        end subroutine setIterDefaultsDC
        !
        !
        subroutine initDiagnosticsForwardSolverIT_DC( self )
            implicit none
            !
            class( ForwardSolverIT_DC_t ), intent( inout ) :: self
            !
            self%n_iter_actual = 0
            self%relResFinal   = R_ZERO
            !
            allocate( self%relResVec( self%max_iter_total ) )
            !
         end subroutine initDiagnosticsForwardSolverIT_DC
         !
         !
         subroutine zeroDiagnosticsForwardSolverIT_DC(self)
            implicit none
            !
            class( ForwardSolverIT_DC_t ), intent( inout ) :: self
            !
            !
            self%relResVec = R_ZERO
            !
            call self%solver%zeroDiagnostics()
            !
         end subroutine zeroDiagnosticsForwardSolverIT_DC
         !
         !
         subroutine getESolutionForwardSolverIT_DC( self, source, e_solution )
            implicit none
            !
            class( ForwardSolverIT_DC_t ), intent( inout ) :: self
            class( Source_t ), intent( in )                :: source
            class( cVector_t ), intent( inout )            :: e_solution
            !
            class( cVector_t ), allocatable :: temp_esol
            class( cScalar_t ), allocatable :: phi0
            integer :: iter
            complex( kind=prec ) :: i_omega_mu
            !
            !
            call self%solver%zeroDiagnostics()
            !
            self%solver%converged = .FALSE.
            self%solver%failed    = .FALSE.
            self%nDivCor = 0
            !
            if( source%non_zero_source ) then
                !
                select type( grid => self%solver%preconditioner%model_operator%metric%grid )
                    class is( Grid3D_SG_t )
                        !
                        allocate( phi0, source = cScalar3D_SG_t( grid, NODE ) )
                        !
                end select
                !
                call self%divergence_correction%rhsDivCor( self%solver%omega, source, phi0 )
                !
                !e_solution = e_solution%Interior()
                !
                !allocate( temp_esol, source = e_solution )
                !
                !self%nDivCor = self%nDivCor + 1
                !call self%divergence_correction%DivCorr( temp_esol, e_solution, phi0 )
                !
                !deallocate( temp_esol )
                !
            endif
            !
            loop: do while ( ( .NOT. self%solver%converged ) .AND. ( .NOT. self%solver%failed ) )
                !
                !
                select type( solver => self%solver )
                    class is( Solver_QMR_t )
                        call solver%solve( source%rhs, e_solution )
                    class default
                        write( *, * ) "ERROR:ForwardSolverIT_DC::getESolutionForwardSolverIT_DC:"
                        stop        "            Unknow solver type."
                end select
                !
                self%solver%converged = self%solver%n_iter .LT. self%solver%max_iter
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
                self%nDivCor = self%nDivCor + 1
                !
                if( self%nDivCor < self%max_div_cor ) then
                    !
                    allocate( temp_esol, source = e_solution )
                    !
                    if( source%non_zero_source ) then
                        !
                        call self%divergence_correction%DivCorr( temp_esol, e_solution, phi0 )
                        !
                    else
                        !
                        call self%divergence_correction%DivCorr( temp_esol, e_solution )
                        !
                    endif
                    !
                    deallocate( temp_esol )
                    !
                else
                    !
                    self%solver%failed = .TRUE.
                    !
                endif
            !
            enddo loop
            !
            !
            if( source%non_zero_source ) deallocate( phi0 )
            !
            self%relResFinal = self%relResVec( self%n_iter_actual )
            !
            if( source%adjt ) then
                !
                select type( model_operator => self%solver%preconditioner%model_operator )
                    class is ( ModelOperator_MF_t )
                        !
                        e_solution = e_solution * model_operator%Metric%Vedge
                        !
                    class default
                        write( *, * ) "ERROR:ForwardSolverIT_DC_t::getESolutionForwardSolverIT_DC:"
                        stop          "    unknow model_operator type"
                end select
                !
            else
                !
                e_solution = e_solution + source%E%Boundary()
                !
            endif
            !
        end subroutine getESolutionForwardSolverIT_DC
        !
end Module ForwardSolverIT_DC
