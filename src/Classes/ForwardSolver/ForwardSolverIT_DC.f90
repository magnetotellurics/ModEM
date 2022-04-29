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
    integer, parameter :: iter_per_div_corDefQMR = 40
    !
    integer, parameter :: iter_per_div_corDefBCG = 80
    !
    integer, parameter :: max_div_corDef = 20
    !
    integer, parameter :: max_iterDivCorDef = 100
    !
    real( kind=prec ), parameter :: tolDivCorDef = 1E-5
    !
    real( kind=prec ), parameter :: tolCurlCurlDef = 1E-7
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
            procedure, public :: setCond         => setCondForwardSolverIT_DC
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
            integer :: max_iter, max_iter_total
            !
            write(*,*) "Constructor ForwardSolverIT_DC_t"
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
                    self%solver = Solver_QMR_t( model_operator )
                    max_iter = iter_per_div_corDefQMR
                case( BiCG )
                    max_iter = iter_per_div_corDefBCG
                    stop "ForwardSolverIT_DC_ctor: Not yet coded for Bi-Conjugate Gradients"
                case default
                    stop "ForwardSolverIT_DC_ctor: Unknow solver"
            end select
            !
            call self%solver%setParameters( max_iter, tolCurlCurlDef )
            !
            call self%setIterDefaultsDC()
            !
            max_iter_total = self%max_div_cor * self%solver%max_iter
            !
            call self%setIterControl( max_iter_total, self%solver%tolerance )
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
            write(*,*) "Destructor ForwardSolverIT_DC_t"
            !
            call self%dealloc()
            !
        end subroutine ForwardSolverIT_DC_dtor
        !
        !
        subroutine setFrequencyForwardSolverIT_DC( self, period )
            implicit none
            !
            class( ForwardSolverIT_DC_t ), intent( inout ) :: self
            real( kind=prec ), intent( in )                :: period
            !
            self%solver%omega = 2.0 * PI / period
            !
            call self%solver%preconditioner%SetPreconditioner( self%solver%omega )
            !
        end subroutine setFrequencyForwardSolverIT_DC
        !
        !
        subroutine setCondForwardSolverIT_DC( self, model_parameter )
            implicit none
            !
            class( ForwardSolverIT_DC_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in )        :: model_parameter
            !
            !
            call self%solver%preconditioner%model_operator%setCond( model_parameter )
            !
            call self%divergence_correction%SetCond()
            !
            call self%solver%preconditioner%SetPreconditioner( self%solver%omega )
            !
        end subroutine setCondForwardSolverIT_DC
        !
        !
        subroutine setIterControlForwardSolverIT_DC( self, maxit, tolerance )
            implicit none
            !
            class( ForwardSolverIT_DC_t ), intent( inout ) :: self
            integer, intent( in )                          :: maxit
            real( kind=prec ), intent( in )                :: tolerance
            !
            integer :: iter_per_dc
            !
            !
            self%tolerance = tolerance
            !
            iter_per_dc = self%solver%max_iter
            !
            self%max_div_cor = maxit / iter_per_dc
            !
            self%max_iter_total = iter_per_dc * self%max_div_cor
            !
            call self%solver%setParameters( iter_per_dc, tolerance )
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
            self%relResFinal = R_ZERO
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
            class( cVector_t ), allocatable :: temp
            class( cScalar_t ), allocatable :: phi0
            integer :: iter
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
            endif
            !
            loop: do while ( ( .NOT. self%solver%converged ) .AND. ( .NOT. self%solver%failed ) )
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
                self%nDivCor = self%nDivCor+1
                !
                if( self%nDivCor < self%max_div_cor ) then
                    !
                    allocate( temp, source = e_solution )
                    !
                    if( source%non_zero_source ) then
                        !
                        call self%divergence_correction%DivCorr( temp, e_solution, phi0 )
                        !
                    else
                        !
                        call self%divergence_correction%DivCorr( temp, e_solution )
                        !
                    endif
                    !
                    deallocate( temp )
                    !
                else
                    !
                    self%solver%failed = .TRUE.
                    !
                endif
            !
            enddo loop
            !
            if( source%non_zero_source ) deallocate( phi0 )
            !
            self%relResFinal = self%relResVec(self%n_iter_actual)
            !
            select type( modOp => self%solver%preconditioner%model_operator )
                class is ( ModelOperator_MF_t )
                    if( source%adjt ) then
                        !
                        e_solution = e_solution * modOp%Metric%Vedge
                        !
                    else
                        !
                        e_solution = e_solution + source%E%Boundary()
                        !
                    endif
                    !
                class default
                    write( *, * ) "ERROR:ForwardSolverIT_DC_t::getESolutionForwardSolverIT_DC:"
                    stop        "    model_operator type unknow"
            end select
            !
        end subroutine getESolutionForwardSolverIT_DC
        !
end Module ForwardSolverIT_DC
