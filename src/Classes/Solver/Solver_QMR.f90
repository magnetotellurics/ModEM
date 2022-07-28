module Solver_QMR
    !
    use Solver
    use ModelOperator
    use PreConditioner_MF_CC
    !
    type, extends( Solver_t ) :: Solver_QMR_t
        !
        ! PROPERTIES HERE
        !
        contains
            !
            final :: Solver_QMR_dtor
            !
            procedure, public :: solve => solveQMR
            procedure, public :: setDefaults => setDefaults_QMR
            !
     end type Solver_QMR_t
     !
     interface Solver_QMR_t
         module procedure Solver_QMR_ctor
     end interface Solver_QMR_t
     !
contains
    !
    function Solver_QMR_ctor( model_operator ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        type( Solver_QMR_t ) :: self
        !
        !write( *, * ) "Constructor Solver_QMR_t"
        !
        call self%init()
        !
        allocate( self%preconditioner, source = PreConditioner_MF_CC_t( model_operator ) )
        !
        call self%setDefaults()
        !
        call self%zeroDiagnostics()
        !
    end function Solver_QMR_ctor
    !
    ! Solver_QMR destructor
    subroutine Solver_QMR_dtor( self )
        implicit none
        !
        type( Solver_QMR_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor Solver_QMR_t"
        !
        call self%dealloc()
        !
    end subroutine Solver_QMR_dtor
    !
    subroutine setDefaults_QMR( self )
        implicit none
        !
        class( Solver_QMR_t ), intent( inout ) :: self
        !
        call self%SetParameters( QMR_iters, tolerance_qmr )
        !
    end subroutine setDefaults_QMR
    !
    !
    subroutine solveQMR( self, b, x )
        implicit none
        !
        class( Solver_QMR_t ), intent( inout ) :: self
        class( Vector_t ), intent( in )        :: b
        class( Vector_t ), intent( inout )     :: x
        !
        class( Vector_t ), allocatable :: R, Y, Z, V, W, YT, ZT, VT, WT, P, Q, PT, D, S
        logical              :: adjoint, ilu_adjt
        complex( kind=prec ) :: ETA, PDE, EPSIL, RDE, BETA, DELTA, RHO, DELTA_EPSIL
        complex( kind=prec ) :: PSI, RHO1, GAMM, GAMM1, THET, THET1, TM2
        complex( kind=prec ) :: bnorm, rnorm
        complex( kind=prec ) :: rhoInv, psiInv
        integer              :: iter
        !
        ! Allocate work CVector objects -- questions as in PCG
        allocate( R, source = x )
        !
        call R%zeros() !  can"t zero x -- if this is to be used as starting guess
                       !  also, never use AX -- which somehow is declared in ModEM!
        !
        allocate( Y, source = R )
        allocate( Z, source = R )
        allocate( V, source = R )
        allocate( W, source = R )
        allocate( YT, source = R )
        allocate( ZT, source = R )
        allocate( VT, source = R )
        allocate( WT, source = R )
        allocate( P, source = R )
        allocate( Q, source = R )
        allocate( PT, source = R )
        allocate( D, source = R )
        allocate( S, source = R )
        !
        ! NOTE: this iterative solver is QMR without look-ahead
        ! patterned after the scheme given on page 24 of Barrett et al.
        ! "Templates for the solution of linear systems of equations:
        ! Building blocks for iterative methods"
        ! Note that there are a couple of small differences, due to
        ! the fact that our system is complex (agrees with
        ! matlab6 version of qmr)
        !
        self%failed = .FALSE.
        adjoint = .FALSE.
        ! R is Ax
        !
        call self%preconditioner%model_operator%Amult( self%omega, x, R, adjoint )
        ! b - Ax, for initial guess x, that has been input to the routine
        call R%linCombS( b, C_MinusOne, C_ONE )
        !
        ! Norm of rhs, residual
        bnorm = CDSQRT( b%dotProd( b ) )
        !
        ! this usually means an inadequate model, in which case Maxwell"s fails
        if( isnan( abs( bnorm ) ) ) then
            stop "Error: b in QMR contains NaNs; exiting..."
        end if
        !
        rnorm = CDSQRT( R%dotProd( R ) )
        !
        ! Initial guess relative error
        iter = 1
        self%relErr( iter ) = real( rnorm / bnorm )
        !
        VT = R 
        ilu_adjt = .FALSE.
        call self%preconditioner%LTsolve( VT, Y, ilu_adjt )
        RHO = CDSQRT( Y%dotProd( Y ) )
        !
        WT = R 
        ilu_adjt = .TRUE.
        call self%preconditioner%UTsolve( WT, Z, ilu_adjt )
        PSI  = CDSQRT( Z%dotProd( Z ) )
        GAMM = C_ONE
        ETA  = C_MinusONE
        !
        ! the do loop goes on while the relative error is greater than the tolerance
        ! and the iterations are less than maxIt
        do while( ( self%relErr( iter ) .GT. self%tolerance ) .AND. ( iter .LT. self%max_iter ) )
            !
            ! Verbosis
            !write( *, * ) "QMR iter, self%relErr( iter )", iter, self%relErr( iter )
            !
            if( ( RHO .EQ. C_ZERO ) .OR. ( PSI .EQ. C_ZERO ) ) then
                !
                self%failed = .TRUE.
                write( *, * ) "QMR FAILED TO CONVERGE : RHO"
                stop "Error: QMR FAILED TO CONVERGE : PSI"
                !
            end if
            !
            rhoInv = ( 1 / RHO )
            psiInv = ( 1 / PSI )
            !
            !    use functions here -- could make subroutines that don"t overwrite
            V = VT
            call V%mult( rhoInv )
            !
            W = WT
            call W%mult( psiInv )
            !
            !  use subroutines here to overwrite with rescalled vectors
            call Y%mult( rhoInv )
            call Z%mult( psiInv )
            !
            DELTA = Z%dotProd( Y )
            if( DELTA .EQ. C_ZERO ) then
                !
                self%failed = .TRUE.
                stop "Error: QMR FAILS TO CONVERGE : DELTA"
                !
            end if
            !
            ilu_adjt = .FALSE.
            call self%preconditioner%UTsolve( Y, YT, ilu_adjt )
            !
            ilu_adjt = .TRUE.
            call self%preconditioner%LTsolve( Z, ZT, ilu_adjt )
            !
            if( iter .EQ. 1 ) then
                !
                P = YT 
                Q = ZT 
                !
            else
                !
                DELTA_EPSIL = DELTA / EPSIL
                PDE = -PSI * DELTA_EPSIL
                RDE = -RHO * CONJG( DELTA_EPSIL )
                !
                call P%linCombS( YT, PDE, C_ONE )
                !
                call Q%linCombS( ZT, RDE, C_ONE )
                !
            end if
            !
            adjoint = .FALSE.
            !
            call PT%Zeros()
            call self%preconditioner%model_operator%Amult( self%omega, P, PT, adjoint )
            EPSIL = Q%dotProd( PT )
            !
            if( EPSIL .EQ. C_ZERO ) then
                self%failed = .TRUE.
                stop "Error: QMR FAILED TO CONVERGE : EPSIL"
            end if
            !
            BETA = EPSIL/DELTA
            if( BETA .EQ. C_ZERO ) then
                self%failed = .TRUE.
                stop "Error: QMR FAILED TO CONVERGE : BETA"
            end if
            !
            VT = PT
            !
            call V%scMultAddS( VT, -BETA ) !  VT = VT - BETA * V
            !
            RHO1 = RHO
            ilu_adjt = .FALSE.
            call self%preconditioner%LTsolve( VT, Y, ilu_adjt )
            RHO = CDSQRT( Y%dotProd( Y ) )
            !
            adjoint = .TRUE.
            !
            call WT%Zeros()
            call self%preconditioner%model_operator%Amult( self%omega, Q, WT, adjoint )
            !
            call W%scMultAddS( WT, -conjg( BETA ) ) !  WT = WT - conjg(BETA) * W
            !
            ilu_adjt = .TRUE.
            call self%preconditioner%UTsolve( WT, Z, ilu_adjt )
            PSI = CDSQRT( Z%dotProd( Z ) )
            !
            if( iter .GT. 1 ) then
                THET1 = THET
            end if
            !
            THET = RHO / ( GAMM * CDABS( BETA ) )
            GAMM1 = GAMM
            GAMM = C_ONE / CDSQRT( C_ONE + THET * THET )
            !
            if( GAMM .EQ. C_ZERO ) then
                self%failed = .TRUE.
                stop "Error: QMR FAILS TO CONVERGE : GAMM"
            end if
            !
            ETA = -ETA * RHO1 * GAMM * GAMM / ( BETA * GAMM1 * GAMM1 )
            !
            if( iter .EQ. 1 ) then
                D = P
                call D%mult( ETA )  !  D = ETA*P
                !
                S = PT
                call S%mult( ETA ) !  S = ETA * PT
            else
                TM2 = THET1 * THET1 * GAMM * GAMM
                call D%linCombS( P, TM2, ETA )  !  D = TM2 * D + ETA * P 
                call S%linCombS( PT, TM2, ETA ) !  S = TM2 * S + ETA * PT 
            end if
            !
            call D%scMultAddS( x, C_ONE )      !  x = x + C_ONE * D
            call S%scMultAddS( R, C_MinusONE ) !  R = R + C_MinusONE * S
            !
            rnorm = CDSQRT( R%dotProd( R ) )
            !
            iter = iter + 1
            !
            self%relErr( iter ) = real( rnorm / bnorm )
            !
        end do
        !
        if( iter .LT. self%max_iter ) then
            write( *, * ) "               Solver QMR converged within ", iter, " : ", self%relErr( iter )
        else
            write( *, * ) "               Solver QMR not converged in ", iter, " : ", self%relErr( iter )
        endif
        !
        deallocate( R )
        deallocate( Y )
        deallocate( Z )
        deallocate( V )
        deallocate( W )
        deallocate( YT )
        deallocate( ZT )
        deallocate( VT )
        deallocate( WT )
        deallocate( P )
        deallocate( Q )
        deallocate( PT )
        deallocate( D )
        deallocate( S )
        !
        self%n_iter = iter
        !
    end subroutine solveQMR
    !
end module Solver_QMR
