module Solver_QMR
    !
    use Solver
    use cVector
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
        write(*,*) "Constructor Solver_QMR_t"
        !
        call self%init()
        !
        self%preconditioner = PreConditioner_MF_CC_t( model_operator )
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
        write(*,*) "Destructor Solver_QMR_t"
        !
        call self%dealloc()
        !
    end subroutine Solver_QMR_dtor
    !
    subroutine setDefaults_QMR( self )
        implicit none
        !
        class( Solver_QMR_t ), intent( inout ) :: self
        !     sets default iteration control parameters for QMR solver
        !     local variables
        integer           :: max_iter
        real( kind=prec ) :: tolerance
        !
        call self%SetParameters( max_iter, tolerance )
        !
    end subroutine setDefaults_QMR
    !
    !
    subroutine solveQMR( self, b, x )
        implicit none
        !
        class( Solver_QMR_t ), intent( inout ) :: self
        class( cVector_t ), intent( in )       :: b
        class( cVector_t ), intent( inout )    :: x
        !
        class( cVector_t ), allocatable :: R, Y, Z, V, W, YT, ZT, VT, WT, P, Q, PT, D, S
        logical              :: adjoint, ilu_adjt
        complex( kind=prec ) :: ETA, PDE, EPSIL, RDE, BETA, DELTA, RHO, DELTA_EPSIL
        complex( kind=prec ) :: PSI, RHO1, GAMM, GAMM1, THET, THET1, TM2
        complex( kind=prec ) :: bnorm,rnorm
        complex( kind=prec ) :: rhoInv,psiInv
        integer              :: iter
        !
        ! Allocate work CVector objects -- questions as in PCG
        allocate( R, source = x )
        !
        call R%zeros() !  can't zero x -- if this is to be used as starting guess
                       !  also, never use AX -- which somehow is declared in ModEM!
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
        ! b - Ax, for inital guess x, that has been input to the routine
        call R%linCombS( b, C_MinusOne, C_ONE )
        !
        ! Norm of rhs, residual
        bnorm = CDSQRT( b%dotProd( b ) )
        rnorm = CDSQRT( R%dotProd( R ) )
        !
        ! this usually means an inadequate model, in which case Maxwell"s fails
        if( isnan( abs( bnorm ) ) ) then
            stop "Error: b in QMR contains NaNs; exiting..."
        end if
        !
        !    iter is iteration counter
        iter = 1
        self%relErr( iter ) = real( rnorm / bnorm )
        !write(*,*) 'in QMR'
        !write(*,*) 'rnorm, bnorm ', rnorm, bnorm
        !write(*,*) 'max_iter, tolerance', self%max_iter,self%tolerance
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
        do while( ( self%relErr( iter ) .gt. self%tolerance ) .AND. ( iter .lt. self%max_iter ) )
            !
            if( ( RHO .eq. C_ZERO ) .or. ( PSI .eq. C_ZERO ) ) then
                !
                self%failed = .TRUE.
                write( *, * ) "QMR FAILED TO CONVERGE : RHO"
                stop "QMR FAILED TO CONVERGE : PSI"
                !
            end if
            !
            rhoInv = ( 1 / RHO )
            psiInv = ( 1 / PSI )
            !
            !    use functions here -- could make subroutines that don"t overwrite
            V = VT%mult( rhoInv )
            W = WT%mult( psiInv )
            !
            !  use subroutines here to overwrite with rescalled vectors
            call Y%multS( rhoInv )
            call Z%multS( psiInv )
            !
            DELTA = Z%dotProd( Y )
            if( DELTA .eq. C_ZERO ) then
                !
                self%failed = .TRUE.
                stop "QMR FAILS TO CONVERGE : DELTA"
                !
            end if
            !
            ilu_adjt = .FALSE.
            call self%preconditioner%UTsolve( Y, YT, ilu_adjt )
            !
            ilu_adjt = .TRUE.
            call self%preconditioner%LTsolve( Z, ZT, ilu_adjt )
            !
            if( iter .eq. 1 ) then
                !
                P = YT 
                Q = ZT 
                !
            else
                ! these calculations are only done when iter > 1
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
            call PT%Zeros()
            call self%preconditioner%model_operator%Amult( self%omega, P, PT, adjoint )
            EPSIL = Q%dotProd( PT )
            !
            if( EPSIL .eq. C_ZERO ) then
                self%failed = .TRUE.
                stop "QMR FAILED TO CONVERGE : EPSIL"
            end if
            !
            BETA = EPSIL/DELTA
            if( BETA .eq. C_ZERO ) then
                self%failed = .TRUE.
                stop "QMR FAILED TO CONVERGE : BETA"
            end if
            !     together these amount to VT = PT-BETA*V
            VT = PT
            !     VT = VT-BETA*V
            call V%scMultAddS( VT, -BETA )
            !
            RHO1 = RHO
            ilu_adjt = .FALSE.
            call self%preconditioner%LTsolve( VT, Y, ilu_adjt )
            RHO = CDSQRT( Y%dotProd( Y ) )
            !
            adjoint = .TRUE.
            call WT%Zeros()
            call self%preconditioner%model_operator%Amult( self%omega, Q, WT, adjoint )
            !
            !     WT = WT - conjg(BETA)*W
            call W%scMultAddS( WT, -conjg( BETA ) )
            !
            ilu_adjt = .TRUE.
            call self%preconditioner%UTsolve( WT, Z, ilu_adjt )
            PSI = CDSQRT( Z%dotProd( Z ) )
            !
            if( iter .gt. 1 ) then
                THET1 = THET
            end if
            !
            THET = RHO / ( GAMM * CDABS( BETA ) )
            GAMM1 = GAMM
            GAMM = C_ONE / CDSQRT( C_ONE + THET * THET )
            !
            if( GAMM .eq. C_ZERO ) then
                self%failed = .TRUE.
                stop "QMR FAILS TO CONVERGE : GAMM"
            end if
            !
            ETA = -ETA * RHO1 * GAMM * GAMM / ( BETA * GAMM1 * GAMM1 )
            !
            if( iter .eq. 1 ) then
                D = P%mult( ETA )     !    using function: D = ETA*P
                S = PT%mult( ETA )     !    using function: S = ETA * PT
            else
                TM2 = THET1 * THET1 * GAMM * GAMM
                call D%linCombS( P, TM2, ETA )     !  D = TM2 * D + ETA * P 
                call S%linCombS( PT, TM2, ETA )    !  S = TM2 * S + ETA * PT 
            end if
            !
            call D%scMultAddS( x, C_ONE )    !  x = x + C_ONE * D
            call S%scMultAddS( R, C_MinusONE )    !  R = R + C_MinusONE * S
            rnorm = CDSQRT( R%dotProd( R ) )
            iter = iter + 1
            !
            ! Keeping track of errors
            ! QMR book-keeping between divergence correction calls
            self%relErr( iter ) = real( rnorm / bnorm )
            !
			!write(*,*) 'iter qmr= ',iter,'    relErr = ', self%relErr(iter)
            !
        end do
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
