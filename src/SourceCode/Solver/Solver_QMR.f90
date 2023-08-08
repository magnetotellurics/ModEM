!
!> Derived class to define a Quasi-Minimal Residue Solver
!
module Solver_QMR
    !
    use Solver_CC
    use ModelOperator_MF_SG
    use ModelOperator_SP
    use PreConditioner_CC_MF
    use PreConditioner_CC_SP
    !
    type, extends( Solver_CC_t ) :: Solver_QMR_t
        !
        !> No derived properties
        !
        contains
            !
            procedure, public :: solve => solve_Solver_QMR
            !
     end type Solver_QMR_t
     !
     interface Solver_QMR_t
         module procedure Solver_QMR_ctor
     end interface Solver_QMR_t
     !
contains
    !
    !> No subroutine briefing
    !
    function Solver_QMR_ctor( model_operator ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        !
        type( Solver_QMR_t ) :: self
        !
        !write( *, * ) "Constructor Solver_QMR_t"
        !
        call self%baseInit
        !
        !> Instantiate the PreConditioner object according to the ModelOperator type
        select type( model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                allocate( self%preconditioner, source = PreConditioner_CC_MF_t( model_operator ) )
                !
            class is( ModelOperator_SP_t )
                !
                allocate( self%preconditioner, source = PreConditioner_CC_SP_t( model_operator ) )
                !
            class default
                call errStop( "Solver_QMR_ctor > Unclassified ModelOperator" )
            !
        end select
        !
        call self%set( max_solver_iters, tolerance_solver )
        !
        call self%zeroDiagnostics
        !
    end function Solver_QMR_ctor
    !
    !> NOTE: this iterative solver is QMR without look-ahead
    !> patterned after the scheme given on page 24 of Barrett et al.
    !> "Templates for the solution of linear systems of equations:
    !> Building blocks for iterative methods"
    !> Note that there are a couple of small differences, due to
    !> the fact that our system is complex (agrees with
    !> matlab6 version of qmr)
    !
    subroutine solve_Solver_QMR( self, b, x )
        implicit none
        !
        class( Solver_QMR_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: b
        class( Vector_t ), intent( inout ) :: x
        !
        class( Vector_t ), allocatable :: R, Y, Z, V, W, YT, ZT, VT, WT, P, Q, PT, D, S
        logical :: adjoint, ilu_adjoint
        complex( kind=prec ) :: ETA, PDE, EPSIL, RDE, BETA, DELTA, RHO, DELTA_EPSIL
        complex( kind=prec ) :: PSI, RHO1, GAMM, GAMM1, THET, THET1, TM2
        complex( kind=prec ) :: bnorm, rnorm
        complex( kind=prec ) :: rhoInv, psiInv
        !
        if( .NOT. x%is_allocated ) then
            call errStop( "solve_Solver_QMR > x not allocated yet" )
        endif
        !
        if( .NOT. b%is_allocated ) then
            call errStop( "solve_Solver_QMR > b not allocated yet" )
        endif
        !
        !> Allocate work Vector objects -- questions as in PCG
        allocate( R, source = x )
        !
        !> can"t zero x -- if this is to be used as starting guess
        !> also, never use AX -- which somehow is declared in ModEM!
        call R%zeros
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
        adjoint = .FALSE.
        !
        !> R is Ax
        call self%preconditioner%model_operator%amult( x, R, self%omega, adjoint )
        !
        !> b - Ax, for initial guess x, that has been input to the routine
        call R%linComb( b, C_MinusOne, C_ONE )
        !
        !> Norm of rhs, residual
        bnorm = SQRT( b%dotProd( b ) )
        !
        !> this usually means an inadequate model, in which case Maxwell"s fails
        if( isnan( real( abs( bnorm ), kind=prec ) ) ) then
            call errStop( "solve_Solver_QMR > b in QMR contains NaNs" )
        endif
        !
        rnorm = SQRT( R%dotProd( R ) )
        !
        !> Initial guess relative error
        self%iter = 1
        self%relErr( self%iter ) = real( rnorm / bnorm )
        !
        VT = R 
        ilu_adjoint = .FALSE.
        call self%preconditioner%LTsolve( VT, Y, ilu_adjoint )
        RHO = SQRT( Y%dotProd( Y ) )
        !
        WT = R 
        ilu_adjoint = .TRUE.
        call self%preconditioner%UTsolve( WT, Z, ilu_adjoint )
        PSI = SQRT( Z%dotProd( Z ) )
        GAMM = C_ONE
        ETA = C_MinusONE
        !
        !> the do loop goes on while the relative error is greater than the tolerance
        !> and the iterations are less than maxIt
        solver_loop: do
            !
            if( ( RHO .EQ. C_ZERO ) .OR. ( PSI .EQ. C_ZERO ) ) then
                call errStop( "solve_Solver_QMR > Failed to converge" )
            endif
            !
            rhoInv = ( 1 / RHO )
            psiInv = ( 1 / PSI )
            !
            !> use functions here -- could make subroutines that don"t overwrite
            V = VT
            call V%mult( rhoInv )
            !
            W = WT
            call W%mult( psiInv )
            !
            !> use subroutines here to overwrite with rescalled vectors
            call Y%mult( rhoInv )
            call Z%mult( psiInv )
            !
            DELTA = Z%dotProd( Y )
            !
            if( DELTA .EQ. C_ZERO ) then
                !
                call errStop( "solve_Solver_QMR > DELTA fails to converge" )
                !
            endif
            !
            ilu_adjoint = .FALSE.
            call self%preconditioner%UTsolve( Y, YT, ilu_adjoint )
            !
            ilu_adjoint = .TRUE.
            call self%preconditioner%LTsolve( Z, ZT, ilu_adjoint )
            !
            if( self%iter .EQ. 1 ) then
                !
                P = YT 
                Q = ZT 
                !
            else
                !
                DELTA_EPSIL = DELTA / EPSIL
                PDE = -PSI * DELTA_EPSIL
                RDE = -RHO * conjg( DELTA_EPSIL )
                !
                call P%linComb( YT, PDE, C_ONE )
                !
                call Q%linComb( ZT, RDE, C_ONE )
                !
            endif
            !
            adjoint = .FALSE.
            !
            call PT%zeros
            call self%preconditioner%model_operator%amult( P, PT, self%omega, adjoint )
            EPSIL = Q%dotProd( PT )
            !
            if( EPSIL .EQ. C_ZERO ) then
                call errStop( "solve_Solver_QMR > EPSIL failed to converge" )
            endif
            !
            BETA = EPSIL/DELTA
            !
            if( BETA .EQ. C_ZERO ) then
                call errStop( "solve_Solver_QMR > BETA failed to converge" )
            endif
            !
            VT = PT
            !
            call VT%multAdd( -BETA, V ) !>  VT = VT - BETA * V
            !
            RHO1 = RHO
            ilu_adjoint = .FALSE.
            call self%preconditioner%LTsolve( VT, Y, ilu_adjoint )
            RHO = SQRT( Y%dotProd( Y ) )
            !
            adjoint = .TRUE.
            !
            call WT%zeros
            call self%preconditioner%model_operator%amult( Q, WT, self%omega, adjoint )
            !
            call WT%multAdd( -conjg( BETA ), W ) !>  WT = WT - conjg(BETA) * W
            !
            ilu_adjoint = .TRUE.
            call self%preconditioner%UTsolve( WT, Z, ilu_adjoint )
            PSI = SQRT( Z%dotProd( Z ) )
            !
            if( self%iter .GT. 1 ) then
                THET1 = THET
            endif
            !
            THET = RHO / ( GAMM * ABS( BETA ) )
            GAMM1 = GAMM
            GAMM = C_ONE / SQRT( C_ONE + THET * THET )
            !
            if( GAMM .EQ. C_ZERO ) then
                call errStop( "solve_Solver_QMR > GAMM fails to converge" )
            endif
            !
            ETA = -ETA * RHO1 * GAMM * GAMM / ( BETA * GAMM1 * GAMM1 )
            !
            if( self%iter .EQ. 1 ) then
                D = P
                call D%mult( ETA ) !> D = ETA*P
                !
                S = PT
                call S%mult( ETA ) !> S = ETA * PT
            else
                TM2 = THET1 * THET1 * GAMM * GAMM
                call D%linComb( P, TM2, ETA )  !> D = TM2 * D + ETA * P 
                call S%linComb( PT, TM2, ETA ) !> S = TM2 * S + ETA * PT 
            endif
            !
            call x%multAdd( C_ONE, D )      !>  x = x + C_ONE * D
            call R%multAdd( C_MinusONE, S ) !>  R = R + C_MinusONE * S
            !
            rnorm = SQRT( R%dotProd( R ) )
            !
            !> Verbose
            !write( *, "( a36, i6, a3, es12.3 )" ) "QMR iter: ", self%iter, " : ", self%relErr( self%iter )
            !
            self%iter = self%iter + 1
            !
            self%relErr( self%iter ) = real( rnorm / bnorm, kind=prec )
            !
            !> Stop Conditions
            if( ( self%relErr( self%iter ) .LE. self%tolerance ) .OR. ( self%iter .GE. self%max_iters ) ) then
                exit
            endif
            !
        enddo solver_loop
        !
        self%converged = self%relErr( self%iter ) .LE. self%tolerance
        !
        if( self%converged ) then
            write( *, "( a51, i6, a7, es12.3 )" ) "->Solver QMR converged within ", self%iter, ": err= ", self%relErr( self%iter )
        else
            write( *, "( a51, i6, a7, es12.3 )" ) "->Solver QMR not converged in ", self%max_iters, ": err= ", self%relErr( self%max_iters )
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
        self%n_iter = self%iter
        !
    end subroutine solve_Solver_QMR
    !
end module Solver_QMR
