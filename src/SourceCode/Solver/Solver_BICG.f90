!
!> Derived class to define a Quasi-Minimal Residue Solver
!
!> Stabilized version of BiConjugate Gradient, set up for solving
!> A x = b using routines in  mult_Aii.
!> solves for the interior (edge) field
!
!> back-ported from the Sparse matrix version, which is modified from my matlab
!> version of BICGstab...
!> so the naming might sound a little different from conventional ones
!
!> interface...........
!> redefining some of the interfaces for our convenience (locally)
!> generic routines for vector operations for edge/ face nodes
!> in a staggered grid
!
!> NOTE: this has not been extensively tested! - I believe it feels a
!> little unstable (despite the name)...
!> if you have time reading this, test it!
!
module Solver_BICG
	!
	use Solver_CC
	use ModelOperator_MF_SG
	use ModelOperator_SP_V1
	use ModelOperator_SP_V2
	use PreConditioner_CC_MF_SG
	use PreConditioner_CC_SP_SG
	use PreConditioner_CC_SP_MR
	!
	type, extends( Solver_CC_t ) :: Solver_BICG_t
	!
	!> No derived properties
	!
	contains
		!
		procedure, public :: solve => solve_Solver_BICG
		!
	end type Solver_BICG_t
	!
	interface Solver_BICG_t
		module procedure Solver_BICG_ctor
	end interface Solver_BICG_t
	!
contains
    !
    !> No subroutine briefing
    !
    function Solver_BICG_ctor( model_operator ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent(in) :: model_operator
        !
        type( Solver_BICG_t ) :: self
        !
        !write( *, * ) "Constructor Solver_BICG_t"
        !
        call self%baseInit
        !
        !> Instantiate the PreConditioner object according to the ModelOperator type
        select type( grid => model_operator%metric%grid )
            !
            class is( Grid3D_SG_t )
                !
                !> Instantiate the PreConditioner object according to the ModelOperator type
                select type( model_operator )
                    !
                    class is( ModelOperator_MF_SG_t )
                        !
                        allocate( self%preconditioner, source = PreConditioner_CC_MF_SG_t( model_operator ) )
                        !
                    class is( ModelOperator_SP_t )
                        !
                        allocate( self%preconditioner, source = PreConditioner_CC_SP_SG_t( model_operator ) )
                        !
                    class default
                        call errStop( "Solver_BICG_ctor > Unclassified SG model_operator" )
                    !
                end select
                !
            class is( Grid3D_MR_t )
                !
                !> Instantiate the PreConditioner object according to the ModelOperator type
                select type( model_operator )
                    !
                    class is( ModelOperator_MF_SG_t )
                        !
                        call errStop( "Solver_BICG_ctor > For MR use model_operator SP" )
                        !
                    class is( ModelOperator_SP_t )
                        !
                        allocate( self%preconditioner, source = PreConditioner_CC_SP_MR_t( model_operator ) )
                        !
                    class default
                        call errStop( "Solver_BICG_ctor > Unclassified MR model_operator" )
                    !
                end select
                !
            class default
                call errStop( "Solver_BICG_ctor > Unclassified grid" )
            !
        end select
        !
        call self%set( max_solver_iters, tolerance_solver )
        !
        call self%zeroDiagnostics
        !
    end function Solver_BICG_ctor
    !
    !> NOTE: this iterative solver is BICG without look-ahead
    !> patterned after the scheme given on page 24 of Barrett et al.
    !> "Templates for the solution of linear systems of equations:
    !> Building blocks for iterative methods"
    !> Note that there are a couple of small differences, due to
    !> the fact that our system is complex (agrees with
    !> matlab6 version of qmr)
    !
    subroutine solve_Solver_BICG( self, b, x )
        implicit none
        !
        class( Solver_BICG_t ), intent(inout) :: self
        class( Vector_t ), intent(in) :: b
        class( Vector_t ), intent(inout) :: x
        !
        class( Vector_t ), allocatable :: R, RT, V, T
        class( Vector_t ), allocatable :: P, PT, PH, S, ST, SH, AX
        class( Vector_t ), allocatable :: xhalf, xmin
        real( kind=prec ) :: rnorm, bnorm, rnormin, btol
        complex( kind=prec ) :: RHO, ALPHA, BETA, OMEGA
        complex( kind=prec ) :: RTV, TT, RHO1
        integer :: iter, imin
        logical :: adjoint, ilu_adjt
        !
        if( .NOT. x%is_allocated ) then
            call errStop( "solve_Solver_BICG > x not allocated yet" )
        endif
        !
        if( .NOT. b%is_allocated ) then
            call errStop( "solve_Solver_BICG > b not allocated yet" )
        endif
        !
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, xhalf )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, xmin)
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, AX )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, R )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, RT )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, P )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, PT )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, PH )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, S )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, ST )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, SH )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, V )
        call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, T )
        !
        self%iter = 1
        !
        !> Norm of rhs, residual
        bnorm = SQRT( b%dotProd( b ) )
        !
        !> this usually means an inadequate model, in which case Maxwell"s fails
        if( isnan( abs( bnorm ) ) ) then
            !
            call errStop( "solve_Solver_BICG > b contains NaNs" )
            !
        elseif( bnorm .EQ. 0.0 ) then ! zero rhs -> zero solution
            !
            call warning( "b in BICG has all zeros, returning zero solution" )
            !
            x = b
            !
            self%iter = 1
            self%n_iter = 1
            self%relErr(1) = 0.0
            !
            return
            !
        endif
        !
        !> now calculate the (original) residual
        adjoint = .FALSE.
        !
        !call self%preconditioned%model_operator%multA_N( x, R, adjoint )
        call self%preconditioner%model_operator%amult( x, R, self%omega, adjoint )
        !
        !> R= b - Ax, for initial guess x, that has been inputted to the routine
        rnorm = CDSQRT( R%dotProd( R ) )
        !
        call R%linComb( b, C_MinusOne, C_ONE )
        !
        !> Norm of residual
        rnorm = CDSQRT( R%dotProd( R ) )
        !
        btol = self%tolerance * bnorm
        !
        if( rnorm .LE. btol ) then ! the first guess is already good enough
            !
            self%n_iter = 1
            !
            call warning( "solve_Solver_BICG > the first guess is already good enough" )
            !
            self%relErr(1) = rnorm / bnorm
            !
            return
            !
        endif
        !
        !> ================= Now start configuring the iteration ===================
        !
        !> the adjoint (shadow) residual
        !
        rnormin = rnorm
        xmin = x
        self%relErr(1) = real( rnormin / bnorm )
        !
        self%converged = .FALSE.
        imin = 1
        RHO = C_ONE
        OMEGA = C_ONE
        RT = R ! use the overloaded =
        !
        !============================== looooops! ================================
        !
        do iter = 2, self%max_iters
            !
            self%iter = iter
            !
            RHO1 = RHO
            RHO = RT%dotProd( R )
            !
            if( RHO .EQ. 0.0 ) then
                !
                call errStop( "solve_Solver_BICG > RHO .EQ. 0.0" )
                !
            endif
            !
            if( self%iter .EQ. 2 ) then
                P = R
            else
                !
                BETA = ( RHO / RHO1 ) * ( ALPHA / OMEGA )
                !
                if( BETA .EQ. 0.0 ) then
                    !
                    call errStop( "solve_Solver_BICG > BETA .EQ. 0.0" )
                    !
                endif
                !
                !> P= R + BETA * (P - OMEGA * V);
                call P%linComb( V, C_One, -OMEGA )
                call P%linComb( R, BETA, C_One )
                !
            endif
            !
            !> L
            ilu_adjt = .FALSE.
            call PT%zeros
            call self%preconditioner%LTsolve( P, PT , ilu_adjt )
            !
            ! U
            ilu_adjt = .FALSE.
            call PH%zeros
            call self%preconditioner%UTsolve( PT, PH , ilu_adjt )
            !
            ! PH = P
            adjoint = .FALSE.
            !
            call V%zeros
            call self%preconditioner%model_operator%amult( PH, V, self%omega, adjoint )
            !
            RTV = RT%dotProd( V )
            !
            if( RTV .EQ. 0.0 ) then
                !
                call errStop( "solve_Solver_BICG > RTV .EQ. 0.0" )
                !
            endif
            !
            ALPHA = RHO / RTV
            !
            if( ALPHA .EQ. 0.0 ) then
                !
                call errStop( "solve_Solver_BICG > ALPHA .EQ. 0.0" )
                !
            endif
            !
            !> xhalf = x + ALPHA*PH ! the first half of iteration
            xhalf = x
            call xhalf%linComb( PH, C_One, ALPHA )
            !
            adjoint = .FALSE.
            !
            call AX%zeros
            call self%preconditioner%model_operator%amult( xhalf, AX, self%omega, adjoint )
            !
            call AX%linComb( b, C_MinusOne, C_One )
            rnorm = CDSQRT( AX%dotProd( AX ) )
            !
            if( rnorm .LT. btol ) then
                !
                x = xhalf
                self%n_iter = self%iter
                self%converged = .TRUE.
                self%relErr( self%iter ) = real( rnorm / bnorm )
                !
                exit
                !
            endif
            !
            if( rnorm .LT. rnormin) then
                !
                rnormin = rnorm
                xmin = xhalf
                imin = self%iter
                !
            endif
            !
            !> S = R - ALPHA*V  !residual for the 0.5 x
            S = R
            call S%linComb( V, C_One, -ALPHA )
            !
            !> L
            ilu_adjt = .FALSE.
            call self%preconditioner%LTsolve( S, ST, ilu_adjt )
            !
            !> U
            ilu_adjt = .FALSE.
            call self%preconditioner%UTsolve( ST, SH, ilu_adjt )
            !
            !> SH = S
            adjoint = .FALSE.
            call T%zeros
            call self%preconditioner%model_operator%amult( SH, T, self%omega, adjoint )
            !
            TT = T%dotProd( T )
            !
            if( TT .EQ. 0.0 ) then
                !
                call errStop( "solve_Solver_BICG > TT .eq. 0.0" )
                !
            endif
            !
            OMEGA = T%dotProd(S) / TT
            !
            if( OMEGA .EQ. 0.0 ) then
                !
                call errStop( "solve_Solver_BICG > OMEGA .eq. 0.0" )
                !
            endif
            !
            !> x = xhalf + OMEGA * SH  ! the second half (shadow) of iteration
            x = xhalf
            call x%linComb( SH, C_One, OMEGA )
            !
            adjoint = .FALSE.
            !
            call AX%zeros
            call self%preconditioner%model_operator%amult( x, AX, self%omega, adjoint )
            !
            call AX%linComb( b, C_MinusOne, C_One )
            !
            rnorm = CDSQRT( AX%dotProd( AX ) )
            !
            self%relErr( self%iter ) = real( rnorm / bnorm )
            !
            if( rnorm .LT. btol ) then
                !
                self%n_iter = self%iter
                self%converged = .TRUE.
                !
                exit
                !
            endif
            !
            if( rnorm .LT. rnormin) then
                !
                rnormin = rnorm
                xmin = x
                imin = self%iter
                !
            endif
            !
            !R = S - OMEGA * T  !residual for the 1.0 x
            R = S
            call R%linComb( T, C_One, -OMEGA )
            !
            !> Verbose
            !write( *, "( a36, i6, a3, es12.3 )" ) "BICG iter: ", self%iter, " : ", self%relErr( self%iter )
            !
        enddo
        !
        if( self%converged ) then
            write( *, "( a52, i6, a7, es12.3 )" ) "->Solver BICG converged within ", self%iter, ": err= ", self%relErr( self%iter )
        else
            write( *, "( a52, i6, a7, es12.3 )" ) "->Solver BICG not converged in ", self%max_iters, ": err= ", self%relErr( self%max_iters )
        endif
        !
        if( .NOT. self%converged ) then
            ! it should be noted that this is the way my matlab version works
            ! the bicg will return the 'best' (smallest residual) iteration
            x = xmin; !comment this line
            self%n_iter = self%max_iters
            self%relErr( self%max_iters ) = self%relErr( imin) ! and this line
            ! to use the last iteration result instead of the 'best'
        endif
        !
        deallocate( xhalf )
        deallocate( xmin)
        deallocate( AX )
        deallocate( R )
        deallocate( RT )
        deallocate( P )
        deallocate( PT )
        deallocate( PH )
        deallocate( S )
        deallocate( ST )
        deallocate( SH )
        deallocate( V )
        deallocate( T )
        !
    end subroutine solve_Solver_BICG
    !
end module Solver_BICG
