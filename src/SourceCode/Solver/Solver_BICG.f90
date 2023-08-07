!
!> Derived class to define a Quasi-Minimal Residue Solver
!
!> Stablized version of BiConjugate Gradient, set up for solving
!> A x = b using routines in  mult_Aii.
!> solves for the interior (edge) field
!
!> backported from the Sparse matrix version, which is modified from my matlab
!> version of BICGstab...
!> so the naming might sound a little different from conventional ones
!
!> interface...........
!> redefining some of the interfaces for our convenience (locally)
!> generic routines for vector operations for edge/ face nodes
!> in a staggered grid
!
!> NOTE: this has not been extensively tested! - I believe it feels a
!> little unstable (dispite the name)...
!> if you have time reading this, test it!
!
module Solver_BICG
    !
    use Solver
    use ModelOperator_MF_SG
    use ModelOperator_SP
    use PreConditioner_CC_MF
    use PreConditioner_CC_SP
    !
    type, extends( Solver_t ) :: Solver_BICG_t
        !
        !> No derived properties
        !
        contains
            !
            procedure, public :: solve => solveBICG
            !
            procedure, public :: setDefaults => setDefaults_BICG
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
        class( ModelOperator_t ), intent( in ) :: model_operator
        !
        type( Solver_BICG_t ) :: self
        !
        !write( *, * ) "Constructor Solver_BICG_t"
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
                call errStop( "Solver_BICG_ctor > Unclassified ModelOperator" )
            !
        end select
        !
        call self%setDefaults
        !
        call self%zeroDiagnostics
        !
    end function Solver_BICG_ctor
    !
    !> No subroutine briefing
    !
    subroutine setDefaults_BICG( self )
        implicit none
        !
        class( Solver_BICG_t ), intent( inout ) :: self
        !
        call self%setParameters( max_solver_iters, tolerance_solver )
        !
    end subroutine setDefaults_BICG
    !
    !> NOTE: this iterative solver is BICG without look-ahead
    !> patterned after the scheme given on page 24 of Barrett et al.
    !> "Templates for the solution of linear systems of equations:
    !> Building blocks for iterative methods"
    !> Note that there are a couple of small differences, due to
    !> the fact that our system is complex (agrees with
    !> matlab6 version of qmr)
    !
    subroutine solveBICG( self, b, x )
        implicit none
        !
        class( Solver_BICG_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: b
        class( Vector_t ), intent( inout ) :: x
        !
        class( Vector_t ), allocatable :: R, RT, V, T
        class( Vector_t ), allocatable :: P, PT, PH, S, ST, SH, AX
        class( Vector_t ), allocatable :: xhalf, xmin
		real( kind=prec ) :: rnorm, bnorm, rnormin, btol
		complex( kind=prec ) :: RHO, ALPHA, BETA, OMEGA
		complex( kind=prec ) :: RTV, TT, RHO1
		integer :: iter, imin
		integer :: maxiter
		logical :: adjoint, ilu_adjt, converged
        !
        if( .NOT. x%is_allocated ) then
            call errStop( "solveBICG > x not allocated yet" )
        endif
        !
        if( .NOT. b%is_allocated ) then
            call errStop( "solveBICG > b not allocated yet" )
        endif
		!
		call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, xhalf )
		call self%preconditioner%model_operator%metric%createVector( complex_t, x%grid_type, xmin )
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
        !> Norm of rhs, residual
        bnorm = SQRT( b%dotProd( b ) )
        !
        !> this usually means an inadequate model, in which case Maxwell"s fails
        if( isnan( real( abs( bnorm ), kind=prec ) ) ) then
            call errStop( "solveQMR > b in QMR contains NaNs" )
		elseif( bnorm .EQ. 0.0 ) then ! zero rhs -> zero solution
			!
			call warning( "b in BICG has all zeros, returning zero solution" )
			x = b
			BICGiter%niter=1
			BICGiter%failed=.false.
			BICGiter%rerr=0.0
			return
			!
		endif
		!
		!> now calculate the (original) residual
		adjoint = .FALSE.
		!
		call A( x, adjoint, R )
		! R= b - Ax, for inital guess x, that has been inputted to the routine
		rnorm = CDSQRT(dotProd(R, R))
		Call linComb(C_ONE,b,C_MinusOne,R,R)
		! Norm of residual
		rnorm = CDSQRT(dotProd(R, R))
		btol = BICGiter%tol * bnorm
		if ( rnorm .le. btol ) then ! the first guess is already good enough
		! returning
		BICGiter%niter=1
		BICGiter%failed=.false.
		BICGiter%rerr(1)=rnorm/bnorm
		return
		end if
		!================= Now start configuring the iteration ===================!
		! the adjoint (shadow) residual
		rnormin = rnorm
		xmin = x
		BICGiter%rerr(1) = real(rnormin/bnorm)
		write(6,*) 'initial residual: ', BICGiter%rerr(1)
		converged = .false.
		maxiter = BICGiter%maxit
		imin = 1
		RHO = C_ONE
		OMEGA = C_ONE
		RT = R ! use the overloaded =
		!============================== looooops! ================================!
		do iter= 2, maxiter
		RHO1 = RHO
		RHO = dotProd(RT, R)
		if (RHO .eq. 0.0) then
		BICGiter%failed = .true.
		exit
		end if
		if (iter .eq. 2) then
		P = R
		else
		BETA = (RHO/RHO1)*(ALPHA/OMEGA)
		if (BETA .eq. 0.0) then
		BICGiter%failed = .true.
		exit
		end if
		! P= R + BETA * (P - OMEGA * V);
		Call linComb(C_One,P,-OMEGA,V,P)
		call linComb(C_One,R,BETA,P,P)
		end if
		! L
		ilu_adjt = .false.
		call M1solve(P,ilu_adjt,PT)
		! U
		ilu_adjt = .false.
		call M2solve(PT,ilu_adjt,PH)
		!      PH = P
		adjoint = .false.
		call zero(V)
		call A(PH,adjoint,V)
		RTV = dotProd(RT, V)
		if (RTV.eq.0.0) then
		BICGiter%failed = .true.
		exit
		end if
		ALPHA = RHO / RTV
		if (ALPHA.eq.0.0) then
		BICGiter%failed = .true.
		exit
		end if
		! xhalf = x + ALPHA*PH ! the first half of iteration
		call linComb(C_One,x,ALPHA,PH,xhalf)
		adjoint = .false.
		call zero(AX)
		call A(xhalf,adjoint,AX)
		call linComb(C_One,b,C_MinusOne,AX,AX)
		rnorm = CDSQRT(dotProd(AX,AX))


		if (rnorm.lt.btol) then
		x = xhalf
		BICGiter%failed = .false.
		BICGiter%niter = iter
		converged = .true.
		BICGiter%rerr(iter)=real(rnorm/bnorm)
		exit
		end if
		if (rnorm .lt. rnormin) then
		rnormin = rnorm
		xmin = xhalf
		imin = iter
		end if
		! S = R - ALPHA*V  !residual for the 0.5 x
		call linComb(C_One,R,-ALPHA,V,S)
		! L
		ilu_adjt = .false.
		call M1solve(S,ilu_adjt,ST)
		! U
		ilu_adjt = .false.
		call M2solve(ST,ilu_adjt,SH)
		!     SH = S
		adjoint = .false.
		call zero(T)
		call A(SH,adjoint,T)
		TT = dotProd(T,T)
		if (TT.eq.0.0) then
		BICGiter%failed = .true.
		exit
		end if
		OMEGA = dotProd(T,S)/TT
		if (OMEGA.eq.0.0) then
		BICGiter%failed = .true.
		exit
		end if
		! x = xhalf + OMEGA * SH  ! the second half (shadow) of iteration
		call linComb(C_One,xhalf,OMEGA,SH,x)
		adjoint = .false.
		call zero(AX)
		call A(x,adjoint,AX)
		call linComb(C_One,b,C_MinusOne,AX,AX)
		rnorm = CDSQRT(dotProd(AX,AX))
		BICGiter%rerr(iter) = real(rnorm / bnorm)
		if (rnorm.lt.btol) then
		BICGiter%failed = .false.
		BICGiter%niter = iter
		converged = .true.
		exit
		end if
		if (rnorm .lt. rnormin) then
		rnormin = rnorm
		xmin = x
		imin = iter
		end if
		!R = S - OMEGA * T  !residual for the 1.0 x
		call linComb(C_One,S,-OMEGA,T,R)
		end do

		if (.not. converged) then
		! it should be noted that this is the way my matlab version works
		! the bicg will return the 'best' (smallest residual) iteration
		x = xmin; !comment this line
		BICGiter%niter=BICGiter%maxit
		BICGiter%rerr(BICGiter%maxit) = BICGiter%rerr(imin) ! and this line
		! to use the last iteration result instead of the 'best'
		end if
		Call deall(xhalf)
		Call deall(xmin)
		Call deall(AX)
		Call deall(R)
		Call deall(RT)
		Call deall(P)
		Call deall(PT)
		Call deall(PH)
		Call deall(S)
		Call deall(ST)
		Call deall(SH)
		Call deall(V)
		Call deall(T)

        self%n_iter = self%iter
        !
    end subroutine solveBICG
    !
end module Solver_BICG
