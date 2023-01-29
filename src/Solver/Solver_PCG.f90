!
!> Derived class to define a Preconditioned Conjugate-Gradient Solver
!
module Solver_PCG
    !
    use Solver
    use PreConditioner_MF_DC
    !
    !> Solver used only for Divergence Correction
    type, extends( Solver_t ) :: Solver_PCG_t
        !
        !> No derived properties
        !
        contains
            !
            procedure, public :: solve => solvePCG
            procedure, public :: setDefaults => setDefaults_PCG
            !
    end type Solver_PCG_t
    !
    interface Solver_PCG_t
        module procedure Solver_PCG_ctor
    end interface Solver_PCG_t
    !
contains
    !
    !> No function briefing
    function Solver_PCG_ctor( model_operator ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        !
        type( Solver_PCG_t ) :: self
        !
        !write( *, * ) "Constructor Solver_PCG_t"
        !
        call self%init()
        !
        allocate( self%preconditioner, source = PreConditioner_MF_DC_t( model_operator ) )
        !
        call self%setDefaults()
        !
        call self%zeroDiagnostics()
        !
    end function Solver_PCG_ctor
    !
    !> No subroutine briefing
    subroutine setDefaults_PCG( self )
        implicit none
        !
        class( Solver_PCG_t ), intent(inout) :: self
        !
        call self%SetParameters( max_divcor_iters, tolerance_divcor )
        !
    end subroutine setDefaults_PCG
    !
    !> No subroutine briefing
    subroutine solvePCG( self, b, x )
        implicit none
        !
        class( Solver_PCG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: b
        class( Scalar_t ), allocatable, intent( inout ) :: x
        !
        !>    these will have to be created in a way to match
        !>     the specific type of the input Scalar_t ...
        !>    Can we just declare these to be of the abstract type?
        class ( Scalar_t ), allocatable :: r, s, p, q
        complex( kind=prec ) :: beta, alpha, delta, deltaOld
        complex( kind=prec ) :: bnorm, rnorm
        integer :: i
        !
        !>  create local cScalar objects -- could we also use modOp%createCScalar?
        allocate( r, source = x )    !> cannot zero x, since it is first guess
        call r%zeros()
        allocate( s, source = r )
        allocate( p, source = r )
        allocate( q, source = r )
        !
        !> just like
        !call self%preconditioner%model_operator%Amult( x, r )
        !>    if we can make AMult generic, with versions that operate on cScalar/cVector
        !>     this could be more generic ...    also could change the name of this operator
        !>      to make this more obvious
        call self%preconditioner%model_operator%divCgrad( x, r )
        !
        !>     r = b-r
        call r%linComb( b, C_MinusOne, C_ONE )
        !
        bnorm = SQRT( b%dotProd(b) )
        rnorm = SQRT( r%dotProd(r) )
        !
        self%relErr(1) = rnorm/bnorm
        !
        i = 0
        !
        loop: do while ( ( self%relErr( i + 1 ) .GT. self%tolerance ).AND.( i + 1 .LT. self%max_inv_iters ) )
            !
            call self%preconditioner%LUsolve( r, s )
            !
            delta = r%dotProd(s)
            if( i .EQ. 0 ) then
                beta = C_ZERO
            else
                beta = delta / deltaOld
            endif
            !
            call p%linComb( s, beta, C_ONE )
            !
            call q%zeros()
            call self%preconditioner%model_operator%divCgrad( p, q )
            !
            alpha = delta / p%dotProd(q)
            !
            call x%multAdd( alpha, p )
            !
            call r%multAdd( -alpha, q )
            !
            deltaOld = delta
            !
            i = i + 1
            !
            rnorm = SQRT( r%dotProd(r) )
            !
            self%relErr( i + 1 ) = rnorm / bnorm
            !
            !write( *, * ) "PCG iter, self%relErr( i + 1 )", i + 1, self%relErr( i + 1 )
            !
        enddo loop
        ! !
        ! if( i + 1 .LT. self%max_inv_iters ) then
            ! write( *, * ) "                    divCorr PCG converged within ", i + 1, " : ", self%relErr( i + 1 )
        ! else
            ! write( *, * ) "                    divCorr PCG not converged in ", i + 1, " : ", self%relErr( i + 1 )
        ! endif
        ! !
        deallocate( r )
        deallocate( s )
        deallocate( p )
        deallocate( q )
        !
        self%n_inv_iter = i
        !
    end subroutine solvePCG !> PCG
    !
end module Solver_PCG
