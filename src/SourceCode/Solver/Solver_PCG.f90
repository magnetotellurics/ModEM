!
!> Derived class to define a Preconditioned Conjugate-Gradient Solver
!
module Solver_PCG
    !
    use Solver
    use ModelOperator_MF_SG
    use ModelOperator_SP
    use PreConditioner_DC_MF
    use PreConditioner_DC_SP
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
    !> No subroutine briefing
    !
    function Solver_PCG_ctor( model_operator ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        !
        type( Solver_PCG_t ) :: self
        !
        !write( *, * ) "Constructor Solver_PCG_t"
        !
        call self%baseInit
        !
        !> Instantiate the PreConditioner object according to the ModelOperator type
        select type( model_operator )
            !
            class is( ModelOperator_MF_SG_t )
                !
                allocate( self%preconditioner, source = PreConditioner_DC_MF_t( model_operator ) )
            !
            class is( ModelOperator_SP_t )
                !
                allocate( self%preconditioner, source = PreConditioner_DC_SP_t( model_operator ) )
                !
            class default
                call errStop( "Solver_PCG_ctor: Unclassified ModelOperator" )
            !
        end select
        !
        call self%setDefaults
        !
        call self%zeroDiagnostics
        !
    end function Solver_PCG_ctor
    !
    !> No subroutine briefing
    subroutine setDefaults_PCG( self )
        implicit none
        !
        class( Solver_PCG_t ), intent( inout ) :: self
        !
        call self%setParameters( max_divcor_iters, tolerance_divcor )
        !
    end subroutine setDefaults_PCG
    !
    !> No subroutine briefing
    subroutine solvePCG( self, b, x )
        implicit none
        !
        class( Solver_PCG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: b
        class( Scalar_t ), intent( inout ) :: x
        !
        !>    these will have to be created in a way to match
        !>     the specific type of the input Scalar_t ...
        !>    Can we just declare these to be of the abstract type?
        class ( Scalar_t ), allocatable :: r, s, p, q
        complex( kind=prec ) :: beta, alpha, delta, deltaOld
        complex( kind=prec ) :: bnorm, rnorm
        !
        if( .NOT. x%is_allocated ) then
            call errStop( "solvePCG > x not allocated yet" )
        endif
        !
        if( .NOT. b%is_allocated ) then
            call errStop( "solvePCG > b not allocated yet" )
        endif
        !
        !>  create local cScalar objects -- could we also use modOp%createCScalar?
        allocate( r, source = x )    !> cannot zero x, since it is first guess
        call r%zeros
        allocate( s, source = r )
        allocate( p, source = r )
        allocate( q, source = r )
        !
        !> just like
        !call self%preconditioner%model_operator%Amult( x, r )
        !>    if we can make AMult generic, with versions that operate on cScalar/cVector
        !>     this could be more generic ...    also could change the name of this operator
        !>      to make this more obvious
        call self%preconditioner%model_operator%divCGrad( x, r )
        !
        !>     r = b-r
        call r%linComb( b, C_MinusOne, C_ONE )
        !
        bnorm = SQRT( b%dotProd(b) )
        rnorm = SQRT( r%dotProd(r) )
        !
        self%iter = 1
        !
        self%relErr( self%iter ) = rnorm / bnorm
        !
        loop: do while ( ( self%relErr( self%iter ) .GT. self%tolerance ).AND.( self%iter .LT. self%max_iters ) )
            !
            call self%preconditioner%LUsolve( r, s )
            !
            delta = r%dotProd(s)
            if( self%iter .EQ. 1 ) then
                beta = C_ZERO
            else
                beta = delta / deltaOld
            endif
            !
            call p%linComb( s, beta, C_ONE )
            !
            call q%zeros
            call self%preconditioner%model_operator%divCGrad( p, q )
            !
            alpha = delta / p%dotProd(q)
            !
            call x%multAdd( alpha, p )
            !
            call r%multAdd( -alpha, q )
            !
            deltaOld = delta
            !
            rnorm = SQRT( r%dotProd(r) )
            !
            !write( *, "( a46, i6, a3, es12.3 )" ) "PCG self%iter, self%relErr( self%iter ):", self%iter, " : ", self%relErr( self%iter )
            !
            self%iter = self%iter + 1
            !
            self%relErr( self%iter ) = rnorm / bnorm
            !
        enddo loop
        !
        if( self%iter .LT. self%max_iters ) then
            write( *, "( a46, i6, a7, es12.3 )" ) "->divCor PCG converged within ", self%iter, ": err= ", self%relErr( self%iter )
        else
            write( *, "( a46, i6, a7, es12.3 )" ) "->divCor PCG not converged in ", self%max_iters, ": err= ", self%relErr( self%max_iters )
        endif
        !
        deallocate( r )
        deallocate( s )
        deallocate( p )
        deallocate( q )
        !
        self%n_iter = self%iter
        !
    end subroutine solvePCG !> PCG
    !
end module Solver_PCG
