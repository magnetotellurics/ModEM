module Solver_PCG
    !
    use Solver
    use cScalar
    use ModelOperator
    use PreConditioner_MF_DC
    !
    !    solver object for PCG -- used only for Divergence Correction
    type, extends( Solver_t ) :: Solver_PCG_t
        !
        ! PROPERTIES HERE
        !
        contains
            !
            final :: Solver_PCG_dtor
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
    function Solver_PCG_ctor( model_operator ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
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
    ! Solver_PCG destructor
    subroutine Solver_PCG_dtor( self )
        implicit none
        !
        type( Solver_PCG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor Solver_PCG_t"
        !
        call self%dealloc()
        !
    end subroutine Solver_PCG_dtor
    !
    subroutine setDefaults_PCG( self )
        implicit none
        !
        class( Solver_PCG_t ), intent(inout) :: self
        !
        call self%SetParameters( max_iterDivCorDef, tolDivCorDef )
        !
    end subroutine setDefaults_PCG
    !
    !************************************************    
    subroutine solvePCG( self, b, x )
        implicit none
        !
        class( Solver_PCG_t ), intent( inout )           :: self
        class( cScalar_t ), intent( in )                 :: b
        class( cScalar_t ), allocatable, intent( inout ) :: x
        ! local variables
        !    these will have to be created in a way to match
        !     the specific type of the input cScalar_t ...
        !    Can we just declare these to be of the abstract type?
        class ( cScalar_t ), allocatable :: r, s, p, q
        complex( kind=prec ) :: beta, alpha, delta, deltaOld
        complex( kind=prec ) :: bnorm, rnorm
        integer              :: i
        !
        !  create local cScalar objects -- could we also use modOp%createCScalar?
        allocate( r, source = x )    ! cannot zero x, since it is first guess
        call r%zeros()
        allocate( s, source = r )
        allocate( p, source = r )
        allocate( q, source = r )
        !
        ! just like
        !call self%preconditioner%model_operator%Amult( x, r )
        !    if we can make AMult generic, with versions that operate on cScalar/cVector
        !     this could be more generic ...    also could change the name of this operator
        !      to make this more obvious
        call self%preconditioner%model_operator%divCgrad( x, r )
        !
        !     r = b-r
        call r%linCombS( b, C_MinusOne, C_ONE )
        !
        bnorm = sqrt(real( b%dotProd(b)))
        rnorm = sqrt(real( r%dotProd(r)))
        !
        self%relErr(1) = rnorm/bnorm
        !
        !write( *, * ) "PCG iter, self%relErr( 1 )", 1, self%relErr( 1 )
        !
        i = 0
        !
        loop: do while ( ( self%relErr( i + 1 ) .GT. self%tolerance ).and.( i + 1 .LT. self%max_iter ) )
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
            call p%linCombS( s, beta, C_ONE )
            !
            call q%zeros()
            call self%preconditioner%model_operator%divCgrad( p, q )
            !
            alpha = delta/p%dotProd(q)
            !
            call p%scMultAddS( x, alpha )
            !
            call q%scMultAddS( r, -alpha )
            !
            deltaOld = delta
            !
            i = i + 1
            !
            rnorm = sqrt( real( r%dotProd(r) ) )
            !
            self%relErr( i + 1 ) = rnorm/bnorm
            !
            !write( *, * ) "PCG iter, self%relErr( i + 1 )", i + 1, self%relErr( i + 1 )
            !
        enddo loop
        !
        !
        if( i + 1 .LT. self%max_iter ) then
            write( *, * ) "               DivCorr PCG converged within ", i + 1, " : ", self%relErr( i + 1 )
        else
            write( *, * ) "               DivCorr PCG not converged in ", i + 1, " : ", self%relErr( i + 1 )
        endif
        !
        deallocate( r )
        deallocate( s )
        deallocate( p )
        deallocate( q )
        !
        self%n_iter = i
        !
    end subroutine solvePCG ! PCG
    !
end module Solver_PCG
