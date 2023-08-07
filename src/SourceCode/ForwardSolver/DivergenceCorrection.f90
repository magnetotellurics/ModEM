!
!> Lone class to define a DivergenceCorrection
!
Module DivergenceCorrection
    !
    use Solver_PCG
    !
    type :: DivergenceCorrection_t
        !
        type( Solver_PCG_t ) :: solver
        !
        real( kind=prec ) :: divJ(2)
        !
        contains
            !
            procedure, public :: setCond => setCond_DivergenceCorrection
            !
            procedure, public :: rhsDivCor => rhsDivCor_DivergenceCorrection
            !
            procedure, public :: divCorr => divCorr_DivergenceCorrection
            !
    end type DivergenceCorrection_t
    !
    interface DivergenceCorrection_t
        module procedure DivergenceCorrection_ctor
    end interface DivergenceCorrection_t
    !
contains
    !
    !> No subroutine briefing
    !
    function DivergenceCorrection_ctor( model_operator ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        !
        type( DivergenceCorrection_t ) :: self
        !
        !write( *, * ) "Constructor DivergenceCorrection_t"
        !
        self%divJ = R_ZERO
        !
        !> Specific Solver PCG
        self%solver = Solver_PCG_t( model_operator )
        !
        !call self%setCond
        !
    end function DivergenceCorrection_ctor
    !
    !> Procedure setCond_DivergenceCorrection
    !> some extra things that need to be done for divergence correction, whenever
    !> conductivity (model parameter) changes
	!
    subroutine setCond_DivergenceCorrection( self, omega )
        implicit none
        !
        class( DivergenceCorrection_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: omega
        !
        self%solver%omega = omega
        !
        call self%solver%preconditioner%setPreconditioner( self%solver%omega )
        !
    end subroutine setCond_DivergenceCorrection
    !
    !> No subroutine briefing
	!
    subroutine rhsDivCor_DivergenceCorrection( self, omega, source_e, phi0 )
        implicit none
        !
        class( DivergenceCorrection_t ), intent( in ) :: self
        real( kind=prec ), intent( in ) :: omega
        class( Vector_t ), intent( in ) :: source_e
        class( Scalar_t ), allocatable, intent( out ) :: phi0
        !
        complex( kind=prec ) :: i_omega_mu, c_factor
        !
        if( .NOT. source_e%is_allocated ) then
            call errStop( "rhsDivCor_DivergenceCorrection > source_e not allocated" )
        endif
        !
        !self%solver%omega = omega
        !
        call self%solver%preconditioner%model_operator%metric%createScalar( complex_t, NODE, phi0 )
        !
        call self%solver%preconditioner%model_operator%div( source_e, phi0 )
        !
        i_omega_mu = cmplx( 0., real( 1.0d0 * isign * mu_0 * omega, kind=prec ), kind=prec )
        !
        c_factor = C_ONE / i_omega_mu
        !
        call phi0%mult( c_factor )
        !
    end subroutine rhsDivCor_DivergenceCorrection
    !
    !> No subroutine briefing
    subroutine divCorr_DivergenceCorrection( self, e_solution, phi0 )
        implicit none
        !
        class( DivergenceCorrection_t ), intent( inout ) :: self
        class( Vector_t ), intent( inout ) :: e_solution
        class( Scalar_t ), intent( in ), optional :: phi0
        !
        class( Vector_t ), allocatable :: temp_e
        class( Scalar_t ), allocatable :: phiRHS, phiSol
        logical :: SourceTerm
        !
        if( .NOT. e_solution%is_allocated ) then
            call errStop( "divCorr_DivergenceCorrection > e_solution not allocated" )
        endif
        !
        SourceTerm = present( phi0 )
        !
        call self%solver%preconditioner%model_operator%metric%createScalar( complex_t, NODE, phiRHS )
        !
        !> compute divergence of currents for input electric field
        call self%solver%preconditioner%model_operator%divC( e_solution, phiRHS )
        !
        !>  If source term is present, subtract from divergence of currents
        !>  probably OK to use function here -- but could replace with subroutine
        if( SourceTerm ) then
            call phiRHS%multAdd( C_MinusOne, phi0 )
        endif
        !
        !> compute the size of current Divergence before (using dot product)
        !>    this will be part of diagnostics
        self%divJ(1) = sqrt( phiRHS%dotProd( phiRHS ) )
        !
        !> point-wise multiplication with volume weights centered on corner nodes
		call phiRHS%mult( self%solver%preconditioner%model_operator%metric%v_node )
		!
        call self%solver%preconditioner%model_operator%metric%createScalar( complex_t, NODE, phiSol )
        !
        !> solve system of equations -- solver will have to know about
        !>     (a) the equations to solve -- the Divergence Correction operator is modOp%divCGrad
        !>     (b) preconditioner: object, and preconditioner matrix
        call self%solver%solve( phiRHS, phiSol )
        !
        !> Temporary Solution - Must be zeroed!????
        allocate( temp_e, source = e_solution )
        call temp_e%zeros
        !
        !> compute gradient of phiSol (Divergence correction for temp_e)
        call self%solver%preconditioner%model_operator%grad( phiSol, temp_e )
        !
        deallocate( phiSol )
        !
        !> subtract Divergence Correction from temp_e
        !>    e_solution = e_solution - temp_e
        !
        call e_solution%sub( temp_e )
        !
        !> divergence of the corrected output electrical field
        call self%solver%preconditioner%model_operator%divC( e_solution, phiRHS )
        !
        !>  If source term is present, subtract from divergence of currents
        if( SourceTerm ) then
            call phiRHS%multAdd( C_MinusOne, phi0 )
        endif
        !
        !> compute the size of current Divergence after
        self%divJ(2) = sqrt( phiRHS%dotProd( phiRHS ) )
        !
        write( *, "( a37, es12.3, a4, es12.3 )" ) "DivJ: ", self%divJ(1), " => ", self%divJ(2)
        !
        deallocate( phiRHS )
        !
    end subroutine divCorr_DivergenceCorrection
    !
end module DivergenceCorrection
!