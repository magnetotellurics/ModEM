!
!> Lone class to define a DivergenceCorrection
!
Module DivergenceCorrection
    !
    use Solver_PCG
    use cScalar3D_SG
    !
    type :: DivergenceCorrection_t
        !
        type( Solver_PCG_t ) :: solver
        !
        real( kind=prec ) :: divJ(2)
        !
        contains
            !
            procedure, public :: setCond => setCondDivergenceCorrection
            !
            procedure, public :: rhsDivCor => rhsDivCorDivergenceCorrection
            !
            procedure, public :: divCorr => divCorrDivergenceCorrection
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
        type( DivergenceCorrection_t ) :: self
        !
        !write( *, * ) "Constructor DivergenceCorrection_t"
        !
        self%divJ = R_ZERO
        !
        !> Specific Solver PCG
        self%solver = Solver_PCG_t( model_operator )
        !
        call self%setCond
        !
    end function DivergenceCorrection_ctor
    !
    !> Procedure setCondDivergenceCorrection
    !> some extra things that need to be done for divergence correction, whenever
    !>      conductivity (model parameter) changes
    subroutine setCondDivergenceCorrection( self )
        implicit none
        !
        class( DivergenceCorrection_t ), intent( inout ) :: self
        !
        !>    set DivCorr arrays in model operator ... 
        call self%solver%preconditioner%model_operator%divCorSetUp
        !>    set preconditioner
        call self%solver%preconditioner%setPreconditioner( self%solver%omega )
        !
    end subroutine setCondDivergenceCorrection
    !
    !> No subroutine briefing
    subroutine rhsDivCorDivergenceCorrection( self, omega, source_e, phi0 )
        implicit none
        !
        class( DivergenceCorrection_t ), intent( in ) :: self
        real( kind=prec ), intent( in ) :: omega
        class( Vector_t ), intent( inout ) :: source_e
        class( Scalar_t ), intent( inout ) :: phi0
        !
        complex( kind=prec ) :: i_omega_mu
        complex( kind=prec ) :: c_factor
        !
        call self%solver%preconditioner%model_operator%Div( source_e, phi0 )
        !
        i_omega_mu = cmplx( 0., real( 1.0d0 * isign * mu_0 * omega, kind=prec ), kind=prec )
        !
        c_factor = C_ONE / i_omega_mu
        !
        call phi0%mult( c_factor )
        !
    end subroutine rhsDivCorDivergenceCorrection
    !
    !> No subroutine briefing
    subroutine divCorrDivergenceCorrection( self, in_e, out_e, phi0 )
        implicit none
        !
        class( DivergenceCorrection_t ), intent( inout ) :: self
        class( Vector_t ), intent( inout ) :: in_e
        class( Vector_t ), intent( inout ) :: out_e
        class( Scalar_t ), intent( in ), optional :: phi0
        !
        class( Scalar_t ), allocatable :: phiSol, phiRHS
        logical :: SourceTerm
        !
        SourceTerm = present( phi0 )
        !
        call self%solver%preconditioner%model_operator%metric%createScalar( complex_t, NODE, phiRHS )
        !
        !> compute divergence of currents for input electric field
        call self%solver%preconditioner%model_operator%DivC( in_e, phiRHS )
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
        !
        !> ???? Interesting point: if changing phiRHS to phiSol, the QMR starts to slowly converge
        call phiRHS%mult( self%solver%preconditioner%model_operator%metric%v_node )
        !
        call self%solver%preconditioner%model_operator%metric%createScalar( complex_t, NODE, phiSol )
        !
        !>    solve system of equations -- solver will have to know about
        !>     (a) the equations to solve -- the divergence correction operator
        !>     is modOp%divCGrad
        !>     (b) preconditioner: object, and preconditioner matrix
        call self%solver%solve( phiRHS, phiSol )
        !
        !>    have to decide how to manage output
        !if(output_level > 2) then
        !write (*,*) "finished divergence correction:", size( self%solver%relErr ), self%solver%n_inv_iter
        !write (*,"(i8, es20.6)") self%solver%n_inv_iter, self%solver%relErr( self%solver%n_inv_iter )
        !endif
        !
        !> compute gradient of phiSol (Divergence correction for in_e)
        call self%solver%preconditioner%model_operator%grad( phiSol, out_e )
        !
        deallocate( phiSol )
        !
        !> subtract Divergence correction from in_e
        !>    out_e = in_e - out_e
        !
        call out_e%linComb( in_e, C_MinusOne, C_ONE )
        !
        !> divergence of the corrected output electrical field
        call self%solver%preconditioner%model_operator%DivC( out_e, phiRHS )
        !
        !>  If source term is present, subtract from divergence of currents
        if( SourceTerm ) then
            call phiRHS%multAdd( C_MinusOne, phi0 )
        endif
        !
        !> compute the size of current Divergence after
        self%divJ(2) = sqrt( phiRHS%dotProd( phiRHS ) )
        !
        !write( *, * ) "                    DivJ: ", self%divJ(1), " => ", self%divJ(2)
        !
        deallocate( phiRHS )
        !
    end subroutine divCorrDivergenceCorrection
    !
end module DivergenceCorrection
