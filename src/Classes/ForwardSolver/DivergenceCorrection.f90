Module DivergenceCorrection
    !
    use Solver_PCG
    use ModelOperator
    use cScalar3D_SG
    use cVector3D_SG
    use Source
    !
    type :: DivergenceCorrection_t
        !
        class( Solver_t ), allocatable :: solver
        !
        ! divergence of currents computed at most recent call to DivCor -- before and after
        real( kind=prec ) :: divJ(2)
        !
        contains
            !
            final :: DivergenceCorrection_dtor
            !
            procedure, public :: setCond   => setCondDivergenceCorrection
            procedure, public :: rhsDivCor => rhsDivCorDivergenceCorrection
            procedure, public :: divCorr   => divCorrDivergenceCorrection
            !
    end type DivergenceCorrection_t
    !
    interface DivergenceCorrection_t
        module procedure DivergenceCorrection_ctor
    end interface DivergenceCorrection_t
    !
contains
    !
    !
    function DivergenceCorrection_ctor( model_operator ) result( self )
        implicit none
        !
        class( ModelOperator_t ), intent( in ) :: model_operator
        type( DivergenceCorrection_t )         :: self
        !
        !write( *, * ) "Constructor DivergenceCorrection_t"
        !
        self%divJ = 0.0
        !
        ! Specific Solver PCG
        allocate( self%solver, source = Solver_PCG_t( model_operator ) )
        !
        call self%setCond()
        !
    end function DivergenceCorrection_ctor
    !
    ! Destructor
    subroutine DivergenceCorrection_dtor( self )
        implicit none
        !
        type( DivergenceCorrection_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor DivergenceCorrection_t"
        !
        deallocate( self%solver )
        !
    end subroutine DivergenceCorrection_dtor
    !
    ! some extra things that need to be done for divergence correction, whenever
    !      conductivity (model parameter) changes
    subroutine setCondDivergenceCorrection( self )
        implicit none
        !
        class( DivergenceCorrection_t ), intent( inout ) :: self
        !
        !    set DivCorr arrays in model operator ... 
        call self%solver%preconditioner%model_operator%divCorSetup
        !    set preconditioner
        call self%solver%preconditioner%setPreconditioner( self%solver%omega )
        !
    end subroutine setCondDivergenceCorrection
    !
    ! 
    subroutine rhsDivCorDivergenceCorrection( self, omega, source, phi0 )
        implicit none
        !
        class( DivergenceCorrection_t ), intent( in ) :: self
        real( kind=prec ), intent( in )               :: omega
        class( Source_t ), intent( in )               :: source
        class( Scalar_t ), intent( inout )            :: phi0
        !
        complex( kind=prec ) :: c_factor
        !
        c_factor = -ONE_I / ( mu_0 * ISIGN * omega ) ! 1 / ( isign * 1i * w * mu )
        !
        call self%solver%preconditioner%model_operator%Div( source%E%interior(), phi0 )
        !
        call phi0%mult( c_factor )
        !
    end subroutine rhsDivCorDivergenceCorrection
    !
    !
    subroutine divCorrDivergenceCorrection( self, inE, outE, phi0 )
        implicit none
        !
        class( DivergenceCorrection_t ), intent( inout ) :: self
        class( Vector_t ), intent( in )                  :: inE
        class( Vector_t ), intent( inout )               :: outE
        class( Scalar_t ), intent( in ), optional        :: phi0
        !
        class( Scalar_t ), allocatable :: phiSol, phiRHS
        logical :: SourceTerm
        !
        !
        SourceTerm = present( phi0 )
        !
        select type( grid => self%solver%preconditioner%model_operator%metric%grid )
            class is( Grid3D_SG_t )
                !
                allocate( phiSol, source = cScalar3D_SG_t( grid, NODE ) )
                !
                call phiSol%zeros()
                !
                allocate( phiRHS, source = cScalar3D_SG_t( grid, NODE ) )
                !
                call phiRHS%zeros()
                !
            class default
                stop "Error: divCorrDivergenceCorrection > unknown grid type"
        end select
        !
        ! compute divergence of currents for input electric field
        call self%solver%preconditioner%model_operator%DivC( inE, phiRHS )
        !
        !  If source term is present, subtract from divergence of currents
        !  probably OK to use function here -- but could replace with subroutine
        if( SourceTerm ) then
            call phi0%scMultAddS( phiRHS, C_MinusOne )
        endif

        ! compute the size of current Divergence before (using dot product)
        !    this will be part of diagnostics
        self%divJ(1) = sqrt( phiRHS .dot. phiRHS )
        !
        ! point-wise multiplication with volume weights centered on corner nodes
        !
        ! ???? Interesting point: if changing phiRHS to phiSol, the QMR starts to slowly converge
        call phiRHS%mult( self%solver%preconditioner%model_operator%metric%Vnode )
        !
        !    solve system of equations -- solver will have to know about
        !     (a) the equations to solve -- the divergence correction operator
        !     is modOp%divCgrad
        !     (b) preconditioner: object, and preconditioner matrix
        !
        select type( solver => self%solver )
            !
            class is( Solver_PCG_t )
                call solver%solve( phiRHS, phiSol )
            !
            class default
                stop "Error: divCorrDivergenceCorrection > Unknown solver type."
            !
        end select
        !
        !    have to decide how to manage output
        !if (output_level > 2) then
        !write (*,*) "finished divergence correction:", size( self%solver%relErr ), self%solver%n_iter
        !write (*,"(i8, es20.6)") self%solver%n_iter, self%solver%relErr( self%solver%n_iter )
        !end if
        !
        ! compute gradient of phiSol (Divergence correction for inE)
        call self%solver%preconditioner%model_operator%grad( phiSol, outE )
        !
        deallocate( phiSol )
        !
        ! subtract Divergence correction from inE
        !    outE = inE - outE
        !
        call outE%linCombS( inE, C_MinusOne, C_ONE )
        !
        ! divergence of the corrected output electrical field
        call self%solver%preconditioner%model_operator%DivC( outE, phiRHS )

        !  If source term is present, subtract from divergence of currents
        if( SourceTerm ) then
            call phi0%scMultAddS( phiRHS, C_MinusOne )
        endif
        !
        ! compute the size of current Divergence after
        self%divJ(2) = sqrt( phiRHS .dot. phiRHS )
        !
        write( *, * ) "               DivJ: ", self%divJ(1), " => ", self%divJ(2)
        !
        deallocate( phiRHS )
        !
    end subroutine divCorrDivergenceCorrection

end module DivergenceCorrection
