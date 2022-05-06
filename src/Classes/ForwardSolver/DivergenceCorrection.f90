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
        real( kind=prec ) :: divJ(2) ! divergence of currents computed at most recent call to DivCor -- before and after
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
        self%solver = Solver_PCG_t( model_operator )
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
    !**********
    !
    subroutine rhsDivCorDivergenceCorrection( self, omega, source, phi0 )
        implicit none
        !
        class( DivergenceCorrection_t ), intent( inout ) :: self
        real ( kind=prec ), intent( in )    :: omega
        class( Source_t ), intent( in )     :: source
        class( cScalar_t ), intent( inout ) :: phi0
        !
        complex( kind=prec ) :: cFactor

        cFactor = -ONE_I/(mu_0*ISIGN*omega)    ! 1/(isign*1i*w*mu)
        !
        !    take divergence of sourceInterior, and return as cScalar of
        !     appropriate explicit type
        call self%solver%preconditioner%model_operator%Div( source%E%interior(), phi0 ) 
        !
        !  multiply result by VNode -- add to rhs of symetrized
        !    current conservation equation
        !
		! GARY'S SUGESTION FOR CSEM: REMOVE NEXT LINE
		!call phi0%mults( self%solver%preconditioner%model_operator%metric%Vnode )
        !
		!  multiply result by cFactor (in place)
        call phi0%mults( cFactor )
        !
    end subroutine rhsDivCorDivergenceCorrection
    !****************************************************************
    subroutine divCorrDivergenceCorrection( self, inE, outE, phi0 )
        implicit none
        ! function to compute divergence correction for input electric
        ! field vector inE, returning result in outE 
        !  Optional argument phi0 is scaled divergence of source term
        !     computed by rhsDivCor
        !
        class( DivergenceCorrection_t ), intent( inout ) :: self
        class( cVector_t ), intent( in )           :: inE
        class( cVector_t ), intent( inout )        :: outE
        class( cScalar_t ), intent( in ), optional :: phi0
        !
        class( cScalar_t ), allocatable :: phiSol, phiRHS
        logical :: SourceTerm
        !
        !
        SourceTerm = present( phi0 )
        !
        ! alocating phiSol, phiRHS  -- these need to be cScalars of explicit
        !    type that matches inE, outE -- phi0 may not be an actual input
        !    so this cannot be used as a prototye
        !    I am writing this under the assumption that there will be
        !     createScalar, createVector in ModelOperator class (should be
        !      declared as procedures in abstract class, implemented to return
        !      cVector or cScalar of appropriate type)

        !    I DO NOT WANT select type at this level -- 
        !      make procedures in ForwardModeling generic, with no reference to
        !      specific classes
        !
        !
        select type( grid => self%solver%preconditioner%model_operator%metric%grid )
            class is( Grid3D_SG_t )
                !
                allocate( phiSol, source = cScalar3D_SG_t( grid, NODE ) )
                allocate( phiRHS, source = cScalar3D_SG_t( grid, NODE ) )
                !
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
		write( *, "(A30, es20.6)" ) "divJ before correction  ", self%divJ(1)
        !
        ! point-wise multiplication with volume weights centered on corner nodes
        !
		call phiRHS%mults( self%solver%preconditioner%model_operator%metric%Vnode )
        !
        !    solve system of equations -- solver will have to know about
        !     (a) the equations to solve -- the divergence correction operator
        !     is modOp%divCgrad
        !     (b) preconditioner: object, and preconditioner matrix
        !
        select type( solver => self%solver )
            class is( Solver_PCG_t )
                call solver%solve( phiRHS, phiSol )
            class default
                write(*, *) "ERROR:DivergenceCorrection::DivCorr:"
                stop          "            Unknow solver type."
        end select
        !
        !    have to decide how to manage output
        !if (output_level > 2) then
        !write (*,*) "finished divergence correction:", size( self%solver%relErr ), self%solver%n_iter
        !write (*,"(i8, es20.6)") self%solver%n_iter, self%solver%relErr( self%solver%n_iter )
        !end if

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
        ! compute the size of current Divergence after
        self%divJ(2) = sqrt( phiRHS .dot. phiRHS )
        !
		write( *, "(A30, es20.6)" ) "divJ after correction  ", self%divJ(2)
        !
        deallocate( phiRHS )

    end subroutine divCorrDivergenceCorrection

end module DivergenceCorrection
