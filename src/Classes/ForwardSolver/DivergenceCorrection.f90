Module DivergenceCorrection
   !
   use Solver_PCG
   use ModelOperator
   use cScalar3D_SG
   use cVector3D_SG
   use Source
   !
   type, public :: DivergenceCorrection_t
      !
      class( Solver_t ), allocatable :: solver ! think this should be abstract class!
         !     want this to work for MR soon enough!
         !   OK -- "Solver" objects will correspond to specific solver methods
         !      solver control parameters, and diagnostics will be propertiess
         !      of these objects.   Right now we have three extensions of the
         !      base class "Solver" -- QMR, BiCG, solver_pcg; one property of the solver
         !      will be a pointer to the "Preconditioner" 
         !    I am assuming here that we only use solver_pcg for divergence correction,
         !      since this is fast and works well for symmetric real problems
         !      But we could easily change this--I think by just using a different
         !      type here?
         !   ALSO: with the scheme I am imagining now, solver_pcg will always require
         !        RHS to be a TScalar, while QMR and BiCG will require a TVector!
      real( kind=prec ) :: divJ(2) = 0.0 !   divergence of currents computed at most
                                   ! recent call to DivCor -- save these in EMsolver
                                   !  object to keep a record of divergence correction
                                   !  For any other diagnostice to collect, add as
                                   ! object properties, and collect after each call to
                                   ! DivCor
   contains
         !
         final :: DivergenceCorrection_dtor
         !
         procedure, public :: rhsDivCor
         procedure, public :: DivCorr
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
         !write(*,*) "Constructor DivergenceCorrection_t"
         !
         ! Specific Solver PCG
         self%solver = Solver_PCG_t( model_operator )
         !
      end function DivergenceCorrection_ctor
      !
      ! Destructor
      subroutine DivergenceCorrection_dtor( self )
        implicit none
        !
        type( DivergenceCorrection_t ), intent( inout ) :: self
        !
        !write(*,*) "Destructor DivergenceCorrection_t"
        !
      end subroutine DivergenceCorrection_dtor
      !
      subroutine rhsDivCor( self, omega, source, phi0 )
      !   this routine will calculate the RHS (phi0) for divergence
      !    correction equations when the source term for the curl-curl
      !    equations has non-zero source
      !
         implicit none
      !
         class( DivergenceCorrection_t ) :: self
         real ( kind=prec ), intent( in ) :: omega
         class( Source_t ), intent( in ) :: source
         class( cScalar_t ), intent(inout)         :: phi0
         !
         ! local variables 
         complex ( kind=prec )  :: cFactor
         !   want this to be abstract also!
         class( cVector_t ), allocatable :: sourceInterior
     
         cFactor = -C_ONE/(mu_0*ISIGN*omega)   ! 1/(isign*1i*w*mu)
         !   I am assuming that in the case of an interior source, this is
         !    stored in the cVector source%E -- i
         !    proedure "interior" zeros the boundary edges
         sourceInterior = source%E%interior()
          
         !   take divergence of sourceInterior, and return as cScalar of
         !    appropriate explicit type
         call self%solver%model_operator%Div( sourceInterior, phi0 ) 
         !  multiply result by cFactor (in place)
         call phi0%mults( cFactor )
         !  multiply result by VNode -- add to rhs of symetrized
         !   current conservation equation
         !
       ! CREATE MULTS TO OVERWRITE RHS -- phi0 ????
       call phi0%mults( self%solver%model_operator%metric%Vnode )
     
      end subroutine rhsDivCor
      !****************************************************************
      subroutine DivCorr( self, inE, outE, phi0 )
        implicit none
        ! function to compute divergence correction for input electric
        ! field vector inE, returning result in outE 
        !  Optional argument phi0 is scaled divergence of source term
        !    computed by rhsDivCor
        !
        class( DivergenceCorrection_t )                         :: self
        class( cVector_t ), allocatable, intent( in )           :: inE
        class( cVector_t ), allocatable, intent(inout)          :: outE
        class( cScalar_t ), allocatable, intent( in ), optional :: phi0
        !
        !  local variables
        class( cScalar_t ), allocatable :: phiSol, phiRHS
        complex( kind=prec)   :: c2
        integer               :: status
        logical               :: SourceTerm
        !
        SourceTerm = present( phi0 )

        ! alocating phiSol, phiRHS  -- these need to be cScalars of explicit
        !   type that matches inE, outE -- phi0 may not be an actual input
        !   so this cannot be used as a prototye
        !   I am writing this under the assumption that there will be
        !    createScalar, createVector in ModelOperator class (should be
        !     declared as procedures in abstract class, implemented to return
        !     cVector or cScalar of appropriate type)

        !   I DO NOT WANT select type at this level -- 
        !     make procedures in ForwardModeling generic, with no reference to
        !     specific classes
      !
        allocate( phiSol, source = self%solver%model_operator%createScalar() )
        allocate( phiRHS, source = self%solver%model_operator%createScalar() )
        !
        ! compute divergence of currents for input electric field
        call self%solver%model_operator%DivC(inE, phiRHS )

        !  If source term is present, subtract from divergence of currents
        !  probably OK to use function here -- but could replace with subroutine
        if( SourceTerm ) then
            phiRHS = phiRHS - phi0
        endif

        ! compute the size of current Divergence before (using dot product)
        !   this will be part of diagnostics -- need to develop before including
        !    anything specific here
        self%divJ(1) = sqrt( phiRHS .dot. phiRHS )

        ! point-wise multiplication with volume weights centered on corner nodes
        !
      ! CREATE MULTS TO OVERWRITE RHS -- phiRHS ????
      call phiRHS%mults( self%solver%model_operator%metric%Vnode )

        !   solve system of equations -- solver will have to know about
        !    (a) the equations to solve -- the divergence correction operator
        !       is modOp%divCgrad
        !    (b) preconditioner: object, and preconditioner matrix
      !
      select type( solver => self%solver )
           class is( Solver_PCG_t )
              call solver%solve( phiRHS, phiSol )
         class default
              write(*, *) 'ERROR:DivergenceCorrection::DivCorr:'
              STOP '         Unknow solver type.'
        end select
        !

        !   have to decide how to manage output
        !if (output_level > 2) then
       !   write (*,"(a12,a32,i5,g15.7)") node_info, &
       !      "finished divergence correction:", solver_pcgiter%niter, solver_pcgiter%rerr(solver_pcgiter%niter)
       !end if

       ! compute gradient of phiSol (Divergence correction for inE)
       call self%solver%model_operator%grad(phiSol,outE)

       ! subtract Divergence correction from inE
       call outE%linCombS(inE,C_MinusOne,C_ONE)
       !   this is outE = inE - outE

       ! divergence of the corrected output electrical field
       call self%solver%model_operator%DivC(outE,phiRHS)

       !  If source term is present, subtract from divergence of currents
       if(SourceTerm) then
          phiRHS = phiRHS - phi0
       endif
       ! compute the size of current Divergence after
       self%divJ(2) = sqrt( phiRHS .dot. phiRHS )

       ! deallocate the temporary work arrays  ... how do we handle deallocation?
       !   do we need to do this????
       deallocate( phiSol )
       deallocate( phiRHS )

   end subroutine DivCorr

end module DivergenceCorrection
