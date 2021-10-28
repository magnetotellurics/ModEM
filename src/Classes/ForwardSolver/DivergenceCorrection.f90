Module DivergenceCorrection
   !
   use ModelOperator
   use ModelOperator_MF
   use Solver_PCG
   use cScalar3D_SG
   use cVector3D_SG
   use Source
   !
   type, public :: DivergenceCorrection_t
      !  need to think about public/private
      class( ModelOperator_t ), pointer :: model_operator !    this should be a pointer
         !  and this is the abstract type ...  when setting up will need to use the
         !   correct model operator sub-type
      type( Solver_PCG_t ) :: PCG
         !   OK -- "Solver" objects will correspond to specific solver methods
         !      solver control parameters, and diagnostics will be properties
         !      of these objects.   Right now we have three extensions of the
         !      base class "Solver" -- QMR, BiCG, PCG; one property of the solver
         !      will be a pointer to the "Preconditioner" 
         !    I am assuming here that we only use PCG for divergence correction,
         !      since this is fast and works well for symmetric real problems
         !      But we could easily change this--I think by just using a different
         !      type here?
         !   ALSO: with the scheme I am imagining now, PCG will always require
         !        RHS to be a TScalar, while QMR and BiCG will require a TVector!
      real( kind=prec ) :: divJ(2) !   divergence of currents computed at most
                                   ! recent call to DivCor -- save these in EMsolver
                                   !  object to keep a record of divergence correction
                                   !  For any other diagnostice to collect, add as
                                   ! object properties, and collect after each call to
                                   ! DivCor
   contains
         !
         final :: DivergenceCorrection_dtor
         !
         !  the usual -- not doing these properly now!
         !procedure, public        :: create     => createDivergenceCorrection
         !procedure, public        :: allocate   => allocateDivergenceCorrection
         !procedure, public        :: deallocate => deallocateDivergenceCorrection
         !    these are the routines with substantive algorithms--note that
         !    all operators are stored and managed by model_operator
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
        !
        type( DivergenceCorrection_t )                 :: self
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        !
        write(*,*) "Constructor DivergenceCorrection_t"
        !
        self%model_operator => model_operator
        !
      end function DivergenceCorrection_ctor
      !
      ! Destructor
      subroutine DivergenceCorrection_dtor( self )
        implicit none
        !
        type( DivergenceCorrection_t ), intent( in out ) :: self
        !
        !write(*,*) "Destructor DivergenceCorrection_t"
        !
      end subroutine DivergenceCorrection_dtor
      !
      function rhsDivCor( self, omega, source ) result( phi0 )
      !   this routine will calculate the RHS (phi0) for divergence
      !    correction equations when the source term for the curl-curl
      !    equations has non-zero source
     !
         implicit none
       !
         class( DivergenceCorrection_t ) :: self
         real ( kind=prec ), intent( in )             :: omega
         class( Source_t ), allocatable, intent( in ) :: source
         class( cScalar_t ), allocatable              :: phi0
         !
         ! local variables 
         complex ( kind=prec )  :: cFactor
         type( cVector3D_SG_t ) :: sourceInterior
     
         cFactor = -C_ONE/(mu_0*ISIGN*omega)   ! 1/(isign*1i*w*mu)
         !   my idea here is that:
         !    a) TVectors (and scalars) will alternatively store
         !         field values in the usual tensor product grids,
         !         or as columnn vectors -- the storage state can
         !         be found, and also switched.   B will be a TVector,
         !         but the state is unknown at compile time
         !    b)  interior is a function (there will also be a "boundary"
         !        function) that returns a TVector (same state as input, B)
         !        
         sourceInterior = source%bdry%interior()
         ! again assumning we have "mult" defined for TVectors, and that for
         !   this function output can overwrite input
         sourceInterior = cFactor * sourceInterior
         !
         select type( modOp => self%model_operator )
         class is (ModelOperator_MF_t)
          !
         ! MODEL OPERATOR DIVIDED????
         !
            !phi0 = modOp / sourceInterior
            !
            !   same assumption for multiplication of TScalars
            phi0 = modOp%metric%VNode * phi0
         class default
            write(*, *) "ERROR:DivergenceCorrection_t::rhsDivCor:"
            STOP        "model_operator type unknow"
         end select
     
      end function rhsDivCor
      !****************************************************************
      function DivCorr( self, inE, phi0 ) result(outE)
        !
        implicit none
        ! function to compute divergence correction for input electric
        ! field vector inE, returning result in outE 
        !  Optional argument phi0 is scaled divergence of source term
        !    computed b rhsDivCor

        !  use sg_scalar ! mult => diagMult_scalar, linComb => linComb_cscalar
        ! rename routines for linear algebra operations; this from old code;
        !    will be different now.  Note: comments like this SHOULD BE DELETED
        !    as soon as they are irrelevant!

        class( DivergenceCorrection_t ) :: self
        type( cVector3D_SG_t ), intent( in )           :: inE
        type( cScalar3D_SG_t ), intent( in ), optional :: phi0
        type( cVector3D_SG_t )                         :: outE
        !
        !  local variables
        type( cScalar3D_SG_t) :: phiSol, phiRHS
        complex( kind=prec)   :: c2
        integer               :: status
        logical               :: SourceTerm
        !
        SourceTerm = present( phi0 )

        ! alocating phiSol, phiRHS
        !   OK -- SOMETHING TO DISCUSS HERE -- want to create TScalar objects
        !    in proper "state" -- MR or CSG ... how to do this automatically?`
        !    Can we have the creation routine just take inE as input (get grid
        !     from this, also get type ... or as I do here use the "copy"
        !     version of creator and model_operator%p -- then do we need to zero?
        !     I guess divergence correction needs to have SP version?
        !       Or just do conversion to TVector in solver???
        !
         select type( modOp => self%model_operator )
         class is (ModelOperator_MF_t)
            !
         ! MODEL_OPERATOR DOESNT HAVE P or PCOPY
         !
         !phiSol =  modOp%pcopy
            !phiRHS =  modOp%p%copy

            ! compute divergence of currents for input electric field
            phiRHS = modOp%DivC( inE )

            !  If source term is present, subtract from divergence of currents
            if(SourceTerm) then
              phiRHS = phiRHS - phi0
            endif

           ! compute the size of current Divergence before (using dot product)
           !   this will be part of diagnostics -- need to develop before including
           !    anything specific here
           self%divJ(1) = sqrt( phiRHS .dot. phiRHS )

           ! point-wise multiplication with volume weights centered on corner nodes
           phiRHS = modOp%metric%Vnode * phiRHS

           ! PCG is a generic pre-conditioned CG algorithm
           !   parameters are set before call to this
           !   We could pass a pointer to the coefficient matrix multiply
           !   but maybe better to set this as a property in the solver
           !    object -- along with pointers to the preconditioners
           !     (also matrix multiply) 
           !Call self%PCG(self%model_operator%divCgrad,phiRHS,phiSol)
           !Call self%PCG%solve( phiRHS, phiSol )

           !   have to decide how to manage output
           !if (output_level > 2) then
           !   write (*,"(a12,a32,i5,g15.7)") node_info, &
           !      "finished divergence correction:", PCGiter%niter, PCGiter%rerr(PCGiter%niter)
           !end if

           ! compute gradient of phiSol (Divergence correction for inE)
           outE = modOp%grad( phiSol )

           ! subtract Divergence correction from inE
           !   note: order here is different than previous use of minus --
           !    still OK to overwrite?
           outE = inE - outE

           ! divergence of the corrected output electrical field
           phiRHS = modOp%DivC(outE)

           !  If source term is present, subtract from divergence of currents
           if(SourceTerm) then
             phiRHS = phiRHS - phi0
           endif

           ! compute the size of current Divergence after
           self%divJ(2) = sqrt( phiRHS .dot. phiRHS )

           ! output level defined in basic file_units module
           !if (output_level > 2) then
           !write(*,"(a12,a47,g15.7)") node_info, "divergence of currents before correction: ", divJ(1, nDivCor)
           !write(*,"(a12,a47,g15.7)") node_info, "divergence of currents  after correction: ", divJ(2, nDivCor)
           !end if

           ! deallocate the temporary work arrays  ... how do we handle deallocation?
           !deallocate( phiSol )
           !deallocate( phiRHS )
         !
         class default
            write(*, *) "ERROR:DivergenceCorrection_t::DivCorr:"
            STOP        "model_operator type unknow"
         end select

   end function DivCorr

end module DivergenceCorrection
