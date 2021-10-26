module DivergenceCorrection
  use Constants
  use ModelOperator
  use Solver_PCG

  type :: DivergenceCorrection_t
     
     type(ModelOperator_t), pointer :: ModOp => null()
     type(Solver_PCG_t) :: PCG
     real(kind = prec)  :: divJ(2) ! Divergence of currents computed at most
                                   ! recent call to DivCor -- save these in EMsolver
                                   ! object to keep a record of divergence correction
                                   ! For any other diagnostice to collect, add as
                                   ! object properties, and collect after each call to
                                   ! DivCor.
   contains
     ! procedure, public        :: Create
     ! procedure, public        :: Allocate
     ! procedure, public        :: DeAllocate

     !         procedure, public	:: rhsDivCor
     !         procedure, public	:: DivCor
  end type DivergenceCorrection_t
  
contains
  
  !  function rhsDivCor(self,omega,src) result(phi0)
  !  !   this routine will calculate the RHS (phi0) for divergence
  !  !    correction equations when the source term for the curl-curl
  !  !    equations has non-zero source
  !     type(DivergeneCorrection_t),intent(in)	:: self
  !     class(Source_t),intent(in)	:: src
  !     class(TScalar_complex),intent(inout)	:: phi0
  !     real (kind=prec)	:: omega
  !     ! local variables 
  !     complex (kind=prec)	:: cFactor
  !     type(TVector3D_SG_complex)	:: srcInterior

  !     cFactor = -C_ONE/(mu_0*ISIGN*omega)   ! 1/(isign*1i*w*mu)
  !     !   my idea here is that:
  !     !    a) TVectors (and scalars) will alternatively store
  !     !         field values in the usual tensor product grids,
  !     !         or as columnn vectors -- the storage state can
  !     !         be found, and also switched.   B will be a TVector,
  !     !         but the state is unknown at compile time
  !     !    b)  interior is a function (there will also be a "boundary"
  !     !        function) that returns a TVector (same state as input, B)
  !     !        
  !     srcInterior = src%B%interior 
  !     ! again assumning we have "mult" defined for TVectors, and that for
  !     !   this function output can overwrite input
  !     srcInterior = srcInterior%mult(cFactor)
  !     phi0 = self%ModOp%Div(srcInterior)
  !     !   same assumption for multiplication of TScalars
  !     phi0 = phi0%mult(self%ModOp%Metric%VNode)

  !  end function rhsDivCor
  !  !****************************************************************
  !  function DivCorr(self,inE,phi0) result(outE)

  !  ! function to compute divergence correction for input electric
  !  ! field vector inE, returning result in outE 
  !  !  Optional argument phi0 is scaled divergence of source term
  !  !    computed b rhsDivCor

  !  !  use sg_scalar ! mult => diagMult_scalar, linComb => linComb_cscalar
  !  ! rename routines for linear algebra operations; this from old code;
  !  !    will be different now.  Note: comments like this SHOULD BE DELETED
  !  !    as soon as they are irrelevant!

  !  type(DivergenceCorrection_t)	:: self
  !  type(TVector3D_SG_complex), intent(in)	:: inE
  !  type(TScalar3D_SG_complex), intent(in), optional	:: phi0
  !  type(TVector3D_SG_complex), intent(inout)	:: outE
  !  !  local variables
  !  implicit none
  !  type (TScalar3D_SG_complex)	        :: phiSol, phiRHS
  !  complex (kind=prec)        	:: c2
  !  integer				:: status
  !  logical				:: SourceTerm

  !  SourceTerm = present(phi0)

  !  ! alocating phiSol, phiRHS
  !  !   OK -- SOMETHING TO DISCUSS HERE -- want to create TScalar objects
  !  !    in proper "state" -- MR or CSG ... how to do this automatically?`
  !  !    Can we have the creation routine just take inE as input (get grid
  !  !     from this, also get type ... or as I do here use the "copy"
  !  !     version of creator and ModOp%p -- then do we need to zero?
  !  !     I guess divergence correction needs to have SP version?
  !  !       Or just do conversion to TVector in solver???
  !  phiSol =  self%ModOp%pcopy
  !  phiRHS =  self%ModOp%p%copy

  !  ! compute divergence of currents for input electric field
  !  phiRHS = self%ModOp%DivC(inE)

  !  !  If source term is present, subtract from divergence of currents
  !  if(SourceTerm) then
  !     phiRHS = phiRHS%minus(phi0)
  !  endif

  ! ! compute the size of current Divergence before (using dot product)
  ! !   this will be part of diagnostics -- need to develop before including
  ! !    anything specific here
  ! self%divJ(1) = sqrt(dotProd(phiRHS,phiRHS))

  ! ! point-wise multiplication with volume weights centered on corner nodes
  ! phiRHS = phiRHS%mult(self%ModOp%Metric%Vnode)

  ! ! PCG is a generic pre-conditioned CG algorithm
  ! !   parameters are set before call to this
  ! !   We could pass a pointer to the coefficient matrix multiply
  ! !   but maybe better to set this as a property in the solver
  ! !    object -- along with pointers to the preconditioners
  ! !     (also matrix multiply) 
  ! !Call self%PCG(self%ModOp%divCgrad,phiRHS,phiSol)
  ! Call self%PCG%solve(phiRHS,phiSol)

  ! !   have to decide how to manage output
  ! !if (output_level > 2) then
  ! !   write (*,'(a12,a32,i5,g15.7)') node_info, &
  ! !      'finished divergence correction:', PCGiter%niter, PCGiter%rerr(PCGiter%niter)
  ! !end if

  ! ! compute gradient of phiSol (Divergence correction for inE)
  ! outE = self%modOP%grad(phiSiol)

  ! ! subtract Divergence correction from inE
  ! !   note: order here is different than previous use of minus --
  ! !    still OK to overwrite?
  ! outE = inE%minus(outE)

  ! ! divergence of the corrected output electrical field
  ! phiRHS = self%ModOp%DivC(outE)

  ! !  If source term is present, subtract from divergence of currents
  ! if(SourceTerm) then
  !    phiRHS = phiRHS%minus(phi0)
  ! endif

  ! ! compute the size of current Divergence after
  ! self%divJ(2) = sqrt(dotProd(phiRHS,phiRHS))

  ! ! output level defined in basic file_units module
  ! !if (output_level > 2) then
  ! !write(*,'(a12,a47,g15.7)') node_info, 'divergence of currents before correction: ', divJ(1, nDivCor)
  ! !write(*,'(a12,a47,g15.7)') node_info, 'divergence of currents  after correction: ', divJ(2, nDivCor)
  ! !end if

  ! ! deallocate the temporary work arrays  ... how do we handle deallocation?
  ! Call deall(phiSol)
  ! Call deall(phiRHS)

end subroutine SdivCorr ! SdivCorr
