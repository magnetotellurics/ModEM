! *****************************************************************************
module SolverSens
!  Module that computes "forcing" for sensitivity calculation
!       (and adjoints).  This module is specific to the numerical
!       implementation of the solver, in this case for global
!       EARTH finite difference modeling.

   use math_constants
   use utilities
   use sg_scalar
   use model_operators
   use modelmap
   use jacobian
   use boundaries
   use output
   use SolnSpace
   use transmitters
   use dataTypes

   implicit none

   public	::  Pmult, PmultT

   Contains

   !**********************************************************************
   subroutine Pmult(h0,m0,m,h)
   !   mapping from model parameter m to source;
   !       P m = h
   !   h0 is input background field solution;
   !    (h is output used for forcing;
   !    The output is allocated before calling this routine)
   !   rho0 is the FIXED background resistivity on the grid
   !     which is used to compute rho from m0;
   !     rho0 NEVER CHANGED during an inversion

   type(solnVector_t), intent(in)       :: h0
   type(modelParam_t), intent(in)	    :: m0 ! used to compute h0
   type(modelParam_t), intent(in)		:: m
   type(rhsVector_t), intent(inout)     :: h
   !  local variables
   type(cvector)                        :: dH,dE,Hj,Ej
   type(rvector)                        :: drhoF
   type(sparsevecc)                     :: Hb
   type(rscalar)                        :: drho,rho
   type(grid_t), pointer                :: grid
   type(rscalar), pointer               :: rho0
   character(80)                        :: cfile

   ! allocate temporary data structures
   Hj = h0%vec
   grid => h0%grid

   ! map from model parameter to faces ... all this is somewhat confusing since
   ! L \rho = l^F \rho^F (S^F)^{-1}. So, L already does multiplication by length elements
   ! division by area elements (parts of the curl on primary and dual grids). Will clean
   ! this up later; for now, keep as is.
   call operatorP(m,drho,grid,m0)
   call operatorL(drho,drhoF,grid)

   ! Insert boundary conditions in Hj
   call createBC(Hb,grid)
   call insertBC(Hb,Hj)

   ! compute Ej = $\bar{C}$ Hj ... we are taking transposes, so don't conjugate Hj
   call operatorD_l_mult(Hj,grid)
   call operatorC(Hj,Ej,grid)
   !call operatorD_Si_divide(Ej,grid)

   ! compute dE = diag($\bar{C}$ Hj) $\delta\rho_F$
   dE = Ej
   call diagMult(Ej,drhoF,dE)

   ! compute dH = $C^\dag$ dE
   !call operatorD_l_mult(dE,grid)
   call operatorCt(dE,dH,grid)
   call operatorD_Si_divide(dH,grid)

   ! Output
   call scMult(MinusONE,dH,h%source)

   ! Output dH/d\rho
   if (output_level > 3) then
       cfile = trim(cUserDef%modelname)//'_'//trim(freqList%info(h%tx)%code)//'.drho'
       call outputModel(cfile,grid,drho%v)
   endif

   ! Clean up
   call deall_cvector(dH)
   call deall_cvector(dE)
   call deall_cvector(Hj)
   call deall_cvector(Ej)
   call deall_rvector(drhoF)
   call deall_rscalar(drho)
   call deall_sparsevecc(Hb)

   end subroutine Pmult

!**********************************************************************
   subroutine PmultT(h0,m0,h,mReal,mImag)
   !   transpose of Pmult, mapping from adjoint soln h to m
   !   mapping from primary field to model parameter;
   !       P^T h = m
   !   h0 is input background field solution;
   !   rho0 is the FIXED background resistivity on the grid
   !     which is used to compute rho from m0;
   !     rho0 NEVER CHANGED during an inversion
   !   NOTE: because the model parameter is real, while h is complex
   !       the adjoint mapping returns separate data structures
   !        for real and imaginary parts; imaginary output is optional ...

   type(solnVector_t), intent(in)                  :: h0
   type(modelParam_t), intent(in)                  :: m0 ! used to compute h0
   type(solnVector_t), intent(in)                  :: h
   type(modelParam_t), intent(inout)               :: mReal
   type(modelParam_t), intent(inout),optional      :: mImag
   !  local variables
   type(cvector)					    :: dH,dE,Hj,Ej
   type(rvector)					    :: dE_real,dE_imag
   type(sparsevecc)                     :: Hb
   type(rscalar)                        :: drho,rho
   type(grid_t), pointer                :: grid
   type(rscalar), pointer               :: rho0
   character(80)                        :: cfile


   !  allocate temporary data structures
   Hj = h0%vec
   grid => h0%grid
   dH = h%vec

   ! compute dE = $(C^\dag)^T$ dH
   call operatorD_Si_divide(dH,grid)
   call operatorC(dH,dE,grid)
   !call operatorD_l_mult(dE,grid)

   ! Insert boundary conditions in Hj
   call createBC(Hb,grid)
   call insertBC(Hb,Hj)

   ! compute Ej = $\bar{C}$ Hj ... we are taking transposes, so don't conjugate Hj
   call operatorD_l_mult(Hj,grid)
   call operatorC(Hj,Ej,grid)
   !call operatorD_Si_divide(Ej,grid)

   ! compute dE = diag($\bar{C}$ Hj) $(C^\dag)^T$ dH
   call diagMult(Ej,dE,dE)

   ! Map from faces back to model parameter space: real part
   ! ... all this is somewhat confusing since
   ! L \rho = l^F \rho^F (S^F)^{-1}. So, L already does multiplication by length elements
   ! division by area elements (parts of the curl on primary and dual grids). Will clean
   ! this up later; for now, keep as is.
   dE_real = real(dE)
   call operatorLt(drho,dE_real,grid)
   call operatorPt(drho,mReal,grid,m0)
   call scMult(MinusONE,mReal,mReal)

   if(present(mImag)) then

	   ! Map from faces back to model parameter space: imag part
	   dE_imag = imag(dE)
	   call operatorLt(drho,dE_imag,grid)
	   call operatorPt(drho,mImag,grid,m0)
       call scMult(MinusONE,mImag,mImag)

   endif

   ! Output dH/d\rho
   if (output_level > 3) then
       cfile = trim(cUserDef%modelname)//'_'//trim(freqList%info(h%tx)%code)//'.drho'
       call outputModel(cfile,grid,drho%v)
   endif

   ! Clean up
   call deall_cvector(dH)
   call deall_cvector(dE)
   call deall_cvector(Hj)
   call deall_cvector(Ej)
   call deall_rvector(dE_real)
   call deall_rvector(dE_imag)
   call deall_rscalar(drho)
   call deall_sparsevecc(Hb)

  end subroutine PmultT

   !**********************************************************************
end module SolverSens
