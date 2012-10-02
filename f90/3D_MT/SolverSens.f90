! *****************************************************************************
!  Module that computes "forcings" for sensitivity calculation
!       (and adjoints).  This module is specific to the numerical
!       implementation of the solver, in this case for 3D
!       MT finite difference modeling.   This module works with the
!       "natural" representations of conductivity: defined on edges
!        of the staggered grid.  Mappings from the potentially
!        more flexible earth conductivity parameter to these fixed,
!        grid-specific representations are to be implemented in module
!	 ModelParam.  This module has no dependence on the specific
!        conductivity parameterization
!
module SolverSens
   use math_constants
   use utilities
   use SolnSpace
   use ModelSpace
   use transmitters
   use datatypes

   implicit none

   !  public routines
   public	::  Pmult, PmultT

   Contains

   !**********************************************************************
   subroutine Pmult(e0,sigma0,dsigma,e)
   !   mapping from modelParam dsigma to source for forward problem
   !    (needed to calculate J*dsigma, where J is sensitivity)
   !   e0 is input background field solution;
   !    e is output ... used for forcing, created before calling
   !    this routine

   type(solnVector_t), intent(in)		    :: e0
   type(modelParam_t), intent(in)	:: sigma0 ! used to compute e0
   type(modelParam_t), intent(in)	:: dsigma
   type(rhsVector_t), intent(inout)		:: e

   !  local variables
   complex(kind=prec)  :: minus_i_omega_mu
   type(rvector_mg)  :: temp
   integer  :: k

   minus_i_omega_mu = cmplx(0.,-ISIGN*MU_0*txDict(e0%tx)%omega,kind=prec)
   call create(e0%grid,temp,EDGE)

   ! map dsigma to edges, storing in array temp
   call dModelParamToEdge(dsigma,temp,sigma0) !multi-grid

   !  multiply temp by i_omeag_mu*e0, put result in e
   do k = 1,e0%nPol
      call diagMult(e0%pol(k),temp,e%b(k)%s)
      call scMult(minus_i_omega_mu,e%b(k)%s,e%b(k)%s)
   enddo

   call deall(temp)
   end subroutine Pmult

   !**********************************************************************
   subroutine PmultT(e0,sigma0,e,dsigmaReal,dsigmaImag)
   !   transpose of Pmult, mapping from adjoint soln e to dsigma
   !        e -> dsigma
   !   e0 is input background field solution
   !   NOTE: because the model parameter is real, while e is complex
   !       the adjoint mapping returns separate data structures
   !        for real and imaginary parts; imaginary output is optional ...

   type(solnVector_t), intent(in)			:: e0
   type(modelParam_t), intent(in)	:: sigma0 ! used to compute e0
   type(solnVector_t), intent(in)			:: e
   type(modelParam_t), intent(inout)		:: dsigmaReal
   type(modelParam_t), intent(inout),optional	:: dsigmaImag

   !  local variables
   complex(kind=prec)			:: minus_i_omega_mu
   type(cvector_mg), pointer    :: Ctemp(:)
   type(rvector_mg)				:: temp

 !  type(rvector_mg)  ::tempMG

   integer					:: k,istat

   minus_i_omega_mu = cmplx(0.,-ISIGN*MU_0*txDict(e0%tx)%omega,kind=prec)
   call create(e0%grid,temp,EDGE) !create_rvector_mg
   allocate(Ctemp(e0%nPol), STAT=istat)
   do k = 1,e0%nPol
      call create(e0%grid,Ctemp(k),EDGE) !create_cvector_mg
   enddo

   ! multiply backward solutions (e) by minus_i_omega_mu * e0
   !   and sum over modes ...
   do k = 1,e0%nPol
      call diagMult(e0%pol(k),e%pol(k),Ctemp(k)) !diagMult_cvector_mg
   enddo
   do k = 2,e0%nPol
      call add(Ctemp(1), Ctemp(k), Ctemp(1))!add_cvector_mg
   enddo
   call scMult(minus_i_omega_mu,Ctemp(1),Ctemp(1)) !scMult_cvector_mg

   ! map real/imag parts onto parameter space
   temp = real(Ctemp(1))
   call dEdgeToModelParam(temp,dsigmaReal,sigma0)

   if(present(dsigmaImag)) then
      ! also compute imaginary part
      temp = imag(Ctemp(1))
      call dEdgeToModelParam(temp,dsigmaImag,sigma0)
   endif

   call deall(temp) !deall_rvector_mg
   do k = 1,e0%nPol
      call deall(Ctemp(k)) !deall_cvector_mg
   enddo
   deallocate(Ctemp, STAT=istat)

   end subroutine PmultT

   !**********************************************************************
end module SolverSens
