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
   complex(kind=prec)		:: minus_i_omega_mu
   type(rvector)			        :: temp
   integer				            :: k

   minus_i_omega_mu = cmplx(0.,-ISIGN*MU_0*txDict(e0%tx)%omega,kind=prec)
   call create_rvector(e0%grid,temp,EDGE)

   ! map dsigma to edges, storing in array temp
   call dModelParamToEdge(dsigma,temp,sigma0)

   !  multiply temp by i_omeag_mu*e0, put result in e
   do k = 1,e0%nPol
      call diagMult_crvector(e0%pol(k),temp,e%b(k)%s)
      call scMult_cvector(minus_i_omega_mu,e%b(k)%s,e%b(k)%s)
   enddo

   call deall_rvector(temp)

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
   type(cvector), pointer		:: Ctemp(:)
   type(rvector)				:: temp
   integer					:: k,istat

   minus_i_omega_mu = cmplx(0.,-ISIGN*MU_0*txDict(e0%tx)%omega,kind=prec)
   call create_rvector(e0%grid,temp,EDGE)
   allocate(Ctemp(e0%nPol), STAT=istat)
   do k = 1,e0%nPol
      call create_cvector(e0%grid,Ctemp(k),EDGE)
   enddo

   ! multiply backward solutions (e) by minus_i_omega_mu * e0
   !   and sum over modes ...
   do k = 1,e0%nPol
      call diagMult_cvector(e0%pol(k),e%pol(k),Ctemp(k))
   enddo
   do k = 2,e0%nPol
      call add_cvector(Ctemp(1), Ctemp(k), Ctemp(1))
   enddo
   call scMult_cvector(minus_i_omega_mu,Ctemp(1),Ctemp(1))

   ! map real/imag parts onto parameter space
   temp = real(Ctemp(1))
   call dEdgeToModelParam(temp,dsigmaReal,sigma0)

   if(present(dsigmaImag)) then
      ! also compute imaginary part
      temp = imag(Ctemp(1))
      call dEdgeToModelParam(temp,dsigmaImag,sigma0)
   endif

   call deall_rvector(temp)
   do k = 1,e0%nPol
      call deall_cvector(Ctemp(k))
   enddo
   deallocate(Ctemp, STAT=istat)

   end subroutine PmultT

   !**********************************************************************
end module SolverSens
