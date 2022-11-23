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

   type(rscalar),pointer,dimension(:),private :: resist,dxx,dxy,dxz,dyy,dyz,dzz  !!!!!
   
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
!
!   minus_i_omega_mu = cmplx(0.,-ISIGN*MU_0*txDict(e0%tx)%omega,kind=prec)
!   call create_rvector(e0%grid,temp,EDGE)
!
!   ! map dsigma to edges, storing in array temp
!   call dModelParamToEdge(dsigma,temp,sigma0)
!
!   !  multiply temp by i_omeag_mu*e0, put result in e
!   do k = 1,e0%nPol
!      call diagMult_crvector(e0%pol(k),temp,e%b(k)%s)
!      call scMult_cvector(minus_i_omega_mu,e%b(k)%s,e%b(k)%s)
!   enddo
!
!   call deall_rvector(temp)

   end subroutine Pmult

   !**********************************************************************
   subroutine PmultT(e0,rho0,e,dsigmaReal,dsigmaImag)
   !   transpose of Pmult, mapping from adjoint soln e to dsigma
   !        e -> dsigma
   !   e0 is input background field solution
   !   NOTE: because the model parameter is real, while e is complex
   !       the adjoint mapping returns separate data structures
   !        for real and imaginary parts; imaginary output is optional ...

   type(solnVector_t), intent(in)			:: e0
   type(modelParam_t), intent(in)	:: rho0 ! used to compute e0
   type(solnVector_t), intent(in)			:: e
   type(modelParam_t), intent(inout)		:: dsigmaReal
   type(modelParam_t), intent(inout),optional	:: dsigmaImag
   character(80)           :: paramType = 'LOGE'
   
   !  local variables
   type(rscalar),allocatable		:: Ctemp(:)
   integer					:: istat,ia
   real(kind=prec) :: omega
   
   omega = txDict(e0%tx)%omega 

   allocate(Ctemp(6))
   do ia = 1,6
     call create_rscalar(e0%grid,Ctemp(ia),CELL_EARTH)
   enddo

   call dCTKmultCvector(rho0,omega,e0%pol(1),e%pol(1),e0%pol(2),e%pol(2),Ctemp)
   
   ! map real/imag parts onto parameter space
   if(present(dsigmaImag)) then
     !暂时什么都不做，因为用不到
   else
     call create_modelParam(e0%grid,paramtype,dsigmaReal,Ctemp)
   endif

   !release memory
   do ia = 1,6
     call deall_rscalar(Ctemp(ia))
   enddo
   deallocate(Ctemp)
   
  end subroutine PmultT
  
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine dCTKmultCvector(sigma0,omega,cInL1,cInR1,cInL2,cInR2,cOut)
     implicit none
     ! calculate h^T*(dK/dm)*e, sigma0 is the true ln(m)
     type(modelParam_t), intent(in)	:: sigma0 
     type(cvector), intent(in)		        :: cInL1,cInR1,cInL2,cInR2
     type(rscalar), intent(inout) :: cOut(6)
     real(kind=prec) :: omega
     ! local variables
     integer i,j,k,ia,km
     integer Nx,Ny,Nz,Nz0,NzAir
     complex(kind=prec) :: dCdm !temporary variable for dC/dm
     ! > common coefficients
     real(kind=prec) :: dxi,dxi2,dxi4,dyj,dyj2,dyj4,dzk,dzk2,dzk4
     real(kind=prec) :: xiyj,xizk,yjxi,yjzk,zkxi,zkyj
     real(kind=prec) :: xi2yjzk,xi2zkyj,xi4yjzk,xi4zkyj
     real(kind=prec) :: yj2xizk,yj2zkxi,yj4xizk,yj4zkxi
     !real(kind=prec) :: zk2xiyj,zk2yjxi,zk4xiyj,zk4yjxi
     complex(kind=prec) sum1(6),sum2(6),up1(6),up2(6),down1(6),down2(6), &
     sumup1,sumdown1,sumup2,sumdown2 !,sum
     type(grid_t)   :: grid 
     character(80)  :: paramType = 'LINEAR'  !!!!!
     complex(kind=prec) :: ciuw

     ciuw = ISIGN * CMPLX(0.0d0,omega*MU_0)
     
     allocate(resist(3),dxx(6),dxy(6),dxz(6),dyy(6),dyz(6),dzz(6))
     
     call getGrid_modelParam(grid,sigma0)
     call getResValue_modelParam(sigma0,paramType,resist)
     call SensModelParamToTensor(sigma0,dxx,dxy,dxz,dyy,dyz,dzz)
     
     Nx = grid%Nx
     Ny = grid%Ny
     Nz = grid%Nz  ! grid%Nz=grid%NzEarth+grid%NzAir
     NzAir = grid%NzAir
     Nz0 = grid%NzAir + 1

     ! In the Earth 
     do k = Nz0, Nz
     
       ! -----------------------
       km = k - NzAir ! index km for sensitivities in Earth
       dzk = grid%Dz(k) !delete
       dzk2 = 0.5d0*grid%Dz(k)
       dzk4 = 0.25d0*grid%Dz(k) 
       ! -----------------------
             
       do j = 1, Ny
       
         ! -----------------------
         dyj = grid%Dy(j) !delete
         dyj2 = 0.5d0*grid%Dy(j)
         dyj4 = 0.25d0*grid%Dy(j)
         yjzk = grid%Dy(j)/grid%Dz(k)
         zkyj = grid%Dz(k)/grid%Dy(j)
         ! -----------------------
                
         do i = 1, Nx
                    
           ! common coefficients       
           dxi = grid%Dx(i) !  delete
           dxi2 = 0.5d0*grid%Dx(i)
           dxi4 = 0.25d0*grid%Dx(i)
           xizk = grid%Dx(i)/grid%Dz(k)
           zkxi = grid%Dz(k)/grid%Dx(i)  
           ! xiyj = grid%Dx(i)/grid%Dy(j)
           ! yjxi = grid%Dy(j)/grid%Dx(i)      

           ! 
           xi2yjzk = dxi2*yjzk
           xi2zkyj = dxi2*zkyj
           xi4yjzk = dxi4*yjzk
           xi4zkyj = dxi4*zkyj
  
           yj2xizk = dyj2*xizk
           yj2zkxi = dyj2*zkxi
           yj4xizk = dyj4*xizk
           yj4zkxi = dyj4*zkxi 
                                             
           ! zk2xiyj = dzk2*xiyj
           ! zk2yjxi = dzk2*yjxi
           ! zk4xiyj = dzk4*xiyj
           ! zk4yjxi = dzk4*yjxi  
           
           ! zero initialization
           do ia = 1,6
             sum1(ia) = C_ZERO 
             sum2(ia) = C_ZERO
             up1(ia) = C_ZERO
             up2(ia) = C_ZERO 
             down1(ia) = C_ZERO
             down2(ia) = C_ZERO                  
           enddo
        
         ! ------------------------
         !
         ! >  dC/dm for Hx equations
         !
         ! ------------------------
  
           ! ================
           ! row index i,j,k -Hx
           ! ================
           
         if(j.ne.1) then
          
 ! ------------------------- X components ------------------------------- 
          
           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
           
             ! ==> dC5 * Hx(i,j,k)
             dCdm=xi4yjzk*dyy(ia)%v(i,j,km)+xi4zkyj*dzz(ia)%v(i,j,km)  &
                  -dxi4*dyz(ia)%v(i,j,km) !1/2*dC5
             sumup1 = sumup1 + dCdm * cInR1%x(i,j,k)
					   sumdown1 = sumdown1 + dCdm * cInL1%x(i,j,k)
					   sumup2 = sumup2 + dCdm * cInR2%x(i,j,k)
					   sumdown2 = sumdown2 + dCdm * cInL2%x(i,j,k)
                     
             ! ==> dC8 * Hx(i,j,k+1)
             dCdm=-xi2yjzk*dyy(ia)%v(i,j,km)
             sumup1 = sumup1 + dCdm * cInR1%x(i,j,k+1)
					   sumdown1 = sumdown1 + dCdm * cInL1%x(i,j,k+1)
					   sumup2 = sumup2 + dCdm * cInR2%x(i,j,k+1)
					   sumdown2 = sumdown2 + dCdm * cInL2%x(i,j,k+1)        

             if(j.ne.Ny) then
               ! ==> dC6 * Hx(i,j+1,k)
               dCdm=-xi2zkyj*dzz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%x(i,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%x(i,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%x(i,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%x(i,j+1,k)
           
               ! ==> dC9 * Hx(i,j+1,k+1)
               dCdm=dxi2*dyz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%x(i,j+1,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%x(i,j+1,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%x(i,j+1,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%x(i,j+1,k+1)      
             endif 
                   
 ! ------------------------- Y components -------------------------------            
                               
             if(i.ne.1) then
               ! ==> dC18 * Hy(i,j,k)
               dCdm=dxi4*dxz(ia)%v(i,j,km)-dzk2*dzz(ia)%v(i,j,km)  &
                    -xi4yjzk*dxy(ia)%v(i,j,km)+dyj4*dyz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%y(i,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k)
           
               ! ==> dC20 * Hy(i,j,k+1)
               dCdm=xi4yjzk*dxy(ia)%v(i,j,km)+dyj4*dyz(ia)%v(i,j,km)  &
                    -dxi4*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%y(i,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k+1)                            
             endif
                   
             if(i.ne.Nx) then
               ! ==> dC19 * Hy(i+1,j,k)
               dCdm=dxi4*dxz(ia)%v(i,j,km)+dzk2*dzz(ia)%v(i,j,km)  &
                    -xi4yjzk*dxy(ia)%v(i,j,km)-dyj4*dyz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k)

               ! ==> dC21 * Hy(i+1,j,k+1)
               dCdm=xi4yjzk*dxy(ia)%v(i,j,km)-dyj4*dyz(ia)%v(i,j,km)  &
                    -dxi4*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k+1)                                     
             endif         
           
 ! ------------------------- Z components -------------------------------             
 
             if(i.ne.1) then
               ! ==> dC30 * Hz(i,j,k)
               dCdm=dzk4*dyz(ia)%v(i,j,km)-xi4zkyj*dxz(ia)%v(i,j,km)  &
                    +dxi4*dxy(ia)%v(i,j,km)-dyj2*dyy(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i,j,k)

               if(j.ne.Ny) then
                 ! ==> dC32 * Hz(i,j+1,k)
                 dCdm=xi4zkyj*dxz(ia)%v(i,j,km)+dzk4*dyz(ia)%v(i,j,km)  &
                      -dxi4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i,j+1,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i,j+1,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i,j+1,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i,j+1,k)  
							 endif 
                    
						 endif
 
             if(i.ne.Nx) then
               ! ==> dC31 * Hz(i+1,j,k)
               dCdm=-dzk4*dyz(ia)%v(i,j,km)-xi4zkyj*dxz(ia)%v(i,j,km)  &
                    +dxi4*dxy(ia)%v(i,j,km)+dyj2*dyy(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i+1,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i+1,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j,k)
					         
               if(j.ne.Ny) then
                 ! ==> dC33 * Hz(i+1,j+1,k)
                 dCdm=xi4zkyj*dxz(ia)%v(i,j,km)-dzk4*dyz(ia)%v(i,j,km)  &
                      -dxi4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
               endif 
                    
             endif 
 
					   ! sum up h^T*e for one model parameter, row Hx(i,j,k)
					   up1(ia) = up1(ia) + sumup1 * cInL1%x(i,j,k)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%x(i,j,k)
					   up2(ia) = up2(ia) + sumup2 * cInL2%x(i,j,k)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%x(i,j,k)          

           end do !ia
               
! -------------------------- End of the First Row i,j,k ---------------------------------

         end if

           ! ================
           ! row index i,j+1,k -Hx
           ! ================
                    
         if(j.ne.Ny) then
           
 ! ------------------------- X components -------------------------------          

           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
					              
             ! ==> dC5 * Hx(i,j+1,k)
             dCdm=dxi4*dyz(ia)%v(i,j,km)+xi4yjzk*dyy(ia)%v(i,j,km)  &
                  +xi4zkyj*dzz(ia)%v(i,j,km)
             sumup1 = sumup1 + dCdm * cInR1%x(i,j+1,k)
					   sumdown1 = sumdown1 + dCdm * cInL1%x(i,j+1,k)
					   sumup2 = sumup2 + dCdm * cInR2%x(i,j+1,k)
					   sumdown2 = sumdown2 + dCdm * cInL2%x(i,j+1,k)

             if(j.ne.1) then
               ! ==> dC7 * Hx(i,j,k+1)
               dCdm=-dxi2*dyz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%x(i,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%x(i,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%x(i,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%x(i,j,k+1)
             endif           
 
             ! ==> dC8 * Hx(i,j+1,k+1)
             dCdm=-xi2yjzk*dyy(ia)%v(i,j,km)
             sumup1 = sumup1 + dCdm * cInR1%x(i,j+1,k+1)
					   sumdown1 = sumdown1 + dCdm * cInL1%x(i,j+1,k+1)
					   sumup2 = sumup2 + dCdm * cInR2%x(i,j+1,k+1)
					   sumdown2 = sumdown2 + dCdm * cInL2%x(i,j+1,k+1)
                    
 ! ------------------------- Y components -------------------------------            
           
             if(i.ne.1) then
               ! ==> dC12 * Hy(i,j,k)
               dCdm=-dxi4*dxz(ia)%v(i,j,km)+dzk2*dzz(ia)%v(i,j,km)  &
                    -xi4yjzk*dxy(ia)%v(i,j,km)+dyj4*dyz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%y(i,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k)
               

               ! ==> dC14 * Hy(i,j,k+1)
               dCdm=xi4yjzk*dxy(ia)%v(i,j,km)+dyj4*dyz(ia)%v(i,j,km)  &
                    +dxi4*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%y(i,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k+1)
             endif
           
             if(i.ne.Nx) then
               ! ==> dC13 * Hy(i+1,j,k)
               dCdm=-dxi4*dxz(ia)%v(i,j,km)-dzk2*dzz(ia)%v(i,j,km)  &
                    -xi4yjzk*dxy(ia)%v(i,j,km)-dyj4*dyz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k)

               ! ==> dC15 * Hy(i+1,j,k+1)
               dCdm=xi4yjzk*dxy(ia)%v(i,j,km)-dyj4*dyz(ia)%v(i,j,km)  &
                    +dxi4*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k+1)
             endif           

 ! ------------------------- Z components ------------------------------- 
         
             if(i.ne.1) then
         
               if(j.ne.1) then
                 ! ==> dC28 * Hz(i,j,k)
                 dCdm=xi4zkyj*dxz(ia)%v(i,j,km)-dzk4*dyz(ia)%v(i,j,km)  &
                      +dxi4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i,j,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i,j,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i,j,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i,j,k)
               endif
           
               ! ==> dC30 * Hz(i,j+1,k)
               dCdm=-dzk4*dyz(ia)%v(i,j,km)-xi4zkyj*dxz(ia)%v(i,j,km)  &
                    -dxi4*dxy(ia)%v(i,j,km)-dyj2*dyy(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i,j+1,k)
                               
             endif
 
             if(i.ne.Nx) then
         
               if(j.ne.1) then
                 ! ==> dC29 * Hz(i+1,j,k)
                 dCdm=xi4zkyj*dxz(ia)%v(i,j,km)+dzk4*dyz(ia)%v(i,j,km)  &
                      +dxi4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i+1,j,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i+1,j,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j,k)
               endif
           
               ! ==> dC31 * Hz(i+1,j+1,k)
               dCdm=dzk4*dyz(ia)%v(i,j,km)-xi4zkyj*dxz(ia)%v(i,j,km)  &
                    -dxi4*dxy(ia)%v(i,j,km)+dyj2*dyy(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
					     
             endif 
 
             ! sum up h^T*e for one model parameter, row Hx(i,j+1,k)
					   up1(ia) = up1(ia) + sumup1 * cInL1%x(i,j+1,k)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%x(i,j+1,k)
					   up2(ia) = up2(ia) + sumup2 * cInL2%x(i,j+1,k)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%x(i,j+1,k)                  

           end do !ia
                      
! -------------------------- End of the Second Row i,j+1,k ---------------------------------
                      
         end if


           ! ================
           ! row index i,j,k+1 -Hx
           ! ================
           
         if(j.ne.1) then
         
 ! ------------------------- X components -------------------------------          

           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
					   
             if(k.ne.Nz) then
               ! ==> dC5 * Hx(i,j,k+1)
               dCdm=dxi4*dyz(ia)%v(i,j,km)+xi4yjzk*dyy(ia)%v(i,j,km)  &
                    +xi4zkyj*dzz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%x(i,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%x(i,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%x(i,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%x(i,j,k+1)
           
               if(j.ne.Ny) then
                 ! ==> dC6 * Hx(i,j+1,k+1)
                 dCdm=-xi2zkyj*dzz(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%x(i,j+1,k+1)
					       sumdown1 = sumdown1 + dCdm * cInL1%x(i,j+1,k+1)
					       sumup2 = sumup2 + dCdm * cInR2%x(i,j+1,k+1)
					       sumdown2 = sumdown2 + dCdm * cInL2%x(i,j+1,k+1)
               endif
                             
             else
               ! ==> dB4 * Hx(i,j,k+1)
               dCdm=dxi4*dyz(ia)%v(i,j,km)+xi4yjzk*dyy(ia)%v(i,j,km) ! -dZyx*dxi2*grid%delY(j) !debug
               sumup1 = sumup1 + dCdm * cInR1%x(i,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%x(i,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%x(i,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%x(i,j,k+1)
             endif
         
 ! ------------------------- Y components -------------------------------          
           
             if(i.ne.1) then
               ! ==> dC16&B9 * Hy(i,j,k)
               dCdm=xi4yjzk*dxy(ia)%v(i,j,km)-dyj4*dyz(ia)%v(i,j,km)  &
                    +dxi4*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%y(i,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k)
					          
               if(k.ne.Nz) then
                 ! ==> dC18 * Hy(i,j,k+1)
                 dCdm=-dxi4*dxz(ia)%v(i,j,km)-dzk2*dzz(ia)%v(i,j,km)  &
                      -xi4yjzk*dxy(ia)%v(i,j,km)-dyj4*dyz(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%y(i,j,k+1)
					       sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k+1)
					       sumup2 = sumup2 + dCdm * cInR2%y(i,j,k+1)
					       sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k+1)
               else
                 ! ==> dB11 * Hy(i,j,k+1)
                 dCdm=-dxi4*dxz(ia)%v(i,j,km)-xi4yjzk*dxy(ia)%v(i,j,km)  &
                      -dyj4*dyz(ia)%v(i,j,km) !-dZyy*dxi2*dyj2 !debug
                 sumup1 = sumup1 + dCdm * cInR1%y(i,j,k+1)
					       sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k+1)
					       sumup2 = sumup2 + dCdm * cInR2%y(i,j,k+1)
					       sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k+1)
               endif
                               
             endif

             if(i.ne.Nx) then
               ! ==> dC17&B10 * Hy(i+1,j,k)
               dCdm=xi4yjzk*dxy(ia)%v(i,j,km)+dyj4*dyz(ia)%v(i,j,km)  &
                    +dxi4*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k)
					     
               if(k.ne.Nz) then
                 ! ==> dC19 * Hy(i+1,j,k+1)
                 dCdm=-dxi4*dxz(ia)%v(i,j,km)+dzk2*dzz(ia)%v(i,j,km)  &
                      -xi4yjzk*dxy(ia)%v(i,j,km)+dyj4*dyz(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k+1)
					       sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k+1)
					       sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k+1)
					       sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k+1)
               else
                 ! ==> dB12 * Hy(i+1,j,k+1)
                 dCdm=-dxi4*dxz(ia)%v(i,j,km)-xi4yjzk*dxy(ia)%v(i,j,km)  &
                      +dyj4*dyz(ia)%v(i,j,km) !-dZyy*dxi2*dyj2 !debug
                 sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k+1)
					       sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k+1)
					       sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k+1)
					       sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k+1)
               endif
           
             endif

 ! ------------------------- Z components -------------------------------             
 
             if(i.ne.1) then
               ! ==> dC24&B15 * Hz(i,j,k)
               dCdm=dzk4*dyz(ia)%v(i,j,km)-xi4zkyj*dxz(ia)%v(i,j,km)  &
                    -dxi4*dxy(ia)%v(i,j,km)+dyj2*dyy(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i,j,k)
					     
               if(j.ne.Ny) then
                 ! ==> dC26&B17 * Hz(i,j+1,k)
                 dCdm=xi4zkyj*dxz(ia)%v(i,j,km)+dzk4*dyz(ia)%v(i,j,km)  &
                      +dxi4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i,j+1,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i,j+1,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i,j+1,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i,j+1,k)
               endif 
           
             endif

             if(i.ne.Nx) then
               ! ==> dC25&B16 * Hz(i+1,j,k)
               dCdm=-dzk4*dyz(ia)%v(i,j,km)-xi4zkyj*dxz(ia)%v(i,j,km)  &
                    -dxi4*dxy(ia)%v(i,j,km)-dyj2*dyy(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i+1,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i+1,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j,k)
					     
               if(j.ne.Ny) then
                 ! ==> dC27&B18 * Hz(i+1,j+1,k)
                 dCdm=xi4zkyj*dxz(ia)%v(i,j,km)-dzk4*dyz(ia)%v(i,j,km)  &
                      +dxi4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
               endif
                     
             endif

             ! sum up h^T*e for one model parameter, row Hx(i,j,k+1)
					   up1(ia) = up1(ia) + sumup1 * cInL1%x(i,j,k+1)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%x(i,j,k+1)
					   up2(ia) = up2(ia) + sumup2 * cInL2%x(i,j,k+1)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%x(i,j,k+1)

           end do !ia 
           
! -------------------------- End of the Third Row i,j,k+1 ---------------------------------                   
           
         end if
         
                            
           ! ================
           ! row index i,j+1,k+1 -Hx
           ! ================

         if(j.ne.Ny) then
                    
 ! ------------------------- X components -------------------------------          
 
           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
					             
             if(k.ne.Nz) then
               ! ==> dC5 * Hx(i,j+1,k+1)
               dCdm=-dxi4*dyz(ia)%v(i,j,km)+xi4yjzk*dyy(ia)%v(i,j,km)  &
                    +xi4zkyj*dzz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%x(i,j+1,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%x(i,j+1,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%x(i,j+1,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%x(i,j+1,k+1)
             else
               ! ==> dB4 * Hx(i,j+1,k+1)
               dCdm=-dxi4*dyz(ia)%v(i,j,km)+xi4yjzk*dyy(ia)%v(i,j,km) !-dZyx*dxi2*grid%delY(j+1) !debug
               sumup1 = sumup1 + dCdm * cInR1%x(i,j+1,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%x(i,j+1,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%x(i,j+1,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%x(i,j+1,k+1)
             endif         
          
 ! ------------------------- Y components -------------------------------   
 
             if(i.ne.1) then                               
               ! ==> dC10&B5 * Hy(i,j,k)
               dCdm=xi4yjzk*dxy(ia)%v(i,j,km)-dyj4*dyz(ia)%v(i,j,km)  &
                    -dxi4*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%y(i,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k)
					     
               if(k.ne.Nz) then
                 ! ==> dC12 * Hy(i,j,k+1)
                 dCdm=dxi4*dxz(ia)%v(i,j,km)+dzk2*dzz(ia)%v(i,j,km)  &
                      -xi4yjzk*dxy(ia)%v(i,j,km)-dyj4*dyz(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%y(i,j,k+1)
					       sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k+1)
					       sumup2 = sumup2 + dCdm * cInR2%y(i,j,k+1)
					       sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k+1)
               else
                 ! ==> dB7 * Hy(i,j,k+1)
                 dCdm=dxi4*dxz(ia)%v(i,j,km)  & ! -dZyy*dxi2*dyj2 !debug
                      -xi4yjzk*dxy(ia)%v(i,j,km)-dyj4*dyz(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%y(i,j,k+1)
					       sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k+1)
					       sumup2 = sumup2 + dCdm * cInR2%y(i,j,k+1)
					       sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k+1)
               endif
         
             endif
         
             if(i.ne.Nx) then                               
               ! ==> dC11&B6 * Hy(i+1,j,k)
               dCdm=xi4yjzk*dxy(ia)%v(i,j,km)+dyj4*dyz(ia)%v(i,j,km)  &
                    -dxi4*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k)
					     
               if(k.ne.Nz) then
                 ! ==> dC13 * Hy(i+1,j,k+1)
                 dCdm=dxi4*dxz(ia)%v(i,j,km)-dzk2*dzz(ia)%v(i,j,km)  &
                      -xi4yjzk*dxy(ia)%v(i,j,km)+dyj4*dyz(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k+1)
					       sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k+1)
					       sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k+1)
					       sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k+1)
               else
                 ! ==> dB8 * Hy(i+1,j,k+1)
                 dCdm=dxi4*dxz(ia)%v(i,j,km)  & ! -dZyy*dxi2*dyj2 !debug
                      -xi4yjzk*dxy(ia)%v(i,j,km)+dyj4*dyz(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k+1)
					       sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k+1)
					       sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k+1)
					       sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k+1)                                 
               endif
                    
             endif         
         
 ! ------------------------- Z components -------------------------------            
         
         
             if(i.ne.1) then
           
               if(j.ne.1) then
                 ! ==> dC22&B13 * Hz(i,j,k)
                 dCdm=xi4zkyj*dxz(ia)%v(i,j,km)-dzk4*dyz(ia)%v(i,j,km)  &
                      -dxi4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i,j,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i,j,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i,j,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i,j,k)
               endif        
           
               ! ==> dC24&B15 * Hz(i,j+1,k)
               dCdm=-dzk4*dyz(ia)%v(i,j,km)-xi4zkyj*dxz(ia)%v(i,j,km)  &
                    +dxi4*dxy(ia)%v(i,j,km)+dyj2*dyy(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i,j+1,k)
					                
             endif
 
             if(i.ne.Nx) then
         
               if(j.ne.1) then
                 ! ==> dC23&B14 * Hz(i+1,j,k)
                 dCdm=xi4zkyj*dxz(ia)%v(i,j,km)+dzk4*dyz(ia)%v(i,j,km)  &
                      -dxi4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i+1,j,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i+1,j,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j,k)
               endif
           
               ! ==> dC25&B16 * Hz(i+1,j+1,k)
               dCdm=dzk4*dyz(ia)%v(i,j,km)-xi4zkyj*dxz(ia)%v(i,j,km)  &
                    +dxi4*dxy(ia)%v(i,j,km)-dyj2*dyy(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
					                    
             endif                
         
             ! sum up h^T*e for one model parameter, row Hx(i,j+1,k+1)
					   up1(ia) = up1(ia) + sumup1 * cInL1%x(i,j+1,k+1)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%x(i,j+1,k+1)
					   up2(ia) = up2(ia) + sumup2 * cInL2%x(i,j+1,k+1)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%x(i,j+1,k+1)
					   
           end do  !ia        

! -------------------------- End of the Fourth Row i,j+1,k+1 --------------------------------- 
          
         end if                

         
         ! ------------------------
         !
         ! >  dC/dm for Hy equations
         !
         ! ------------------------   

           ! ================
           ! row index i,j,k -Hy
           ! ================
           
         if(i.ne.1) then
          
 ! ------------------------- Y components ------------------------------- 

           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
					    
             ! ==> dC17 * Hy(i,j,k)
             dCdm=-dyj4*dxz(ia)%v(i,j,km)+yj4xizk*dxx(ia)%v(i,j,km)  &
                  +yj4zkxi*dzz(ia)%v(i,j,km)
             sumup1 = sumup1 + dCdm * cInR1%y(i,j,k)
					   sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k)
					   sumup2 = sumup2 + dCdm * cInR2%y(i,j,k)
					   sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k)
					   
             ! ==> dC20 * Hy(i,j,k+1)
             dCdm=-yj2xizk*dxx(ia)%v(i,j,km)
             sumup1 = sumup1 + dCdm * cInR1%y(i,j,k+1)
					   sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k+1)
					   sumup2 = sumup2 + dCdm * cInR2%y(i,j,k+1)
					   sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k+1)
         
             if(i.ne.Nx) then
               ! ==> dC18 * Hy(i+1,j,k)
               dCdm=-yj2zkxi*dzz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k)
					     
               ! ==> dC21 * Hy(i+1,j,k+1)
               dCdm=dyj2*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k+1)
             endif 
                  
 ! ------------------------- Z components -------------------------------         
 
             if(j.ne.1) then
               ! ==> dC29 * Hz(i,j,k)
               dCdm=dyj4*dxy(ia)%v(i,j,km)-dxi2*dxx(ia)%v(i,j,km)  &
                    -yj4zkxi*dyz(ia)%v(i,j,km)+dzk4*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i,j,k)
               
               if(i.ne.Nx) then
                 ! ==> dC30 * Hz(i+1,j,k)
                 dCdm=dzk4*dxz(ia)%v(i,j,km)+yj4zkxi*dyz(ia)%v(i,j,km)  &
                      -dyj4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i+1,j,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i+1,j,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j,k)
               endif
                    
             endif
           
             if(j.ne.Ny) then
               ! ==> dC32 * Hz(i,j+1,k)
               dCdm=dyj4*dxy(ia)%v(i,j,km)+dxi2*dxx(ia)%v(i,j,km)  &
                    -yj4zkxi*dyz(ia)%v(i,j,km)-dzk4*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i,j+1,k)
               
               if(i.ne.Nx) then
                 ! ==> dC33 * Hz(i+1,j+1,k)
                 dCdm=-dzk4*dxz(ia)%v(i,j,km)+yj4zkxi*dyz(ia)%v(i,j,km)  &
                      -dyj4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
               endif
                    
             endif

             ! sum up h^T*e for one model parameter, row Hy(i,j,k)
					   up1(ia) = up1(ia) + sumup1 * cInL1%y(i,j,k)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%y(i,j,k)
					   up2(ia) = up2(ia) + sumup2 * cInL2%y(i,j,k)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%y(i,j,k)
					   
           enddo !ia                  
         
! -------------------------- End of the First Row i,j,k --------------------------------- 
         
         endif  
         
           ! ================
           ! row index i+1,j,k -Hy
           ! ================
           
         if(i.ne.Nx) then
         
 ! ------------------------- Y components ------------------------------- 

           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
					    
             ! ==> dC17 * Hy(i+1,j,k)
             dCdm=dyj4*dxz(ia)%v(i,j,km)+yj4xizk*dxx(ia)%v(i,j,km)  &
                  +yj4zkxi*dzz(ia)%v(i,j,km)
             sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k)
					   sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k)
					   sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k)
					   sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k)
         
             if(i.ne.1) then
               ! ==> dC19 * Hy(i,j,k+1)
               dCdm=-dyj2*dxz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%y(i,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k+1)
             endif

             ! ==> dC20 * Hy(i+1,j,k+1)
             dCdm=-yj2xizk*dxx(ia)%v(i,j,km)
             sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k+1)
					   sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k+1)
					   sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k+1)
					   sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k+1)
         
 ! ------------------------- Z components ------------------------------- 
 
             if(j.ne.1) then
         
               if(i.ne.1) then
                 ! ==> dC28 * Hz(i,j,k)
                 dCdm=yj4zkxi*dyz(ia)%v(i,j,km)-dzk4*dxz(ia)%v(i,j,km)  &
                      +dyj4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i,j,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i,j,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i,j,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i,j,k)
               endif
           
               ! ==> dC29 * Hz(i+1,j,k)
               dCdm=-yj4zkxi*dyz(ia)%v(i,j,km)-dzk4*dxz(ia)%v(i,j,km)  &
                    -dyj4*dxy(ia)%v(i,j,km)-dxi2*dxx(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i+1,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i+1,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j,k)
               
             endif

             if(j.ne.Ny) then
         
               if(i.ne.1) then
                 ! ==> dC31 * Hz(i,j+1,k)
                 dCdm=yj4zkxi*dyz(ia)%v(i,j,km)+dzk4*dxz(ia)%v(i,j,km)  &
                      +dyj4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i,j+1,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i,j+1,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i,j+1,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i,j+1,k)
               endif
           
               ! ==> dC32 * Hz(i+1,j+1,k)
               dCdm=-yj4zkxi*dyz(ia)%v(i,j,km)+dzk4*dxz(ia)%v(i,j,km)  &
                    -dyj4*dxy(ia)%v(i,j,km)+dxi2*dxx(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
               
             endif

             ! sum up h^T*e for one model parameter, row Hy(i+1,j,k)
					   up1(ia) = up1(ia) + sumup1 * cInL1%y(i+1,j,k)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%y(i+1,j,k)
					   up2(ia) = up2(ia) + sumup2 * cInL2%y(i+1,j,k)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%y(i+1,j,k)
					 
					 end do !ia                 
         
! -------------------------- End of the Second Row i+1,j,k --------------------------------- 

         endif            
           
           ! ================
           ! row index i,j,k+1 -Hy
           ! ================
           
         if(i.ne.1) then

 ! ------------------------- Y components -------------------------------
 
           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
					               
             if(k.ne.Nz) then
               ! ==> dC17 * Hy(i,j,k+1)
               dCdm=dyj4*dxz(ia)%v(i,j,km)+yj4xizk*dxx(ia)%v(i,j,km)  &
                    +yj4zkxi*dzz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%y(i,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k+1)
               
               if(i.ne.Nx) then
                 ! ==> dC18 * Hy(i+1,j,k+1)
                 dCdm=-yj2zkxi*dzz(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k+1)
					       sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k+1)
					       sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k+1)
					       sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k+1)
               endif
                  
             else
               ! ==> dB12 * Hy(i,j,k+1)
               dCdm=dyj4*dxz(ia)%v(i,j,km)+yj4xizk*dxx(ia)%v(i,j,km) ! +dZxy*grid%delX(i)*dyj2 !debug 
               sumup1 = sumup1 + dCdm * cInR1%y(i,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%y(i,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i,j,k+1)
             endif

 ! ------------------------- Z components -------------------------------

             if(j.ne.1) then
               ! ==> dC23&B14 * Hz(i,j,k)
               dCdm=-dyj4*dxy(ia)%v(i,j,km)+dxi2*dxx(ia)%v(i,j,km)  &
                    +dzk4*dxz(ia)%v(i,j,km)-yj4zkxi*dyz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i,j,k)
					     
               if(i.ne.Nx) then
                 ! ==> dC24&B15 * Hz(i+1,j,k)
                 dCdm=dyj4*dxy(ia)%v(i,j,km)  &
                      +dzk4*dxz(ia)%v(i,j,km)+yj4zkxi*dyz(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i+1,j,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i+1,j,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j,k)
               endif
                   
             endif

             if(j.ne.Ny) then
               ! ==> dC26&B17 * Hz(i,j+1,k)
               dCdm=-dyj4*dxy(ia)%v(i,j,km)-dxi2*dxx(ia)%v(i,j,km)  &
                    -dzk4*dxz(ia)%v(i,j,km)-yj4zkxi*dyz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i,j+1,k)
               
               if(i.ne.Nx) then
                 ! ==> dC27&B18 * Hz(i+1,j+1,k)
                 dCdm=dyj4*dxy(ia)%v(i,j,km)  &
                      -dzk4*dxz(ia)%v(i,j,km)+yj4zkxi*dyz(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
               endif
                   
             endif
     
             ! sum up h^T*e for one model parameter, row Hy(i,j,k+1)
					   up1(ia) = up1(ia) + sumup1 * cInL1%y(i,j,k+1)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%y(i,j,k+1)
					   up2(ia) = up2(ia) + sumup2 * cInL2%y(i,j,k+1)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%y(i,j,k+1)
					    
           end do !ia
         
! -------------------------- End of the Third Row i,j,k+1 --------------------------------- 

         endif            
           
           ! ================
           ! row index i+1,j,k+1 -Hy
           ! ================
           
         if(i.ne.Nx) then
         
 ! ------------------------- Y components -------------------------------

           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
					   
             if(k.ne.Nz) then
               ! ==> dC17 * Hy(i+1,j,k+1)
               dCdm=-dyj4*dxz(ia)%v(i,j,km)+yj4xizk*dxx(ia)%v(i,j,km)  &
                    +yj4zkxi*dzz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k+1)
             else
               ! ==> dB12 * Hy(i+1,j,k+1)
               dCdm=-dyj4*dxz(ia)%v(i,j,km)+yj4xizk*dxx(ia)%v(i,j,km) ! +dZxy*grid%delX(i)*dyj2 !debug 
               sumup1 = sumup1 + dCdm * cInR1%y(i+1,j,k+1)
					     sumdown1 = sumdown1 + dCdm * cInL1%y(i+1,j,k+1)
					     sumup2 = sumup2 + dCdm * cInR2%y(i+1,j,k+1)
					     sumdown2 = sumdown2 + dCdm * cInL2%y(i+1,j,k+1)
             endif

 ! ------------------------- Z components -------------------------------

             if(j.ne.1) then
               
               if(i.ne.1) then
                 ! ==> dC22&B13 * Hz(i,j,k)
                 dCdm=yj4zkxi*dyz(ia)%v(i,j,km)-dzk4*dxz(ia)%v(i,j,km)  &
                      -dyj4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i,j,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i,j,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i,j,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i,j,k)
               endif
         
               ! ==> dC23&B14 * Hz(i+1,j,k)
               dCdm=dyj4*dxy(ia)%v(i,j,km)+dxi2*dxx(ia)%v(i,j,km)  &
                    -dzk4*dxz(ia)%v(i,j,km)-yj4zkxi*dyz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i+1,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i+1,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j,k)
             endif


             if(j.ne.Ny) then
             
               if(i.ne.1) then
                 ! ==> dC25&B16 * Hz(i,j+1,k)
                 dCdm=yj4zkxi*dyz(ia)%v(i,j,km)+dzk4*dxz(ia)%v(i,j,km)  &
                      -dyj4*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i,j+1,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i,j+1,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i,j+1,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i,j+1,k)
               endif
         
               ! ==> dC26&B17 * Hz(i+1,j+1,k)
               dCdm=dyj4*dxy(ia)%v(i,j,km)-dxi2*dxx(ia)%v(i,j,km)  &
                    +dzk4*dxz(ia)%v(i,j,km)-yj4zkxi*dyz(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
             endif

             ! sum up h^T*e for one model parameter, row Hy(i+1,j,k+1)
					   up1(ia) = up1(ia) + sumup1 * cInL1%y(i+1,j,k+1)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%y(i+1,j,k+1)
					   up2(ia) = up2(ia) + sumup2 * cInL2%y(i+1,j,k+1)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%y(i+1,j,k+1)
					   
           end do !ia
         
! -------------------------- End of the Fourth Row i+1,j,k+1 --------------------------------- 

         endif                                       

         
         ! ------------------------
         !
         ! >  dC/dm for Hz equations
         !
         ! ------------------------  
         
           ! ================
           ! row index i,j,k -Hz
           ! ================
           
         if(i.ne.1.and.j.ne.1) then
          
 ! ------------------------- Z components -------------------------------

           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
					    
             ! dC29 * Hz(i,j,k)
             dCdm=-dzk4*dxy(ia)%v(i,j,km)+yj4zkxi*dyy(ia)%v(i,j,km)  &
                  +xi4zkyj*dxx(ia)%v(i,j,km)
             sumup1 = sumup1 + dCdm * cInR1%z(i,j,k)
					   sumdown1 = sumdown1 + dCdm * cInL1%z(i,j,k)
					   sumup2 = sumup2 + dCdm * cInR2%z(i,j,k)
					   sumdown2 = sumdown2 + dCdm * cInL2%z(i,j,k)
         
             if(i.ne.Nx) then 
               ! dC30 * Hz(i+1,j,k)
               dCdm=-yj2zkxi*dyy(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i+1,j,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i+1,j,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j,k)
             endif
         
             if(j.ne.Ny) then
             
               ! dC32 * Hz(i,j+1,k)
               dCdm=-xi2zkyj*dxx(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i,j+1,k)
               
               if(i.ne.Nx) then
                 ! dC33 * Hz(i+1,j+1,k)
                 dCdm=dzk2*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
               endif
         
             endif
         
             ! sum up h^T*e for one model parameter, row Hz(i,j,k)
					   up1(ia) = up1(ia) + sumup1 * cInL1%z(i,j,k)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%z(i,j,k)
					   up2(ia) = up2(ia) + sumup2 * cInL2%z(i,j,k)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%z(i,j,k)
					   
           end do !ia                
         
! -------------------------- End of the First Row i,j,k --------------------------------- 
                   
         endif
         
         
           ! ================
           ! row index i+1,j,k -Hz
           ! ================ 
                   
         if(i.ne.Nx.and.j.ne.1) then
          
 ! ------------------------- Z components -------------------------------          

           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
					            
             ! dC29 * Hz(i+1,j,k)
             dCdm=dzk4*dxy(ia)%v(i,j,km)+yj4zkxi*dyy(ia)%v(i,j,km)  &
                  +xi4zkyj*dxx(ia)%v(i,j,km)
             sumup1 = sumup1 + dCdm * cInR1%z(i+1,j,k)
					   sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j,k)
					   sumup2 = sumup2 + dCdm * cInR2%z(i+1,j,k)
					   sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j,k)
             
             if(j.ne.Ny) then
         
               if(i.ne.1) then
                 ! dC31 * Hz(i,j+1,k)
                 dCdm=-dzk2*dxy(ia)%v(i,j,km)
                 sumup1 = sumup1 + dCdm * cInR1%z(i,j+1,k)
					       sumdown1 = sumdown1 + dCdm * cInL1%z(i,j+1,k)
					       sumup2 = sumup2 + dCdm * cInR2%z(i,j+1,k)
					       sumdown2 = sumdown2 + dCdm * cInL2%z(i,j+1,k)
               endif

               ! dC32 * Hz(i+1,j+1,k)
               dCdm=-xi2zkyj*dxx(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
					       
             endif
                
             ! sum up h^T*e for one model parameter, row Hz(i+1,j,k)
					   up1(ia) = up1(ia) + sumup1 * cInL1%z(i+1,j,k)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%z(i+1,j,k)
					   up2(ia) = up2(ia) + sumup2 * cInL2%z(i+1,j,k)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%z(i+1,j,k)
					   
           end do !ia                 
         
! -------------------------- End of the Second Row i+1,j,k --------------------------------- 
                           
         endif         
         
           ! ================
           ! row index i,j+1,k -Hz
           ! ================         
          
         if(i.ne.1.and.j.ne.Ny) then
          
 ! ------------------------- Z components -------------------------------          

           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
					            
             ! dC29 * Hz(i,j+1,k)
             dCdm=dzk4*dxy(ia)%v(i,j,km)+yj4zkxi*dyy(ia)%v(i,j,km)  &
                  +xi4zkyj*dxx(ia)%v(i,j,km)
             sumup1 = sumup1 + dCdm * cInR1%z(i,j+1,k)
					   sumdown1 = sumdown1 + dCdm * cInL1%z(i,j+1,k)
					   sumup2 = sumup2 + dCdm * cInR2%z(i,j+1,k)
					   sumdown2 = sumdown2 + dCdm * cInL2%z(i,j+1,k)

             if(i.ne.Nx) then
               ! dC30 * Hz(i+1,j+1,k)
               dCdm=-yj2zkxi*dyy(ia)%v(i,j,km)
               sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					     sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					     sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					     sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
             endif 
                   
             ! sum up h^T*e for one model parameter, row Hz(i,j+1,k)
					   up1(ia) = up1(ia) + sumup1 * cInL1%z(i,j+1,k)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%z(i,j+1,k)
					   up2(ia) = up2(ia) + sumup2 * cInL2%z(i,j+1,k)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%z(i,j+1,k)
					   
           end do !ia                        
         
! -------------------------- End of the Third Row i,j+1,k ---------------------------------
                   
         endif
         
           ! ================
           ! row index i+1,j+1,k -Hz
           ! ================   
           
         if(i.ne.Nx.and.j.ne.Ny) then
          
 ! ------------------------- Z components -------------------------------          

           do ia = 1, 6

					   sumup1 = C_ZERO
					   sumup2 = C_ZERO
					   sumdown1 = C_ZERO
					   sumdown2 = C_ZERO
					            
             ! dC29 * Hz(i+1,j+1,k)
             dCdm=-dzk4*dxy(ia)%v(i,j,km)+yj4zkxi*dyy(ia)%v(i,j,km)  &
                  +xi4zkyj*dxx(ia)%v(i,j,km)
             sumup1 = sumup1 + dCdm * cInR1%z(i+1,j+1,k)
					   sumdown1 = sumdown1 + dCdm * cInL1%z(i+1,j+1,k)
					   sumup2 = sumup2 + dCdm * cInR2%z(i+1,j+1,k)
					   sumdown2 = sumdown2 + dCdm * cInL2%z(i+1,j+1,k)
					   
             ! sum up h^T*e for one model parameter, row Hz(i+1,j+1,k)
					   up1(ia) = up1(ia) + sumup1 * cInL1%z(i+1,j+1,k)
					   down1(ia) = down1(ia) + sumdown1 * cInR1%z(i+1,j+1,k)
					   up2(ia) = up2(ia) + sumup2 * cInL2%z(i+1,j+1,k)
					   down2(ia) = down2(ia) + sumdown2 * cInR2%z(i+1,j+1,k)
					   
           end do !ia                  
         
! -------------------------- End of the Fourth Row i+1,j+1,k ---------------------------------
                   
         endif           
                  
           ! final assignment
           do ia = 1, 6
             sum1(ia) = up1(ia) + down1(ia)
             sum2(ia) = up2(ia) + down2(ia)
           enddo
           
           do ia = 1, 3
             cOut(ia)%v(i,j,km) = -real(sum1(ia) + sum2(ia)) * resist(ia)%v(i,j,km)
           enddo
           do ia = 4, 6
             cOut(ia)%v(i,j,km) = -real(sum1(ia) + sum2(ia)) * D2R ! supposed to delete * D2R, but retaining it can suppress artifacts in the inversion model and has obtained better convergence for the field data we tested
           enddo
                           
         end do         
       end do
     end do
  
     ! release memory
     call deall_grid(grid)
     do ia = 1,3
       call deall_rscalar(resist(ia))
     enddo
     do ia = 1,6
       call deall_rscalar(dxx(ia))
       call deall_rscalar(dxy(ia))
       call deall_rscalar(dxz(ia))
       call deall_rscalar(dyy(ia))
       call deall_rscalar(dyz(ia))
       call deall_rscalar(dzz(ia))
     enddo
     
     deallocate(resist,dxx,dxy,dxz,dyy,dyz,dzz)
     
   end subroutine dCTKmultCvector


end module SolverSens
