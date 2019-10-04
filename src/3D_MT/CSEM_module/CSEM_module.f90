module CSEM_module
! A lower lavel module which contains all routines and functions 
! required to construct the right hand side (b0%s) for CSEM case
! Created by Naser Meqbel 04.10.2019

use math_constants
use datafunc
use dataspace
use solnspace
use emsolve3d
use transmitters
!This module uses Kerry Key's 1D code required to construct the b0%s in the secondary field formulation
use Dipole1D

implicit none
   type(cvector), save, public  :: E_p         ! Primary field ---> make it public to check when initlizing the FWD: 
                                               ! if it is already allocated DON'T RECOMPUTE for each iteration in the inversion.
   type(rvector), save, private :: condNomaly  ! Nomalous conductivity (on edge) -->  the 3D model constructed from the 1D model
   type(rvector), save, private :: condAnomaly ! Anomalous conductivity (on edge) --> the 3D model which contains the differance: cond-condNomaly; cond is the current 3D conductivity model    
 Contains
 
subroutine get_source_for_csem(sigma,grid,iTx,source)
 
 type(modelParam_t),intent(in)		:: sigma
 type(grid_t), intent(in)        :: grid 
 integer, intent(in)                        :: iTx
 type (cvector), intent(inout)              :: source
 !Local
 real(kind=prec)	    :: period, omega
 complex(kind=prec)	    :: i_omega_mu
 
 
 !b0%s=i_omega_mu*(sigma-sigma1d)*Ep
 call initilize_1d_vectors(grid)
 call set1DModel(sigma)
 call setAnomConductivity(sigma)
 call comp_dipole1D                  ! Calculate E-Field by Key's code
 call create_Ep(grid)
 
   omega = txDict(iTx)%omega
   period = txDict(iTx)%period
   i_omega_mu = cmplx(0.,-1.0d0*ISIGN*MU_0*omega,kind=prec)
   
   call diagMult(condAnomaly,E_P,source)
   call scMult(i_omega_mu,source,source)
 
 
 
 
end subroutine get_source_for_csem
!#############################################
subroutine create_Ep(grid)

type(grid_t), intent(in)        :: grid  
!Local 
integer ix,iy,iz,counter



    call create_cvector(grid,E_p,EDGE)
 	  counter = 1
	  ! E-field corresponing to these nodes is Ex
	   Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny+1 !Edge Y
			  Do ix = 1,grid%Nx !Center X	  		  	  
				  E_p%x(ix,iy,iz) = ex1D(counter)				  
				  counter = counter + 1
			  End Do
		  End Do
	  End Do
		  
	  ! E-field corresponing to these nodes is Ey
	  Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny !Center y
			  Do ix = 1,grid%Nx+1 !Edge x	  				  
				  E_p%y(ix,iy,iz) = ey1D(counter)
				  counter = counter + 1
			  End Do
		  End Do
	  End Do
	  
	  ! E-field corresponing to these nodes is Ez
	  Do iz = 1,grid%Nz !Center Z
		  Do iy = 1,grid%Ny+1 !Edge y
			  Do ix = 1,grid%Nx+1 !Edge x
				  E_p%z(ix,iy,iz) = jz1D(counter)
				  counter = counter + 1
			  End Do
		  End Do
	  End Do
		
end subroutine create_Ep  
!#############################################
subroutine initilize_1d_vectors(grid)
 type(grid_t), intent(in)        :: grid 
!Local
integer n1D,counter,ix,iy,iz


	  n1D = (grid%Nx)*(grid%Ny+1)*(grid%Nz+1)
	  n1D = n1D + (grid%Nx+1)*(grid%Ny)*(grid%Nz+1)
	  n1D = n1D + (grid%Nx+1)*(grid%Ny+1)*(grid%Nz)

	  if (allocated (x1D)) then  
		 Deallocate(x1D, y1D, z1D)
		 Deallocate(ex1D,ey1D,jz1D)
		 Deallocate(bx1D,by1D,bz1D)
	  end if
	  
	  allocate (x1D(n1D), y1D(n1D), z1D(n1D))
	  allocate (ex1D(n1D),ey1D(n1D),jz1D(n1D))
	  allocate (bx1D(n1D),by1D(n1D),bz1D(n1D))
	
	 
	  
	
	  !====================================================================
	  ! Create position vector that the primary field has to be calculated
	  !====================================================================
	  counter = 1
	  ! E-field corresponing to these nodes is Ex
	  Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny+1 !Edge Y
			  Do ix = 1,grid%Nx !Center X
				x1D(counter) = grid%xCenter(ix)
				y1D(counter) = grid%yEdge(iy)
				z1D(counter) = grid%zEdge(iz)
				counter = counter + 1
				End Do
		  End Do
	  End Do
	  
	  ! E-field corresponing to these nodes is Ey
	  Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny !Center y
			  Do ix = 1,grid%Nx+1 !Edge x
				  x1D(counter) = grid%xEdge(ix)
				  y1D(counter) = grid%yCenter(iy)
				  z1D(counter) = grid%zEdge(iz)
				  counter = counter + 1
			  End Do
		  End Do
	  End Do
	  
	  ! E-field corresponing to these nodes is Ez
	  Do iz = 1,grid%Nz !Center Z
		  Do iy = 1,grid%Ny+1 !Edge y
			  Do ix = 1,grid%Nx+1 !Edge x
				  x1D(counter) = grid%xEdge(ix)
				  y1D(counter) = grid%yEdge(iy)
				  z1D(counter) = grid%zCenter(iz)
				  counter = counter + 1
			  End Do
		  End Do
	  End Do
end subroutine 	initilize_1d_vectors  
!#############################################
subroutine set1DModel(sigma)

   !   this is a private routine, used to extract layer averages from
   !   a 3D conductivity parameter (sigma) and set up
   !   (1) nlay1D    ! Number of layers
   !   (2) sig1D => ! (S/m) Layer conductivities 
   !   (3) zlay1D => ! (m)   Depth to top of each layer, first layer ignored the 1D Model  z_P, sigma_P
  
type(modelParam_t),intent(in)		:: sigma 


   !   local variables ... this is an easy, but not necessarily most efficient
   !   way to get an average background layered conductivity ...
   !    could add routines to modelParameter module to do this more directly
   type(rscalar)	::	 sigmaCell 
   character(len=80)    ::       paramtype,modelFile
   type(rscalar)        ::       model
   type(modelParam_t)   ::       aModel,Anomalous_model,sigma_temp,sigma_temp1
   type(rvector)::  cond


   integer	:: nzEarth,Nz,nzAir,i,j,k,ixTx,iyTx
   real(kind=prec)	:: wt,vAir,asigma
   character(len=256)   ::       PrimaryFile
   


	!   first define conductivity on cells  
   	!   (extract into variable which is public)
   		call modelParamToCell(sigma, sigmaCell, paramtype)
   		nlay1D = sigmaCell%nz
   		nzEarth = sigmaCell%grid%nzEarth
   		nzAir = sigmaCell%grid%nzAir
          
		 ixTx= minNode(xTx1D, sigmaCell%grid%xEdge)  
		 iyTx= minNode(yTx1D, sigmaCell%grid%yEdge)  

   		!   for layer boundaries use z-edges of 3D grid
   		  	 if(allocated(zlay1D)) then
			    Deallocate(zlay1D, sig1D)
		     end if
			 
			 allocate(zlay1D(nlay1D))
			 allocate(sig1D(nlay1D))
			  do k=1,nlay1D
			  zlay1D(k) = sigmaCell%grid%zEdge(k)
			  end do
			 



   		! For create sig1D, we divide this process into two parts (1) for air layers and 
   		!    (2) for earth layers
   		! For air layer, sig1D equal to air layer conductivity
   		! For earth layer, The Geometric mean is be used to create sig1D

   		sig1D(1:nzAir) = sigmaCell%v(1,1,1:nzAir)
   		do k = nzAir+1,nlay1D
 		    sig1D(k) = R_ZERO
 			wt = R_ZERO
       			do i = 1,sigmaCell%grid%Nx
          			do j = 1,sigmaCell%grid%Ny
              				wt = wt + sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
              				sig1D(k) = sig1D(k) + log(sigmaCell%v(i,j,k))* &
 						    sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
          			end do
       			end do           
       			sig1D(k) = sig1D(k)/wt
      			sig1D(k) = exp(sig1D(k)) !sigmaCell%v(1,1,k) !exp(sig1D(k)) 
				write(22,*)k,1.0/sig1D(k),1.0/sigmaCell%v(1,1,k)
   		end do
   		call getValue_modelParam(sigma,paramType,model,vAir)



   
	
   ! Put the background (Primary) "condNomaly" conductivities in ModEM model format
   model%v=R_ZERO
   do k = 1,nzEarth
      asigma = sig1D(k+nzAir)
      if( trim(ParamType) == LOGE) asigma = log(asigma)
      do i = 1,sigmaCell%grid%Nx
         do j = 1,sigmaCell%grid%Ny	    
            model%v(i,j,k) = asigma
         end do
     end do
   end do   
   
 
   
   call copy_modelParam(amodel,sigma)   
   call setType_modelParam(amodel,paramType)
   call setValue_modelParam(amodel,paramType,model,vAir)   
   call ModelParamToEdge(amodel,condNomaly)
	  
  ! clean up
   call deall_modelParam(amodel)
   call deall_rscalar(model)
   call deall_rscalar(sigmaCell)
end subroutine set1DModel
!#########################################################################

subroutine setAnomConductivity(sigma)
   !   This is a private routine that sets anomalous conductivity
   !    in module variable condAnomaly using input model parameter sigma,
   !     and layered background conductivity (already set in module 
   !     variables z_P and sigma_P by a call to setPrimaryCond)

   type(modelParam_t),intent(in), target		:: sigma
   
   type(rvector)::  cond
   real (kind = kind(0.0d0)) :: temp
   type(modelParam_t)                   		:: out_sigma,sigma_temp
   character(len=80)    ::       paramtype,modelFile
   type(rscalar)	::	 sigmaCell 
   integer   k,i,j,ix,iy,iz
   real(kind=prec)	:: vAir

   
  
	  
   ! map conductivity onto edges
   Call ModelParamToEdge(sigma,cond)
   CondAnomaly = subtract_rvector_f(cond , condNomaly )     
   call deall_rvector(cond)
   
end subroutine setAnomConductivity


end module CSEM_module