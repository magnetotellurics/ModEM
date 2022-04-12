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
!use ioascii
!This module uses Kerry Key's 1D code required to construct the b0%s in the secondary field formulation
use Dipole1D
!This module uses Rita Streich's 1D code required to construct the b0%s in the secondary field formulation
use EM1D


implicit none
   type(cvector), save, public  :: E_p         ! Primary field ---> make it public to check when initlizing the FWD: 
                                               ! if it is already allocated DON'T RECOMPUTE for each iteration in the inversion.
   ! Local variables used within this module
   type(rvector), save, private :: condNomaly  ! Nomalous conductivity (on edge) -->  the 3D model constructed from the 1D model
   type(rvector), save, private :: condAnomaly ! Anomalous conductivity (on edge) --> the 3D model which contains the differance: cond-condNomaly; cond is the current 3D conductivity model    
  
  ! The local 1D model which is used in any 1D code
   integer, private                             :: nlay1D_temp    ! Number of layers
   real(8), dimension(:), allocatable, private  :: sig1D_temp     ! (S/m) Layer conductivities 
   real(8), dimension(:), allocatable, private  :: zlay1D_temp    ! (m)   Depth to top of each layer, first layer ignored 
 
   real(8), dimension(:), allocatable, private  :: sig1D_temp_h     ! (S/m) Layer conductivities in the horizontal direction used in VTI
   real(8), dimension(:), allocatable, private  :: sig1D_temp_v     ! (S/m) Layer conductivities in the vertical direction used in VTI   
   type(rvector), save, private :: condNomaly_h  ! Nomalous conductivity (on edge) -->  the 3D model constructed from the 1D model    
   type(rvector), save, private :: condNomaly_v  ! Nomalous conductivity (on edge) -->  the 3D model constructed from the 1D model     
 
   type(rvector), save, private :: condAnomaly_h ! Anomalous conductivity (on edge) --> the 3D model which contains the differance: cond-condNomaly; cond is the current 3D conductivity model    
   type(rvector), save, private :: condAnomaly_v ! Anomalous conductivity (on edge) --> the 3D model which contains the differance: cond-condNomaly; cond is the current 3D conductivity model     
    Contains


    
subroutine get_source_for_csem(sigma,grid,iTx,source)
 	 type(modelParam_t),intent(in)		:: sigma
	 type(grid_t), intent(in)           :: grid 
	 integer, intent(in)                :: iTx
	 type (cvector), intent(inout)      :: source

        !call get_vti(sigma)
        !write(*,*) 'In get 1D'

    write(*,*) trim(node_info), 'Start using ', trim(compute_1D_from),' to solve the 1D FWD problem'
    if (trim(compute_1D_from)=="Dipole1D") then
       call  get_source_for_csem_Dipole1D(sigma,grid,iTx,source)
    elseif (trim(compute_1D_from)=="EM1D") then
       call get_source_for_csem_EM1D(sigma,grid,iTx,source)
    end if
	
    write(*,*) trim(node_info), 'finish using ', trim(compute_1D_from),' to solve the 1D FWD problem'
    
  
    
end subroutine

!#########################################################################
subroutine get_source_for_csem_EM1D(sigma,grid,iTx,source)
	 type(modelParam_t),intent(in)		:: sigma
	 type(grid_t), intent(in)           :: grid 
	 integer, intent(in)                :: iTx
     type (cvector), intent(inout)      :: source

	  !Local
    type(backgrounddata)  :: bgdat      !model description, coordinates, data
    type(sorec)                       :: src    !source specification
    type(sorec),dimension(:),pointer  :: receivers  !receiver specification
    type(freqdata)                    :: freqdat    !frequency dependent specifications
    type(refl_struct)                 :: refl_var   !all variables that have to be remembered while computing 1D fields
  
    integer	:: ifreq,icur,comm,ix,iy,iz
    real(kind=prec)	        :: omega
    complex(kind=prec)	    :: i_omega_mu 
    type(solnVectorMTX_t)      	:: eAll_temp
    type (solnVector_t)     	:: e_temp
    integer status
 
     
 call set1DModel_VTI(sigma,txDict(iTx)%xyzTx(1),txDict(iTx)%xyzTx(2))
 call setAnomConductivity_VTI(sigma)
 
 bgdat%omega = txDict(iTx)%omega 
 bgdat%dowhat = 1
 call create_background_data(grid,bgdat)    
 call create_source_data(iTx,src,freqdat)    
     
  ifreq=1
  icur=1
  comm=1
  refl_var%nzrecHxy=0
  call reflectivity_unified(src,bgdat,refl_var,ifreq,icur,comm) ! Output fiel will be saved in bgdat
  call create_Ep_from_EM1D(grid,bgdat) ! Put the 1D field into E_p
  
  
   omega = txDict(iTx)%omega
   i_omega_mu = cmplx(0.,-1.0d0*ISIGN*MU_0*omega,kind=prec)
   !call diagMult(condAnomaly,E_P,source)
            ! write out EM solutions
   
            call create_solnVectorMTX(1,eAll_temp)
            call create_solnVector(grid,1,e_temp)
            e_temp%pol(1)=E_P
            call copy_solnVector(eAll_temp%solns(1),e_temp) 
         write(6,*) 'Saving the EM solution...'
        ! call write_solnVectorMTX(eAll_temp,'E_P.sol')
 
         
             source%x = condAnomaly_h%x * E_P%x
             source%y = condAnomaly_h%y * E_P%y
             source%z = condAnomaly_v%z * E_P%z
   
            call scMult(i_omega_mu,source,source)    
			
			
! Clean stuff			
call deall(eAll_temp) 


deallocate(bgdat%sigv, STAT=status)
deallocate(bgdat%sigh, STAT=status)
deallocate(bgdat%epsrv, STAT=status)
deallocate(bgdat%epsrh, STAT=status)
deallocate(bgdat%zbound, STAT=status)
! deallocate(bgdat%Expos, STAT=status)
! deallocate(bgdat%Eypos, STAT=status)
! deallocate(bgdat%Ezpos, STAT=status)
!deallocate(bgdat%Ex, STAT=status)
!deallocate(bgdat%Ey, STAT=status)
!deallocate(bgdat%Ez, STAT=status)
deallocate(src%nelem, STAT=status)
deallocate(src%pos, STAT=status)
deallocate(src%ljx, STAT=status)
deallocate(src%ljy, STAT=status)
deallocate(src%ljz, STAT=status)
deallocate(src%akx, STAT=status)
deallocate(src%aky, STAT=status)
deallocate(src%akz, STAT=status)

!clean
 call deall_rvector(CondAnomaly_h)
 call deall_rvector(condAnomaly_v)

end subroutine get_source_for_csem_EM1D
!#########################################################################
subroutine create_Ep_from_EM1D(grid,bgdat)

type(grid_t), intent(in)                :: grid  
type(backgrounddata), intent(in)        :: bgdat 
!Local 
integer ix,iy,iz,counter
write(51,*)grid%nx,grid%ny,grid%nz
      call create_cvector(grid,E_p,EDGE)
 	  counter = 1
	  ! E-field corresponing to these nodes is Ex
	   Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny+1 !Edge Y
			  Do ix = 1,grid%Nx !Center X	  		  	  
				  E_p%x(ix,iy,iz) = bgdat%Ex(counter)
                 			  
				  counter = counter + 1
                 
			  End Do
		  End Do
	  End Do
	 counter = 1	  
	  ! E-field corresponing to these nodes is Ey
	  Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny !Center y
			  Do ix = 1,grid%Nx+1 !Edge x	  				  
				  E_p%y(ix,iy,iz) = bgdat%Ey(counter)				   
				  counter = counter + 1
			  End Do
		  End Do
	  End Do
	   counter = 1
	  ! E-field corresponing to these nodes is Ez
	  Do iz = 1,grid%Nz !Center Z
		  Do iy = 1,grid%Ny+1 !Edge y
			  Do ix = 1,grid%Nx+1 !Edge x
				  E_p%z(ix,iy,iz) = bgdat%Ez(counter)
				  counter = counter + 1
			  End Do
		  End Do
      End Do
      

		
end subroutine create_Ep_from_EM1D  
!#############################################

subroutine create_source_data(iTx,src,freqdat) 
	integer, intent(in)             :: iTx
	type(sorec), intent(inout)      :: src    !source specification
    type(freqdata),intent(inout)     :: freqdat    !frequency dependent specifications
! Local
    integer nelem
    
    
! fill the Freq object with the required information        
   freqdat%nfreq = 1
   !allocate frequency vector
   allocate(freqdat%omega(1), stat=ierr)
   freqdat%omega(1)=txDict(iTx)%omega    !angular frequencies (2*pi*f, f in Hertz)
   
! fill the sources object with the required information    
   src%type=1   !receiver = 0, dipole = 1, wire = 2, star = 3
   src%srcname="Tx"  ! dummy name

  nelem=1
  allocate(src%nelem(nelem), stat=ierr)
  allocate(src%pos(3,nelem), stat=ierr)
  allocate(src%ljx(nelem), stat=ierr)
  allocate(src%ljy(nelem), stat=ierr)
  allocate(src%ljz(nelem), stat=ierr)
 
  allocate(src%akx(nelem), stat=ierr)
  allocate(src%aky(nelem), stat=ierr)
  allocate(src%akz(nelem), stat=ierr)
  

  src%nelem(nelem)= nelem
  src%pos(1,nelem)=txDict(iTx)%xyzTx(1)
  src%pos(2,nelem)=txDict(iTx)%xyzTx(2)
  src%pos(3,nelem)=txDict(iTx)%xyzTx(3)
    
  src%ljx(nelem)=cos(D2R*(txDict(iTx)%azimuthTx))
  src%ljy(nelem)=sin(D2R*(txDict(iTx)%azimuthTx))
  src%ljz(nelem)=0.0
 
  src%akx(nelem)=0.0
  src%aky(nelem)=0.0
  src%akz(nelem)=0.0
  
  src%elsrc = .true. 
  
  allocate(src%cur(src%nelem(1),1),stat=ierr)
  src%cur(1,1)=cmplx(1.0,0.0)
  
  
end subroutine create_source_data
!#########################################################################

subroutine create_background_data(grid,bgdat)  
 type(grid_t), intent(in)        :: grid 
 type(backgrounddata), intent(inout)        :: bgdat 

 !Local
 integer	:: counter,ilay,ix,iy,iz
 integer(kind=int32)             :: nx1,ny1,nz1          !nr of points in my domain for which fields are computed
 
     bgdat%nlay= nlay1D_temp
    !allocate vectors for medium properties
    allocate(bgdat%sigv(bgdat%nlay),bgdat%sigh(bgdat%nlay),bgdat%epsrv(bgdat%nlay),bgdat%epsrh(bgdat%nlay), stat=ierr)
    !if (ierr.ne.0) call alloc_error(pid,'readinput','sig, epsr',ierr)
    !allocate vector for layer boundary depths: 1 element less than nr of layers
    allocate(bgdat%zbound(bgdat%nlay-1),stat=ierr)
    !if (ierr.ne.0) call alloc_error(pid,'readinput','zbound',ierr)
  
   bgdat%rsplmin=50.0
   !bgdat%aniso = vti
   bgdat%aniso = vti !default vti, change to iso if any layer is isotropic!!! OR keep it VTI in all cases: if isotropic case, sig1D_temp_h==sig1D_temp_v 
   
    !write(*,*) 'bgdat%nlay', bgdat%nlay
   do ilay=1,bgdat%nlay
			bgdat%sigh(ilay) = sig1D_temp_h(ilay)
			bgdat%sigv(ilay) = sig1D_temp_v(ilay)
			bgdat%epsrh(ilay) = 1.0
			bgdat%epsrv(ilay) = 1.0
            write(65,*) ilay,sig1D_temp_h(ilay),sig1D_temp_v(ilay)
   end do  
   do ilay=1,bgdat%nlay-1
        bgdat%zbound(ilay)=zlay1D_temp(ilay)
   end do
   
   
   
   	  nx1 = (grid%Nx)*(grid%Ny+1)*(grid%Nz+1)
	  ny1 = (grid%Nx+1)*(grid%Ny)*(grid%Nz+1)
	  nz1 = (grid%Nx+1)*(grid%Ny+1)*(grid%Nz)
    
    bgdat%nExy=	0
    bgdat%nEx = nx1
    bgdat%nEy = ny1
    bgdat%nEz = nz1
	
	bgdat%nHxy=0    
    bgdat%nHx = 0
    bgdat%nHy = 0      
    bgdat%nHz = 0 

    bgdat%allcomp_samecoord = .false.
    bgdat%allzrec_samecoord = .true.
   
  allocate(bgdat%Expos(bgdat%nEx,3),bgdat%Eypos(bgdat%nEy,3),bgdat%Ezpos(bgdat%nEz,3), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'In-backgroundfield','Epos',ierr)
 
 ! allocate(bgdat%Hxpos(bgdat%nHx,3),bgdat%Hypos(bgdat%nHy,3),bgdat%Hzpos(bgdat%nHz,3), stat=ierr)
 ! if (ierr.ne.0) call alloc_error(pid,'In-backgroundfield','Hpos',ierr)
  
  
  allocate(bgdat%Ex(bgdat%nEx),bgdat%Ey(bgdat%nEy),bgdat%Ez(bgdat%nEz), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'Out-backgroundfield','E fields',ierr)    

 ! allocate(bgdat%Hx(nxyz),bgdat%Hy(nxyz),bgdat%Hz(nxyz), stat=ierr)
 ! if (ierr.ne.0) call alloc_error(pid,'Out-backgroundfield','H fields',ierr)  
  
  bgdat%Expos=0
  bgdat%Eypos=0
  bgdat%Ezpos=0
  
    bgdat%Ex = 0._real64
    bgdat%Ey = 0._real64
    bgdat%Ez = 0._real64
 	  !====================================================================
	  ! Create position vector that the primary field has to be calculated
	  !====================================================================
	  counter = 1
	  ! E-field corresponding to these nodes is Ex
	  Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny+1 !Edge Y
			  Do ix = 1,grid%Nx !Center X
				bgdat%Expos(counter,1) = grid%xCenter(ix)
				bgdat%Expos(counter,2) = grid%yEdge(iy)
				bgdat%Expos(counter,3) = grid%zEdge(iz)
				counter = counter + 1
				End Do
		  End Do
      End Do
	  
       counter = 1
	  ! E-field corresponing to these nodes is Ey
	  Do iz = 1,grid%Nz+1 !Edge Z
		  Do iy = 1,grid%Ny !Center y
			  Do ix = 1,grid%Nx+1 !Edge x
				  bgdat%Eypos(counter,1) = grid%xEdge(ix)
                  bgdat%Eypos(counter,2) = grid%yCenter(iy)
                  bgdat%Eypos(counter,3) = grid%zEdge(iz)
				  counter = counter + 1
			  End Do
		  End Do
      End Do
	  
      counter = 1
	  ! E-field corresponing to these nodes is Ez
	  Do iz = 1,grid%Nz !Center Z
		  Do iy = 1,grid%Ny+1 !Edge y
			  Do ix = 1,grid%Nx+1 !Edge x
                  bgdat%Ezpos(counter,1)= grid%xEdge(ix)
                  bgdat%Ezpos(counter,2) = grid%yEdge(iy)
				  bgdat%Ezpos(counter,3) = grid%zCenter(iz)
				  counter = counter + 1
			  End Do
		  End Do
	  End Do 
 
 

end subroutine create_background_data
!#########################################################################

subroutine get_source_for_csem_Dipole1D(sigma,grid,iTx,source)
 
 type(modelParam_t),intent(in)		:: sigma
 type(grid_t), intent(in)        :: grid 
 integer, intent(in)                        :: iTx
 type (cvector), intent(inout)      :: source

!Local
     real(kind=prec)        :: omega
     complex(kind=prec)	    :: i_omega_mu 
 
 ! Get the Transmitter setting:
		xTx1D = txDict(iTx)%xyzTx(1)
		yTx1D = txDict(iTx)%xyzTx(2)
		zTx1D = txDict(iTx)%xyzTx(3)
		ftx1D = 1.0d0/txDict(iTx)%PERIOD
		sdm1D = txDict(iTx)%Moment           ! (Am), dipole moment. Normalize to unit source moment
		azimuthTx1D = txDict(iTx)%azimuthTx ! (degrees) 
		dipTx1D     = txDict(iTx)%dipTx
		write(*,*)trim(node_info), " Tx Azi ", azimuthTx1D
		HTmethod1D      = 'kk_ht_201'    ! Use 201 point HT digital filters.
		outputdomain1D  = 'spatial'      ! Assume spatial domain comps
		lbcomp          = .false.        ! This is changed to true if magnetics in data file
		lUseSpline1D    = .true.         ! Use spline interpolation for faster 1D computations
		linversion      = .false.        ! Compute derivatives with respect to sigma(layers)
		
		phaseConvention = 'lag'          ! The usual default is lag, where phase becomes larger 
										 !    positive values with increasing range.
		lenTx1D         = 00.d0        ! (m) Dipole length 0 = point dipole
		numIntegPts     = 0             ! Number of points to use for Gauss quadrature integration for finite dipole
 
 
 call set1DModel(sigma,xTx1D,yTx1D)
 
 ! Put the privte sig1D_temp and zlay1D_temp into public sig1D and zlay1D required in Dipole1D
	 nlay1D=nlay1D_temp
	 if(allocated(zlay1D)) then
		 Deallocate(zlay1D, sig1D)
	 end if
	allocate(zlay1D(nlay1D),sig1D(nlay1D))
	zlay1D=zlay1D_temp
    sig1D=sig1D_temp
    
 call setAnomConductivity(sigma)
 call initilize_1d_vectors(grid)     ! Initilaize the 1D vectors where to compupte the E field
 call comp_dipole1D                  ! Calculate E-Field by Key's code
 call create_Ep_from_Dipole1D(grid)

    omega = txDict(iTx)%omega
   i_omega_mu = cmplx(0.,-1.0d0*ISIGN*MU_0*omega,kind=prec)
   call diagMult(CondAnomaly_h,E_P,source)
   call scMult(i_omega_mu,source,source)   
   

!clean
 call deall_rvector(CondAnomaly_h)
 
end subroutine get_source_for_csem_Dipole1D
!#############################################


subroutine create_Ep_from_Dipole1D(grid)

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
	  
		 Deallocate(x1D, y1D, z1D)
		 Deallocate(ex1D,ey1D,jz1D)
		 Deallocate(bx1D,by1D,bz1D)
		 
end subroutine create_Ep_from_Dipole1D  
!#############################################
subroutine initilize_1d_vectors(grid)
 type(grid_t), intent(in)        :: grid 
!Local
integer counter,ix,iy,iz


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
subroutine set1DModel(sigma,xTx1D,yTx1D,FromFile)

   !   this is a private routine, used to extract layer averages from
   !   a 3D conductivity parameter (sigma) and set up
   !   (1) nlay1D    ! Number of layers
   !   (2) sig1D => ! (S/m) Layer conductivities 
   !   (3) zlay1D => ! (m)   Depth to top of each layer, first layer ignored the 1D Model  z_P, sigma_P
  
type(modelParam_t),intent(in)		:: sigma 
real(kind=prec),intent(in)                           :: xTx1D,yTx1D 
logical, intent(in), optional                        :: FromFile

   !   local variables ... this is an easy, but not necessarily most efficient
   !   way to get an average background layered conductivity ...
   !    could add routines to modelParameter module to do this more directly
   type(rscalar)	::	 sigmaCell 
   character(len=80)    ::       paramtype,modelFile
   type(rscalar)        ::       model
   type(modelParam_t)   ::       aModel,Anomalous_model
   type(rvector)::  cond


   integer	:: nzEarth,Nz,nzAir,i,j,k,ixTx,iyTx,izTx,counter
   real(kind=prec)	:: wt,vAir,asigma,temp_sigma_value
   character(len=256)   ::       PrimaryFile
   !character(len=20)    ::       get_1D_from
   


	!   first define conductivity on cells  
   	!   (extract into variable which is public)
   		call modelParamToCell(sigma, sigmaCell, paramtype)
   		nlay1D_temp = sigmaCell%nz
   		nzEarth = sigmaCell%grid%nzEarth
   		nzAir = sigmaCell%grid%nzAir
          
		 ixTx= minNode(xTx1D, sigmaCell%grid%xEdge)  
		 iyTx= minNode(yTx1D, sigmaCell%grid%yEdge)
         izTx= minNode(zTx1D, sigmaCell%grid%zEdge)  		 
         !get_1D_from= Trim(solverControl%get_1D_from)
   		!   for layer boundaries use z-edges of 3D grid
   		  	 if(allocated(zlay1D_temp)) then
			    Deallocate(zlay1D_temp, sig1D_temp)
		     end if
			 
			 allocate(zlay1D_temp(nlay1D_temp))
			 allocate(sig1D_temp(nlay1D_temp))
			  do k=1,nlay1D_temp
			  zlay1D_temp(k) = sigmaCell%grid%zEdge(k)
			  end do
			 
   		! For create sig1D, we divide this process into two parts (1) for air layers and 
   		!    (2) for earth layers
   		! For air layer, sig1D equal to air layer conductivity
   		! For earth layer, The Geometric mean is be used to create sig1D

   		sig1D_temp(1:nzAir) = sigmaCell%v(1,1,1:nzAir)

			     if (trim(get_1D_from) =="Geometric_mean") then
				       do k = nzAir+1,nlay1D_temp
						wt = R_ZERO
						temp_sigma_value=R_ZERO
							do i = 1,sigmaCell%grid%Nx
								do j = 1,sigmaCell%grid%Ny
									if (log(sigmaCell%v(i,j,k)) .gt. -20.0) then
										wt = wt + sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
										temp_sigma_value = temp_sigma_value + log(sigmaCell%v(i,j,k))* &
										sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
									end if
										
								end do
                            end do
							   sig1D_temp(k) = exp(temp_sigma_value/wt)
					       write(220,*)k,zlay1D_temp(k),1.0/sig1D_temp(k),sig1D_temp(k),get_1d_from
   		               end do
					elseif (trim(get_1D_from) =="At_Tx_Position") then
					   do k = nzAir+1,nlay1D_temp
					       sig1D_temp(k)=sigmaCell%v(ixTx,iyTx,k)
					       !write(230,*)k,zlay1D(k),1.0/sig1D(k),sig1D(k),get_1d_from
   		               end do						
					elseif (trim(get_1d_from)=="Geometric_mean_around_Tx") then
					    do k = nzAir+1,nlay1D_temp
						  wt = R_ZERO
							do i = ixTx-5,ixTx+5
								do j = iyTx-5,iyTx+5
								    if (log(sigmaCell%v(i,j,k)) .gt. -20.0) then
										wt = wt + sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
										sig1D_temp(k) = sig1D_temp(k) + log(sigmaCell%v(i,j,k))* &
										sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
									 end if
								end do
							end do           
							sig1D_temp(k) = exp(sig1D_temp(k)/wt)	
						  !write(240,*)k,zlay1D(k),1.0/sig1D(k),sig1D(k),get_1d_from
   		                end do
					elseif (trim(get_1d_from)=="Full_Geometric_mean") then
					    wt = R_ZERO
						temp_sigma_value=R_ZERO
						counter=0
						do k = nzAir+1,nlay1D_temp
							do i = 1,sigmaCell%grid%Nx
								do j = 1,sigmaCell%grid%Ny
									if (log(sigmaCell%v(i,j,k)) .gt. -20.0  ) then
								        counter=counter+1
										wt = wt + sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)*sigmaCell%grid%dz(k)
										temp_sigma_value = temp_sigma_value + log(sigmaCell%v(i,j,k))
									end if
								end do
							end do           
   		                end do
                        do k = nzAir+1,nlay1D_temp
						  sig1D_temp(k) = exp(temp_sigma_value/counter)	
						  !write(250,*)k,zlay1D(k),1.0/sig1D(k),sig1D(k),get_1d_from
                        end do	
                    elseif (trim(get_1d_from)=="Fixed_Value") then
					    temp_sigma_value=sigmaCell%v(ixTx,iyTx,nzAir+1) !the value exactly below the Tx
	                     do k = nzAir+1,nlay1D_temp
						  sig1D_temp(k) = temp_sigma_value
						  !write(260,*)k,zlay1D(k),1.0/sig1D(k),sig1D(k),get_1d_from
                        end do					
                    end if	
   		call getValue_modelParam(sigma,paramType,model,vAir)



   
	
   ! Put the background (Primary) "condNomaly" conductivities in ModEM model format
   model%v=R_ZERO
   do k = 1,nzEarth
      asigma = sig1D_temp(k+nzAir)
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
   call ModelParamToEdge(amodel,condNomaly_h)
	  
  ! clean up
   call deall_modelParam(amodel)
   call deall_rscalar(model)
   call deall_rscalar(sigmaCell)
end subroutine set1DModel
!#########################################################################
subroutine set1DModel_VTI(sigma,xTx1D,yTx1D,FromFile)

   !   this is a private routine, used to extract layer averages from
   !   a 3D conductivity parameter (sigma) and set up
   !   (1) nlay1D    ! Number of layers
   !   (2) sig1D => ! (S/m) Layer conductivities 
   !   (3) zlay1D => ! (m)   Depth to top of each layer, first layer ignored the 1D Model  z_P, sigma_P
  
type(modelParam_t),intent(in)		:: sigma 
real(kind=prec),intent(in)          :: xTx1D,yTx1D 
logical, intent(in), optional       :: FromFile

   !   local variables ... this is an easy, but not necessarily most efficient
   !   way to get an average background layered conductivity ...
   !    could add routines to modelParameter module to do this more directly
   type(rscalar)		::	     sigmaCell_h, sigmaCell_v 
   character(len=80)    ::       paramtype,modelFile
   type(rscalar)        ::       model
   type(modelParam_t)   ::       aModel,Anomalous_model
   type(rvector)::  cond


   integer	:: nzEarth,i,j,k,ixTx,iyTx,counter, Nx,Ny,Nz,nzAir
   real(kind=prec)	:: wt,vAir,asigma,temp_sigma_value


   


	    ! first define conductivity on cells: in h and v
        !call get_vti(sigma)
       ! write(*,*) 'In get 1D'
   		call modelParamToCell(sigma, sigmaCell_h, paramtype,cCond_v=sigmaCell_v)
   		
   		! Get the grid spec from either h or v sigmaCell
        Nx=sigmaCell_h%grid%Nx
        Ny=sigmaCell_h%grid%Ny
        Nz=sigmaCell_h%grid%Nz
       
        nzEarth = sigmaCell_h%grid%nzEarth
   		nzAir = sigmaCell_h%grid%nzAir
       !  write(55,*)Nx,Ny,Nz,nzAir
        nlay1D_temp = Nz
        
        ! Get the Tx position (cell #) in X and Y
		 ixTx= minNode(xTx1D, sigmaCell_h%grid%xEdge)  
		 iyTx= minNode(yTx1D, sigmaCell_h%grid%yEdge)  

         
   		! for layer boundaries use z-edges of 3D grid
   		  	 if(allocated(zlay1D_temp)) then
			    Deallocate(zlay1D_temp, sig1D_temp_h,sig1D_temp_v)
		     end if
			 
			 allocate(zlay1D_temp(nlay1D_temp))
			 allocate(sig1D_temp_h(nlay1D_temp))
             allocate(sig1D_temp_v(nlay1D_temp))
			  do k=1,nlay1D_temp
			  zlay1D_temp(k) = sigmaCell_h%grid%zEdge(k)
			  end do
			 
   		! For create sig1D, we divide this process into two parts (1) for air layers and 
   		!    (2) for earth layers
   		! For air layer, sig1D equal to air layer conductivity
   		! For earth layer, The Geometric mean is be used to create sig1D

   		sig1D_temp_h(1:nzAir) = sigmaCell_h%v(1,1,1:nzAir)
        sig1D_temp_v(1:nzAir) = sigmaCell_v%v(1,1,1:nzAir)
        write(99,*)'get_1d_from= ',get_1d_from
			     if (trim(get_1D_from) =="Geometric_mean") then
				       do k = nzAir+1,nlay1D_temp
						wt = R_ZERO
						temp_sigma_value=R_ZERO
							do i = 1,Nx
								do j = 1,Ny
								  if (log(sigmaCell_h%v(i,j,k)) .gt. -20.0 ) then
										wt = wt + sigmaCell_h%grid%dx(i)*sigmaCell_h%grid%dy(j)
										temp_sigma_value = temp_sigma_value + log(sigmaCell_h%v(i,j,k))* &
										sigmaCell_h%grid%dx(i)*sigmaCell_h%grid%dy(j)
								  end if
								end do
                            end do
							   sig1D_temp_h(k) = exp(temp_sigma_value/wt)
                       end do
                       
				       do k = nzAir+1,nlay1D_temp
						wt = R_ZERO
						temp_sigma_value=R_ZERO
							do i = 1,Nx
								do j = 1,Ny
								 if (log(sigmaCell_v%v(i,j,k)) .gt. -20.0  ) then
										wt = wt + sigmaCell_v%grid%dx(i)*sigmaCell_v%grid%dy(j)
										temp_sigma_value = temp_sigma_value + log(sigmaCell_v%v(i,j,k))* &
										sigmaCell_v%grid%dx(i)*sigmaCell_v%grid%dy(j)
								  end if
								end do
                            end do
							   sig1D_temp_v(k) = exp(temp_sigma_value/wt)
                       end do                       
                       
					elseif (trim(get_1D_from) =="At_Tx_Position") then
					   do k = nzAir+1,nlay1D_temp
					       sig1D_temp_h(k)=sigmaCell_h%v(ixTx,iyTx,k)
                       end do
 					   do k = nzAir+1,nlay1D_temp
					       sig1D_temp_v(k)=sigmaCell_v%v(ixTx,iyTx,k)
                       end do  
                        do k = nzAir+1,nlay1D_temp
                            write(70,*) k, sig1D_temp_h(k),sig1D_temp_v(k)
                       end do
					elseif (trim(get_1d_from)=="Geometric_mean_around_Tx") then
					    do k = nzAir+1,nlay1D_temp
						  wt = R_ZERO
							do i = ixTx-5,ixTx+5
								do j = iyTx-5,iyTx+5
									if (log(sigmaCell_h%v(i,j,k)) .gt. -20.0  ) then
											wt = wt + sigmaCell_h%grid%dx(i)*sigmaCell_h%grid%dy(j)
											sig1D_temp_h(k) = sig1D_temp_h(k) + log(sigmaCell_h%v(i,j,k))* &
											sigmaCell_h%grid%dx(i)*sigmaCell_h%grid%dy(j)
									end if
								end do
							end do           
							sig1D_temp_h(k) = exp(sig1D_temp_h(k)/wt)	
                        end do

 					    do k = nzAir+1,nlay1D_temp
						  wt = R_ZERO
							do i = ixTx-5,ixTx+5
								do j = iyTx-5,iyTx+5
									if (log(sigmaCell_v%v(i,j,k)) .gt. -20.0  ) then
											wt = wt + sigmaCell_v%grid%dx(i)*sigmaCell_v%grid%dy(j)
											sig1D_temp_v(k) = sig1D_temp_v(k) + log(sigmaCell_v%v(i,j,k))* &
											sigmaCell_v%grid%dx(i)*sigmaCell_v%grid%dy(j)
									end if
								end do
							end do           
							sig1D_temp_v(k) = exp(sig1D_temp_v(k)/wt)	
                        end do                       
					elseif (trim(get_1d_from)=="Full_Geometric_mean") then
					    wt = R_ZERO
						temp_sigma_value=R_ZERO
						counter=0
						do k = nzAir+1,nlay1D_temp
							do i = 1,Nx
								do j = 1,Ny
									if (log(sigmaCell_h%v(i,j,k)) .gt. -20.0  ) then
											counter=counter+1
											wt = wt + sigmaCell_h%grid%dx(i)*sigmaCell_h%grid%dy(j)*sigmaCell_h%grid%dz(k)
											temp_sigma_value = temp_sigma_value + log(sigmaCell_h%v(i,j,k))
									end if
								end do
							end do           
   		                end do
                        do k = nzAir+1,nlay1D_temp
						  sig1D_temp_h(k) = exp(temp_sigma_value/counter)	
                        end do	

					    wt = R_ZERO
						temp_sigma_value=R_ZERO
						counter=0
						do k = nzAir+1,nlay1D_temp
							do i = 1,Nx
								do j = 1,Ny
									if (log(sigmaCell_v%v(i,j,k)) .gt. -20.0  ) then
											counter=counter+1
											wt = wt + sigmaCell_v%grid%dx(i)*sigmaCell_v%grid%dy(j)*sigmaCell_v%grid%dz(k)
											temp_sigma_value = temp_sigma_value + log(sigmaCell_v%v(i,j,k))
									end if
								end do
							end do           
   		                end do
                        do k = nzAir+1,nlay1D_temp
						  sig1D_temp_v(k) = exp(temp_sigma_value/counter)	
                        end do	                        
                        
                        
                    elseif (trim(get_1d_from)=="Fixed_Value") then
					    temp_sigma_value=sigmaCell_h%v(ixTx,iyTx,nzAir+1) !the value exactly below the Tx
	                     do k = nzAir+1,nlay1D_temp
						  sig1D_temp_h(k) = temp_sigma_value
                         end do	
                         
					    temp_sigma_value=sigmaCell_v%v(ixTx,iyTx,nzAir+1) !the value exactly below the Tx
	                     do k = nzAir+1,nlay1D_temp
						  sig1D_temp_v(k) = temp_sigma_value
                        end do	                         
                    end if	
   		call getValue_modelParam(sigma,paramType,model,vAir)  !It is just to get model structure (place holder)



   
	
   ! Put the background (Primary) h "condNomaly" conductivities in ModEM model format
   model%v=R_ZERO
   do k = 1,nzEarth
      asigma = sig1D_temp_h(k+nzAir)
      if( trim(ParamType) == LOGE) asigma = log(asigma)
      do i = 1,Nx
         do j = 1,Ny	    
            model%v(i,j,k) = asigma
         end do
     end do
   end do   

   call copy_modelParam(amodel,sigma)   
   call setType_modelParam(amodel,paramType)
   call setValue_modelParam(amodel,paramType,model,vAir)   
   call ModelParamToEdge(amodel,condNomaly_h)

   ! Put the background (Primary) v "condNomaly" conductivities in ModEM model format
   model%v=R_ZERO
   do k = 1,nzEarth
      asigma = sig1D_temp_v(k+nzAir)
      if( trim(ParamType) == LOGE) asigma = log(asigma)
      do i = 1,Nx
         do j = 1,Ny	    
            model%v(i,j,k) = asigma
            
         end do
     end do
   end do   

   call copy_modelParam(amodel,sigma)   
   call setType_modelParam(amodel,paramType)
   call setValue_modelParam(amodel,paramType,model,vAir)   
   call ModelParamToEdge(amodel,condNomaly_v)
   
  ! clean up
   call deall_modelParam(amodel)
   call deall_rscalar(model)
   call deall_rscalar(sigmaCell_h)
   call deall_rscalar(sigmaCell_v)
end subroutine set1DModel_VTI
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
   CondAnomaly_h = subtract_rvector_f(cond , condNomaly_h ) 
   call deall_rvector(cond)
   
end subroutine setAnomConductivity
!#########################################################################

subroutine setAnomConductivity_VTI(sigma)
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
   CondAnomaly_h = subtract_rvector_f(cond , condNomaly_h ) 
   CondAnomaly_v = subtract_rvector_f(cond , condNomaly_v )  
   call deall_rvector(cond)
   
end subroutine setAnomConductivity_VTI
!#########################################################################
subroutine get_source_for_csem_EM1D_test(sigma,grid,iTx,source)
	 type(modelParam_t),intent(in)		:: sigma
	 type(grid_t), intent(in)        :: grid 
	 integer, intent(in)                        :: iTx
	 type (cvector), intent(inout)              :: source

	  !Local
  type(backgrounddata)  :: bgdat      !model description, coordinates, data
  type(sorec),dimension(:),pointer  :: src    !source specification
  type(sorec),dimension(:),pointer  :: receivers  !receiver specification
  type(freqdata)                    :: freqdat    !frequency dependent specifications
  type(refl_struct)                 :: refl_var   !all variables that have to be remembered while computing 1D fields
  
    integer	:: nzEarth,Nz,nzAir,i,j,k,ixTx,iyTx,counter,ilay,nelem,nsorec,ifreq,icur,comm,ix,iy,iz,nrec,nfreq,ii
  
  nsorec=1
  allocate(src(nsorec), stat=ierr)
  src(nsorec)%type=1 
  src(nsorec)%srcname='trx'
  
  nelem=1
  allocate(src(nsorec)%nelem(1), stat=ierr)
  src(nsorec)%nelem(1) = nelem
  
  allocate(src(nsorec)%pos(3,nelem), stat=ierr)
  allocate(src(nsorec)%ljx(nelem),src(nsorec)%ljy(nelem),src(nsorec)%ljz(nelem),src(nsorec)%akx(nelem),src(nsorec)%aky(nelem),src(nsorec)%akz(nelem), stat=ierr)
  
  src(nsorec)%pos(1,nelem)=0.0
  src(nsorec)%pos(2,nelem)=0.0
  src(nsorec)%pos(3,nelem)=-0.1

  src(nsorec)%ljx(nelem)=1.0
  src(nsorec)%ljy(nelem)=0.0
  src(nsorec)%ljz(nelem)=0.0
  
  src(nsorec)%akx(nelem)=0.0
  src(nsorec)%aky(nelem)=0.0
  src(nsorec)%akz(nelem)=0.0  
  
  
  src(nsorec)%ncur = 1
  src(nsorec)%elsrc = .true. 
  !src(nsorec)%shiftx = sum(src(nsorec)%pos(1,:)) / nelem 
  
  
!################################### 1D model
     bgdat%nlay=5
    !allocate vectors for medium properties
    allocate(bgdat%sigv(bgdat%nlay),bgdat%sigh(bgdat%nlay),bgdat%epsrv(bgdat%nlay),bgdat%epsrh(bgdat%nlay), stat=ierr)
    if (ierr.ne.0) call alloc_error(pid,'readinput','sig, epsr',ierr)
    !allocate vector for layer boundary depths: 1 element less than nr of layers
    allocate(bgdat%zbound(bgdat%nlay-1),stat=ierr)
    if (ierr.ne.0) call alloc_error(pid,'readinput','zbound',ierr)
  
   bgdat%sigv(1)=0.0
   bgdat%sigv(2)=3.6
   bgdat%sigv(3)=0.333333333333333
   bgdat%sigv(4)=0.01
   bgdat%sigv(5)=1.0
   
   bgdat%sigh(1)=0.0
   bgdat%sigh(2)=3.6
   bgdat%sigh(3)=0.333333333333333
   bgdat%sigh(4)=0.01
   bgdat%sigh(5)=1.0  
  
   bgdat%epsrv(1)=1.0
   bgdat%epsrv(2)=1.0
   bgdat%epsrv(3)=1.0
   bgdat%epsrv(4)=1.0
   bgdat%epsrv(5)=1.0
   
   bgdat%epsrh(1)=1.0
   bgdat%epsrh(2)=1.0
   bgdat%epsrh(3)=1.0
   bgdat%epsrh(4)=1.0
   bgdat%epsrh(5)=1.0
   
   bgdat%aniso = vti
   
   bgdat%zbound(1)=0.0
   bgdat%zbound(2)=300.0
   bgdat%zbound(3)=1100.0
   bgdat%zbound(4)=1200.0
   
  !################################### Receivers 
  nrec=1
  allocate(receivers(nrec), stat=ierr)
  receivers(nrec)%type=0 
  receivers(nrec)%srcname='recgrp'
  
   nelem=1
   allocate(receivers(nrec)%nelem(1), stat=ierr)
   allocate(receivers(nrec)%recnames(nelem), stat=ierr)
   receivers(nrec)%nelem(1) = nelem
   allocate(receivers(nrec)%pos(3,nelem), stat=ierr)
  
   receivers(nrec)%pos(1,nelem)=0.000000 
   receivers(nrec)%pos(2,nelem)=0.000000 
   receivers(nrec)%pos(3,nelem)=299.9
   receivers(nrec)%recnames(nelem)='rr'
   
  !#######################################
   bgdat%rsplmin=50.0
   
   
  !####################################### Freq   
    nfreq=1
    freqdat%nfreq = nfreq
    allocate(freqdat%omega(nfreq), stat=ierr)
    freqdat%omega(nfreq)=0.5 * dtwopi  ! Freq= 0.5 Hz
    
     !####################################### current  
    allocate(src(nsorec)%cur(src(nsorec)%nelem(1),nfreq),stat=ierr)
    src(nsorec)%cur(nelem,nfreq)=cmplx(1.000,0.000)
    

    
  nrec=1
  allocate(bgdat%Exypos(nrec,3), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'init_dat_pos','bgdat%Exypos',ierr)
  allocate(bgdat%Expos(0,3), bgdat%Eypos(0,3), bgdat%Hxpos(0,3), bgdat%Hypos(0,3), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'init_dat_pos','bgdat%Expos etc.',ierr)

  bgdat%Hxypos => bgdat%Exypos
  bgdat%Ezpos => bgdat%Exypos
  bgdat%Hzpos => bgdat%Exypos

  bgdat%Exypos = 0._real64
    
    allocate(bgdat%Ex(nrec),bgdat%Ey(nrec),bgdat%Ez(nrec),bgdat%Hx(nrec),bgdat%Hy(nrec),bgdat%Hz(nrec), stat=ierr)
    if (ierr.ne.0) call alloc_error(pid,'init_dat_pos','bgdat%Ex etc.',ierr)
    bgdat%Ex = 0._real64
    bgdat%Ey = 0._real64
    bgdat%Ez = 0._real64
    bgdat%Hx = 0._real64
    bgdat%Hy = 0._real64
    bgdat%Hz = 0._real64   
    
    bgdat%allcomp_samecoord = .true.
    bgdat%allzrec_samecoord = .false.


  !set positions for first (and in most cases the only) receiver group
  nrec = receivers(1)%nelem(1)
  do ii = 1,nrec
    bgdat%Exypos(ii,:) = receivers(1)%pos(:,ii)
  enddo
  bgdat%nExy = nrec
  bgdat%nEx = 0
  bgdat%nEy = 0
  bgdat%nEz = nrec
  bgdat%nHxy = nrec
  bgdat%nHx = 0
  bgdat%nHy = 0
  bgdat%nHz = nrec  
  
  ifreq=1
  icur=1
  comm=1
write(*,*)' nrec',nrec, bgdat%Exypos(1,1), bgdat%Exypos(1,2) , bgdat%Exypos(1,3) ,freqdat%omega(ifreq)

  
  

     bgdat%omega = freqdat%omega(ifreq)
     bgdat%dowhat = 1
        call reflectivity_unified(src(nsorec),bgdat,refl_var,ifreq,icur,comm)

        !write field data for this current
        WRITE(6, *) ' bgdat%Ex(1) ', bgdat%Ex(1)
  
end subroutine get_source_for_csem_EM1D_test


end module CSEM_module