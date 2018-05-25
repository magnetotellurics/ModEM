Module DC_GeoElec
! The following cod is orginally written by Klaus Spitzer (Freiberg Uni.; see Klaus's following comments)
!I have reorganized the code to fit within the framwork of ModEM (Naser Meqbel, 14.14.2012) 
!=============================================================================c
!                                                                             c
!   A 3D FINITE DIFFERENCE DC RESISTIVITY FORWARD MODELING ALGORITHM FOR      c
!      ARBITRARY SURFACE AND SUBSURFACE SOURCE AND RECEIVER POSITIONS         c
!                                                                             c
!=============================================================================c
!                                                                             c
! Author: Klaus Spitzer                                                       c
!                                                                             c
!         Ecole Polytechnique                                                 c
!         P.O. Box 6079                                                       c
!         Succ. Centre-Ville                                                  c
!         Montreal H3C 3A7                                                    c
!         Canada                                                              c
!                                                                             c
!         Tel.: +1 514 340 4563, Fax.: +1 514 340 3970                        c
!         E-Mail: spitzer@geo.polymtl.ca                                      c
!                                                                             c
!=============================================================================c
! Version: 6.62                                           DATE: APR 6, 2001   c
!                                                                             c
! including updates for single pole and dipole sources at the surface or      c
! subsurface, respectively. A singularity removal procedure is applied for    c
! irregular grid spacings (version 3.0). Modified mixed boundary conditions   c
! are used for all boundaries except the surface boundary reducing the        c
! discretization error considerably. This version is especially designed for  c
! crosshole applications with nsource (= number of sources) sources, e.g.,    c
! along a transmitter borehole. Thus, nsource forward runs are performed.     c
! Nevertheless, it is applicable to any kind of survey using any electrode    c
! configuration. Respectively, the total or anomalous potential is written    c
! to standard output along definable planes (version 6.6). A set of           c
! potentials along predefined receiver positions is additionally written to   c
! file. Versions 6.3 and higher are equipped with the grid-independent        c
! electrode positioning technique. Sources are detached from grid nodes       c
! through the singularity removal procedure, receiver locations are arbitrary c
! through trilinear interpolation of the ANOMALOUS field. Enhanced            c
! flexibility is obtained in conjunction with CREATE version 6.6 that allows  c
! to build a 3D model from 300 sizeable blocks using 62 different             c
! resistivities.                                                              c
!                                                                             c
!=============================================================================c
!                                                                             c
! INPUT FILE:                                                                 c
!                                                                             c
! iprint           : PRINT FLAG FOR OUTPUT OF THE MODEL  0 = NO OUTPUT        c
!                                                        1 = OUTPUT           c
! hv               : PRINT PARAMETER FOR OUTPUT OF THE POTENTIAL ALONG A      c
!                    VERTICAL SECTION THROUGH SOURCE NO. 1 PARALLEL TO THE    c
!                    Y-AXIS (v) OR HORIZONTALLY ALONG THE SURFACE (h)         c
! iflags           : SPECIFIES SOURCE CONFIGURATION, 1 = POLE, 2 = DIPOLE     c 
! nsource          : NUMBER OF SOURCES                                        c
! rec              : RECEIVER POSITIONS                                       c
! im, jm, km       : NUMBER OF NODES IN X-, Y-, AND Z-DIRECTION               c
! comment          : COMMENT ON THE MODEL, MAXIMUM 67 CHARACTERS              c
! rs1, rs2         : SOURCE POSITIONS (POLES OR DIPOLES)                      c
! x(i), y(j), z(k) : COORDINATES OF GRID NODES IN X-, Y-, AND Z-DIRECTION     c
! rhow(i)          : MODEL RESISTIVITIES                                      c
! irhokenn(i)      : RESISTIVITY CODE FOR rhow(i)                             c
!                                                                             c
!                    RESISTIVITIES ASSIGNED TO GRID CELLS                     c
!                    ARE DESCRIBED IN FORM OF A RESISTIVITY CODE,             c
!                    FOR EXAMPLE: rhow(1) =  10                               c
!                                 rhow(2) = 100  .                            c
!                    ONE HORIZONTAL MODEL LAYER IS, E.G.,   1111              c
!                                                           1221              c
!                                                           1221              c
!                                                           1111              c
!                                                                             c
!                    DENOTING A CENTRAL 100 OHM*M BLOCK SURROUNDED BY         c
!                    10 OHM*M CELLS. HORIZONTAL LAYERS ARE ARRANGED           c
!                    CONSECUTIVELY.                                           c
! OUTPUT FILE:                                                                c
! ALL DATA ARE PRINTED TO THE STANDARD OUTPUT. USUALLY THEY ARE REDIRECTED    c
! INTO A LOG FILE *.LOG BY ANY KIND OF COMMAND PROCEDURE. THE POTENTIAL ALONG c
! THE SPECIFIED BOREHOLE (I.E., AT THE DEFINED RECEIVER LOCATIONS) IS         c
! ADDITIONALLY PRINTED ON inputfilename.vbh.                                  c
!                                                                             c
!=============================================================================c  
!
   implicit none
!---------------------------------------------------------------------
! Public data:
!---------------------------------------------------------------------

! Transmitter paramters
	 real(8), dimension(1,3), public    :: rs1,rs2
! Model	 
     integer, public                                 :: imax,jmax,kmax,iflags
	 real(8), dimension(:), allocatable, public      :: x,y,z
	 real(8), dimension(:,:,:), allocatable, public  :: sigmac,v,b_RHS_DC 		 
     character,        public                              :: AdjFWD_DC*80
!---------------------------------------------------------------------
! Private data:
!--------------------------------------------------------------------- 	 
     real(8), dimension(:,:,:), allocatable, private  :: c1,c2,c3,c4,c5,c6,c0,rhoa   
     real(8), dimension(:), allocatable, private      :: bnum,c1num,c2num,c3num,c4num,c5num,c6num,c0num,vnum,d
     real(8), dimension(:), allocatable, private      :: f,g,h
	 real(8), dimension(:), allocatable, private      :: r1,r2,rhow
	 real(8), dimension(:,:), allocatable, private      :: rec
     integer(4), dimension(:), allocatable, private      :: iposc0,iposc1,iposc2,iposc3,iposc4,iposc5,iposc6	 
     integer(4),private                                 :: ima,jma,kma,isurf1,isurf2,k,i,j,im,jm,km,imjmkm,iprint,iz,ipout,nsource
	  real(8), private       :: squer,sigmac1,sigmacn
      character, private  :: comment*80,infile*40,axumfile*40,irhokenn(62)*1,hv*2,totanom*1
        real(8), private :: xstart 
Contains
!############################################################################# Create the RHS for FWD problem
subroutine RHS_DC_FWD

       do 1005,i=1,3
 1005  r1(i)=rs1(1,i)

! CALCULATION OF THE SINGULARITY COEFFICIENTS

      !write(6,'(1x,a)') 'NEUMANN BOUNDARY CONDITIONS AT THE SURFACE'
      do 1000,k=1,kmax
      do 1000,j=1,jmax
      do 1000,i=1,imax
      c1(i,j,k)=-1/f(i-1)/8.*((sigmac(i-1,j-1,k-1)-squer)*g(j-1)*h(k-1) +(sigmac(i-1,j,k-1)-squer)*g(j)*h(k-1)+(sigmac(i-1,j-1,k)-squer) *g(j-1)*h(k)+(sigmac(i-1,j,k)-squer)*g(j)*h(k))
      c2(i,j,k)=-1/f(i)/8.*((sigmac(i,j-1,k-1)-squer)*g(j-1)*h(k-1)+(sigmac(i,j,k-1)-squer)*g(j)*h(k-1)+(sigmac(i,j-1,k)-squer) *g(j-1)*h(k)+(sigmac(i,j,k)-squer)*g(j)*h(k))
      c3(i,j,k)=-1/g(j-1)/8.*((sigmac(i-1,j-1,k-1)-squer)*f(i-1)*h(k-1) +(sigmac(i,j-1,k-1)-squer)*f(i)*h(k-1)+(sigmac(i-1,j-1,k)-squer)  *f(i-1)*h(k)+(sigmac(i,j-1,k)-squer)*f(i)*h(k))
      c4(i,j,k)=-1/g(j)/8.*((sigmac(i-1,j,k-1)-squer)*f(i-1)*h(k-1) +(sigmac(i,j,k-1)-squer)*f(i)*h(k-1)+(sigmac(i-1,j,k)-squer) *f(i-1)*h(k)+(sigmac(i,j,k)-squer)*f(i)*h(k))
      c5(i,j,k)=-1/h(k-1)/8.*((sigmac(i-1,j,k-1)-squer)*f(i-1)*g(j) +(sigmac(i,j,k-1)-squer)*f(i)*g(j)+(sigmac(i-1,j-1,k-1)-squer) *f(i-1)*g(j-1)+(sigmac(i,j-1,k-1)-squer)*f(i)*g(j-1))
      c6(i,j,k)=-1/h(k)/8.*((sigmac(i-1,j,k)-squer)*f(i-1)*g(j) +(sigmac(i,j,k)-squer)*f(i)*g(j)+(sigmac(i-1,j-1,k)-squer) *f(i-1)*g(j-1)+(sigmac(i,j-1,k)-squer)*f(i)*g(j-1))
 1000 c0(i,j,k)=-(c1(i,j,k)+c2(i,j,k)+c3(i,j,k)+c4(i,j,k)+c5(i,j,k)+c6(i,j,k))

!write(6,*)'squer #####', sigmac(1,1,1),squer
!  NEUMANN BOUNDARY CONDITIONS AT THE SURFACE (k=1)
 

      do 11001,j=1,jmax
      do 11001,i=1,imax
      c6(i,j,1)=c5(i,j,1)+c6(i,j,1)
11001 c5(i,j,1)=0.

! MIXED BOUNDARY CONDITIONS AT ALL OTHER BOUNDARIES

1401 call mixedb(c1,c2,c3,c4,c5,c6,c0,imax,jmax,kmax,b_RHS_DC,f,g,h)   
   
     
   !open(10,file=AdjFWD_DC)
      do 1015,k=1,kmax
      do 1015,j=1,jmax
      do 1015,i=1,imax
1015  b_RHS_DC(i,j,k)=-(c1(i,j,k)*vacalc3(x(i-1),y(j),z(k),r1,r2,squer,iflags,isurf1,isurf2) &
                      +c2(i,j,k)*vacalc3(x(i+1),y(j),z(k),r1,r2,squer,iflags,isurf1,isurf2) &
                      +c3(i,j,k)*vacalc3(x(i),y(j-1),z(k),r1,r2,squer,iflags,isurf1,isurf2) &
                      +c4(i,j,k)*vacalc3(x(i),y(j+1),z(k),r1,r2,squer,iflags,isurf1,isurf2) &
                      +c5(i,j,k)*vacalc3(x(i),y(j),z(k-1),r1,r2,squer,iflags,isurf1,isurf2) &
                      +c6(i,j,k)*vacalc3(x(i),y(j),z(k+1),r1,r2,squer,iflags,isurf1,isurf2) &
                      +c0(i,j,k)*vacalc3(x(i),y(j),z(k),r1,r2,squer,iflags,isurf1,isurf2))
     !write(10,*) k,j,i,b_RHS_DC(i,j,k)
    !  close(10)
end subroutine RHS_DC_FWD
!########################################################################
subroutine de_ini_private_data_DC
 if (allocated (c1)) then
	     deallocate (c1,c2,c3,c4,c5,c6, c0,b_RHS_DC,bnum)
         deallocate (d,c1num,c2num,c3num,c4num,c5num,c6num,c0num,rhoa)
         deallocate (v,vnum,rhow,f,g,h,r1,r2,rec)
		 deallocate (iposc0,iposc1,iposc2, iposc3,iposc4,iposc5,iposc6)	  
 end if
 
end subroutine de_ini_private_data_DC
!############################################################################# Rotine to initialize data used within this module
subroutine ini_private_data_DC
	   if (allocated (c1)) then
         deallocate (c1,c2,c3,c4,c5,c6, c0,b_RHS_DC,bnum)
         deallocate (d,c1num,c2num,c3num,c4num,c5num,c6num,c0num,rhoa)
         deallocate (v,vnum,rhow,f,g,h,r1,r2,rec)
		 deallocate (iposc0,iposc1,iposc2, iposc3,iposc4,iposc5,iposc6)	
       end if
       
 	     allocate (c1(imax,jmax,kmax),c2(imax,jmax,kmax),c3(imax,jmax,kmax),c4(imax,jmax,kmax),c5(imax,jmax,kmax),c6(imax,jmax,kmax), c0(imax,jmax,kmax),b_RHS_DC(imax,jmax,kmax),bnum(0:imax*jmax*kmax))
         allocate (d(0:imax*jmax*kmax),c1num(0:imax*jmax*kmax),c2num(0:imax*jmax*kmax),c3num(0:imax*jmax*kmax),c4num(0:imax*jmax*kmax),c5num(0:imax*jmax*kmax),c6num(0:imax*jmax*kmax),c0num(0:imax*jmax*kmax),rhoa(0:imax+1,0:jmax+1,0:kmax+1))
         allocate (v(0:imax+1,0:jmax+1,0:kmax+1),vnum(0:imax*jmax*kmax),rhow(62),f(0:imax),g(0:jmax),h(0:kmax),r1(3),r2(3),rec(1,3))
		 allocate (iposc0(0:imax*jmax*kmax),iposc1(0:imax*jmax*kmax),iposc2(0:imax*jmax*kmax), iposc3(0:imax*jmax*kmax),iposc4(0:imax*jmax*kmax),iposc5(0:imax*jmax*kmax),iposc6(0:imax*jmax*kmax))	    

         im=imax
		 jm=jmax
		 km=kmax
		 
 !################# few flags        
	  iprint=0	 
	 !Pole source
	  iflags=1
	  
	 ! Total potential 
	  totanom='t'
	  
	 ! Plan for output Potential 
	  hv='xy'
	  
	  ! Index of the output slice
	  ipout=1
	  
	  ! Number of sources
	  nsource=1
	  
	   imjmkm=im*jm*km
 
      
!############################# Source stuff ###############     
       do 1005,i=1,3
 1005  r1(i)=rs1(1,i)
       
           if(rs1(1,3).gt.0.0) then 
                      isurf1=1
					  isurf2=1
                      !print*
                     ! print '(1x,a,i4,a)', '* SOURCE SUBSURFACE'
          else if(rs1(1,3).eq.0.0) then
                      !print*
                      !print'(1x,a,i4,a)','* SOURCE AT THE SURFACE'
                      isurf1=0
          else
            stop 'SOURCE CONFIGURATION ERROR'
          endif      
  
!############################# Sigma stuff ###############           
! CALCULATION OF CELL RESISTIVITIES SIGMAC FOR GRID CELLS ALONG  
! THE OUTER BOUNDARIES                                

      do 240,k=1,km-1
      do 240,j=1,jm-1
      sigmac(0,j,k)=sigmac(1,j,k)
  240 sigmac(im,j,k)=sigmac(im-1,j,k)
      do 230,k=1,km-1
      do 230,i=0,im
      sigmac(i,0,k)=sigmac(i,1,k)
  230 sigmac(i,jm,k)=sigmac(i,jm-1,k)
      do 220,j=0,jm
      do 220,i=0,im
      sigmac(i,j,0)=sigmac(i,j,1)
  220 sigmac(i,j,km)=sigmac(i,j,km-1)


! CONDUCTIVITY ASSIGNMENT TEST

      do 260,k=1,km
      do 260,j=1,jm
      do 260,i=1,im
      if(sigmac(i,j,k).eq.0) then
      write(55,*) i,j,k,sigmac(i,j,k) 	  
	  stop 'READ ERROR FOR SIGMA! CHECK GRID AND CONDUCTIVITY ASSIGNMENT!'
	  end if
 260  continue   

      sigmac1=sigmac(1,1,1)
      sigmacn=sigmac(im,jm,km)     
!################################## anomalous conductivity
    ! Compute the anomalous conductivity     
     call calcsquer(imax,jmax,kmax,squer)   
	 !squer=1.0/500.0
!###################################### Grid stuff ######################
  ! CALCULATION OF GRID SPACINGS IN X-DIRECTION (f(i)), Y-DIRECTION (g(j)) AND
! Z-DIRECTION (h(k))
! GRID TEST
      do 401,i=1,im-1
  401 if(x(i+1).le.x(i)) stop 'GRID ERROR IN X'
      do 402,j=1,jm-1
  402 if(y(j+1).le.y(j)) stop 'GRID ERROR IN Y'
      do 403,k=1,km-1
  403 if(z(k+1).le.z(k)) stop 'GRID ERROR IN Z'

      do 400,i=1,im-1
  400 f(i)=dabs(x(i)-x(i+1))
      f(0)=f(1)
      f(im)=f(im-1)
      do 500,j=1,jm-1
  500 g(j)=dabs(y(j)-y(j+1))
      g(0)=g(1)
      g(jm)=g(jm-1)
      do 600,k=1,km-1
  600 h(k)=dabs(z(k)-z(k+1))
      h(0)=h(1)
      h(km)=h(km-1)


! CALCULATION OF x(0,im+1) ,y(0,jm+1) AND z(0,km+1)

      x(0)=x(1)-f(1)
      x(im+1)=x(im)+f(im)
      y(0)=y(1)-g(1)
      y(jm+1)=y(jm)+g(jm)
      z(0)=z(1)-h(1)
      z(km+1)=z(km)+h(km)

! DETERMINATION OF THE POSITIONING VECTORS FOR COMPACT STORAGE

! DETERMINATION OF THE POSITIONING VECTORS FOR COMPACT STORAGE
iposc1=0
iposc2=0
iposc3=0
iposc4=0
iposc5=0
iposc6=0
iposc0=0
      do 2000,iz=1,imjmkm
 2000 iposc0(iz)=iz
      
      do 2100,i=1,jm*km
      do 2100,iz=(i-1)*im+1,i*im-1
      iposc1(iz+1)=iz
 2100 iposc2(iz)=iz+1
      
      do 2200,iz=1,(km-1)*im*jm
 2200 iposc6(iz)=iz+im*jm
      
      do 2300,iz=1+im*jm,imjmkm
 2300 iposc5(iz)=iz-im*jm
      
      do 2400,k=1,km
      do 2400,iz=(k-1)*im*jm+1,(k-1)*im*jm+(jm-1)*im
      iposc4(iz)=iz+im
 2400 iposc3(iz+im)=iz


!################# Start value for the solver
      xstart=0.1
      do 261,k=1,km
      do 261,j=1,jm
      do 261,i=1,im
      if((sigmac(i,j,k)-squer).ge.10e-5) xstart=0.0
 261  continue 

     
end subroutine ini_private_data_DC
!#############################################################################	   
subroutine FWDsolve3D_DC
	   implicit real*8 (a-h,o-z) 






      ima=imax
      jma=jmax
      kma=kmax

      pi=3.141592654d0
	  
	    
  

         im=imax
		 jm=jmax
		 km=kmax
		 




       do 1005,i=1,3
 1005  r1(i)=rs1(1,i)
      

 




!write(6,*)trim(AdjFWD_DC)



 
 do 10002,k=1,km
      do 10002,j=1,jm
      do 10002,i=1,im
      c1(i,j,k)=-1/f(i-1)/8.*(sigmac(i-1,j-1,k-1)*g(j-1)*h(k-1)+sigmac(i-1,j,k-1)*g(j)*h(k-1)+sigmac(i-1,j-1,k)*g(j-1)*h(k)+sigmac(i-1,j,k)*g(j)*h(k))
      c2(i,j,k)=-1/f(i)/8.*(sigmac(i,j-1,k-1)*g(j-1)*h(k-1)+sigmac(i,j,k-1)*g(j)*h(k-1)+sigmac(i,j-1,k)*g(j-1)*h(k) +sigmac(i,j,k)*g(j)*h(k))
      c3(i,j,k)=-1/g(j-1)/8.*(sigmac(i-1,j-1,k-1)*f(i-1)*h(k-1) +sigmac(i,j-1,k-1)*f(i)*h(k-1)+sigmac(i-1,j-1,k)*f(i-1)*h(k)+sigmac(i,j-1,k)*f(i)*h(k))
      c4(i,j,k)=-1/g(j)/8.*(sigmac(i-1,j,k-1)*f(i-1)*h(k-1)+sigmac(i,j,k-1)*f(i)*h(k-1)+sigmac(i-1,j,k)*f(i-1)*h(k)+sigmac(i,j,k)*f(i)*h(k))
      c5(i,j,k)=-1/h(k-1)/8.*(sigmac(i-1,j,k-1)*f(i-1)*g(j)+sigmac(i,j,k-1)*f(i)*g(j)+sigmac(i-1,j-1,k-1)*f(i-1)*g(j-1)+sigmac(i,j-1,k-1)*f(i)*g(j-1))
      c6(i,j,k)=-1/h(k)/8.*(sigmac(i-1,j,k)*f(i-1)*g(j)+sigmac(i,j,k)*f(i)*g(j)+sigmac(i-1,j-1,k)*f(i-1)*g(j-1)+sigmac(i,j-1,k)*f(i)*g(j-1))
10002 c0(i,j,k)=-(c1(i,j,k)+c2(i,j,k)+c3(i,j,k)+c4(i,j,k)+c5(i,j,k)+c6(i,j,k))


      do 1100,j=1,jm
      do 1100,i=1,im
      c6(i,j,1)=c5(i,j,1)+c6(i,j,1)
 1100 c5(i,j,1)=0.

! MIXED BOUNDARY CONDITIONS AT ALL OTHER BOUNDARIES
 1400 call mixedb(c1,c2,c3,c4,c5,c6,c0,im,jm,km,b_RHS_DC,f,g,h)
! RENUMBERING OF v,b AND c

      do 1800,k=1,km
      do 1800,j=1,jm
      do 1800,i=1,im
      c0num((j-1)*im+i+(k-1)*im*jm)=c0(i,j,k)
      c1num((j-1)*im+i+(k-1)*im*jm)=c1(i,j,k)
      c2num((j-1)*im+i+(k-1)*im*jm)=c2(i,j,k)
      c3num((j-1)*im+i+(k-1)*im*jm)=c3(i,j,k)
      c4num((j-1)*im+i+(k-1)*im*jm)=c4(i,j,k)
      c5num((j-1)*im+i+(k-1)*im*jm)=c5(i,j,k)
      c6num((j-1)*im+i+(k-1)*im*jm)=c6(i,j,k)
 1800 bnum((j-1)*im+i+(k-1)*im*jm)=b_RHS_DC(i,j,k)
      
! SYMMETRIZATION OF THE MATRIX
 2750 call sym(c0num,c1num,c2num,c3num,c4num,c5num,c6num,bnum,im,jm,km)
      call cgpc(c0num,c1num,c2num,c3num,c4num,c5num,c6num,iposc0,iposc1,iposc2,iposc3,iposc4,iposc5,iposc6,bnum,vnum,im,jm,km,xstart)

 ! ---- 3D ARRAY NUMBERING OF V

      do 2800,k=1,km
      do 2800,j=1,jm
      do 2800,i=1,im
 2800 v(i,j,k)=vnum(i+(j-1)*im+(k-1)*im*jm)
     
! convert anomalous v to total v
if (trim(AdjFWD_DC)=='FWD_DC') then
      do 2900,k=1,km
      do 2900,j=1,jm
      do 2900,i=1,im
 2900 v(i,j,k)=v(i,j,k)+vacalc3(x(i),y(j),z(k),r1,r2,squer,iflags,isurf1,isurf2)
end if

end subroutine FWDsolve3D_DC

!####################################################33  
subroutine compute_Pri_potential

      do 2900,k=1,km
      do 2900,j=1,jm
      do 2900,i=1,im
 2900 v(i,j,k)= vacalc3(x(i),y(j),z(k),r1,r2,squer,iflags,isurf1,isurf2)

end  subroutine compute_Pri_potential
!#########################################

!=============================================================================c
!                                                                             c
!   SUBROUTINE CALCSQUER CALCULATES A REFERENCE CONDUCTIVITY FOR THE          c 
!   SINGULARITY REMOVAL PROCEDURE IN THE FORM OF AN AVERAGE OVER CELLS        c
!   NEIGHBORING THE SOURCES
!                                                                             c
!=============================================================================c
!                                                                             c
! Author: Klaus Spitzer                                                       c
!                                                                             c
!         Ecole Polytechnique                                                 c
!         P.O. Box 6079                                                       c
!         Succ. Centre-Ville                                                  c
!         Montreal H3C 3A7                                                    c
!         Canada                                                              c
!                                                                             c
!         Tel.: +1 514 340 4563, Fax.: +1 514 340 3970                        c
!         E-Mail: spitzer@geo.polymtl.ca                                      c
!                                                                             c
!=============================================================================c
! Version: *.5                                          DATE:  JUN 16, 1999   c
!=============================================================================c
      subroutine calcsquer(im,jm,km,squer)
      implicit real*8 (a-h,o-z)   

      integer        ::ixtrans,iint,jint,kint,im,jm,km,jytrans,kztrans	
 
     

! Calculation of squer

        do 2201,ixtrans=1,im
        if(x(ixtrans).ge.rs1(1,1)) then
         iint=ixtrans
         goto 2202
        endif
 2201   continue
 2202   continue
        do 2203,jytrans=1,jm
        if(y(jytrans).ge.rs1(1,2)) then
         jint=jytrans
         goto 2204
        endif
 2203   continue
 2204   continue
        do 2205,kztrans=1,km
        if(z(kztrans).ge.rs1(1,3)) then
         kint=kztrans
         goto 2206
        endif
 2205   continue
 2206   continue
        squer=(sigmac(iint,jint,kint)+sigmac(iint,jint-1,kint)+ sigmac(iint,jint,kint-1)+sigmac(iint,jint-1,kint-1)+ sigmac(iint-1,jint,kint)+sigmac(iint-1,jint-1,kint)+ sigmac(iint-1,jint,kint-1)+sigmac(iint-1,jint-1,kint-1))/8.


      return
end subroutine calcsquer
!##################################################################################

subroutine mixedb(c1,c2,c3,c4,c5,c6,c0,im,jm,km,b_RHS_DC,f,g,h)
      implicit real*8 (a-h,o-z)

      !common /dimmax/ imax,jmax,kmax
!      parameter(imax=100,jmax=100,kmax=100)
      integer        ::im,jm,km,k,i,j
      dimension c1(imax,jmax,kmax),c2(imax,jmax,kmax),c3(imax,jmax,kmax),c4(imax,jmax,kmax),c5(imax,jmax,kmax), c6(imax,jmax,kmax),c0(imax,jmax,kmax),b_RHS_DC(imax,jmax,kmax)
      dimension f(0:imax),g(0:jmax),h(0:kmax)
      dimension evec(3),evec1(3),evec2(3)
      real(8), dimension (3):: rs11,rs22


      eps=1.0000001

	  do i=1,3
	  rs11(i)=rs1(1,i)
	  rs22(i)=rs2(1,i)
	  end do
	  
!------------------------------
! BOUNDARY Y BACK (j=1) 
!------------------------------

      j=1

! Flache innen, Kante oben innen

      evec(1)=0.
      evec(2)=-1.
      evec(3)=0.

      do 100,k=1,km-1
      do 100,i=2,im-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*calpha1*(g(1)+g(0))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(g(1)+g(0))
      endif
  100 c4(i,j,k)=c3(i,j,k)+c4(i,j,k)
  

      evec1(1)=-1.
      evec1(2)=0.
      evec1(3)=0.

      i=1
      do 120,k=1,km-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      if(iflags.eq.1) then
       c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*calpha1*(g(1)+g(0))/r1abs-c1(i,j,k)*calpha3*(f(1)+f(0))/r1abs
      else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
       c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*(calpha1/(r1abs**2*sum1r)- calpha2/(r2abs**2*sum1r))*(g(1)+g(0))-c1(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(f(1)+f(0))
      endif
      c4(i,j,k)=c3(i,j,k)+c4(i,j,k)
  120 c2(i,j,k)=c1(i,j,k)+c2(i,j,k)

  
  
  
    evec1(1)=1.
      evec1(2)=0.
      evec1(3)=0.

      i=im
      do 130,k=1,km-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*calpha1*(g(1)+g(0))/r1abs-c2(i,j,k)*calpha3*(f(im)+f(im-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(g(1)+g(0))-c2(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(f(im)+f(im-1))
      endif
      c4(i,j,k)=c3(i,j,k)+c4(i,j,k)
  130 c1(i,j,k)=c1(i,j,k)+c2(i,j,k)  

      evec1(1)=0.
      evec1(2)=0.
      evec1(3)=1.

      k=km
      do 140,i=2,im-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*calpha1*(g(1)+g(0))/r1abs-c6(i,j,k)*calpha3*(h(km)+h(km-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(g(1)+g(0))-c6(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(h(km)+h(km-1))
      endif
      c4(i,j,k)=c3(i,j,k)+c4(i,j,k)
  140 c5(i,j,k)=c5(i,j,k)+c6(i,j,k)


  
  
      evec1(1)=-1.
      evec1(2)=0.
      evec1(3)=0.

      evec2(1)=0.
      evec2(2)=0.
      evec2(3)=1.

      i=1
      k=km    
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha5,calpha6,evec2,r1abs,r2abs)

      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*calpha1*(g(1)+g(0))/r1abs- c1(i,j,k)*calpha3*(f(1)+f(0))/r1abs-c6(i,j,k)*calpha5*(h(km)+h(km-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(g(1)+g(0))-c1(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(f(1)+f(0))-c6(i,j,k)*(calpha5/(r1abs**2*sum1r)-calpha6/(r2abs**2*sum1r))*(h(km)+h(km-1))
      endif
      c4(i,j,k)=c3(i,j,k)+c4(i,j,k)
      c2(i,j,k)=c1(i,j,k)+c2(i,j,k)
      c5(i,j,k)=c5(i,j,k)+c6(i,j,k)


      evec1(1)=1.
      evec1(2)=0.
      evec1(3)=0.

      evec2(1)=0.
      evec2(2)=0.
      evec2(3)=1.

      i=im
      k=km    
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha5,calpha6,evec2,r1abs,r2abs)

      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*calpha1*(g(1)+g(0))/r1abs-c2(i,j,k)*calpha3*(f(im)+f(im-1))/r1abs-c6(i,j,k)*calpha5*(h(km)+h(km-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c3(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(g(1)+g(0))-c2(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(f(im)+f(im-1))-c6(i,j,k)*(calpha5/(r1abs**2*sum1r)-calpha6/(r2abs**2*sum1r))*(h(km)+h(km-1))
      endif
      c4(i,j,k)=c3(i,j,k)+c4(i,j,k)
      c1(i,j,k)=c1(i,j,k)+c2(i,j,k)
      c5(i,j,k)=c5(i,j,k)+c6(i,j,k)

      i=1


      evec(1)=-1.
      evec(2)=0.
      evec(3)=0.

      do 200,k=1,km-1               
      do 200,j=2,jm-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c1(i,j,k)*calpha1*(f(1)+f(0))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c1(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(f(1)+f(0))
      endif
  200 c2(i,j,k)=c1(i,j,k)+c2(i,j,k)


      evec1(1)=0.
      evec1(2)=1.
      evec1(3)=0.

      j=jm
      do 220,k=1,km-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c1(i,j,k)*calpha1*(f(1)+f(0))/r1abs-c4(i,j,k)*calpha3*(g(jm)+g(jm-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c1(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(f(1)+f(0))-c4(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(g(jm)+g(jm-1))
      endif
      c2(i,j,k)=c1(i,j,k)+c2(i,j,k)
  220 c3(i,j,k)=c3(i,j,k)+c4(i,j,k)



! Kante unten innen
!      print*,'       Kante unten innen'

      evec1(1)=0.
      evec1(2)=0.
      evec1(3)=1.

      k=km
      do 230,j=2,jm-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c1(i,j,k)*calpha1*(f(1)+f(0))/r1abs-c6(i,j,k)*calpha3*(h(km)+h(km-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c1(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(f(1)+f(0))-c6(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(h(km)+h(km-1))
      endif
      c2(i,j,k)=c1(i,j,k)+c2(i,j,k)
  230 c5(i,j,k)=c5(i,j,k)+c6(i,j,k)


! Ecke unten vorn


      evec1(1)=0.
      evec1(2)=1.
      evec1(3)=0.

      evec2(1)=0.
      evec2(2)=0.
      evec2(3)=1.

      j=jm
      k=km
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha5,calpha6,evec2,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c1(i,j,k)*calpha1*(f(1)+f(0))/r1abs-c4(i,j,k)*calpha3*(g(jm)+g(jm-1))/r1abs-c6(i,j,k)*calpha5*(h(km)+h(km-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c1(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(f(1)+f(0))-c4(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(g(jm)+g(jm-1))-c6(i,j,k)*(calpha5/(r1abs**2*sum1r)-calpha6/(r2abs**2*sum1r))*(h(km)+h(km-1))
      endif
      c2(i,j,k)=c1(i,j,k)+c2(i,j,k)
      c3(i,j,k)=c3(i,j,k)+c4(i,j,k)
      c5(i,j,k)=c5(i,j,k)+c6(i,j,k)


!----------------------------------
! Rand z unten (k=km)
!---------------------------------- 


      k=km

! Fl�che unten innen

      evec(1)=0.
      evec(2)=0.
      evec(3)=1.

      do 300,j=2,jm-1
      do 300,i=2,im-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c6(i,j,k)*calpha1*(h(km)+h(km-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c6(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(h(km)+h(km-1))
      endif
  300 c5(i,j,k)=c5(i,j,k)+c6(i,j,k)


! Kante vorn innen


      evec1(1)=0.
      evec1(2)=1.
      evec1(3)=0.

      j=jm
      do 320,i=2,im-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c6(i,j,k)*calpha1*(h(km)+h(km-1))/r1abs-c4(i,j,k)*calpha3*(g(jm)+g(jm-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c6(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(h(km)+h(km-1))-c4(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(g(jm)+g(jm-1))
      endif
      c5(i,j,k)=c5(i,j,k)+c6(i,j,k)
  320 c3(i,j,k)=c3(i,j,k)+c4(i,j,k)


      
! Kante rechts innen


      evec1(1)=1.
      evec1(2)=0.
      evec1(3)=0.

      i=im
      do 330,j=2,jm-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c6(i,j,k)*calpha1*(h(km)+h(km-1))/r1abs-c2(i,j,k)*calpha3*(f(im)+f(im-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c6(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(h(km)+h(km-1))-c2(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(f(im)+f(im-1))
      endif
      c5(i,j,k)=c5(i,j,k)+c6(i,j,k)
  330 c1(i,j,k)=c1(i,j,k)+c2(i,j,k)


      
! Ecke rechts vorn


      evec1(1)=0.
      evec1(2)=1.
      evec1(3)=0.

      evec2(1)=1.
      evec2(2)=0.
      evec2(3)=0.

      j=jm
      i=im
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha5,calpha6,evec2,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c6(i,j,k)*calpha1*(h(km)+h(km-1))/r1abs-c4(i,j,k)*calpha3*(g(jm)+g(jm-1))/r1abs-c2(i,j,k)*calpha5*(f(im)+f(im-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c6(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(h(km)+h(km-1))-c4(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(g(jm)+g(jm-1))-c2(i,j,k)*(calpha5/(r1abs**2*sum1r)-calpha6/(r2abs**2*sum1r))*(f(im)+f(im-1))
      endif
      c5(i,j,k)=c5(i,j,k)+c6(i,j,k)
      c3(i,j,k)=c3(i,j,k)+c4(i,j,k)
      c1(i,j,k)=c1(i,j,k)+c2(i,j,k)

!-------------------------------
! Rand y vorn (j=jm)
!-------------------------------


      j=jm

! Fl�che innen und Kante oben innen


      evec(1)=0.
      evec(2)=1.
      evec(3)=0.

      do 400,k=1,km-1
      do 400,i=2,im-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c4(i,j,k)*calpha1*(g(jm)+g(jm-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c4(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(g(jm)+g(jm-1))
      endif
  400 c3(i,j,k)=c3(i,j,k)+c4(i,j,k)


  
      
! Kante rechts innen, Ecke oben


      evec1(1)=1.
      evec1(2)=0.
      evec1(3)=0.

      i=im
      do 440,k=1,km-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha3,calpha4,evec1,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c4(i,j,k)*calpha1*(g(jm)+g(jm-1))/r1abs-c2(i,j,k)*calpha3*(f(im)+f(im-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c4(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(g(jm)+g(jm-1))-c2(i,j,k)*(calpha3/(r1abs**2*sum1r)-calpha4/(r2abs**2*sum1r))*(f(im)+f(im-1))
      endif
      c3(i,j,k)=c3(i,j,k)+c4(i,j,k)
  440 c1(i,j,k)=c1(i,j,k)+c2(i,j,k)


!---------------------------------------
! Rand x rechts (i=im)
!---------------------------------------


      i=im

! Flche innen


      evec(1)=1.
      evec(2)=0.
      evec(3)=0.

      do 500,k=1,km-1
      do 500,j=2,jm-1
      call alfacalc(x(i),y(j),z(k),rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      if(iflags.eq.1) then
      c0(i,j,k)=c0(i,j,k)-c2(i,j,k)*calpha1*(f(im)+f(im-1))/r1abs
       else
       if(r1abs.eq.r2abs) r1abs=r1abs*eps
       sum1r=1./r1abs-1./r2abs
      c0(i,j,k)=c0(i,j,k)-c2(i,j,k)*(calpha1/(r1abs**2*sum1r)-calpha2/(r2abs**2*sum1r))*(f(im)+f(im-1))
      endif
  500 c1(i,j,k)=c1(i,j,k)+c2(i,j,k)



      j=1
      do 610,k=1,km
      do 610,i=1,im
  610 c3(i,j,k)=0.
      i=1
      do 710,k=1,km
      do 710,j=1,jm
  710 c1(i,j,k)=0.
      k=km
      do 810,j=1,jm
      do 810,i=1,im
  810 c6(i,j,k)=0.
      j=jm
      do 910,k=1,km
      do 910,i=1,im
  910 c4(i,j,k)=0.
      i=im
      do 1110,k=1,km
      do 1110,j=1,jm
 1110 c2(i,j,k)=0.
    

      return  
  
  
end subroutine mixedb	  
!=============================================================================c
!                                                                             c
!   SUBROUTINE ALFACALC CALCULATES THE ANGLE BETWEEN THE RADIAL DISTANCE TO   c
!   THE SOURCE AND THE OUTWARD NORMAL                                         c
!                                                                             c
!=============================================================================c
!                                                                             c
! Author: Klaus Spitzer                                                       c
!                                                                             c
!         Ecole Polytechnique                                                 c
!         P.O. Box 6079                                                       c
!         Succ. Centre-Ville                                                  c
!         Montreal H3C 3A7                                                    c
!         Canada                                                              c
!                                                                             c
!         Tel.: +1 514 340 4563, Fax.: +1 514 340 3970                        c
!         E-Mail: spitzer@geo.polymtl.ca                                      c
!                                                                             c
!=============================================================================c
! Version: 1.0                                          DATE:  NOV 06, 1996   c
!=============================================================================c

      subroutine alfacalc (rv1,rv2,rv3,rs11,rs22,calpha1,calpha2,evec,r1abs,r2abs)
      implicit real*8 (a-h,o-z)
      dimension r(3),rs11(3),rs22(3),r1(3),r2(3),evec(3)
      integer l
      r(1)=rv1
      r(2)=rv2
      r(3)=rv3


      do  110,l=1,3
      r1(l)=r(l)-rs11(l)
  110 r2(l)=r(l)-rs22(l)
  
      r1abs=dsqrt(r1(1)**2+r1(2)**2+r1(3)**2)
      r2abs=dsqrt(r2(1)**2+r2(2)**2+r2(3)**2)


      calpha1=(r1(1)*evec(1)+r1(2)*evec(2)+r1(3)*evec(3))/r1abs
      calpha2=(r2(1)*evec(1)+r2(2)*evec(2)+r2(3)*evec(3))/r2abs


      return
      end subroutine alfacalc
!=============================================================================c
!                                                                             c
! FUNCTION VACALC3 FOR CALCULATING THE ANALYTICAL POTENTIAL VALUES FOR THE    c
! HOMOGENEOUS HALFSPACE AND SINGLE POLE/DIPOLE SOURCES AT THE SURFACE         c
! AND SUBSURFACE                                                              c
!                                                                             c
!=============================================================================c
!                                                                             c
! Author: Klaus Spitzer                                                       c
!                                                                             c
!         Ecole Polytechnique                                                 c
!         P.O. Box 6079                                                       c
!         Succ. Centre-Ville                                                  c
!         Montreal H3C 3A7                                                    c
!         Canada                                                              c
!                                                                             c
!         Tel.: +1 514 340 4563, Fax.: +1 514 340 3970                        c
!         E-Mail: spitzer@geo.polymtl.ca                                      c
!                                                                             c
!=============================================================================c
! Version: 4.3                                          DATE:  DEC 08, 1998   c
!=============================================================================c
!                                                                             c
! SWITCHES:                                                                   c
!                                                                             c
! iflags   = 1  single pole source                                            c
!          = 2  dipole source                                                 c
!                                                                             c
! isurf1   = 0  source 1 at the surface                                       c
!          = 1  source 1 subsurface                                           c
!                                                                             c
! isurf2   = 0  source 2 at the surface                                       c
!          = 1  source 2 subsurface                                           c
!                                                                             c
!=============================================================================c

      real*8 function vacalc3(rv1,rv2,rv3,rs11,rs22,sigmaq,iflags,isurf1,isurf2)
      implicit real*8 (a-h,o-z)


      integer  ::iflags,isurf1,isurf2
      dimension r(3),rs11(3),rs22(3)

      pi=asin(cos(0.))*2.


      r(1)=rv1
      r(2)=rv2
      r(3)=rv3


      r1abs=dsqrt((r(1)-rs11(1))**2+(r(2)-rs11(2))**2+(r(3)-rs11(3))**2)
      if(iflags.eq.2) r2abs=dsqrt((r(1)-rs22(1))**2+(r(2)-rs22(2))**2+(r(3)-rs22(3))**2)

      r1strabs=dsqrt((r(1)-rs11(1))**2+(r(2)-rs11(2))**2+(-r(3)-rs11(3))**2)
      if(iflags.eq.2) r2strabs=dsqrt((r(1)-rs22(1))**2+(r(2)-rs22(2))**2+(-r(3)-rs22(3))**2)


      if(iflags.eq.1) then
       if(r1abs.le.1d-6.or.r(3).lt.0) then 
                 vacalc3=30.
                 return
       endif
      endif

      if(iflags.eq.2) then
       if(r1abs.le.1d-6.or.r2abs.le.1d-6.or.r(3).lt.0) then 
                 vacalc3=30.
                 return
       endif
      endif


! SINGLE POLE SOURCE

      if(iflags.eq.1) then

! ------- SURFACE

        if(isurf1.eq.0) vacalc3=1./(2.*pi*sigmaq)*(1./r1abs)

! ------- SUBSURFACE

        if(isurf1.eq.1) vacalc3=1./(4.*pi*sigmaq)*(1./r1abs+1./r1strabs)

!		DIPOLE SOURCES

      else if(iflags.eq.2) then

! ------- BOTH SOURCES AT THE SURFACE

      if(isurf1.eq.0.and.isurf2.eq.0) vacalc3=1./(2.*pi*sigmaq)*(1./r1abs-1./r2abs)

! ------- BOTH SOURCES SUBSURFACE

      if(isurf1.eq.1.and.isurf2.eq.1) vacalc3=1./(4.*pi*sigmaq)*(1./r1abs+1./r1strabs-1./r2abs-1./r2strabs)

! ------- SOURCE 1 AT THE SURFACE, SOURCE 2 SUBSURFACE

      if(isurf1.eq.0.and.isurf2.eq.1) vacalc3=1./(2.*pi*sigmaq)*(1./r1abs-(1./r2abs+1./r2strabs)/2.)

! ------- SOURCE 1 SUBSURFACE, SOURCE 2 AT THE SURFACE

      if(isurf2.eq.0.and.isurf1.eq.1) vacalc3=1./(2.*pi*sigmaq)*((1./r1abs+1./r1strabs)/2.-1./r2abs)

      else 
       stop 'SOURCE CONFIGURATION ERROR IN VACALC3'
      endif


      return
      end  function vacalc3
!#########################################################################################################3
!=============================================================================c
!                                                                             c
! SUBROUTINE SYM FOR SYMMETRIZING THE COEFFICIENT MATRIX                      c
!                                                                             c
!=============================================================================c
!                                                                             c
! Author: Klaus Spitzer                                                       c
!                                                                             c
!         Niedersaechsisches Landesamt fuer Bodenforschung                    c
!         - Geowissenschaftliche Gemeinschaftsaufgaben -                      c
!         Stilleweg 2                                                         c
!         30655 Hannover                                                      c
!         Germany                                                             c
!                                                                             c
!         Tel.: +49 511 643 3531, Fax.: +49 511 643 2304                      c
!         E-Mail: SPITZER@GATE1.BGR.D400.DE                                   c
!                                                                             c
!=============================================================================c
! Version: 1.0                                          DATE:  05.09.94       c
!=============================================================================c

      subroutine sym(c0num,c1num,c2num,c3num,c4num,c5num,c6num,bnum,im,jm,km)
      implicit real*8 (a-h,o-z)

      !common /dimmax/ imax,jmax,kmax
      integer                 :: im,jm,km,i
      dimension bnum(0:imax*jmax*kmax),c1num(0:imax*jmax*kmax),c2num(0:imax*jmax*kmax),c3num(0:imax*jmax*kmax),c4num(0:imax*jmax*kmax),c5num(0:imax*jmax*kmax),c6num(0:imax*jmax*kmax),c0num(0:imax*jmax*kmax)

! SYMMETRIZING THE ROWS 2 TO im

      do 2510,i=2,im
      symfakt=c2num(i-1)/c1num(i)
      c0num(i)=c0num(i)*symfakt
      c1num(i)=c1num(i)*symfakt
      c2num(i)=c2num(i)*symfakt
      c3num(i)=c3num(i)*symfakt
      c4num(i)=c4num(i)*symfakt
      c5num(i)=c5num(i)*symfakt
      c6num(i)=c6num(i)*symfakt
 2510 bnum(i)=bnum(i)*symfakt


! SYMMETRIZING THE ROWS im+1 TO im*jm

      do 2520,i=im+1,im*jm
      symfakt=c4num(i-im)/c3num(i)
      c0num(i)=c0num(i)*symfakt
      c1num(i)=c1num(i)*symfakt
      c2num(i)=c2num(i)*symfakt
      c3num(i)=c3num(i)*symfakt
      c4num(i)=c4num(i)*symfakt
      c5num(i)=c5num(i)*symfakt
      c6num(i)=c6num(i)*symfakt
 2520 bnum(i)=bnum(i)*symfakt


! SYMMETRIZING THE ROWS im*jm+1 TO im*jm*km

      do 2530,i=im*jm+1,im*jm*km
      symfakt=c6num(i-im*jm)/c5num(i)
      c0num(i)=c0num(i)*symfakt
      c1num(i)=c1num(i)*symfakt
      c2num(i)=c2num(i)*symfakt
      c3num(i)=c3num(i)*symfakt
      c4num(i)=c4num(i)*symfakt
      c5num(i)=c5num(i)*symfakt
      c6num(i)=c6num(i)*symfakt
 2530 bnum(i)=bnum(i)*symfakt

      return 
      end  subroutine sym
!###################################################################################
!=============================================================================c
!                                                                             c
! SUBROUTINE CGPC FOR SOLVING THE LINEAR SET OF EQUATIONS BY A                c
! PRECONDITIONED METHOD OF CONJUGATE GRADIENTS  ACCORDING TO SCHWARZ (1991)   c
!                                                                             c
!=============================================================================c
!                                                                             c
! Author: Klaus Spitzer                                                       c
!                                                                             c
!         Niedersaechsisches Landesamt fuer Bodenforschung                    c
!         - Geowissenschaftliche Gemeinschaftsaufgaben -                      c
!         Stilleweg 2                                                         c
!         30655 Hannover                                                      c
!         Germany                                                             c
!                                                                             c
!         Tel.: +49 511 643 3531, Fax.: +49 511 643 2304                      c
!         E-Mail: SPITZER@GATE1.BGR.D400.DE                                   c
!                                                                             c
!=============================================================================c
! Version: 1.0               UPDATE                     DATE:  16.03.95       c
!=============================================================================c
                                                                               
                                                                               
      subroutine cgpc(c0num,c1num,c2num,c3num,c4num,c5num,c6num,iposc0,iposc1,iposc2,iposc3,iposc4,iposc5,iposc6,bnum,xn,im,jm,km,xstart)               
                                                                               
      implicit real*8 (a-h,o-z)                                                


      !parameter (imax=100,jmax=100,kmax=100)
       integer  :: im,jm,km,It,i,imjmkm
      dimension xn(0:imax*jmax*kmax),bnum(0:imax*jmax*kmax)
      dimension  rn(0:imax*jmax*kmax),rnp1(imax*jmax*kmax),pn(0:imax*jmax*kmax)
      dimension  pnm1(0:imax*jmax*kmax),amalx(0:imax*jmax*kmax)
      dimension  amalpn(0:imax*jmax*kmax),xnp1(0:imax*jmax*kmax)
      dimension c0num(0:imax*jmax*kmax),c1num(0:imax*jmax*kmax)
      dimension  c2num(0:imax*jmax*kmax),c3num(0:imax*jmax*kmax)
      dimension  c4num(0:imax*jmax*kmax),c5num(0:imax*jmax*kmax)
      dimension  c6num(0:imax*jmax*kmax),d(0:imax*jmax*kmax),ron(0:imax*jmax*kmax)
    
      integer*4 iposc0(0:imax*jmax*kmax)
      integer*4 iposc1(0:imax*jmax*kmax),iposc2(0:imax*jmax*kmax)
      integer*4 iposc3(0:imax*jmax*kmax),iposc4(0:imax*jmax*kmax)
      integer*4 iposc5(0:imax*jmax*kmax),iposc6(0:imax*jmax*kmax)

      !equivalence (pn(1),pnm1(1)),(rn(1),rnp1(1)),(amalx(1),amalpn(1))
      imjmkm=im*jm*km

      !write(6,'(//1x,a/1x,a/1x,a/1x,a)')'**************************************************************', &
      !'     THE SET OF EQUATIONS IS SOLVED BY THE PRECONDITIONED','               CONJUGATE GRADIENT METHOD CGPC', &
      !          '**************************************************************'

! DETERMINATION OF OMEGA

      omega=1.4

      !print '(/1x,a,f4.2)','SUCCESSIVE OVERRELAXATION FACTOR OMEGA   : ',omega
      !print '(1x,a,f4.2//)','STARTING VALUE FOR THE ITERATION PROCESS : ',xstart

! SCALING OF THE COEFFICIENT MATRIX 
      call skal(c0num,c1num,c2num,c3num,c4num,c5num,c6num,iposc0,iposc1,iposc2,iposc3,iposc4,iposc5,iposc6,bnum,imjmkm,d)
! TEST ON SYMMETRY

!      do 2600,i =1,im*jm*km-1
! 2600 if(dabs(c2num(i)-c1num(i+1)).ge.1.0d-6) print '(1x,a,1x,2d10.5)','CAUTION!! MATRIX NOT SYMMETRIC!',c2num(i), c1num(i+1)
!      do 2700,i =1,im*jm*km-im
! 2700 if(dabs(c4num(i)-c3num(i+im)).ge.1.0d-6) print '(1x,a,1x,2d10.5)','CAUTION!! MATRIX NOT SYMMETRIC!',c4num(i), c3num(i+im)
!      do 2800,i =1,im*jm*km-im*jm
! 2800 if(dabs(c6num(i)-c5num(i+im*jm)).ge.1.0d-6) print '(1x,a,1x,2d10.5)','CAUTION!! MATRIX NOT SYMMETRIC!',c6num(i), c5num(i+im*jm)
 
! FIRST STEP OF THE ITERATION

      do 50,i=1,imjmkm
   50 xn(i)=xstart
      call ab(c0num,c1num,c2num,c3num,c4num,c5num,c6num,iposc0,iposc1,iposc2,iposc3,iposc4,iposc5,iposc6,xn,amalx,imjmkm)
      do 100,i=1,imjmkm
  100 rn(i)=bnum(i)-amalx(i)
      do 125,i=1,imjmkm
  125 pnm1(i)=0
      rhonm1=1.
! BEGIN ITERATION LOOP

      !print '(/1x,a,a)','WORKING ON ...', trim(AdjFWD_DC)

      do 1000,it=1,25000

      do 150,i=1,imjmkm
  150 ron(i)=rn(i)
      call eqsolv(c1num,c2num,c3num,c4num,c5num,c6num,iposc1,iposc2,iposc3,iposc4,iposc5,iposc6,ron,imjmkm,omega)
      rhon=0.
      do 300,i=1,imjmkm
  300 rhon=rhon+rn(i)*ron(i)

      betan=rhon/rhonm1
      do 400,i=1,imjmkm
  400 pn(i)=ron(i)+betan*pnm1(i)
      call ab(c0num,c1num,c2num,c3num,c4num,c5num,c6num,iposc0,iposc1,iposc2,iposc3,iposc4,iposc5,iposc6,pn,amalpn,imjmkm)
      sigman=0.
      do 500,i=1,imjmkm
  500 sigman=sigman+pn(i)*amalpn(i)
      alphan=rhon/sigman

      do 600,i=1,imjmkm
      rnp1(i)=rn(i)-alphan*amalpn(i)
  600 xnp1(i)=xn(i)+alphan*pn(i)

      do 700,i=1,imjmkm
      rhonm1=rhon
      rn(i)=rnp1(i)
      pnm1(i)=pn(i)
  700 xn(i)=xnp1(i)
      rnp1abs=0

! STOPPING CRITERION

      do 800,i=1,imjmkm
  800 rnp1abs=rnp1abs+rnp1(i)*rnp1(i)
      rnp1abs=dsqrt(rnp1abs)
      !if(mod(it,100).eq.0) write(6,'(1x,a3,i5,1x,a8,d11.2,1x,a5,d11.2,1x,a6,d11.2,1x,a7,d11.2,1x,a7,d11.2)') 'it:',it,'rnp1abs:',rnp1abs,'rhon:',rhon,'betan:',betan,'sigman:',sigman,'alphan:',alphan

      if(rnp1abs.lt.1d-10) goto 1200

! END OF ITERATION LOOP

 1000 continue

! RESCALING

 1200 continue
      do 1500,i=1,imjmkm
 1500 xn(i)=xn(i)*d(i)

! OUTPUT OF FINAL RESIDUAL

      !write(6,'(1x,a3,i5,1x,a8,d11.2,1x,a5,d11.2,1x,a6,d11.2,1x,a7,d11.2,1x,a7,d11.2)') 'it:',it,'rnp1abs:',rnp1abs,'rhon:',rhon,'betan:',betan,'sigman:',sigman,'alphan:',alphan

 2000 return
      end subroutine cgpc
!##############################################################################################
!=============================================================================c
!                                                                             c
! SUBROUTINE SKAL FOR SCALING THE COEFFICIENT MATRIX                          c
!                                                                             c
!=============================================================================c
!                                                                             c
! Author: Klaus Spitzer                                                       c
!                                                                             c
!         Niedersaechsisches Landesamt fuer Bodenforschung                    c
!         - Geowissenschaftliche Gemeinschaftsaufgaben -                      c
!         Stilleweg 2                                                         c
!         30655 Hannover                                                      c
!         Germany                                                             c
!                                                                             c
!         Tel.: +49 511 643 3531, Fax.: +49 511 643 2304                      c
!         E-Mail: SPITZER@GATE1.BGR.D400.DE                                   c
!                                                                             c
!=============================================================================c
! Version: 1.0                                          DATE:  05.09.94       c
!=============================================================================c

      subroutine skal(c0num,c1num,c2num,c3num,c4num,c5num,c6num,iposc0,iposc1,iposc2,iposc3,iposc4,iposc5,iposc6,bnum,imjmkm,d)

      implicit real*8 (a-h,o-z)

      !parameter (imax=100,jmax=100,kmax=100)
      integer   :: imjmkm,i
      dimension c0num(0:imax*jmax*kmax),c1num(0:imax*jmax*kmax)
     dimension c2num(0:imax*jmax*kmax),c3num(0:imax*jmax*kmax)
     dimension c4num(0:imax*jmax*kmax),c5num(0:imax*jmax*kmax)
     dimension c6num(0:imax*jmax*kmax),d(0:imjmkm)
     dimension bnum(0:imax*jmax*kmax)
    
      integer*4 iposc0(0:imax*jmax*kmax)
     integer*4 iposc1(0:imax*jmax*kmax),iposc2(0:imax*jmax*kmax)
     integer*4 iposc3(0:imax*jmax*kmax),iposc4(0:imax*jmax*kmax)
     integer*4 iposc5(0:imax*jmax*kmax),iposc6(0:imax*jmax*kmax)




! SCALING OF ALL NON ZERO ELEMENTS

      do  i=1,imjmkm
       if(c0num(i).lt.0.) then
         !print'(/1x,a,i6,a,g15.5/)','NEGATIVE C0 CORRECTED AT NODE NO. ',i,', VALUE: ',c0num(i)
         c0num(i)=dabs(c0num(i))
       endif
       d(i)=1/dsqrt(c0num(i))
      end do
      
       do i=1,imjmkm
        bnum(i)=d(i)*bnum(i)
       c0num(i)=d(i)*c0num(i)*d(iposc0(i))
       c1num(i)=d(i)*c1num(i)*d(iposc1(i))
       c2num(i)=d(i)*c2num(i)*d(iposc2(i))
       c3num(i)=d(i)*c3num(i)*d(iposc3(i))
       c4num(i)=d(i)*c4num(i)*d(iposc4(i))
       c5num(i)=d(i)*c5num(i)*d(iposc5(i))
      c6num(i)=d(i)*c6num(i)*d(iposc6(i))
       end do

             return
 end subroutine skal
!=============================================================================c
!                                                                             c
!    SUBROUTINE AB FOR MATRIX - VECTOR MULTIPLICATION UNDER COMPACT STORAGE   c
!                                                                             c
!=============================================================================c
!                                                                             c
! Author: Klaus Spitzer                                                       c
!                                                                             c
!         Niedersaechsisches Landesamt fuer Bodenforschung                    c
!         - Geowissenschaftliche Gemeinschaftsaufgaben -                      c
!         Stilleweg 2                                                         c
!         30655 Hannover                                                      c
!         Germany                                                             c
!                                                                             c
!         Tel.: +49 511 643 3531, Fax.: +49 511 643 2304                      c
!         E-Mail: SPITZER@GATE1.BGR.D400.DE                                   c
!                                                                             c
!=============================================================================c
! Version: 1.0                                          DATE:  05.09.94       c
!=============================================================================c
                                                                        
     subroutine ab(c0num,c1num,c2num,c3num,c4num,c5num,c6num,iposc0,iposc1,iposc2,iposc3,iposc4,iposc5,iposc6,x,amalx,imjmkm)   

      implicit real*8 (a-h,o-z)                                               
                                                                              
                                               
     integer   imjmkm,i                                                                        
      dimension c0num(0:imax*jmax*kmax),c1num(0:imax*jmax*kmax)
      dimension c2num(0:imax*jmax*kmax),c3num(0:imax*jmax*kmax)
     dimension c4num(0:imax*jmax*kmax),c5num(0:imax*jmax*kmax)
      dimension c6num(0:imax*jmax*kmax)
      dimension x(0:imax*jmax*kmax),amalx(0:imax*jmax*kmax)
      integer*4 iposc0(0:imax*jmax*kmax)
     integer*4  iposc1(0:imax*jmax*kmax),iposc2(0:imax*jmax*kmax)
     integer*4  iposc3(0:imax*jmax*kmax),iposc4(0:imax*jmax*kmax)
     integer*4  iposc5(0:imax*jmax*kmax),iposc6(0:imax*jmax*kmax)
     
! MULTIPLICATION OF THE MATRIX A BY THE VECTOR X 
    do 10,i=1,imjmkm
   10  amalx(i)=c5num(i)*x(iposc5(i))+c3num(i)*x(iposc3(i))+ c1num(i)*x(iposc1(i))+c0num(i)*x(iposc0(i))+c2num(i)*x(iposc2(i))+c4num(i)*x(iposc4(i))+c6num(i)*x(iposc6(i))
       return
      end subroutine ab
!=============================================================================c
!                                                                             c
! SUBROUTINE EQSOLV FOR SOLVING A LINEAR SET OF EQUATIONS USING THE           c
! GAUSS METHOD (in connection with preconditioning)                           c
!                                                                             c
!=============================================================================c
!                                                                             c
! Author: Klaus Spitzer                                                       c
!                                                                             c
!         Niedersaechsisches Landesamt fuer Bodenforschung                    c
!         - Geowissenschaftliche Gemeinschaftsaufgaben -                      c
!         Stilleweg 2                                                         c
!         30655 Hannover                                                      c
!         Germany                                                             c
!                                                                             c
!         Tel.: +49 511 643 3531, Fax.: +49 511 643 2304                      c
!         E-Mail: SPITZER@GATE1.BGR.D400.DE                                   c
!                                                                             c
!=============================================================================c
! Version: 1.0                                          DATE:  05.09.94       c
!=============================================================================c

      subroutine eqsolv(c1num,c2num,c3num,c4num,c5num,c6num,iposc1,iposc2,iposc3,iposc4,iposc5,iposc6,rn,imjmkm,omega)

      implicit real*8 (a-h,o-z)

      !parameter (imax=100,jmax=100,kmax=100)

      dimension rn(0:imax*jmax*kmax),c1num(0:imax*jmax*kmax),c2num(0:imax*jmax*kmax),c3num(0:imax*jmax*kmax)
     dimension c4num(0:imax*jmax*kmax),c5num(0:imax*jmax*kmax)
     dimension c6num(0:imax*jmax*kmax)
     integer imjmkm, I

     integer*4 iposc1(0:imax*jmax*kmax),iposc2(0:imax*jmax*kmax)
     integer*4 iposc3(0:imax*jmax*kmax),iposc4(0:imax*jmax*kmax)
     integer*4 iposc5(0:imax*jmax*kmax),iposc6(0:imax*jmax*kmax)

! FORWARD SUBSTITUTION

      do 100,i=2,imjmkm
      sum=omega*(c1num(i)*rn(iposc1(i))+c3num(i)*rn(iposc3(i))+c5num(i)*rn(iposc5(i)))
  100 rn(i)=(rn(i)-sum)

! BACKSUBSTITUTION

      do 200,i=imjmkm-1,1,-1
      sum=omega*(c2num(i)*rn(iposc2(i))+c4num(i)*rn(iposc4(i))+c6num(i)*rn(iposc6(i)))
  200 rn(i)=(rn(i)-sum)

      return
      end subroutine eqsolv


	  
end module DC_GeoElec
