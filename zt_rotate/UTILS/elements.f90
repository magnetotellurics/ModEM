! *****************************************************************************
module elements
	! This module is designed to contain the basic *3D spherical staggered-grid*
	! specific subroutines calculating the line and area elements (edge lengths
	! and face areas for staggered-grid cells, respectively), "order numbers" 
	! for EM-related vectors, etc

  implicit none


Contains

  ! ***************************************************************************
  ! * following are basic subroutines previously found in program earth, used
  ! * by the forward modelling subroutines


      subroutine volume_vijk(i,j,k,ph,th,r,vijk)

!-------------------------------------------------------------
! subroutine for calculating primary cell volume element Vijk
!-------------------------------------------------------------
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!   -----------------------------------------------------

      implicit none

      integer                        :: i,j,k
      real(8),dimension(:)           :: ph
      real(8),dimension(:)           :: th
      real(8),dimension(:)           :: r
      real(8)                        :: vijk,r1,th1

      th1=(th(j  )+th(j+1))/2.d0
       r1=( r(k  )+ r(k+1))/2.d0

      vijk=r1*dsin( th1 )*abs( (r(k)**2)-(r(k+1)**2) )*abs( th(j+1)-th(j) )*ph(i)/2.d0

      return
      end subroutine volume_vijk



      subroutine volume_vijk2(i1,i2,j,k,ph,th,r,vijk2)

!-------------------------------------------------------------
! subroutine for calculating dual cell volume element Vijk_2
!-------------------------------------------------------------
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!   -----------------------------------------------------

      implicit none

      integer                        :: i1,i2,j,k
      real(8),dimension(:)           :: ph
      real(8),dimension(:)           :: th
      real(8),dimension(:)           :: r
      real(8)                        :: vijk2,th1,th2,r1,r2,dph

      dph=(ph(i1 )+ph(i2 ))/2.d0
      th1=(th(j  )+th(j+1))/2.d0
      th2=(th(j+1)+th(j+2))/2.d0
       r1=( r(k  )+ r(k+1))/2.d0
       r2=( r(k+1)+ r(k+2))/2.d0

      vijk2=r(k+1)*dsin( th(j+1) )*abs( (r1**2)-(r2**2) )*abs( th2-th1 )*dph/2.d0

      return
      end subroutine volume_vijk2


      subroutine area_sijk(j,k,th,r,sijk)
 
!-------------------------------------------------------------
!      subroutine for calculating area element S(i)jk
!-------------------------------------------------------------
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!   -----------------------------------------------------

      implicit none

      integer                        :: j,k
      real(8),dimension(:)	         :: th
      real(8),dimension(:)	         :: r
      real(8)                        :: sijk
 
      sijk=abs( (r(k)**2)-(r(k+1)**2) )*abs( th(j+1)-th(j) )/2.d0
 
      return
      end subroutine area_sijk



      subroutine area_sijk2(j,k,th,r,sijk2)
 
!-------------------------------------------------------------
!      subroutine for calculating area element S(i)jk_2
!
!-------------------------------------------------------------
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!   -----------------------------------------------------

      implicit none

      integer						 :: j,k
      real(8),dimension(:)			 :: th
      real(8),dimension(:)			 :: r
      real(8)						 :: sijk2,th1,th2,r1,r2
 
      th1=(th(j  )+th(j+1))/2.d0
      th2=(th(j+1)+th(j+2))/2.d0
       r1=( r(k  )+ r(k+1))/2.d0
       r2=( r(k+1)+ r(k+2))/2.d0
      sijk2=abs( (r1**2)-(r2**2) )*abs( th2-th1 )/2.d0
 
      return
      end subroutine area_sijk2



      subroutine area_sjki(i,j,k,ph,th,r,sjki)
 
!-------------------------------------------------------------
!      subroutine for calculating area element Sjki
!-------------------------------------------------------------
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!   -----------------------------------------------------

      implicit none

      integer                      :: i,j,k
      real(8),dimension(:)		   :: ph
      real(8),dimension(:)         :: th
      real(8),dimension(:)         :: r
      real(8) sjki
 
      sjki=abs( (r(k)**2)-(r(k+1)**2) )* dsin( th(j) )*ph(i)/2.0d0
 
      return
      end subroutine area_sjki



      subroutine area_sjki2(i1,i2,j,k,ph,th,r,sjki2)
 
!-------------------------------------------------------------
!      subroutine for calculating area element Sjki_2
!-------------------------------------------------------------
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!   -----------------------------------------------------

      implicit none

      integer                          :: i1,i2,j,k
      real(8),dimension(:)			   :: ph
      real(8),dimension(:)             :: th
      real(8),dimension(:)             :: r
      real(8)                          :: sjki2
      real(8)                          :: dph,th1,r1,r2
 
      dph = ( ph(i1 )+ph(i2 ) )/2.0d0
      th1 = ( th(j  )+th(j+1) )/2.0d0
       r1 = (  r(k  )+ r(k+1) )/2.0d0
       r2 = (  r(k+1)+ r(k+2) )/2.0d0
      sjki2=abs( (r1**2)-(r2**2) )* dsin( th1 )*dph/2.0d0
 
      return
      end subroutine area_sjki2



      subroutine area_skij(i,j,k,ph,th,r,skij)
 
!-------------------------------------------------------------
!      subroutine for calculating area element Skij
!-------------------------------------------------------------
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!   -----------------------------------------------------

      implicit none

      integer						:: i,j,k
      real(8),dimension(:)			:: ph
      real(8),dimension(:)			:: th
      real(8),dimension(:)			:: r
      real(8)						:: skij
 
      skij=(r(k)**2)*abs( dcos( th(j) )- dcos( th(j+1)) )*ph(i)
 
      return
      end subroutine area_skij



      subroutine area_skij2(i1,i2,j,k,ph,th,r,skij2)
 
!-------------------------------------------------------------
!      subroutine for calculating area element Skij_2
!-------------------------------------------------------------
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!   -----------------------------------------------------

      implicit none

      integer                       :: i1,i2,j,k
      real(8),dimension(:)			:: ph
      real(8),dimension(:)			:: th
      real(8),dimension(:)			:: r
      real(8)                       :: skij2
      real(8)                       :: r1,th1,th2,dph
 
      dph = ( ph(i1)+ph(i2) )/2.d0
      th1 = ( th(j  )+th(j+1) )/2.d0
      th2 = ( th(j+1)+th(j+2) )/2.d0
       r1 = (  r(k  )+ r(k+1) )/2.d0
      skij2=(r1**2)*abs( dcos( th1 )- dcos( th2 ) )*dph
 
      return
      end subroutine area_skij2

! *****************************************************************************

      subroutine leng_xijk(i,j,k,ph,th,r,xijk)

!-------------------------------------------------------------
!      subroutine for calculating line element Xijk
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!-------------------------------------------------------------

      implicit none

      integer                         :: i,j,k
      real(8),dimension(:)			  :: ph
      real(8),dimension(:)			  :: th
      real(8),dimension(:)			  :: r
      real(8)                         :: xijk
 
      xijk=r(k)*dsin( th(j) )*ph(i)
 
      return
      end subroutine leng_xijk


      subroutine leng_xijk2(i1,i2,j,k,ph,th,r,xijk2)

!-------------------------------------------------------------
!      subroutine for calculating line element Xijk2
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!-------------------------------------------------------------

      implicit none

      integer                         :: i1,i2,j,k
      real(8),dimension(:)            :: ph
      real(8),dimension(:)            :: th
      real(8),dimension(:)            :: r
      real(8)                         :: xijk2,r1,th1,dph

      dph = ( ph(i1 )+ph(i2 ) )/2.0d0
      th1 = ( th(j  )+th(j+1) )/2.0d0
       r1 = (  r(k  )+ r(k+1) )/2.0d0

      xijk2=r1*dsin( th1 )*dph

      return
      end subroutine leng_xijk2

 
      subroutine leng_xijk2_shifted(i,j,k,ph,th,r,xijk2)

!-------------------------------------------------------------
!      subroutine for calculating line element Xijk2  - A.K.
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!-------------------------------------------------------------

      implicit none

      integer                         :: i,j,k
      real(8),dimension(:)			  :: ph
      real(8),dimension(:)			  :: th
      real(8),dimension(:)			  :: r
      real(8)                         :: xijk2,r1,th1
 
	  th1=( th(j  )+th(j+1) )/2.d0
       r1=(  r(k  )+ r(k+1) )/2.d0

      xijk2=r1*dsin( th1 )*ph(i)
 
      return
      end subroutine leng_xijk2_shifted


      subroutine leng_yijk(j,k,th,r,yijk)

!-------------------------------------------------------------
!      subroutine for calculating line element Y(i)jk
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!-------------------------------------------------------------

      implicit none

      integer								:: j,k
      real(8),dimension(:)					:: th
      real(8),dimension(:)					:: r
      real(8)								:: yijk
 
      yijk=r(k)*abs( th(j+1)-th(j) )
 
      return
      end subroutine leng_yijk


      subroutine leng_yijk2(j,k,th,r,yijk2)

!-------------------------------------------------------------
!      subroutine for calculating line element Yijk2
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!-------------------------------------------------------------

      implicit none

      integer                               :: j,k
      real(8),dimension(:)                  :: th
      real(8),dimension(:)                  :: r
      real(8)                               :: yijk2,r1,th1,th2

      th1 = ( th(j  )+th(j+1) )/2.0d0
      th2 = ( th(j+1)+th(j+2) )/2.0d0
       r1 = (  r(k  )+ r(k+1) )/2.0d0

      yijk2=r1*abs( th2-th1 )

      return
      end subroutine leng_yijk2



      subroutine leng_yijk2_shifted(j,k,th,r,yijk2)

!-------------------------------------------------------------
!      subroutine for calculating line element Y(i)jk2 - A.K.
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!-------------------------------------------------------------

      implicit none

      integer								:: j,k
      real(8),dimension(:)					:: th
      real(8),dimension(:)					:: r
      real(8)								:: yijk2,r1
 
      r1 = ( r(k )+ r(k+1) )/2.d0

      yijk2=r1*abs( th(j+1)-th(j) )
 
      return
      end subroutine leng_yijk2_shifted



      subroutine leng_zijk(k,r,zijk)

!-------------------------------------------------------------
!      subroutine for calculating line element Z(ij)k
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!-------------------------------------------------------------

      implicit none

      integer							:: k
      real(8),dimension(:)				:: r
      real(8)							:: zijk
 
      zijk=abs( r(k)-r(k+1) )
 
      return
      end subroutine leng_zijk


      subroutine leng_zijk2(k,r,zijk2)

!-------------------------------------------------------------
!      subroutine for calculating line element Zijk2
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!-------------------------------------------------------------

      implicit none

      integer                           :: k
      real(8),dimension(:)              :: r
      real(8)                           :: zijk2,r1,r2

      r1 = (  r(k  )+ r(k+1) )/2.0d0
      r2 = (  r(k+1)+ r(k+2) )/2.0d0

      zijk2=abs( r1-r2 )

      return
      end subroutine leng_zijk2


      subroutine leng_zijk2_shifted(k,r,zijk2)

!-------------------------------------------------------------
!      subroutine for calculating line element Z(ij)k2	- A.K.
!
!   -----------------------------------------------------
!    ------ phai and theta in radian
!    ------ phai denotes interval
!    ------ theta denotes angle from the north pole
!    ------ r is measured from the center of the earth
!-------------------------------------------------------------

      implicit none

      integer							:: k
      real(8),dimension(:)				:: r
      real(8)							:: zijk2
 
      zijk2=abs( r(k)-r(k+1) )
 
      return
      end subroutine leng_zijk2_shifted


! *****************************************************************************

      subroutine n_allhxijk(l,m,n,i,j,k,ii)
!-------------------------------------------------------------
!       subroutine for calculating order number of Hx-vector for
!             index (i,j,k)
!       including both variables and constants
!        i:(1,L) j:(2,M) k:(1,N+1)
!-------------------------------------------------------------
 
      implicit none

      integer                  :: l,m,n,i,j,k,ii

      ii=(m-1)*(n+1)*(i-1)+(m-1)*(k-1)+(j-1)

      return
      end subroutine n_allhxijk



      subroutine n_allhyijk(l,m,n,i,j,k,ii)
!-------------------------------------------------------------
!       subroutine for calculating order number of Hy-vector for
!             index (i,j,k)
!       including both variables and constants
!             i:(1,L) j:(1,M) k:(1,N+1)
!-------------------------------------------------------------

      implicit none

      integer                       :: l,m,n,i,j,k,ii

      ii=l*(n+1)*(j-1)+l*(k-1)+i
      
      return
      end subroutine n_allhyijk



      subroutine n_allhzijk(l,m,n,i,j,k,ii)
 
!-------------------------------------------------------------
!       subroutine for calculating order number of Hz-vector for
!             index (i,j,k)
!       including both variables and constants
!             i:(1,L) j:(1,M+1) k:(1,N)
!                      <----- no i dependency for j=1 or M+1
!-------------------------------------------------------------

      implicit none

      integer                         :: l,m,n,i,j,k,ii

      if (j > 1.and.j < m+1) then
         ii=(l*(m-1)+2)*(k-1)+1+l*(j-2)+i
      else if (j == 1) then
         ii=(l*(m-1)+2)*(k-1)+1
      else if (j == m+1) then
         ii=(l*(m-1)+2)*(k-1)+1+l*(j-2)+1
      end if

      return
      end subroutine n_allhzijk


! *****************************************************************************


	  subroutine n_allhijk(l,m,n,i,j,k,xyz,ii)
!-------------------------------------------------------------
!       subroutine for calculating order number of h-vector for
!             index (i,j,k) and xyz=1,2,3 (x,y,z respectively)
!       A.K.
!       Hx: i:(1,L) j:(2,M) k:(1,N+1)
!		Hy: i:(1,L) j:(1,M) k:(1,N+1)
!		Hz: i:(1,L) j:(1,M+1) k:(1,N)
!       Hz:   <----- no i dependency for j=1 or M+1
!-------------------------------------------------------------
 
      implicit none

      integer					:: l,m,n,i,j,k,xyz,ii
	  integer					:: allx, ally

	  allx = (m-1)*(n-1)*(l-1)+(m-1)*(n-2)+(m-1)
	  ally = l*(n-1)*(m-1)+l*(n-2)+l

	  if (xyz == 1) then
		ii=(m-1)*(n-1)*(i-1)+(m-1)*(k-2)+(j-1)
	  else if (xyz == 2) then
		ii=allx+l*(n-1)*(j-1)+l*(k-2)+i
	  else if (xyz == 3) then
		if (j > 1.and.j < m+1) then
		   ii=allx+ally+(l*(m-1)+2)*(k-2)+1+l*(j-2)+i
		else if (j == 1) then
		   ii=allx+ally+(l*(m-1)+2)*(k-2)+1
		else if (j == m+1) then
		   ii=allx+ally+(l*(m-1)+2)*(k-2)+1+l*(j-2)+1
		end if
	  else
		write(0,*) 'Error: (n_allhijk) wrong usage: xyz =',xyz
		stop
	  end if
	  
      return
      end subroutine n_allhijk


! *****************************************************************************

      subroutine n_allpotijk(l,m,n,i,j,k,ii)

!-------------------------------------------------------------
!      to calculate vector number for pot1(i,j,k)
!-------------------------------------------------------------

      integer              :: l,m,n,i,j,k,ii

      if (j > 1.and.j < m+1) then
         ii=(k-1)*( (m-1)*l+2 )+1+(j-2)*l+i
      else if (j == 1) then
         ii=(k-1)*( (m-1)*l+2 )+1
      else if (j == m+1) then
         ii=k*( (m-1)*l+2 )
      end if

      return
      end subroutine n_allpotijk



      subroutine n_dvpotijkBC3(l,m,n,i,j,k,ii)

!-------------------------------------------------------------
!      to calculate vector number for div(i,j,k) and pot0(i,j,k)
!toh:29/FEB/2000
!-------------------------------------------------------------

      integer                     :: l,m,n,i,j,k,ii

      if (j > 1 .and. j < m+1) then
         ii=(k-2)*( (m-1)*l+2 )+1+(j-2)*l+i
      else if (j == 1) then
         ii=(k-2)*( (m-1)*l+2 )+1
      else if (j == m+1) then
         ii=(k-1)*( (m-1)*l+2 )
      end if

      return
      end subroutine n_dvpotijkBC3

! *****************************************************************************

      subroutine n_hxijk(l,m,n,i,j,k,ii)

!------------------------------------------------------------- 
!       subroutine for calculating order number of Hx-vector for
!             index (i,j,k)
!-------------------------------------------------------------

      implicit none

      integer                    :: l,m,n,i,j,k,ii

      ii=(m-1)*(n-1)*(i-1)+(m-1)*(k-2)+(j-1)

      return
      end subroutine n_hxijk



      subroutine n_hyijk(nhx,l,m,n,i,j,k,ii)

!-------------------------------------------------------------
!       subroutine for calculating order number of Hy-vector for
!             index (i,j,k)
!-------------------------------------------------------------

      implicit none

      integer nhx,l,m,n,i,j,k,ii

      ii=nhx+l*(n-1)*(j-1)+l*(k-2)+i

      return
      end subroutine n_hyijk



      subroutine n_hzijk(nhx,nhy,l,m,n,i,j,k,ii)
 
!-------------------------------------------------------------
!       subroutine for calculating order number of Hz-vector for
!             index (i,j,k)
!-------------------------------------------------------------

      implicit none

      integer                  :: nhx,nhy,l,m,n,i,j,k,ii

      ii=nhx+nhy
      if (j > 1.and.j < m+1) then
         ii=ii+(l*(m-1)+2)*(k-2)+1+l*(j-2)+i
      else if (j == 1) then
         ii=ii+(l*(m-1)+2)*(k-2)+1
      else if (j == m+1) then
         ii=ii+(l*(m-1)+2)*(k-1)
      end if

      return
      end subroutine n_hzijk


  ! ***************************************************************************
  ! * following are basic subroutines previously found in program earth. None
  ! * but volume_sph are currently used by the forward modelling subroutines

      subroutine volume_sph(k1,k2,z,dim,vols)
!-------------------------------------------------------------
!      to calculate volume of a portion of sphere
!         between z=z(k1) and z(k2)
!-------------------------------------------------------------

      implicit none
      integer                         :: k1,k2,dim
      real(8),dimension(dim)          :: z
      real(8)                         :: vols,pi
  
      pi = 4.0D0*datan(1.0D0)
      vols=4.0d0/3.0d0*pi*abs( (z(k1)**3)-(z(k2)**3) )
 
      return
      end subroutine volume_sph


! *****************************************************************************

      subroutine cart2sphere(x,y,z,rlat,rlon,radius)

!----------------------------------------------------------------
! transforms cartesian coordinates (x,y,z) 
! into spherical coordinates on the unit sphere
! last mod: 25 July 1995
! by A. Schultz, Cambridge University
!----------------------------------------------------------------
	  implicit none

      real(8)                      :: x,y,z,degrad,pi
      real(8)                      :: rlat,rlon,radius

      pi = 4.0D0*datan(1.0D0)
      degrad = pi/180.0D0

      radius = dsqrt(x**2 + y**2 + z**2)
      rlat   = pi/2.0d0 - dacos(z/radius)

      if (pi/2.0d0-dabs(rlat) < 1.D-8) then                          ! Pole
            rlon = 0.0D0
      else if (y >= 0.0D0) then                                   ! Quads I and II
            rlon = dacos(x/dsqrt(x**2 + y**2))
      else                                                        ! Quads III and IV
            rlon = 2.0d0*pi-dacos(x/dsqrt(x**2 + y**2))
      end if

      rlat = rlat/degrad
      rlon = rlon/degrad

      return
      end subroutine cart2sphere


! *****************************************************************************

      real(8) function arclen (p,q)

      implicit real(8) (a-h,o-z)

      real(8),dimension(3)    :: p, q
      real(8)                 :: d
      integer                 :: i
!
!------------------------------------------------------------
!
!                                               adam schultz
!                                               Univ Cantab
!
! arc-length (radians) between 2 points on unit sphere.
!
! input - p,q - vectors of length 3 containing
!               the x-, y-, and z-coordinates (in
!               that order) of points on the unit
!               sphere.
!
! output - arclen - angle in radians between the
!               unit vectors p and q.  0 .le.
!               arclen .le. pi.
!
!-----------------------------------------------------------
!
      d = 0.0d0
      do i = 1,3
        d = d + (p(i) + q(i))**2
      end do

      if (d == 0.0d0) then
         arclen = 4.0d0*datan(1.0d0)
      else if (d >= 4.0d0) then
         arclen = 0.0d0
      else
         arclen = 2.0d0*datan(dsqrt((4.0d0-d)/d))
      end if

      return
      end function arclen



end module elements
