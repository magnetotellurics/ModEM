! *****************************************************************************
module ringcurrent
	! Module containing the subroutines required to set the boundary conditions
	! for the global 3D forward solver. The older version is set_initial_p10.
	! We currently use set_initial_ringcurrent (written by Ikuko Fujii) instead,
	! since this is modelling the reality with greater accuracy.

  use polpak, only: legendre_associated
  use math_constants
  use elements
  implicit none

  private						:: qsimp,trapzd,arth,br,btheta,cel


Contains

      subroutine set_initial_p10(l,m,n,x,y,z,hx,hy,hz,dim)

!-------------------------------------------------------------
!     to set boundary values of P10 at source
!            and zero at the bottom
! 27-05-99 modified by H. Toh for hyt & hzt to be tied up with mm
!-------------------------------------------------------------
 
      integer                   :: i,j,k,l,m,n,ii,dim

      complex(8),dimension(dim) :: hx,hy,hz
      real(8),dimension(l)	    :: x
      real(8),dimension(m+1)    :: y
      real(8),dimension(n+1)    :: z

!toh:
      real(8),dimension(m+1)     :: hyt
      real(8),dimension(m+1)     :: hzt
!      real(8),dimension(180)     :: hyt
!      real(8),dimension(181)     :: hzt
!toh:
      real(8)                   :: hyr,hzr,z1,z2,zz
!toh:4/MAR/2000
      integer                   :: im
      real(8)                   :: sr0, srp, stp, htp, stm, htm
      real(8)                   :: thm, thp, dph, sumsh
!     pi = 4.0d0*datan(1.0d0) ! pi is declared in module statement.
!toh:

      z1=z(1)
      z2=z(n+1)

!
!  First Hy...
!
      k=1
        do j=1,m
          hyt(j)=( dcos(y(j))-dcos(y(j+1)) )/( y(j+1)-y(j) )
          do i=1,l
            call n_allhyijk(l,m,n,i,j,k,ii)
            hy(ii)=dcmplx(hyt(j),0.0d0)
          end do
        end do
      do k=2,n+1
        zz=z(k)
        do j=1,m
          hyr=hyt(j)*(zz-z2)/(z1-z2)
          do i=1,l
            call n_allhyijk(l,m,n,i,j,k,ii)
            hy(ii)=dcmplx(hyr,0.0d0)
          end do
        end do
      end do

!
!  Hx
!
      do k=1,n+1
        do j=2,m
          do i=1,l
            call n_allhxijk(l,m,n,i,j,k,ii)
            hx(ii)=dcmplx(0.0d0,0.0d0)
          end do
        end do
      end do

      z1=( z(1)+z(2) )/2.0d0
      z2=( z(n)+z(n+1) )/2.0d0
!
!  Hz
!
      k=1
        j=1
!toh
!toh: factor of 2 was removed from hzt on 26/JUN/99.
!toh: Makoto mixed up Laplace and Helmholtz representation of the P10 field
!toh: and wrongly multiplied the source Hz by 2.
!toh
!toh:4/MAR/2000 Modified to force SDC at nodes on the top shell.
!toh:           Note, however, that the modified formulae for Hr at k=1
!toh:           are still valid only for zonal current soruces, i.e.,
!toh:           those without Hphi. If there's further necessity to include
!toh:           3D sources other than simple zonal ones, terms associated
!toh:           with Hphi should be added to RHS of the Hr formula.
!toh:
!          hzt(j)=dcos( y(j) )
          thp = ( y(j) + y(j+1) )/2.0d0
          sr0 = 2.0d0*pi*z(1)*z(1)*( 1.0d0 - dcos(thp) )
          srp = 2.0d0*pi*z1*z1*( 1.0d0 - dcos(thp) )
          sumsh = 0.0d0
          do i=1,l
             if ( i == 1 ) then
                im = l
             else
                im = i - 1
             end if
             dph = ( x(im) + x(i) )/2.0d0
             Stp = ( z(1)*z(1) - z1*z1 )*dsin(thp)*dph/2.0d0
             call n_allhyijk(l,m,n,i,1,1,ii)
             sumsh = sumsh + stp*dreal(hy(ii))
          end do
          hzt(j)=( Sr0*dcos( y(j) ) - Sumsh )/Srp
!toh
          call n_allhzijk(l,m,n,i,j,k,ii) ! ii is independent from i if j=1
          hz(ii)=dcmplx(hzt(j),0.0d0)
        do j=2,m
!toh:4/MAR/2000
!          hzt(j)=dcos( y(j) )
!toh
          do i=1,l
!toh:4/MAR/2000
             thp = ( y(j) + y(j+1) )/2.0d0
             thm = ( y(j-1) + y(j) )/2.0d0
             if ( i == 1 ) then
                im = l
             else
                im = i - 1
             end if
             dph = ( x(im) + x(i) )/2.0d0
             Sr0 = z(1)*z(1)*( dcos(thm) - dcos(thp) )*dph
             Stp = ( z(1)*z(1) - z1*z1 )*dsin(thp)*dph/2.0d0
             call n_allhyijk(l,m,n,i,j,1,ii)
             htp=dreal(hy(ii))
             Stm = ( z(1)*z(1) - z1*z1 )*dsin(thm)*dph/2.0d0
             call n_allhyijk(l,m,n,i,j-1,1,ii)
             htm=dreal(hy(ii))
             Srp = z1*z1*( dcos(thm) - dcos(thp) )*dph
             hzt(j)=( Sr0*dcos( y(j) ) - Stp*htp + Stm*htm )/Srp
!toh:
            call n_allhzijk(l,m,n,i,j,k,ii)
            hz(ii)=dcmplx(hzt(j),0.0d0)
          end do
        end do
        j=m+1
!toh:4/MAR/2000
!          hzt(j)=dcos( y(j) )
!toh:
          thp = ( y(j) + y(j-1) )/2.0d0
          sr0 = 2.0d0*pi*z(1)*z(1)*( 1.0d0 + dcos(thp) )
          srp = 2.0d0*pi*z1*z1*( 1.0d0 + dcos(thp) )
          sumsh = 0.0d0
          do i=1,l
             if ( i == 1 ) then
                im = l
             else
                im = i - 1
             end if
             dph = ( x(im) + x(i) )/2.0d0
             Stp = ( z(1)*z(1) - z1*z1 )*dsin(thp)*dph/2.0d0
             call n_allhyijk(l,m,n,i,m,1,ii)
             sumsh = sumsh + stp*dreal(hy(ii))
          end do
          hzt(j)=( Sr0*dcos( y(j) ) - Sumsh )/Srp
!toh
          call n_allhzijk(l,m,n,i,j,k,ii) ! ii is independent from i if j=m+1
          hz(ii)=dcmplx(hzt(j),0.0d0)
      do k=2,n
        zz=( z(k)+z(k+1) )/2.0d0
        do j=1,m+1
          hzr=hzt(j)*(zz-z2)/(z1-z2)
          do i=1,l
            call n_allhzijk(l,m,n,i,j,k,ii)
            hz(ii)=dcmplx(hzr,0.0d0)
          end do
        end do
      end do
!toh

      ! (AK) print *,'end of set-p01'
      return
      end subroutine set_initial_p10


! Start of Ikuko Fujii subroutines
      subroutine set_initial_ringcurrent(l,m,n,x,y,z,hx,hy,hz,dim)

!-------------------------------------------------------------
!     to set boundary values of a ring current at source
!            and zero at the bottom
!-------------------------------------------------------------
 
      integer                   :: i,j,k,l,m,n,ii,dim

      complex(8),dimension(dim) :: hx,hy,hz
      real(8),dimension(l)	    :: x
      real(8),dimension(m+1)    :: y
      real(8),dimension(n+1)    :: z

      real(8),dimension(m)		:: hyt
      real(8),dimension(m+1)    :: hzt
!      real(8),dimension(180)    :: hyt
!      real(8),dimension(181)    :: hzt
      real(8)                   :: hyr,hzr,z1,z2,zz

!!! fujii !!! Variables for ringcurrent systems
      integer,parameter        :: nring = 3    ! # of ringcurrents
      real(8),dimension(nring) :: colatring,height,crrnt 
      !Data colatring /90.D0,23.D0,157.D0/ ! colatitudes (degree) of ringcurrents
      !Data height /19013.D0,50.D0,50.D0/  ! heights (km) of ringcurrents from the boundary
      !Data crrnt /1.D9,1.D7,1.D7/         ! currents (A) of ringcurrents
      real(8)           :: re             ! radius of the boundary (km)
      real(8)           :: br1, btheta1	  ! Br and integrated B_theta
      real(8)           :: scalefactor    ! scaling factor
      integer           :: ir
!Interface
!  Function btheta(colat,colatring,height,re)
!    Implicit none
!    Real(8), Dimension(:), Intent(IN) :: colat
!    Real(8), Intent(IN) :: colatring, height, re
!    Real(8), Dimension(size(colat)) :: btheta
!  End Function btheta
!End Interface
!Interface
!  Function br(colat,colatring,height,re)
!    Implicit none
!    Real(8), Intent(IN) :: colat
!    Real(8), Intent(IN) :: colatring, height, re
!    Real(8) :: br
!  End Function br
!End Interface

!!! A.K. !!! These lines seem to be an artifact of the constant array dimensions.
!!! A.K. !!! If something goes wrong, check if there are any more constant
!!! A.K. !!! dimensions left in the code. Should no longer be the case.
!If ( m .gt. 180 ) then
!  Print *, 'm is too large, ',m
!  Stop
!End if
!If ( n .gt. 180 ) then
!  Print *, 'n is too large, ',n
!  Stop
!End if

! A.K. Corrections to Ikuko Fujii's code:
! Using the DATA statement to initialize colatring and other variables
! was causing problems. In particular, colatring is modified later in the code.
! Running this subroutine the second time in the program did not initialize
! colatring through a data statement, but rather used the vector computed the
! first time the subroutine was run. The incorrect initial value of colatring
! creates wrong values of btheta1 through the local variable btheta_intg, and
! eventually causes a program exception through stack overflow.
! The solution is to initialize this variable in the first line of the code
! rather than through the obsolete DATA statement.
! Also, we now initialize br1 and btheta1 with zeros.

	  btheta1 = 0.0d0
	  br1 = 0.0d0
      colatring = (/ 90.D0,23.D0,157.D0 /) ! colatitudes (degree) of ringcurrents
      height = (/ 19013.D0,50.D0,50.D0 /)  ! heights (km) of ringcurrents from the boundary
      crrnt = (/ 1.D9,1.D7,1.D7 /)         ! currents (A) of ringcurrents
	   
      re = abs(z(1))/1000.D0
      Print *, 'boundary radius (km) = ', re

      colatring = colatring * asin(1.D0) / 90.D0
      crrnt = crrnt * 1.D-7

      z1=z(1)
      z2=z(n+1)
 
!
!  First Hy...
!

      k=1
        do j=1,m
!!! fujii !!! original (P_1_0)  
!!!       hyt(j)=(dcos(y(j))-dcos(y(j+1)))/( y(j+1)-y(j) )
!!! fujii !!! P_n_0             
!!!       hyt(j) = -(plgndr(n,0,dcos(y(j+1)))-plgndr(n,0,dcos(y(j))))/(y(j+1)-y(j)) 
!!! fujii !!! ringcurrents at 'height's            
          hyt(j) = 0.D0
          Do ir = 1, nring
!!! fujii !!! Numerical evaluation of the theta-integral of btheta
			!write(0,'(a10,2i4,5g15.7)') 'j,a,b,..= ',j,ir,colatring(ir)
            Call qsimp(btheta,y(j),y(j+1),colatring(ir),height(ir),re,btheta1)
            hyt(j) = hyt(j) + crrnt(ir)*btheta1
          End do
          hyt(j) = hyt(j) / ( y(j+1)-y(j) )
          !print *,hyt(j)
        end do
        scalefactor = abs(maxval(hyt(1:m)))
        hyt(1:m) = hyt(1:m)/scalefactor

        do j=1,m
          do i=1,l
            call n_allhyijk(l,m,n,i,j,k,ii)
            hy(ii)=dcmplx(hyt(j),0.0d0)
          end do
        end do
      do k=2,n+1
        zz=z(k)
        do j=1,m
          hyr=hyt(j)*(zz-z2)/(z1-z2)
          do i=1,l
            call n_allhyijk(l,m,n,i,j,k,ii)
            hy(ii)=dcmplx(hyr,0.0d0)
          end do
        end do
      end do

!
!  Hx
!
      do k=1,n+1
        do j=2,m
          do i=1,l
            call n_allhxijk(l,m,n,i,j,k,ii)
            hx(ii)=dcmplx(0.0d0,0.0d0)
          end do
        end do
      end do

      z1=( z(1)+z(2) )/2.0d0
      z2=( z(n)+z(n+1) )/2.0d0
!
!  Hz
!
      k=1
        j=1
!!! fujii !!! original (P_1_0)  
!         hzt(j)=dcos( y(j) )
!!! fujii !!! P_n_0             
!         hzt(j) = dble(n) * plgndr(n,0,dcos(y(j)))
!!! fujii !!! ringcurrents at 'height'             
          hzt(j) = 0.D0
          Do ir = 1, nring
            br1 = br(y(j),colatring(ir),height(ir),re)
            hzt(j) = hzt(j) + crrnt(ir)*br1/scalefactor
          End do
          call n_allhzijk(l,m,n,i,j,k,ii)
          hz(ii)=dcmplx(hzt(j),0.0d0)
        do j=2,m
!!! fujii !!! original (P_1_0)         
!         hzt(j)=dcos( y(j) )
!!! fujii !!! P_n_0             
!         hzt(j) = dble(n) * plgndr(n,0,dcos(y(j)))
!!! fujii !!! a ringcurrent at 'height'             
          hzt(j) = 0.D0
          Do ir = 1, nring
            br1 = br(y(j),colatring(ir),height(ir),re)
            hzt(j) = hzt(j)+crrnt(ir)*br1/scalefactor
          End do
          do i=1,l
            call n_allhzijk(l,m,n,i,j,k,ii)
            hz(ii)=dcmplx(hzt(j),0.0d0)
          end do
        end do
        j=m+1
!!! fujii !!! original (P_1_0)          
!         hzt(j)=dcos( y(j) )
!!! fujii !!! P_n_0 
!         hzt(j) = dble(n) * plgndr(n,0,dcos(y(j)))
!!! fujii !!! a ringcurrent at 'height' 
          hzt(j) = 0.D0
          Do ir = 1, nring
            br1 = br(y(j),colatring(ir),height(ir),re)
            hzt(j) = hzt(j)+crrnt(ir)*br1/scalefactor
          End do
          call n_allhzijk(l,m,n,i,j,k,ii)
          hz(ii)=dcmplx(hzt(j),0.0d0)
      do k=2,n
        zz=( z(k)+z(k+1) )/2.0d0
        do j=1,m+1
          hzr=hzt(j)*(zz-z2)/(z1-z2)
          do i=1,l
            call n_allhzijk(l,m,n,i,j,k,ii)
            hz(ii)=dcmplx(hzr,0.0d0)
          end do
        end do
      end do

      print *,'end of set-ringcurrent'
      return
      end subroutine set_initial_ringcurrent



Subroutine qsimp(btheta,a,b,colatring,height,re,btheta_intg)

Implicit none

Real(8), Intent(IN) :: a,b,colatring,height,re
Real(8), Intent(OUT) :: btheta_intg
Interface
  Function btheta(colat,colatring,height,re)
    Implicit none
    Real(8), Dimension(:), Intent(IN) :: colat
    Real(8), Intent(IN) :: colatring, height, re
    Real(8), Dimension(size(colat)) :: btheta
  End Function btheta
End Interface
Integer, parameter :: jmax = 20 
Real(8), parameter :: eps=1.D-6, unlikely = -1.D30
Integer :: j
Real(8) :: os,ost,st

ost = unlikely
os  = unlikely
Do j = 1, jmax
  Call trapzd(btheta,a,b,colatring,height,re,st,j)
  !write(0,'(a10,i4,5g15.7)') 'j,a,b,..= ',j,a,b,colatring,height,re
  !write(0,'(a20,i4,3g15.7)') 'j,abs(),eps*abs(os)= ',j,btheta_intg,os,st
  btheta_intg = (4.D0*st - ost) / 3.D0
  If ( abs(btheta_intg-os) < eps*abs(os) ) Return
  If ( btheta_intg==0.D0 .and. os ==0.D0 .and. j>6 ) Return
  os = btheta_intg
  ost = st
End do

End Subroutine qsimp

Subroutine trapzd(btheta,a,b,colatring,height,re,s,n)

Implicit none

Real(8), Intent(IN) :: a,b,colatring,height,re
Real(8), Intent(INOUT) :: s
Integer, Intent(IN) :: n
Interface 
  Function btheta(colat,colatring,height,re)
    Implicit none
    Real(8), Dimension(:), Intent(IN) :: colat
    Real(8), Intent(IN) :: colatring, height, re
    Real(8), Dimension(size(colat)) :: btheta
  End Function btheta
End Interface 
!Interface
!  Pure Function arth(first,increment,n)
!    Implicit none
!    Real(8), intent(in) :: first, increment
!    Integer, intent(in) :: n
!    Real(8), dimension(n) :: arth
!  End Function arth
!End Interface
Real(8) :: del, fsum
Integer :: it

If ( n==1 ) then
  s = 0.5D0*(b-a)*sum(btheta((/a,b/),colatring,height,re))
Else
  it = 2**(n-2)
  del = (b-a)/it
  !print *, 'size(colat) = ',size(btheta(arth(a+0.5D0*del,del,it),colatring,height,re))
  ! NB : size(arth(a+0.5D0*del,del,it)) = it
  !print *, 'n,arth(a+...),del = ',n,size(arth(a+0.5D0*del,del,it)),del
  fsum = sum(btheta(arth(a+0.5D0*del,del,it),colatring,height,re))
  s = 0.5D0*(s+del*fsum)
End if

End Subroutine trapzd
 
! End of Ikuko Fujii's internal subroutines

Pure Function arth(first,increment,n)

  Implicit none

  Integer, Parameter :: NPAR_ARTH=16, NPAR2_ARTH=8

  Real(8), intent(in) :: first, increment
  Integer, intent(in) :: n
  Real(8), dimension(n) :: arth
  Integer :: k, k2
  Real(8) :: temp

  If ( n > 0 ) arth(1) = first
  If ( n <= NPAR_ARTH ) then
    Do k = 2, n
      arth(k) = arth(k-1) + increment 
    End do
  Else
    Do k = 2, NPAR2_ARTH
      arth(k) = arth(k-1) + increment 
    End do
    temp = increment * NPAR2_ARTH
    k = NPAR2_ARTH
    Do 
       If ( k >= n ) exit
       k2 = k + k
       arth(k+1:min(k2,n)) = temp + arth(1:min(k,n-k))
       temp = temp + temp
       k = k2
    End do
  End if

End Function arth

Function br(colat,colatring,height,re)

! ******************************************************************************
! *   Function to calculate the vertical component of the geomagnetic 
! *   field on the surface of a sphere (radius = re ) due to a ring 
! *   current at given height and colatitude.     *
! *                                                                            *
! *   Complete ellipse integration in Numerical Recipes is used.               *
! *   (subroutine "cel.f90")                                                   *
! *                                                                            *
! *   See the formulation in Fujii & Schultz (2001)                           *
! *                                                                            *
! * -------------------------------------------------------------------------- *
! *   By : Ikuko Fujii (Kakioka Magnetic Observatory)         *
! *   Last modified : August. 17, 2001                                            *
! * -------------------------------------------------------------------------- *
! *                                                                            *
! *   Input;                                                                   *
! *     colat = colatitude(radian) where the geomagnetic field is requested.   *
! *     colatring = colatitude(radian) at the centre of the ring current       *
! *     height = height of the ring current (km) at colatitude=colatring       *
! *     re = radius of a sphere on which obsevation is made (km)               *
! *   Variable;                                                                *
! *     a = radius of the ring current (km)... re+height                       *
! *   Output;                                                                  *
! *     br = r component of the geoamgnetic field                           *
! ******************************************************************************

Implicit none

! ******* Input 
Real(8), intent(IN) :: colat
Real(8), intent(IN) :: colatring,height,re
! ******* Variable 
Integer i
Real(8) a, c0, c, s, s0, k, k2, kc, factor, func1, func2 
! ******* Output 
Real(8) :: br
!Interface
!  Real(8) function cel(qqc,pp,aa,bb)
!    Implicit none
!    Real(8), intent(IN) :: qqc,pp,aa,bb
!    Real(8), parameter :: ca=.0003
!    Real(8) qc,a,b,p,q,e,em,f,g,pio2
!  End function cel
!End Interface
! ===========================================================
       
! -------- Set constant parameters
      a = re + height
      c0 = cos(colatring)
      s0 = sin(colatring)
    
!      Print *, colatring
!      Print *, colat

! -------- Set parameters
  
     If ( colat <= 0.D0+EPS_GRID ) then
        c = cos(colat+0.0017453292)
        s = sin(colat+0.0017453292)
     ElseIf ( colat >= 2.D0*asin(1.D0)-EPS_GRID ) then
        c = cos(colat-0.0017453292)
        s = sin(colat-0.0017453292)
     Else
        c = cos(colat)
        s = sin(colat)
     End if
        k2 = 4.D0*re*a*s*s0/(a*a+re*re-2.D0*a*re*cos(colat+colatring))
        k = sqrt(k2)
        kc = sqrt(1.D0-k2)
        factor=1.D0/re/s/sqrt(a*a+re*re-2.D0*a*re*cos(colat+colatring))
        func1 = (2.D0-k2)*cel(kc,1.D0,1.D0,1.D0) &
                       - 2.D0*cel(kc,1.D0,1.D0,kc*kc)
        func2 = -k2*cel(kc,1.D0,1.D0,1.D0) &
                  - 2.D0*cel(kc,1.D0,1.D0,kc*kc) &
                    + (2.D0-k2)*cel(kc,-k2+1.D0,1.D0,1.D0)

! -------- Compute Bx and Br at colat
        br = factor/re*( re*a*sin(colat+colatring)*func1 + &
            0.5D0/s*((a*a+re*re)*c-2.D0*re*a*c0)*func2 )

End Function br

Function btheta(colat,colatring,height,re)

! ******************************************************************************
! *   Function to calculate the theta component of the geomagnetic field 
! *   on the surface of a sphere (radius = re ) due to a ring current 
! *   at given height and colatitude.     *
! *                                                                            *
! *   Complete ellipse integration in Numerical Recipes is used.               *
! *   (subroutine "cel.f90")                                                   *
! *                                                                            *
! *   See the formulation in Fujii & Schultz (2001)                           *
! *                                                                            *
! * -------------------------------------------------------------------------- *
! *   By : Ikuko Fujii (Kakioka Magnetic Observatory)         *
! *   Last modified : Aug. 17, 2001                                            *
! * -------------------------------------------------------------------------- *
! *                                                                            *
! *   Input;                                                                   *
! *     colat = colatitude(radian) where the geomagnetic field is requested.   *
! *     colatring = colatitude(radian) at the centre of the ring current       *
! *     height = height of the ring current (km) at colatitude=colatring       *
! *     re = radius of a sphere on which obsevation is made (km)               *
! *   Variable;                                                                *
! *     a = radius of the ring current (km)... re+height      
! *   Output;                                                                  *
! *     btheta(:) = theta component of the geoamgnetic field                   *
! ******************************************************************************

Implicit none

! ******* Input 
Real(8), Dimension(:), Intent(IN) :: colat
Real(8), intent(IN) :: colatring,height,re
! ******* Variable 
Integer i,n
Real(8) a, c0, c, s, s0, k, k2, kc, factor, func1, func2 
! ******* Output 
Real(8), Dimension(size(colat)) :: btheta
! ******* Interface
!Interface
!  Real(8) function cel(qqc,pp,aa,bb)
!    Implicit none
!    Real(8), intent(IN) :: qqc,pp,aa,bb
!    Real(8), parameter :: ca=.0003
!    Real(8) qc,a,b,p,q,e,em,f,g,pio2
!  End function cel
!End Interface
! ===========================================================
       
! -------- Set constant parameters
      a = re + height
      c0 = cos(colatring)
      s0 = sin(colatring)
      n = size(colat)
    
  Do i = 1, n
! -------- Set parameters
  
     If ( colat(i) .eq. 0.0D0 ) then
        c = cos(colat(i)+0.0017453292)
        s = sin(colat(i)+0.0017453292)
     ElseIf ( colat(i) .eq. 2.0D0*asin(1.0D0) ) then
        c = cos(colat(i)-0.0017453292)
        s = sin(colat(i)-0.0017453292)
     Else
        c = cos(colat(i))
        s = sin(colat(i))
     End if
        k2 = 4.0D0*re*a*s*s0/(a*a+re*re-2.0D0*a*re*cos(colat(i)+colatring))
        k = sqrt(k2)
        kc = sqrt(1.D0-k2)
        factor=1.0D0/re/s/sqrt(a*a+re*re-2.0D0*a*re*cos(colat(i)+colatring))
        func1 = (2.0D0-k2)*cel(kc,1.0D0,1.0D0,1.0D0) &
                       - 2.D0*cel(kc,1.0D0,1.0D0,kc*kc)
        func2 = -k2*cel(kc,1.0D0,1.0D0,1.0D0) &
                  - 2.0D0*cel(kc,1.0D0,1.0D0,kc*kc) &
                    + (2.0D0-k2)*cel(kc,-k2+1.0D0,1.0D0,1.0D0)

! -------- Compute Btheta at colat
        btheta(i) = -factor*( (re-a*cos(colat(i)+colatring))*func1 + &
            0.5D0*(a*a-re*re)/re*func2 )

  End do

End Function btheta

Real(8) function cel(qqc,pp,aa,bb)

! Compute complete elliptic integrals 
! by Numerical Recipes (first edition)
 
Implicit none

Real(8), intent(IN) :: qqc,pp,aa,bb
Real(8), parameter :: ca=.0003
Real(8) qc,a,b,p,q,e,em,f,g,pio2

pio2=asin(1.D0)

If(qqc.eq.0.) then
  Print *, 'failure in cel'
  Stop
End if

      qc=abs(qqc)
      a=aa
      b=bb
      p=pp
      e=qc
      em=1.
      if(p.gt.0.)then
        p=sqrt(p)
        b=b/p
      else
        f=qc*qc
        q=1.-f
        g=1.-p
        f=f-p
        q=q*(b-a*p)
        p=sqrt(f/g)
        a=(a-b)/g
        b=-q/(g*g*p)+a*p
      endif
      Do
        f=a
        a=a+b/p
        g=e/p
        b=b+f*g
        b=b+b
        p=g+p
        g=em
        em=qc+em
        if(abs(g-qc).le.g*ca)exit
        qc=sqrt(e)
        qc=qc+qc
        e=qc*em
      End do
 
      cel=pio2*(b+a*em)/(em*(em+p))
      
end function cel


! End of Ikuko Fujii's external subroutines

! ***************************************************************************
! * computes the value for the (l,m) legendre polynomial; terribly inefficient
! * but this isn't used anywhere outside of the ringcurrent module so not
! * a concern. For efficient calculations, use paramfunc module (SphHarm).
function plgndr(l,m,x)

    integer, intent(in)     :: l,m
    real(8), intent(in)     :: x
    real(8)                 :: plgndr
    real(8), dimension(0:l) :: P_lm

    call legendre_associated(l,m,x,P_lm)

    plgndr = P_lm(l)

end function plgndr


end module ringcurrent
