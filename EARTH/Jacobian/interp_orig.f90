! *****************************************************************************
module interp_orig
	! Basic module containing utility subroutines from the original earth code

  implicit none


Contains

  ! ***************************************************************************
  ! * following are basic spline interpolation and other utility subroutines 
  !	* previously found in program earth. They are not currently used by the
  !	* forward modelling modules, nor have they been used in the original code

      subroutine polint(xa,ya,n,x,y,dy)

      integer                  :: n,m,i,ns
      real(8), dimension(n)    :: xa,ya
      real(8), dimension(10)   :: c,d 
      real(8)                  :: x,y,dy,dif,dift,den,ho,hp,w

      ns=1
      dif=dabs(x-xa(1))

      do i=1,n 
        dift=dabs(x-xa(i))
        if (dift < dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
      end do

      y=ya(ns)
      ns=ns-1

      do m=1,n-1

        do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den == 0.0d0)stop 1
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        end do

        if (2*ns < n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy

      end do

      return
      end subroutine polint

! *****************************************************************************

      subroutine ratint(xa,ya,n,x,y,dy)
      
	  implicit none
      integer                     :: m,n,i,ns
      real(8), dimension(n)       :: xa,ya
      real(8), dimension(10)      :: c,d
      real(8)                     :: w,h,hh,x,y,dy,t,dd,tiny

      parameter (tiny=1.0d-25)

      ns=1
      hh=dabs(x-xa(1))
      do 11 i=1,n
        h=dabs(x-xa(i))
        if (h == 0.0d0)then
          y=ya(i)
          dy=0.0d0
          return
        else if (h < hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+tiny
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd == 0.0d0)stop 2
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
        if (2*ns < n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end subroutine ratint

! *****************************************************************************

      subroutine spline(x,y,n,yp1,ypn,y2)

	  implicit none

      integer                  :: n,i,k
      real(8), dimension(3)    :: x,y,y2
      real(8), dimension(100)  :: u
      real(8)                  :: yp1,ypn,sig,p,qn,un



      if (yp1 > 0.99d30) then
        y2(1) = 0.0d0
        u(1)  = 0.0d0
      else
        y2(1) = -0.50d0
        u(1)  = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p     = sig*y2(i-1)+2.0d0
        y2(i) = (sig-1.0d0)/p
        u(i)  = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/ &
                (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do

      if (ypn > 0.99d30) then
        qn = 0.0d0
        un = 0.0d0
      else
        qn = 0.5d0
        un = (3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)

      do k=n-1,1,-1
        y2(k) = y2(k)*y2(k+1)+u(k)
      end do

      return
      end subroutine spline

! *****************************************************************************

      subroutine splint(xa,ya,y2a,n,x,y)

	  implicit none

      integer                :: n,k,klo,khi
      real(8),dimension(3)   :: xa,ya,y2a
      real(8)                :: x,y,h,a,b

      klo=1
      khi=n

1     if (khi-klo > 1) then
        k=(khi+klo)/2
        if (xa(k) > x) then
          khi=k
        else
          klo=k
        endif
        goto 1
      endif

      h = xa(khi)-xa(klo)

      if (h == 0.0d0) stop 7

      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)+b*ya(khi)+ &
            ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0

      return
      end subroutine splint


end module interp_orig !	interp