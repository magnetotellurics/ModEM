!******************************************************************************
!f90 adaptation of numerical recipes routine spline
!
!given arrays x(1:n) and y(1:n) containing a tabulated function, 
!i.e. yi = f(xi), with x1<x2<...<xn,
!and given values yp1 and ypn for the first derivative of the interpolating
!function at points 1 and n, respectively, this routine returns an array y2(1:n)
!that contains the second derivatives of the interpolating function at the 
!tabulated points xi. If yp1 and/or ypn are equal to 1x10^30 or larger, the 
!routine is signaled to set the corresponding boundary condition for a natural 
!spline, with zero second derivative on that boundary.
!
!******************************************************************************

SUBROUTINE spline(x,y,n,yp1,ypn,y2)
 use, intrinsic :: iso_fortran_env
       INTEGER(kind=int32) :: n
       REAL(kind=real64) :: yp1,ypn,x(n),y(n),y2(n)  
       INTEGER(kind=int32) :: i,k  
       REAL(kind=real64) :: p,qn,sig,un
       real(kind=real64),dimension(:),allocatable :: u
       integer(kind=int32) :: ierr


       allocate(u(n), stat=ierr)
       if(ierr .NE. 0) call alloc_error(-1,'spline','u vector',ierr)

 
       if(yp1 .gt. .99d30) then  
         y2(1)=0.  
         u(1)=0.  
 
       else  
         y2(1)=-0.5  
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)  
 
       endif  
 
       do 11 i=2,n-1  
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))  
         p=sig*y2(i-1)+2.  
         y2(i)=(sig-1.)/p  
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p  
 
 11    continue  
 
       if(ypn.gt..99d30) then  
         qn=0.  
         un=0.  
 
       else  
         qn=0.5  
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))  
 
       endif  
 
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)  

       do 12 k=n-1,1,-1  
         y2(k)=y2(k)*y2(k+1)+u(k)  
 12    continue  
 

       deallocate(u, stat=ierr)

       return

endsubroutine spline



!******************************************************************************
!f90 adaptation of numerical recipes routine splint
!
!Given the arrays xa(1:n) and ya(1:n), which tabulate a function (with the 
!xai's in order), and given the array y2a(1:n), which is the output from 
!spline, and given a value of x, this routine returns a cubic-spline 
!interpolated value y.
!******************************************************************************
SUBROUTINE splint(xa,ya,y2a,n,x,y)  
    USE, INTRINSIC :: ISO_FORTRAN_ENV
       INTEGER(kind=int32) :: n  
       REAL(kind=real64) :: x,y,xa(n),y2a(n),ya(n)  
       INTEGER(kind=int32) :: k,khi,klo
       REAL(kind=real64) :: a,b,h
#ifdef USE_MPI
       integer(kind=int32) :: ierr
#endif
 
       klo=1  
       khi=n  
 
 1     if(khi-klo.gt.1) then  
         k=(khi+klo)/2  
         if(xa(k).gt.x)then  
           khi=k  
         else  
           klo=k  
         endif  
       goto 1  
       endif  
 
       h=xa(khi)-xa(klo)  
       if(h.EQ.0.) then
         write(*,'(a)') 'ERROR in spline interpolation: bad xa input in splint'
#ifdef USE_MPI
         call MPI_Abort(MPI_COMM_WORLD,-49,ierr)
#else
         stop 49
#endif
       endif

       a=(xa(khi)-x)/h  
       b=(x-xa(klo))/h  
       y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.  
 
       return

endsubroutine splint
 

!f90 adaptation of numerical recipes routine splie2
SUBROUTINE splie2(x2a,ya,m,n,y2a)  
    USE, INTRINSIC :: ISO_FORTRAN_ENV
       INTEGER(kind=int32) :: m,n
       REAL(kind=real64) :: x2a(n),y2a(m,n),ya(m,n)  
 
!>      USES spline  
 
       INTEGER(kind=int32) :: j,k
       REAL(kind=real64),dimension(:),allocatable :: y2tmp,ytmp
       integer(kind=int32) :: ierr
 

       allocate(y2tmp(n),ytmp(n), stat=ierr)
       if(ierr.NE.0) call alloc_error(-1,'splie2','ytmp, y2tmp',ierr)

       do 13 j=1,m  
         do 11 k=1,n  
           ytmp(k)=ya(j,k)  
 11      continue  
         call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)  
 
         do 12 k=1,n  
           y2a(j,k)=y2tmp(k)  
 12      continue  
 
 13    continue  
 
       deallocate(ytmp, y2tmp, stat=ierr)

endsubroutine splie2
 


!f90 adaptation of numerical recipes routine splin2
SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)  
    USE, INTRINSIC :: ISO_FORTRAN_ENV
       INTEGER(kind=int32) :: m,n
       REAL(kind=real64) :: x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)  
 
!>      USES spline,splint  
 
       INTEGER(kind=int32) :: j,k  
       REAL(kind=real64),dimension(:),allocatable :: y2tmp,ytmp,yytmp
       integer(kind=int32) :: ierr
 

       allocate(y2tmp(n),ytmp(n),yytmp(n), stat=ierr)
       if(ierr.NE.0) call alloc_error(-1,'splin2','ytmp, y2tmp, yytmp',ierr)


       do 12 j=1,m  
         do 11 k=1,n  
           ytmp(k)=ya(j,k)  
           y2tmp(k)=y2a(j,k)  
 11      continue  
 
         call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))  
 12    continue  
 
       call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)  
       call splint(x1a,yytmp,y2tmp,m,x1,y)  
 
       deallocate(ytmp,y2tmp,yytmp, stat=ierr)

endsubroutine splin2
 
