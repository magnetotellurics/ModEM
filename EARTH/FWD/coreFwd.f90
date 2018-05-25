! *****************************************************************************
module coreFwd
  ! Main module in the project earth3d, containing all the necessary routines
  ! to solve numerically the integral form of the quasi-static approximation
  ! to Maxwell's equations on the given grid

  use math_constants
  use elements
  use dimensions
  use iospec
  implicit none


Contains


  ! ***************************************************************************
  ! * subroutines required to run the forward solver

      subroutine asolvec(b,x,sb,ija,np)
!-------------------------------------------------------------
!    Solve Ax=b, A preconditioner, determined by incomplete
!    Cholesky decomposition, stored in sb and ija arrays.
!-------------------------------------------------------------

      complex(8),dimension(np2)      :: b,x
!toh
      complex(8)                     :: sum, xi
      complex(8),dimension(np4)      :: sb
!toh
      integer,dimension(np4)         :: ija
      integer                        :: i,m,k,np

      x(1:np) = dcmplx(0.0d0,0.0d0)

!------forward substitution, solve L(t)*y = b, storing y in x
!toh
      x(1)=b(1)/sb(1)
!toh
      do i=1, np-1
!toh
          xi = x(i)
!toh
          do m=ija(i), ija(i+1)-1

             x(ija(m))=x(ija(m))+sb(m)*xi

          end do
!toh
          x(i+1)=(b(i+1)-x(i+1))/sb(i+1)
!toh
      end do

! -----backward substitution, solve L*x = y -----------------
!toh
      x(np)=x(np)/sb(np)
!toh
      do i=np-1, 1, -1

          sum=dcmplx(0.0d0,0.0d0)

          do k=ija(i), ija(i+1)-1
             sum=sum+sb(k)*x(ija(k))
          end do
!toh
          x(i)=(x(i)-sum)/sb(i)
!toh
      end do

      return
      end subroutine asolvec


! *****************************************************************************



      subroutine asolver(n,b,x,sb,ija)

!-------------------------------------------------------------
!    Solve Ax=b, A preconditioner, determined by incomplete
!    Cholesky decomposition, stored in sb and ija arrays.
!-------------------------------------------------------------

      real(8),dimension(np1)       :: b,x
      real(8)                      :: sum
      real(8),dimension(np6)       :: sb
      integer,dimension(np6)       :: ija
      integer                      :: i,k,m,n

      x(1:n)=0.0d0

!------forward substitution, solve L(t)*y = b, storing y in x

      x(1)=b(1)/sb(1)
      do i=1, n-1
          do m=ija(i), ija(i+1)-1
             x(ija(m))=x(ija(m))+sb(m)*x(i)
          end do
          x(i+1)=(b(i+1)-x(i+1))/sb(i+1)
      end do

! -----backward substitution, solve L*x = y -----------------

      x(n)=x(n)/sb(n)
      do i=n-1, 1, -1
          sum=0.0d0
          do k=ija(i), ija(i+1)-1
             sum=sum+sb(k)*x(ija(k))
          end do
          x(i)=(x(i)-sum)/sb(i)
      end do

      return
      end subroutine asolver


! *****************************************************************************



      subroutine atimes_s(x,r,sa,ija,dp,np)

!-------------------------------------------------------------
!
!   Multiply a matrix in row-index sparse storage arrays sa and ija by a
!   vector x(1:n), giving a vector r(1:n).   However, note, sa only
!   contains upper triangular part, because A is symmetric.
!   By: H. Toh & A. Schultz
!   last mod: 18 Feb 99 by A.S.
!-------------------------------------------------------------


      integer                           :: np

      complex(8), dimension(np2)        :: r
      real(8),    dimension(np+1:np5)   :: sa               ! sparse matrix offdiagonal elements
      complex(8), dimension(np)         :: dp               ! sparse matrix diagonal elements
      complex(8), dimension(np2)        :: x
      integer,    dimension(np5)        :: ija
      integer                           :: i,k,istart,iend


      if(ija(1).ne.np+2) then
        write(0,*) 'Error: (atimes_s) mismatched vector and matrix'
      end if

      r(1:np) = dcmplx(0.0d0,0.0d0)
      do i=1,np
         istart = ija(i)
         iend   = ija(i+1)-1
         r(i) = r(i)+dp(i)*x(i)
         do k=istart,iend
            r(i) = r(i)+sa(k)*x(ija(k))
            r(ija(k)) = r(ija(k))+sa(k)*x(i)
         end do
      end do

!      r(1:np) = dp(1:np)*x(1:np)

!      do i=1,np
!         istart = ija(i)
!         iend   = ija(i+1)-1
! 		 r(i) = r(i) + dot_product(sa(istart:iend),x(ija(istart:iend)))
! 		 r(ija(istart:iend)) = r(ija(istart:iend)) + sa(istart:iend)*x(i)
!      end do

      return
      end subroutine atimes_s

! *****************************************************************************

      subroutine atimes_mod(x,r,sa,ija,dp,np,nmax)

!-------------------------------------------------------------
!
!   Multiply a matrix in row-index sparse storage arrays sa and ija by a
!   vector x(1:n), giving a vector r(1:n).   However, note, sa only
!   contains upper triangular part, because A is symmetric.
!   By: H. Toh & A. Schultz
!   last mod by A. Schultz: 18 Feb 1999
!   Storing the the local ija indices in temporary arrays for efficiency;
!   however, this doesn't seem to give much (if any) improvement.
!   last mod by A. Kelbert: 21 Oct 2010
!-------------------------------------------------------------


      integer                           :: np,nmax

      complex(8), dimension(np2)        :: r
      real(8),    dimension(np+1:np5)   :: sa               ! sparse matrix offdiagonal elements
      real(8),    dimension(nmax)       :: sa_row
      complex(8), dimension(np)         :: dp               ! sparse matrix diagonal elements
      complex(8), dimension(np2)        :: x
      integer,    dimension(np5)        :: ija
      integer,    dimension(nmax)       :: ija_row
      integer                           :: i,k,istart,iend,n


      if(ija(1).ne.np+2) then
        write(0,*) 'Error: (atimes) mismatched vector and matrix'
      end if

      ! Initialize temporary variables for efficiency
      ija_row(1:nmax) = 0
      sa_row(1:nmax) = 0.0d0

      r(1:np) = dp(1:np)*x(1:np)

      do i=1,np

         istart = ija(i)
         iend   = ija(i+1)-1

         n = ija(i+1) - ija(i)
         ija_row(1:n) = ija(istart:iend)
         sa_row(1:n)  = sa(istart:iend)

         r(i) = r(i) + sum(sa_row(1:n) * x(ija_row(1:n)))
 	     r(ija_row(1:n)) = r(ija_row(1:n)) + sa_row(1:n)*x(i)

      end do

      return
      end subroutine atimes_mod


 ! *****************************************************************************

      subroutine atimes(x,r,sa,ija,dp,np,nmax)

!-------------------------------------------------------------
!
!   Multiply a matrix in row-index sparse storage arrays sa and ija by a
!   vector x(1:n), giving a vector r(1:n).   However, note, sa only
!   contains upper triangular part, because A is symmetric.
!   By: H. Toh & A. Schultz
!   last mod: 18 Feb 99 by A.S.
!-------------------------------------------------------------

      integer                           :: np
      integer                           :: nmax             ! for temporary arrays; not used

      complex(8), dimension(np2)        :: r
      real(8),    dimension(np+1:np5)   :: sa               ! sparse matrix offdiagonal elements
      complex(8), dimension(np)         :: dp               ! sparse matrix diagonal elements
      complex(8), dimension(np2)        :: x
      integer,    dimension(np5)        :: ija
      integer                           :: i,k,istart,iend


      if(ija(1).ne.np+2) then
        print*,'mismatched vector and matrix'
      end if

      r(1:np) = dp(1:np)*x(1:np)

      do i=1,np
         istart = ija(i)
         iend   = ija(i+1)-1
         r(i) = r(i) + sum(sa(istart:iend) * x(ija(istart:iend)))
         r(ija(istart:iend)) = r(ija(istart:iend)) + sa(istart:iend)*x(i)
      end do


      return
      end subroutine atimes

! *****************************************************************************


      subroutine atimesr(n,x,r,sa,ija)

!-------------------------------------------------------------
!   Multiply a matrix in row-index sparse storage arrays sa and ija by a
!   vector x(1:n), giving a vector r(1:n).   However, note, sa only
!   contains upper triangular part, because A is symmetric.
!-------------------------------------------------------------

      real(8),dimension(np1)  :: r,x
      real(8),dimension(np6)  :: sa
      integer,dimension(np6)  :: ija
      integer                 :: i,n,k

      if(ija(1).ne.n+2) then
        write(0,*) 'Error: (atimesr) mismatched vector and matrix'
      end if

      r(1:n)=0.0d0

      do i=1,n
         r(i)=r(i)+sa(i)*x(i)
         do k=ija(i),ija(i+1)-1
            r(i)=r(i)+sa(k)*x(ija(k))
            r(ija(k))=r(ija(k))+sa(k)*x(i)
         end do
      end do

      return
      end subroutine atimesr

! *****************************************************************************



      subroutine calc_vdivBC3(l,m,n,hx,hy,hz,divr,divi,x,y,z)

!-------------------------------------------------------------
! subroutine to calculate divH at internal nodes
! assume Hz=0 at k=n and above air layers
!
! modified for boundary values at the top of the air layers
! to be fixed at their starting values.
!
! modified to use real space values
!
! A.K.: note that the only reason the boundary conditions are
! required in this subroutine is that they are currently not
! only defined by H%x,H%y at k=1,n+1, but also include H%z at k=1.
! This is strictly incorrect, since these values should be
! computed along with the other interior field components.
! If H%z at k=1 were included in the interior domain, the
! divergence correction would not need to use the boundary
! conditions at all. Divergence is only defined in the interior.
! Date of comment: Nov 5, 2005
!-------------------------------------------------------------

	  integer							:: l,m,n
      complex(8),dimension(np1)         :: hx,hy,hz
      complex(8)                        :: sumdp
      real(8),dimension(np1)            :: divr,divi
      real(8)                           :: dx,dy,dz
      real(8),dimension(l)	            :: x
      real(8),dimension(m+1)            :: y
      real(8),dimension(n+1)            :: z
      real(8)                           :: zm,zp,yp,ym
      real(8)                           :: sijk2,sjki2,skij2
      integer                           :: i,ii,j,k,a,d

!     pi=4.d0*datan(1.d0)
!toh:29/FEB/2000
      do k=2,n ! do not use H%x,H%y at the boundary layers k=1, k=n+1
        zm=( z(k)+z(k-1) )/2.d0
        zp=( z(k)+z(k+1) )/2.d0
        j=1
          yp=( y(j)+y(j+1) )/2.d0
          skij2=2.d0*pi*(zm**2)*(1.d0-dcos(yp))
          call n_allhzijk(l,m,n,1,j,k-1,ii) ! using H%z(1,1,1)
          sumdp =      -hz(ii)*skij2
          skij2=2.d0*pi*(zp**2)*(1.d0-dcos(yp))
          call n_allhzijk(l,m,n,1,j,k,ii)
          sumdp = sumdp+hz(ii)*skij2
          do i=1,l
            if (i.ne.1) then
               a=i-1
            else
               a=l
            end if
            if (i.ne.l) then
               d=i+1
            else
               d=1
            end if
            call area_sjki2(a,i,j,k-1,x,y,z,sjki2)
            call n_allhyijk(l,m,n,i,j,k,ii) ! does not use b.c.
            sumdp = sumdp+hy(ii)*sjki2
          end do
          call n_dvpotijkBC3(l,m,n,1,j,k,ii)
          divr(ii)=dble(sumdp)
          divi(ii)=dimag(sumdp)
!       end j=1
        do j=2,m
          do i=1,l
            if (i.ne.1) then
               a=i-1
            else
               a=l
            end if
            if (i.ne.l) then
               d=i+1
            else
               d=1
            end if
            call area_sijk2(j-1,k-1,y,z,sijk2)
            call n_allhxijk(l,m,n,a,j,k,ii)
            sumdp =      -hx(ii)*sijk2
            call n_allhxijk(l,m,n,i,j,k,ii)
            sumdp = sumdp+hx(ii)*sijk2
            call area_sjki2(a,i,j-1,k-1,x,y,z,sjki2)
            call n_allhyijk(l,m,n,i,j-1,k,ii)
            sumdp = sumdp-hy(ii)*sjki2
            call area_sjki2(a,i,j,k-1,x,y,z,sjki2)
            call n_allhyijk(l,m,n,i,j,k,ii)
            sumdp = sumdp+hy(ii)*sjki2
            call area_skij2(a,i,j-1,k-1,x,y,z,skij2)
            call n_allhzijk(l,m,n,i,j,k-1,ii) ! using H%z(i,j,1)
            sumdp = sumdp-hz(ii)*skij2
            call area_skij2(a,i,j-1,k,x,y,z,skij2)
            call n_allhzijk(l,m,n,i,j,k,ii)
            sumdp = sumdp+hz(ii)*skij2
            call n_dvpotijkBC3(l,m,n,i,j,k,ii)
            divr(ii)=dble(sumdp)
            divi(ii)=dimag(sumdp)
          end do
        end do
        j=m+1
          ym=( y(j)+y(j-1) )/2.d0
          skij2=2.d0*pi*(zm**2)*(dcos(ym)+1.d0)
          call n_allhzijk(l,m,n,1,j,k-1,ii) ! using H%z(1,m+1,1)
          sumdp =      -hz(ii)*skij2
          skij2=2.d0*pi*(zp**2)*(dcos(ym)+1.d0)
          call n_allhzijk(l,m,n,1,j,k,ii)
          sumdp = sumdp+hz(ii)*skij2
          do i=1,l
            if (i.ne.1) then
               a=i-1
            else
               a=l
            end if
            if (i.ne.l) then
               d=i+1
            else
               d=1
            end if
            call area_sjki2(a,i,j-1,k-1,x,y,z,sjki2)
            call n_allhyijk(l,m,n,i,j-1,k,ii)
            sumdp = sumdp-hy(ii)*sjki2
          end do
          call n_dvpotijkBC3(l,m,n,1,j,k,ii)
          divr(ii)=dble(sumdp)
          divi(ii)=dimag(sumdp)
!       end j=M+1
      end do

      return
      end subroutine calc_vdivBC3

! *****************************************************************************



       subroutine calcb(l,m,n,hx,hy,hz,bvec,resist,x,y,z)

!-------------------------------------------------------------
! subroutine to compute right hand side of equation for given
! input Hx,Hy,Hz
!
! modified for spherical problem (26/11/97)
!-------------------------------------------------------------


      implicit none

      integer                     :: l,m,n,i,j,k,ii,jj
      integer                     :: a,b,c,d,e,f

      real(8),dimension(l,m,n)	  :: resist

      real(8),dimension(l)	      :: x
      real(8),dimension(m+1)      :: y
      real(8),dimension(n+1)      :: z
      complex(8),dimension(np1)   :: hx,hy,hz
      complex(8),dimension(np2)   :: bvec
      complex(8)                  :: sum

      real(8)                     :: rijcy,rijky,rijkz,ribkz,rijkx
      real(8)                     :: rijcx,rajkz,rajky,ribkx
      real(8)                     :: xx,yy,zz
      real(8)                     :: ym,yp,zm,zp
      real(8)                     :: zc,zk
      real(8)                     :: xipm,xapm,xipp,xapp,ybm,yjm
      real(8)                     :: ybp,yjp,ximp,xamp
      real(8)                     :: skib,skij,sjci,sjki
      real(8)                     :: skaj,sijc,sijk
      real(8)                     :: sibk,sjka

!
! First, Hx...
!
      jj=0
!     -----------
      do 10 i=1,l
!     -----------
        if (i == 1) then
           a=l
        else
           a=i-1
        end if
        if (i == l) then
           d=1
        else
           d=i+1
        end if
!     -------------
        do 10 k=2,n
!     -------------
          c=k-1
          f=k+1
          zm=(z(k)+z(k-1))/2.d0
          zp=(z(k)+z(k+1))/2.d0
          call leng_zijk(c,z,zc)
          call leng_zijk(k,z,zk)
!     ---------------
          do 10 j=2,m
!     ---------------
            b=j-1
            e=j+1
            ym=(y(j)+y(j-1))/2.d0
            yp=(y(j)+y(j+1))/2.d0
            ybm=zm*( y(b+1)-y(b) )
            yjm=zm*( y(j+1)-y(j) )
            ybp=zp*( y(b+1)-y(b) )
            yjp=zp*( y(j+1)-y(j) )
!   -----------------------------------------------------
!             set of line elements
!   -----------------------------------------------------
            call area_skij(i,b,k,x,y,z,skib)
            call area_skij(i,j,k,x,y,z,skij)
            call area_sjki(i,j,c,x,y,z,sjci)
            call area_sjki(i,j,k,x,y,z,sjki)
!   -----------------------------------------------------
!             set of surface elements 1
!   -----------------------------------------------------
            rijcy=(resist(i,j,c)*yjm+resist(i,b,c)*ybm)/2.d0/sjci
            rijky=(resist(i,j,k)*yjp+resist(i,b,k)*ybp)/2.d0/sjki
            rijkz=(resist(i,j,k)*zk+resist(i,j,c)*zc)/2.d0/skij
            ribkz=(resist(i,b,k)*zk+resist(i,b,c)*zc)/2.d0/skib
!----------------------------------------------------------------
            jj=jj+1
            sum=0.0d0
!----------------------------------------------------------------
! 1 Hxijc
            if (k == 2) then
!           ----------> top Hx
              call leng_xijk(i,j,c,x,y,z,xx)
              call n_allhxijk(L,M,N,i,j,c,ii)
              sum=sum+rijcy*xx*hx(ii)
            end if
! 2 Hxijf
            if (k == n) then
!           ----------> bottom Hx
              call leng_xijk(i,j,f,x,y,z,xx)
              call n_allhxijk(L,M,N,i,j,f,ii)
              sum=sum+rijky*xx*hx(ii)
            end if
! 3 Hzijc
           if (k == 2) then
!           ----------> top Hz
              call n_allhzijk(L,M,N,i,j,c,ii)
              sum=sum-rijcy*zc*hz(ii)
            end if
! 4 Hzdjc
           if (k == 2) then
!           ----------> top Hz
              call n_allhzijk(L,M,N,d,j,c,ii)
              sum=sum+rijcy*zc*hz(ii)
            end if

            bvec(jj)=sum
10    continue

!
! Now for Hy...
!
!     -----------
      do 20 j=1,m
!     -----------
        b=j-1
        e=j+1
        yp=(y(j)+y(j+1))/2.d0
!     -------------
        do 20 k=2,n
!     -------------
          c=k-1
          f=k+1
          zm=(z(k)+z(k-1))/2.d0
          zp=(z(k)+z(k+1))/2.d0
          call leng_zijk(c,z,zc)
          call leng_zijk(k,z,zk)
          call area_sijk(j,c,y,z,sijc)
          call area_sijk(j,k,y,z,sijk)
!     ---------------
          do 20 i=1,l
!     ---------------
            if (i == 1) then
               a=l
            else
               a=i-1
            end if
            if (i == l) then
               d=1
            else
               d=i+1
            end if

            xipm=zm*dsin(yp)*x(i)
            xapm=zm*dsin(yp)*x(a)
            xipp=zp*dsin(yp)*x(i)
            xapp=zp*dsin(yp)*x(a)
!   -----------------------------------------------------
!             set of line elements
!   -----------------------------------------------------
            call area_skij(a,j,k,x,y,z,skaj)
            call area_skij(i,j,k,x,y,z,skij)
!   -----------------------------------------------------
!             set of surface elements 1
!   -----------------------------------------------------
            rijkz=(resist(i,j,k)*zk+resist(i,j,c)*zc)/2.0d0/skij
            rijcx=(resist(i,j,c)*xipm+resist(a,j,c)*xapm)/2.0d0/sijc
            rajkz=(resist(a,j,k)*zk+resist(a,j,c)*zc)/2.0d0/skaj
            rijkx=(resist(i,j,k)*xipp+resist(a,j,k)*xapp)/2.0d0/sijk
!   -----------------------------------------------------
            jj=jj+1
            sum=0.d0
!----------------------------------------------------------------
! 1 Hyijc
            if (k == 2) then
!           ----------> top Hy
              call leng_yijk(j,c,y,z,yy)
              call n_allhyijk(l,m,n,i,j,c,ii)
              sum=sum+rijcx*yy*hy(ii)
            end if
! 2 Hyijf
            if (k == n) then
!           ----------> bottom Hy
              call leng_yijk(j,f,y,z,yy)
              call n_allhyijk(l,m,n,i,j,f,ii)
              sum=sum+rijkx*yy*hy(ii)
            end if
! 3 Hzijc
            if (k == 2) then
!           ----------> top Hz
              call n_allhzijk(l,m,n,i,j,c,ii)
              sum=sum-rijcx*zc*hz(ii)
            end if
! 4 Hziec
            if (k == 2) then
!           ----------> top Hz
              call n_allhzijk(l,m,n,i,e,c,ii)
              sum=sum+rijcx*zc*hz(ii)
            end if

            bvec(jj)=sum
20    continue

!
! Now for Hz...
!
!     -------------
      do 30 k=2,n
!     -------------
        c=k-1
        f=k+1
        call leng_zijk(k,z,zz)
        zm=(z(k)+z(k-1))/2.0d0
        zp=(z(k)+z(k+1))/2.0d0
!     -------------
        do 30 j=1,m+1
!     -------------
          b=j-1
          e=j+1
!         ++++++++++++++++
          if (j == 1) then
!         < N-pole >
!         ++++++++++++++++
            call area_sijk(j,k,y,z,sijk)
            call leng_yijk(j,f,y,z,yy)
!   -----------------------------------------------------
            jj=jj+1
            sum=0.d0
!----------------------------------------------------------------
! 1 Hyijf
            if (k == n) then
!           ----------> bottom Hy
!             -----------
              do i=1,l
!             -----------
                if (i == 1) then
                  a=l
                else
                  a=i-1
                end if
                if (i == l) then
                  d=1
                else
                  d=i+1
                end if
                xipp=zp*dsin(yp)*x(i)
                xapp=zp*dsin(yp)*x(a)
                rijkx=(resist(i,j,k)*xipp+resist(a,j,k)*xapp)/2.0d0/sijk
                call n_allhyijk(l,m,n,i,j,f,ii)
                sum=sum-rijkx*yy*hy(ii)
!             -----------
              end do
!             -----------
            end if
            bvec(jj)=sum

!         +++++++++++++++++++++++
          else if (j == m+1) then
!         < S-pole >
!         +++++++++++++++++++++++
            call area_sijk(b,k,y,z,sibk)
            call leng_yijk(b,f,y,z,yy)
!   -----------------------------------------------------
            jj=jj+1
            sum=0.d0
!----------------------------------------------------------------
! 1 Hyibf
           if (k == n) then
!           ----------> bottom Hy
!             -----------
              do i=1,l
!             -----------
                if (i == 1) then
                  a=l
                else
                  a=i-1
                end if
                if (i == l) then
                  d=1
                else
                  d=i+1
                end if
                ximp=zp*dsin(ym)*x(i)
                xamp=zp*dsin(ym)*x(a)
                ribkx=(resist(i,b,k)*ximp+resist(a,b,k)*xamp)/2.0d0/sibk
                call n_allhyijk(l,m,n,i,b,f,ii)
                sum=sum+ribkx*yy*hy(ii)
!             -----------
              end do
!             -----------
            end if
            bvec(jj)=sum

!         ++++++++++++++++
          else
!         < neither N-pole nor S-pole >
!         ++++++++++++++++
            yp=(y(j)+y(j+1))/2.0d0
            ym=(y(j)+y(j-1))/2.0d0
            ybp=zp*( y(b+1)-y(b) )
            yjp=zp*( y(j+1)-y(j) )
            call area_sijk(j,k,y,z,sijk)
            call area_sijk(b,k,y,z,sibk)
!     ---------------
          do i=1,l
!     ---------------

            if (i == 1) then
               a=l
            else
               a=i-1
            end if
            if (i == l) then
               d=1
            else
               d=i+1
            end if
            ximp=zp*dsin(ym)*x(i)
            xamp=zp*dsin(ym)*x(a)
            xipp=zp*dsin(yp)*x(i)
            xapp=zp*dsin(yp)*x(a)
!   -----------------------------------------------------
!             set of line elements
!   -----------------------------------------------------
            call area_sjki(i,j,k,x,y,z,sjki)
            call area_sjki(a,j,k,x,y,z,sjka)
!   -----------------------------------------------------
!             set of surface elements 1
!   -----------------------------------------------------
            rijky=(resist(i,j,k)*yjp+resist(i,b,k)*ybp)/2.0d0/sjki
            rijkx=(resist(i,j,k)*xipp+resist(a,j,k)*xapp)/2.0d0/sijk
            rajky=(resist(a,j,k)*yjp+resist(a,b,k)*ybp)/2.0d0/sjka
            ribkx=(resist(i,b,k)*ximp+resist(a,b,k)*xamp)/2.0d0/sibk
!   -----------------------------------------------------
            jj=jj+1
            sum=0.d0
!----------------------------------------------------------------
! 1 Hxijf
            if (k == n) then
!           ----------> bottom Hx
              call leng_xijk(i,j,f,x,y,z,xx)
              call n_allhxijk(l,m,n,i,j,f,ii)
              sum=sum-rijky*xx*hx(ii)
            end if
! 2 Hxajf
            if (k == n) then
!           ----------> bottom Hx
              call leng_xijk(a,j,f,x,y,z,xx)
              call n_allhxijk(l,m,n,a,j,f,ii)
              sum=sum+rajky*xx*hx(ii)
            end if
! 3 Hyijf
            if (k == n) then
!           ----------> bottom Hy
              call leng_yijk(j,f,y,z,yy)
              call n_allhyijk(l,m,n,i,j,f,ii)
              sum=sum-rijkx*yy*hy(ii)
            end if
! 4 Hyibf
            if (k == n) then
!           ----------> bottom Hy
              call leng_yijk(b,f,y,z,yy)
              call n_allhyijk(l,m,n,i,b,f,ii)
              sum=sum+ribkx*yy*hy(ii)
            end if

            bvec(jj)=sum
!           -----------
            end do
!           -----------
!         ++++++++++++++++
          end if
!         ++++++++++++++++
30    continue

      return
      end subroutine calcb


! *****************************************************************************



      subroutine copyd1_d3_a(l,m,n,rx,ry,rz,rvec,x,y,z)
!-------------------------------------------------------------
!
! subroutine to copy 1 vector to 3 vectors
!
! modified to assume boundary values at top of air layers are
! fixed to their starting values
!-------------------------------------------------------------

	  integer					:: l,m,n
      complex(8),dimension(np1) :: rx,ry,rz
      real(8),dimension(l)	    :: x
      real(8),dimension(m+1)    :: y
      real(8),dimension(n+1)    :: z
      complex(8),dimension(np2) :: rvec
      real(8)                   :: xlen,ylen,zlen
      integer                   :: ic,i,j,k,ii

      ic=0
      do i=1,l
        do k=2,n
          do j=2,m
            ic=ic+1
            call leng_xijk(i,j,k,x,y,z,xlen)
            call n_allhxijk(l,m,n,i,j,k,ii)
            rx(ii)=rvec(ic)/xlen
			!print *, 'd1->d3: x ',i,j,k,ii,ic,rvec(ic)
!           rx(ii)=rvec(ic)
          end do
        end do
      end do

      do j=1,m
        do k=2,n
          call leng_yijk(j,k,y,z,ylen)
          do i=1,l
            ic=ic+1
            call n_allhyijk(l,m,n,i,j,k,ii)
            ry(ii)=rvec(ic)/ylen
			!print *, 'd1->d3: y ',i,j,k,ii,ic,rvec(ic)
!           ry(ii)=rvec(ic)
          end do
        end do
      end do

	  ! Values corresponding to k=1 are stored in the first l*(m-1)+2
	  ! entries of rz, and are never changed from their initial values
      do k=2,n
        ic=ic+1
        call leng_zijk(k,z,zlen)
        call n_allhzijk(l,m,n,1,1,k,ii)
        rz(ii)=rvec(ic)/zlen
!       rz(ii)=rvec(ic)
        do j=2,m
          do i=1,l
            ic=ic+1
            call n_allhzijk(l,m,n,i,j,k,ii)
            rz(ii)=rvec(ic)/zlen
			!print *, 'd1->d3: z ',i,j,k,ii,ic,rvec(ic)
!           rz(ii)=rvec(ic)
          end do
        end do
        ic=ic+1
        call n_allhzijk(l,m,n,1,m+1,k,ii)
        rz(ii)=rvec(ic)/zlen
!       rz(ii)=rvec(ic)
      end do

      return
      end subroutine copyd1_d3_a


! *****************************************************************************


      subroutine copyd1_pot_bBC3(L,M,N,vi,vo)

!-------------------------------------------------------------
! subroutine to copy 1 vector to potential vector
!
! modified to assume the boundary values at the top of the
! air layers are fixed to their starting values.
!toh:29/FEB/2000
!-------------------------------------------------------------


      real(8),dimension(np1) :: vi,vo
      integer                :: ic,k,l,m,n,i,ii,j

      ic=0
      do k=2,n
!toh:29/FEB/2000
        j=1
          ic=ic+1
          call n_allpotijk(l,m,n,0,j,k,ii)
          vo(ii)=vi(ic)
        do j=2,m
          do i=1,l
            ic=ic+1
            call n_allpotijk(l,m,n,i,j,k,ii)
            vo(ii)=vi(ic)
          end do
        end do
        j=m+1
          ic=ic+1
          call n_allpotijk(l,m,n,0,j,k,ii)
          vo(ii)=vi(ic)
      end do

      return
      end subroutine copyd1_pot_bBC3


! *****************************************************************************


      subroutine copyd3_d1_a(l,m,n,rx,ry,rz,rvec,x,y,z)

!-------------------------------------------------------------
! subroutine to copy 3 vectors to 1 vector
!-------------------------------------------------------------

	  integer					 :: l,m,n
      complex(8),dimension(np1)  :: rx,ry,rz
      real(8),dimension(l)	     :: x
      real(8),dimension(m+1)     :: y
      real(8),dimension(n+1)     :: z
      complex(8),dimension(np2)  :: rvec
      real(8)                    :: xlen,ylen,zlen
      integer                    :: ic,i,j,k,ii

      ic=0
      do i=1,l
        do k=2,n
          do j=2,m
            ic=ic+1
            call leng_xijk(i,j,k,x,y,z,xlen)
            call n_allhxijk(l,m,n,i,j,k,ii)
            rvec(ic)=rx(ii)*xlen
			!print *, 'd3->d1: x ',i,j,k,ii,ic,rvec(ic)
!           rvec(ic)=rx(ii)
          end do
        end do
      end do

      do j=1,m
        do k=2,n
          call leng_yijk(j,k,y,z,ylen)
          do i=1,l
            ic=ic+1
            call n_allhyijk(l,m,n,i,j,k,ii)
            rvec(ic)=ry(ii)*ylen
			!print *, 'd3->d1: y ',i,j,k,ii,ic,rvec(ic)
!           rvec(ic)=ry(ii)
          end do
        end do
      end do

	  ! Values corresponding to k=1 are stored in the first l*(m-1)+2
	  ! entries of rz, and are never changed from their initial values
      do k=2,n
        ic=ic+1
        call leng_zijk(k,z,zlen)
        call n_allhzijk(l,m,n,1,1,k,ii)
        rvec(ic)=rz(ii)*zlen
		!i=1;j=1
		!print *, 'd3->d1: z ',i,j,k,ii,ic,rvec(ic)
!       rvec(ic)=rz(ii)
        do j=2,m
          do i=1,l
            ic=ic+1
            call n_allhzijk(l,m,n,i,j,k,ii)
            rvec(ic)=rz(ii)*zlen
			!print *, 'd3->d1: z ',i,j,k,ii,ic,rvec(ic)
!           rvec(ic)=rz(ii)
          end do
        end do
        ic=ic+1
        call n_allhzijk(l,m,n,1,m+1,k,ii)
        rvec(ic)=rz(ii)*zlen
		!i=1;j=m+1
		!print *, 'd3->d1: z ',i,j,k,ii,ic,rvec(ic)
!       rvec(ic)=rz(ii)
      end do

      return
      end subroutine copyd3_d1_a


! *****************************************************************************



      subroutine d3residpot3BC3(l,m,n,div,r)

!-------------------------------------------------------------
! subroutine to calculate the inital residual needed in the
! divH=0 correction
!
! the initial residual is actually only the divergence of the
! H field since initially the potential is set equal to zero.
!
! modified to assume boundary values at top of air layers are
! fixed to their starting values
!
! modified for real space values
!toh:29/FEB/2000
!-------------------------------------------------------------


      real(8),dimension(np1)       :: div,r
      integer                      :: ic,i,j,k,l,m,n,ii,kk

      ic=0
      do k=2,n
!toh:29/FEB/2000
        j=1
          ic=ic+1
          r(ic)=div(ic)

          call n_dvpotijkBC3(l,m,n,0,j,k,ii)
          if (ii.ne.ic) then
			inquire(ioERR,OPENED=opened)
			if (opened) then
			  write(ioERR,*)'warning.. ii vs ic',ii,ic
			end if
          end if

        do j=2,m
          do i=1,l
            ic=ic+1
            r(ic)=div(ic)

          call n_dvpotijkBC3(l,m,n,i,j,k,ii)
          if (ii.ne.ic) then
			inquire(ioERR,OPENED=opened)
			if (opened) then
			  write(ioERR,*)'warning.. ii vs ic',ii,ic
			end if
          end if

          end do
        end do
        j=m+1
          ic=ic+1
          r(ic)=div(ic)

          call n_dvpotijkBC3(l,m,n,0,j,k,ii)
          if (ii.ne.ic) then
			inquire(ioERR,OPENED=opened)
			if (opened) then
			  write(ioERR,*)'warning.. ii vs ic',ii,ic
			end if
          end if

      end do

      return
      end subroutine d3residpot3BC3


! *****************************************************************************



      subroutine zdot(v1,v2,zresult,np)
!-------------------------------------------------------------
! subroutine to compute the dot product of two double precision
! complex vectors  - not conjugated
! by: A. Schultz, 18.6.98, ITG Cantab
!-------------------------------------------------------------


      complex(8),dimension(np2)   :: v1,v2
      complex(8)                  :: zresult,zsum
      integer                     :: i,np

      zsum=dcmplx(0.0d0,0.0d0)

      do i=1,np
        zsum=zsum+v1(i)*v2(i)
      end do

      zresult = zsum

      return
      end subroutine zdot

! *****************************************************************************



      subroutine zerr(v1,a2,v2,err,np)
!-------------------------------------------------------------
! subroutine to compute error in most recent iteration
! the error is defined as absolute value of change in field
! divided by the field.
!-------------------------------------------------------------


      complex(8),dimension(np2)      :: v1,v2
      complex(8)                     :: a2,err,sum1,sum2,dum
      integer                        :: i,np

      sum1=0.0d0
      do i=1,np
        dum=a2*v2(i)
        sum1=sum1+dconjg(dum)*dum
      end do

      sum2=0.0d0
      do i=1,np
        sum2=sum2+dconjg(v1(i))*v1(i)
      end do

      err=cdsqrt(sum1/sum2)

      return

      end subroutine zerr

! *****************************************************************************



      subroutine newzerr(v1,a2,v2,err,np,alpha,w,x,y,z,l,m,n,nhx,nhy,nhz,iflg)
!-------------------------------------------------------------------
!toh: This new subroutine is designed for computation of a weighted
!toh: sum, rather than a former simple sum, of the numerator (the
!toh: change in field) alone to put higher priority on either surface
!toh: or CMB (for future application of the code to screening effects)
!toh: field. Also, penalty on roughness can be imposed, if activated,
!toh: on the field components (the denominator) with smaller weights
!toh: to suppress their ill behavior due to lack of constraints. The
!toh: trade-off between the weighted sum (the numerator) and the penalty
!toh: (the denominator) can be an arbitrary choice (unity etc.) at this
!toh: moment since we just can not afford time to calculate it (by ABIC
!toh: and so on). Life's too short for decent science but a night's too
!toh: long for a bottle of SAKI !
!toh:
!toh: By H Toh, ITG, University of Cambridge on 27/FEB/99
!toh:
!toh: Calls 'calrough' to compute roughness if alpha > 0.
!toh:
!toh: alpha : trade-off between the weighted sum and the penalty.
!toh:         0 to disable this feature
!toh:         positive to activate penalty on roughness
!toh:
!toh: w : the weight function (the soul of this routine) which is NOT
!toh:     calculated here. It's a waste of time to compute this in the
!toh:     relaxation loop. The function will be calculated only once
!toh:     in 'sph3d' by calling 'calwh' before relaxation starts.
!-------------------------------------------------------------------

	  integer						 :: l,m,n
      complex(8),dimension(np2)      :: v1,v2
      real(8),   dimension(np2)      :: w
      real(8),dimension(l)	         :: x
      real(8),dimension(m+1)         :: y
      real(8),dimension(n+1)         :: z
      complex(8)                     :: a2,err,dum
      real(8)                        :: sum,rough,alpha,r
      integer                        :: i,np,nhx,nhy,nhz,iflg

      rough = 0.0d0

      sum = 0.0d0

      if ( iflg == 0 ) then

         do i = 1, np

            dum = a2*v2(i)
            sum = sum+dconjg(dum)*dum

         end do

      else

         do i = 1, np

            dum = a2*w(i)*v2(i)
            sum = sum+dconjg(dum)*dum

         end do

      end if

!toh: It's possible here before calling calrough to apply w to v1...

      if ( alpha > 0.0d0 ) call calrough(v1,np,rough,x,y,z,l,m,n,nhx,nhy,nhz)

      sum = sum/dot_product(v1(1:np),v1(1:np)) + alpha*rough

!      print *, 'rough in newzerr =', rough

      err=sqrt(sum)

      return

      end subroutine newzerr

! *****************************************************************************

      subroutine calwh(w,np,iflg,omega,z,k0,l,m,n,nhx,nhy,nhz)

!toh: iflg = 0 : equivalent to the original 'zerr'    (put equal weights)
!toh:      >   : put higher priority to SURFACE field (positive iflg)
!toh:      <   : put higher priority to CMB field     (negative iflg)
!toh:
!toh: rho : sdpt is a skin depth in a homogeneous media of resistivity, rho.
!toh:
!toh: calls n_h?ijk to idetify where the magnetic component lies.
!toh:
!toh: H Toh @ Univ. Gegege on 3/MAR/99 (ITG, Univ. Cambridge)


	  integer						 :: l,m,n
      real(8),dimension(np2)         :: w
      real(8),dimension(n+1)         :: z
      real(8)                        :: omega,sdpt,rho,r !AK: sum
      integer                        :: np,k0,k1,nhx,nhy,nhz,iflg
      integer                        :: ii,i,j,k,js,je

!     pi = 4.0d0*datan(1.0d0)
      rho = 100.d0
      sdpt = dsqrt(10.0d0*rho/omega)*1.d+3 ! skin depth in m

      if ( iflg > 0 ) then

!CASE I: Now assign weights as a function of radius. First for Hx...

         do i = 1, l

            do j = 2, m

               do k = 2, n

                  call n_hxijk(l,m,n,i,j,k,ii)
                  r = z(k) - z(k0)

                  if ( k >= k0 ) then

                     w(ii) = exp(r/sdpt)  ! Decay inwards from SRF

                  else

                     w(ii) = exp(-r/sdpt) ! Decay outwards from SRF

                  end if

               end do

            end do

         end do

!Then, Hy...

         do i = 1, l

            do j = 1, m

               do k = 2, n

                  call n_hyijk(nhx,l,m,n,i,j,k,ii)
                  r = z(k) - z(k0)

                  if ( k >= k0 ) then

                     w(ii) = exp(r/sdpt)  ! Decay inwards from SRF

                  else

                     w(ii) = exp(-r/sdpt) ! Decay outwards from SRF

                  end if

               end do

            end do

         end do

!Finally, Hz...

         do i = 1, l

            if ( i == 1 ) then

               js = 1
               je = m+1

            else

               js = 2
               je = m

            end if

            do j = js, je

               do k = 2, n

                  call n_hzijk(nhx,nhy,l,m,n,i,j,k,ii)
                  r = z(k) - z(k0)

                  if ( k >= k0 ) then

                     w(ii) = exp(r/sdpt)  ! Decay inwards from SRF

                  else

                     w(ii) = exp(-r/sdpt) ! Decay outwards from SRF

                  end if

               end do

            end do

         end do

      else if ( iflg < 0 ) then

!CASE II: Assign weights as a function of radius. First for Hx...

         k1 = ( n - k0 + 1 )/2 + k0 - 1

         do i = 1, l

            do j = 2, m

               do k = 2, n

                  call n_hxijk(l,m,n,i,j,k,ii)

                  if ( k >= k0 .and. k <= k1 ) then

                     w(ii) = exp((z(k)-z(k0))/sdpt) ! Decay inwards from SRF

                  else if ( k > k1 ) then

                     w(ii) = exp((z(n)-z(k))/sdpt)  ! Decay outwards from CMB

                  else

                     w(ii) = exp((z(k0)-z(k))/sdpt) ! Decay outwards from SRF

                  end if

               end do

            end do

         end do

!Then, Hy...

         do i = 1, l

            do j = 1, m

               do k = 2, n

                  call n_hyijk(nhx,l,m,n,i,j,k,ii)

                  if ( k >= k0 .and. k <= k1 ) then

                     w(ii) = exp((z(k)-z(k0))/sdpt) ! Decay inwards from SRF

                  else if ( k > k1 ) then

                     w(ii) = exp((z(n)-z(k))/sdpt)  ! Decay outwards from CMB

                  else

                     w(ii) = exp((z(k0)-z(k))/sdpt) ! Decay outwards from SRF

                  end if

               end do

            end do

         end do

!Finally, Hz...

         do i = 1, l

            if ( i == 1 ) then

               js = 1
               je = m+1

            else

               js = 2
               je = m

            end if

            do j = js, je

               do k = 2, n

                  call n_hzijk(nhx,nhy,l,m,n,i,j,k,ii)

                  if ( k >= k0 .and. k <= k1 ) then

                     w(ii) = exp((z(k)-z(k0))/sdpt) ! Decay inwards from SRF

                  else if ( k > k1 ) then

                     w(ii) = exp((z(n)-z(k))/sdpt)  ! Decay outwards from CMB

                  else

                     w(ii) = exp((z(k0)-z(k))/sdpt) ! Decay outwards from SRF

                  end if

               end do

            end do

         end do

      end if

! Normalize weights so that sum of them is equal to np.

      if ( iflg /= 0 ) w(1:np) = w(1:np)/sum(w(1:np))*dfloat(np)

      return

      end subroutine calwh

! *****************************************************************************

      subroutine calrough(v1,np,rough,x,y,z,l,m,n,nhx,nhy,nhz)

!toh: In this routine, roughness is defined as sum of 2nd order difference
!toh: of Hx, Hy & Hz at each shell, viz.,
!toh:
!toh: roughness = [sum over shells](Hx-roughness + Hy-roughness + Hz-roughness)
!toh:
!toh: At a (i,j) node, the 2nd order difference of a field component is
!toh: given by:
!toh:           f(i-1,j,k) - 2f(i,j,k) + f(i+1,j,k)
!toh:         + f(i,j-1,k) - 2f(i,j,k) + f(i,j+1,k)
!toh:         + f(i,j,k-1) - 2f(i,j,k) + f(i,j,k+1)
!toh:
!toh: (Here,I dropped geometric factors for simplicity.)
!toh: However, thinking of the rather expensive nature of this computation,
!toh: it is advisable to users of this code that first disable this function by
!toh: setting alpha = 0, then gradually enlarge the value of alpha to suppress
!toh: very rough behaviors of any components.
!toh:
!toh: Calls n_h?ijk to P/U adjacent field values from hvec.
!toh:
!toh: H Toh @ Univ. Gegege on 3/MAR/99 (ITG, Univ. Cambridge)


	  integer						 :: l,m,n
      complex(8),dimension(np2)      :: v1
      complex(8)                     :: fx1,fx2,fx3,fy1,fy2,fy3,fz1,fz2,fz3
      complex(8)                     :: rx, ry, rz
      real(8),dimension(l)	         :: x
      real(8),dimension(m+1)         :: y
      real(8),dimension(n+1)         :: z
      real(8)                        :: rough, z2
      integer                        :: np,nhx,nhy,nhz,ii,i1,i2,i3,i4
	  integer						 :: i,j,k

! Now cal roughness first for Hx...

      do i = 1, l

         if ( i == 1 ) then

            i1 = l

         else

            i1 = i-1

         end if

         if ( i == l ) then

            i2 = 1

         else

            i2 = i+1

         end if

         do j = 3, m-1

            do k = 3, n-1

               call n_hxijk(l,m,n,i1,j,k,ii)
               fx1 = v1(ii)
               call n_hxijk(l,m,n,i,j,k,ii)
               fx2 = v1(ii)
               call n_hxijk(l,m,n,i2,j,k,ii)
               fx3 = v1(ii)

               call n_hxijk(l,m,n,i,j-1,k,ii)
               fy1 = v1(ii)
               call n_hxijk(l,m,n,i,j,k,ii)
               fy2 = v1(ii)
               call n_hxijk(l,m,n,i,j+1,k,ii)
               fy3 = v1(ii)

               call n_hxijk(l,m,n,i,j,k-1,ii)
               fz3 = v1(ii)
               call n_hxijk(l,m,n,i,j,k,ii)
               fz2 = v1(ii)
               call n_hxijk(l,m,n,i,j,k+1,ii)
               fz1 = v1(ii)

               z2 = ( z(k) + z(k-1) )/2.0d0

               rx = ( (fx3-fx2)/(x(i)*z2*sin(y(j))) &
                    - (fx2-fx1)/(x(i1)*z2*sin(y(j))) ) &
                    /(z2*sin(y(j))*( x(i) + x(i1) )/2.0d0) &
                    + ( (fy3-fy2)/((y(j+1)-y(j))*z2) &
                      - (fy2-fy1)/((y(j)-y(j-1))*z2) ) &
                      /(z2*( y(j+1) - y(j-1) )/2.0d0) &
                    + ( (fz3-fz2)/(z(k-1)-z(k)) - (fz2-fz1)/(z(k)-z(k-1)) ) &
                      /(( z(k-1) - z(k+1) )/2.0d0)

            end do

         end do

      end do

!Then, Hy...

      do i = 1, l

         if ( i == 1 ) then

            i1 = l

         else

            i1 = i-1

         end if

         if ( i == l ) then

            i2 = 1

         else

            i2 = i+1

         end if

         do j = 2, m-1

            do k = 3, n-1

               call n_hyijk(nhx,l,m,n,i1,j,k,ii)
               fx1 = v1(ii)
               call n_hyijk(nhx,l,m,n,i,j,k,ii)
               fx2 = v1(ii)
               call n_hyijk(nhx,l,m,n,i2,j,k,ii)
               fx3 = v1(ii)

               call n_hyijk(nhx,l,m,n,i,j-1,k,ii)
               fy1 = v1(ii)
               call n_hyijk(nhx,l,m,n,i,j,k,ii)
               fy2 = v1(ii)
               call n_hyijk(nhx,l,m,n,i,j+1,k,ii)
               fy3 = v1(ii)

               call n_hyijk(nhx,l,m,n,i,j,k-1,ii)
               fz3 = v1(ii)
               call n_hyijk(nhx,l,m,n,i,j,k,ii)
               fz2 = v1(ii)
               call n_hyijk(nhx,l,m,n,i,j,k+1,ii)
               fz1 = v1(ii)

               z2 = ( z(k) + z(k-1) )/2.0d0

               ry = ( (fx3-fx2)/(x(i)*z2*sin(y(j))) &
                     -(fx2-fx1)/(x(i1)*z2*sin(y(j))) ) &
                    /(z2*sin(y(j))*( x(i) + x(i1) )/2.0d0) &
                    + ( (fy3-fy2)/((y(j+1)-y(j))*z2) &
                      - (fy2-fy1)/((y(j)-y(j-1))*z2) ) &
                      /(z2*( y(j+1) - y(j-1) )/2.0d0) &
                    + ( (fz3-fz2)/(z(k-1)-z(k)) - (fz2-fz1)/(z(k)-z(k-1)) ) &
                      /(( z(k-1) - z(k+1) )/2.0d0)

            end do

         end do

      end do

!Finally, Hz...

      do i = 1, l

         if ( i == 1 ) then

            i1 = l

         else

            i1 = i-1

         end if

         if ( i == l ) then

            i2 = 1

         else

            i2 = i+1

         end if

         do j = 2, m

            if ( j == 2 ) then

               i3 = 1

            else

               i3 = i

            end if

            if ( j == m ) then

               i4 = 1

            else

               i4 = i

            end if

            do k = 3, n-2 ! 1:n+1 for grid => 1:n for Hz => 3:n-2 for roughness

               call n_hzijk(nhx,nhy,l,m,n,i1,j,k,ii)
               fx1 = v1(ii)
               call n_hzijk(nhx,nhy,l,m,n,i,j,k,ii)
               fx2 = v1(ii)
               call n_hzijk(nhx,nhy,l,m,n,i2,j,k,ii)
               fx3 = v1(ii)

               call n_hzijk(nhx,nhy,l,m,n,i3,j-1,k,ii)
               fy1 = v1(ii)
               call n_hzijk(nhx,nhy,l,m,n,i,j,k,ii)
               fy2 = v1(ii)
               call n_hzijk(nhx,nhy,l,m,n,i4,j+1,k,ii)
               fy3 = v1(ii)

               call n_hzijk(nhx,nhy,l,m,n,i,j,k-1,ii)
               fz3 = v1(ii)
               call n_hzijk(nhx,nhy,l,m,n,i,j,k,ii)
               fz2 = v1(ii)
               call n_hzijk(nhx,nhy,l,m,n,i,j,k+1,ii)
               fz1 = v1(ii)

               z2 = ( z(k) + z(k-1) )/2.0d0

               rz = ( (fx3-fx2)/(x(i)*z2*sin(y(j))) &
                     -(fx2-fx1)/(x(i1)*z2*sin(y(j)))) &
                    /(z2*sin(y(j))*( x(i) + x(i1) )/2.0d0) &
                    + ( (fy3-fy2)/((y(j+1)-y(j))*z2) &
                      - (fy2-fy1)/((y(j)-y(j-1))*z2) ) &
                      /(z2*( y(j+1) - y(j-1) )/2.0d0) &
                    + ( (fz3-fz2)/(z(k-1)-z(k)) - (fz2-fz1)/(z(k)-z(k-1)) ) &
                      /(( z(k-1) - z(k+1) )/2.0d0)

            end do

         end do

      end do

! Take sum over components and the absolute value.

      rough = cdabs( rx + ry + rz )

      return

      end subroutine calrough

! *****************************************************************************

!toh
      subroutine icdecompc(ija,sb,np)
!toh
!-------------------------------------------------------------
! subroutine to do incomplete M=H*H(t) decomposition with
! row-indexed matrix.  First, make a copy of sa to sb, because
! sa need to be reserved, used in atimes routine. This routine
! returns H(t) (upper triangular matrix) in sb, which will be
! used in asolver.
!-------------------------------------------------------------


!toh
      complex(8),dimension(np4)    :: sb
!toh
      integer,dimension(np4)       :: ija
      integer                      :: i,ii,iii,n,k,m,md,np

      do k=1, np
!toh
         sb(k)=cdsqrt(sb(k))
!toh
         do m=ija(k), ija(k+1)-1

            if (sb(k) == 0.0d0) then

               ! (AK) print *, 'k=',k
               ! (AK) print *, 'm=',m

            end if

            sb(m)=sb(m)/sb(k)

          end do

          md=ija(k+1)-ija(k)
          ii=ija(k)

          if (md >= 1) then

             do iii=ii,ii+md-1

                sb(ija(iii))=sb(ija(iii))-sb(iii)**2

             end do

          end if

      end do

      ! (AK) print *,'end of icdecomp......'

      return
      end subroutine icdecompc


! *****************************************************************************


      subroutine icdecompr(n,ija,sb,l)
!-------------------------------------------------------------
! subroutine to do incomplete M=H*H(t) decomposition with
! row-indexed matrix.  First, make a copy of sa to sb, because
! sa need to be reserved, used in atimes routine. This routine
! returns H(t) (upper triangular matrix) in sb, which will be
! used in asolver.
!-------------------------------------------------------------

      integer,dimension(np6)             :: ija
      real(8),dimension(np6)             :: sb
      integer                            :: md,i,j,k,l,m,n,ii,iii

      do k=1, n
         if (sb(k) <= 0.0d0) then
            print *,'cholesky decomposition fails at...'
            print *,'k,sbk = ',k,sb(k)
            stop
         end if
         sb(k)=dsqrt(sb(k))
         do m=ija(k), ija(k+1)-1
            sb(m)=sb(m)/sb(k)
            if (m == 1844) then
               ! (AK) print *,'1844 was changed!!!'
            end if
         end do

         md=ija(k+1)-ija(k)
         ii=ija(k)
         if(md >= 1) then
            do iii=ii,ii+md-1
              sb(ija(iii))=sb(ija(iii))-sb(iii)**2
            end do
         end if

      end do

      return
      end subroutine icdecompr

! *****************************************************************************

	         subroutine initpot3uBC3(l,m,n,np,sp,sp1,ija,x,y,z)

!-------------------------------------------------------------
! subroutine to initiate the incomplete LDL decomposition of the
! operator used to correct the H field by a divH=0 computation.
! This subroutine determines the locations of the
! non-zero elements of the operator and stores them in a form that
! is easy for the decomposition algorithm to access. Remember that
! the operator is symmetric, so that only half of the non-zero
! elements are stored.
!
! modified to assume boundary values at top of air layers are fixed
! to their starting values
!
! modified to use real space values
! (13-01-98) modified to fix the pot at the first shell as zero
!-------------------------------------------------------------

      implicit none

      integer                       :: l,m,n,np
      integer,dimension(np6)        :: ija
      real(8),dimension(np6)        :: sp,sp1
      real(8),dimension(l)	        :: x
      real(8),dimension(m+1)        :: y
      real(8),dimension(n+1)        :: z

      real(8)                       :: xijk,xajk,yijk,yibk,zijk,zijc
      real(8)                       :: sipmm,sjpmm,skpmm
      real(8)                       :: simmm,sjmmm,skmmm
      real(8)                       :: zp,zm,yp,ym,skpcap,skmcap,sum
      integer                       :: a,d
      integer                       :: i,j,k,ii,jj,kk

      integer                       :: ii1,ii2

!     pi=4.d0*datan(1.d0)
      jj=0
      ija(1)=np+2
      kk=np+1
!toh:29/FEB/2000
      do k=2,n
!          -------> variables pot are started from the 2nd layer..

        j=1
          call leng_yijk(j,k,y,z,yijk)
          call leng_zijk(k  ,z,zijk)
          call leng_zijk(k-1,z,zijc)
          zp=(z(k)+z(k+1))/2.d0
          zm=(z(k)+z(k-1))/2.d0
          yp=(y(j)+y(j+1))/2.d0
          skpcap=2.d0*pi*(zp**2)*(1.d0-dcos(yp))
          skmcap=2.d0*pi*(zm**2)*(1.d0-dcos(yp))

! Phi(i,j,k)

          jj=jj+1
          sp(jj)=skpcap/zijk+skmcap/zijc
          do i=1,l
            if (i.ne.1) then
              a=i-1
            else
              a=l
            end if
            call area_sjki2(a,i,j,k-1,x,y,z,sjpmm)
            sp(jj)=sp(jj)+sjpmm/yijk

! Phi(i,j+1,k)

            kk=kk+1
            call n_dvpotijkBC3(l,m,n,i,j+1,k,ii)
            sp(kk)=-sjpmm/yijk
            ija(kk)=ii
          end do

! Phi(i,j,k+1)

          if (k.ne.n) then
            kk=kk+1
            call n_dvpotijkBC3(l,m,n,i,j,k+1,ii)
            sp(kk)=-skpcap/zijk
            ija(kk)=ii
          end if
          ija(jj+1)=kk+1
!       end j=1

        do j=2,m
          do i=1,l
            if (i.ne.1) then
              a=i-1
            else
              a=l
            end if
            if (i.ne.l) then
              d=i+1
            else
              d=1
            end if

            call leng_xijk(i,j,k,x,y,z,xijk)
            call leng_xijk(a,j,k,x,y,z,xajk)
            call leng_yijk(j,k,y,z,yijk)
            call leng_yijk(j-1,k,y,z,yibk)
            call leng_zijk(k,z,zijk)
            call leng_zijk(k-1,z,zijc)
            call area_sijk2(j-1,k-1,y,z,sipmm)
            call area_sijk2(j-1,k-1,y,z,simmm)
            call area_sjki2(a,i,j  ,k-1,x,y,z,sjpmm)
            call area_sjki2(a,i,j-1,k-1,x,y,z,sjmmm)
            call area_skij2(a,i,j-1,k  ,x,y,z,skpmm)
            call area_skij2(a,i,j-1,k-1,x,y,z,skmmm)

! Phi(i,j,k)

            jj=jj+1
            sp(jj)=sipmm/xijk+sjpmm/yijk+skpmm/zijk+ &
                   simmm/xajk+sjmmm/yibk+skmmm/zijc

! Phi(i+1,j,k)

           if (d.ne.1) then
            kk=kk+1
            call n_dvpotijkBC3(l,m,n,d,j,k,ii)
            sp(kk)=-sipmm/xijk
            ija(kk)=ii
           end if

! Phi(i-1,j,k)

           if (a == L) then
            kk=kk+1
            call n_dvpotijkBC3(l,m,n,a,j,k,ii)
            sp(kk)=-simmm/xajk
            ija(kk)=ii
           end if

! Phi(i,j+1,k)

            kk=kk+1
            call n_dvpotijkBC3(l,m,n,i,j+1,k,ii)
            sp(kk)=-sjpmm/yijk
            ija(kk)=ii

! Phi(i,j,k+1)

            if(k.ne.N) then
              kk=kk+1
              call n_dvpotijkBC3(l,m,n,i,j,k+1,ii)
              sp(kk)=-skpmm/zijk
              ija(kk)=ii
            end if
            ija(jj+1)=kk+1
          end do
        end do

        j=m+1
          call leng_yijk(j-1,k,y,z,yibk)
          call leng_zijk(k  ,z,zijk)
          call leng_zijk(k-1,z,zijc)
          zp=(z(k)+z(k+1))/2.d0
          zm=(z(k)+z(k-1))/2.d0
          ym=(y(j)+y(j-1))/2.d0
          skpcap=2.d0*pi*(zp**2)*(1.d0+dcos(ym))
          skmcap=2.d0*pi*(zm**2)*(1.d0+dcos(ym))

! Phi(i,j,k)

          jj=jj+1
          sp(jj)=skpcap/zijk+skmcap/zijc
          do i=1,l
            if (i.ne.1) then
              a=i-1
            else
              a=l
            end if
            call area_sjki2(a,i,j-1,k-1,x,y,z,sjmmm)
            sp(jj)=sp(jj)+sjmmm/yijk
          end do

! Phi(i,j,k+1)

          if (k.ne.n) then
            kk=kk+1
            call n_dvpotijkBC3(l,m,n,i,j,k+1,ii)
            sp(kk)=-skpcap/zijk
            ija(kk)=ii
          end if
          ija(jj+1)=kk+1
!---- end j=M+1

      end do

! copy sp to sp1, which is used in later calculations


      sp1(1:kk)=sp(1:kk)


      ! (AK) print *,'kk=',kk
      ! (AK) print *,'end initpot....'

      return
      end subroutine initpot3uBC3

! *****************************************************************************


      subroutine relaxpot_chol3BC3(l,m,n,pot,div,r,z,q,p,s,d,sp, &
                               sp1,ija,wk,nrel)
!-------------------------------------------------------------
! subroutine to do min rTr relaxation of divH=0 problem
! the results are used to correct the H field in the
! real relaxation problem
!
! an incomplete ll(t) decomposition is used to accelerate the
! relaxation
!
! modified to assume boundary values on top of air layers are
! fixed to their starting values
!-------------------------------------------------------------

      integer                  :: l,m,n,ic,nrel,nclr,npl

      real(8),dimension(np1)   :: div,pot
      real(8),dimension(np1)   :: r,z,q,wk
      real(8),dimension(np1)   :: p,s,d
      real(8)                  :: c1,b1,e1,cone,d1,l1,cmone

      real(8),dimension(np6)   :: sp,sp1
      integer,dimension(np6)   :: ija



! clear vectors

      nclr=(2+l*(m-1))*(n+1)
      pot(1:nclr) = 0.0d0
      wk(1:nclr)  = 0.0d0
      q(1:nclr)   = 0.0d0
      r(1:nclr)   = 0.0d0
      z(1:nclr)   = 0.0d0
      p(1:nclr)   = 0.0d0
      s(1:nclr)   = 0.0d0
      d(1:nclr)   = 0.0d0

! calculate initial residual

      npl=(2+l*(m-1))*(n-1)
!toh:29/FEB/2000
      cone=1.0d0
      cmone=-1.0d0
      call d3residpot3BC3(l,m,n,div,r)

! begin relaxation
! set e1=0.0d0 (this was omitted & appears to be a bug  -  A.S. 23.2.99
!

      e1 = 0.0d0      ! A.S.


      do ic=1,nrel
        call asolver(npl,r,z,sp,ija)

! update vectors

        call atimesr(npl,z,q,sp1,ija)

        c1 = dot_product(q(1:npl),s(1:npl))

        if((dabs(e1) == 0.0d0).or.(ic == 1)) then
          b1=0.0d0
        else
          b1=-c1/e1
        end if
        d(1:npl) = cmone*z(1:npl) + b1*d(1:npl)
        p(1:npl) =  cone*q(1:npl) + b1*p(1:npl)
        call asolver(npl,p,s,sp,ija)

        e1 = dot_product(p(1:npl),s(1:npl))
        d1 = dot_product(z(1:npl),p(1:npl))

        if(dabs(e1) == 0.0d0) then
          l1=0.0d0
        else
          l1=-d1/e1
        end if
        wk(1:npl) = cone*wk(1:npl) + l1*d(1:npl)
        r(1:npl)  = cone*r(1:npl)  + l1*p(1:npl)
      end do

      call copyd1_pot_bBC3(l,m,n,wk,pot)

      return
      end subroutine relaxpot_chol3BC3


! *****************************************************************************

!toh
       subroutine set_op3u(l,m,n,sp,ija,sp1,ija1,resist,x,y,z, &
                          omega,dp,np,rowlen)
!toh
!-------------------------------------------------------------
! subroutine to initialize the sparse operator that we are trying
! to solve. This subroutine determines the locations of the
! non-zero elements of the operator and stores them in a form that
! is easy for the decomposition algorithm to access. Remember that
! the operator is symmetric, so that only half of the non-zero
! elements are stored.
!
! also stores the diagonal sub-blocks used in the incomplete
! cholesky decomposition
!
! modified so that the boundary values for k=1 are fixed to their
! starting values
!
! modified so that the boundary values for k=n+1 are fixed to their
! starting values
!
! modified for real space values
!
! modified for spherical problem (25/11/97)
! modified for upper triangular storage into sp (4/12/97)
!
! modified by Anna Kelbert (20/10/2010) to output the maximum
! "row length" of sp, computed as max(ija(i+1)-ija(i)).
! Use this in atimes to preallocate temporary arrays for efficiency.
!-------------------------------------------------------------

      implicit none

      integer                            :: np,nhx,nhy,nhz,l,m,n,i,j
      integer                            :: k,ii,jj,jj1,kk,kk1
      integer,dimension(np5)             :: ija
      integer,dimension(np4)             :: ija1
      integer                            :: a,b,c,d,e,f
      integer                            :: ijalen,rowlen

      real(8),dimension(l,m,n)	         :: resist
      real(8)                            :: omega

      real(8),dimension(l)	             :: x
      real(8),dimension(m+1)             :: y
      real(8),dimension(n+1)             :: z

      real(8)                            :: rijcy,rijky,rijkz,ribkz
      real(8)                            :: rijkx,rijcx,rajkz,rajky,ribkx
      real(8)                            :: uw
      real(8)                            :: xx,yy,zz
      real(8)                            :: xm,xp,ym,yp,zm,zp
      real(8)                            :: zc,zk
      real(8)                            :: xipm,xapm,xipp,xapp,ybm
      real(8)                            :: yjm,ybp,yjp,ximp,xamp
      real(8)                            :: skib,skij,sjci,sjki,sijk2
      real(8)                            :: skaj,sijc,sijk,sjki2
      real(8)                            :: sibk,sjka,skij2

!toh
      complex(8),dimension(np)           :: dp
      real(8),dimension(np+1:np5)        :: sp
      complex(8),dimension(np4)          :: sp1
!toh
      complex(8)                         :: ciuw

      nhx=l*(m-1)*(n-1)
      nhy=l*m*(n-1)

!     nhz=( l*(m-1)+2 )*(n-1)
!       nhx <-- number of variables to be determined as H(phai)
!       nhy <-- number of variables to be determined as H(theta)
!       nhz <-- number of variables to be determined as H(z)

      uw=omega * MU_0				! =	  omega*mu0
      ciuw=dcmplx(0.0d0,uw)			! =-i*omega*mu0

!     uw=4.0d-7*pi*2.0d0*pi/dreal(period)
!     ciuw=dcmplx(0.0,uw)


! First, Hx...

      jj=0
      ija(1)=np+2
      ija1(1)=np+2
      kk=np+1
      kk1=np+1

!     -----------
      do 10 i=1,l
!     -----------
        if (i == 1) then
           a=L
        else
           a=i-1
        end if
        if (i == l) then
           d=1
        else
           d=i+1
        end if
!     -------------
        do 10 k=2,n
!     -------------
          c=k-1
          f=k+1
          zm=(z(k)+z(k-1))/2.0d0
          zp=(z(k)+z(k+1))/2.0d0
          call leng_zijk(c,z,zc)
          call leng_zijk(k,z,zk)
!     ---------------
          do 10 j=2,m
!     ---------------
            b=j-1
            e=j+1
            call leng_xijk(i,j,k,x,y,z,xx)
!           ym=(y(j)+y(j-1))/2.d0
!           yp=(y(j)+y(j+1))/2.d0
            ybm=zm*( y(b+1)-y(b) )
            yjm=zm*( y(j+1)-y(j) )
            ybp=zp*( y(b+1)-y(b) )
            yjp=zp*( y(j+1)-y(j) )
!   -----------------------------------------------------
!             set of line elements
!   -----------------------------------------------------
            call area_skij(i,b,k,x,y,z,skib)
            call area_skij(i,j,k,x,y,z,skij)
            call area_sjki(i,j,c,x,y,z,sjci)
            call area_sjki(i,j,k,x,y,z,sjki)
!   -----------------------------------------------------
!             set of surface elements 1
!   -----------------------------------------------------
            call area_sijk2(b,c,y,z,sijk2)
!           sijk2=((zm**2)-(zp**2))*(yp-ym)/2.d0
!           < S(i+1/2)(j-1/2)(k-1/2) >
!   -----------------------------------------------------
!             set of surface element  2
!   -----------------------------------------------------

            rijcy=(resist(i,j,c)*yjm+resist(i,b,c)*ybm)/2.0d0/sjci
            rijky=(resist(i,j,k)*yjp+resist(i,b,k)*ybp)/2.0d0/sjki
            rijkz=(resist(i,j,k)*zk+resist(i,j,c)*zc)/2.0d0/skij
            ribkz=(resist(i,b,k)*zk+resist(i,b,c)*zc)/2.0d0/skib

!----------------------------------------------------------------
! 0 Hxijk
            jj=jj+1
!toh
            dp(jj)=ciuw*sijk2/xx+rijcy+rijky+rijkz+ribkz
            sp1(jj)=dp(jj)
!toh
!----------------------------------------------------------------
! 1 Hxiek
            if (e.ne.m+1) then
!           ----------> there are no Hx at S-pole
              kk=kk+1
              kk1=kk1+1
              call n_hxijk(l,m,n,i,e,k,ii)
              sp(kk)=-rijkz
!toh
              sp1(kk1)=dcmplx(sp(kk),0.d0)
!toh
              ija(kk)=ii
              ija1(kk1)=ii
            end if
! 2 Hxijf
            if(f.ne.n+1) then
!           ----------> bottom Hx is fixed
              kk=kk+1
              kk1=kk1+1
              call n_hxijk(l,m,n,i,j,f,ii)
              sp(kk)=-rijky
!toh
              sp1(kk1)=dcmplx(sp(kk),0.d0)
!toh
              ija(kk)=ii
              ija1(kk1)=ii
            end if
! 3 Hyibk
              kk=kk+1
              call n_hyijk(nhx,l,m,n,i,b,k,ii)
              sp(kk)=ribkz
              ija(kk)=ii
! 4 Hydbk
              kk=kk+1
              call n_hyijk(nhx,l,m,n,d,b,k,ii)
              sp(kk)=-ribkz
              ija(kk)=ii
! 5 Hyijk
              kk=kk+1
              call n_hyijk(nhx,l,m,n,i,j,k,ii)
              sp(kk)=-rijkz
              ija(kk)=ii
! 6 Hydjk
              kk=kk+1
              call n_hyijk(nhx,l,m,n,d,j,k,ii)
              sp(kk)=rijkz
              ija(kk)=ii
! 7 Hzijc
            if (c.ne.1) then
!           ----------> top Hz is fixed
              kk=kk+1
              call n_hzijk(nhx,nhy,l,m,n,i,j,c,ii)
              sp(kk)=rijcy
              ija(kk)=ii
            end if
! 8 Hzdjc
            if (c.ne.1) then
!           ----------> top Hz is fixed
              kk=kk+1
              call n_hzijk(nhx,nhy,l,m,n,d,j,c,ii)
              sp(kk)=-rijcy
              ija(kk)=ii
            end if
! Hzijk
              kk=kk+1
              call n_hzijk(nhx,nhy,l,m,n,i,j,k,ii)
              sp(kk)=-rijky
              ija(kk)=ii
! Hzdjk
              kk=kk+1
              call n_hzijk(nhx,nhy,l,m,n,d,j,k,ii)
              sp(kk)=rijky
              ija(kk)=ii

            ija(jj+1)=kk+1
            ija1(jj+1)=kk1+1
10    continue


! Now for Hy...

!     -----------
      do 20 j=1,m
!     -----------
        b=j-1
        e=j+1
        yp=(y(j)+y(j+1))/2.0d0
!     -------------
        do 20 k=2,n
!     -------------
          c=k-1
          f=k+1
          call leng_yijk(j,k,y,z,yy)
          zm=(z(k)+z(k-1))/2.0d0
          zp=(z(k)+z(k+1))/2.0d0
          call leng_zijk(c,z,zc)
          call leng_zijk(k,z,zk)
          call area_sijk(j,c,y,z,sijc)
          call area_sijk(j,k,y,z,sijk)
!     ---------------
          do 20 i=1,l
!     ---------------
            if (i == 1) then
               a=l
            else
               a=i-1
            end if
            if (i == l) then
               d=1
            else
               d=i+1
            end if
!           xm=x(a)/2.0d0
!           xp=x(j)/2.0d0
            xipm=zm*dsin(yp)*x(i)
            xapm=zm*dsin(yp)*x(a)
            xipp=zp*dsin(yp)*x(i)
            xapp=zp*dsin(yp)*x(a)
!   -----------------------------------------------------
!             set of line elements
!   -----------------------------------------------------
            call area_skij(a,j,k,x,y,z,skaj)
            call area_skij(i,j,k,x,y,z,skij)
!   -----------------------------------------------------
!             set of surface elements 1
!   -----------------------------------------------------
            call area_sjki2(a,i,j,c,x,y,z,sjki2)
!           sjki2=((zm**2)-(zp**2))*sin(yp)*(xp+xm)/2.d0
!           < S(j+1/2)(k-1/2)(i-1/2) >
!   -----------------------------------------------------
!             set of surface elements 2
!   -----------------------------------------------------

            rijkz=(resist(i,j,k)*zk+resist(i,j,c)*zc)/2.0d0/skij
            rijcx=(resist(i,j,c)*xipm+resist(a,j,c)*xapm)/2.0d0/sijc
            rajkz=(resist(a,j,k)*zk+resist(a,j,c)*zc)/2.0d0/skaj
            rijkx=(resist(i,j,k)*xipp+resist(a,j,k)*xapp)/2.0d0/sijk

! 0 Hyijk
            jj=jj+1
!toh
            dp(jj)=ciuw*sjki2/yy+rijkx+rijcx+rajkz+rijkz
            sp1(jj)=dp(jj)
!toh
! 1 Hydjk
           if (d.ne.1) then
            kk=kk+1
            kk1=kk1+1
            call n_hyijk(nhx,l,m,n,d,j,k,ii)
            sp(kk)=-rijkz
!toh
              sp1(kk1)=dcmplx(sp(kk),0.d0)
!toh
            ija(kk)=ii
            ija1(kk1)=ii
           end if
! 1.5 Hyajk
           if (a == l) then
            kk=kk+1
            kk1=kk1+1
            call n_hyijk(nhx,l,m,n,a,j,k,ii)
            sp(kk)=-rajkz
!toh
              sp1(kk1)=dcmplx(sp(kk),0.d0)
!toh
            ija(kk)=ii
            ija1(kk1)=ii
           end if
! 2 Hyijf
            if(f.ne.n+1) then
!           ----------> bottom Hy is fixed
              kk=kk+1
              kk1=kk1+1
              call n_hyijk(nhx,l,m,n,i,j,f,ii)
              sp(kk)=-rijkx
!toh
              sp1(kk1)=dcmplx(sp(kk),0.d0)
!toh
              ija(kk)=ii
              ija1(kk1)=ii
            end if
! 3 Hzijc
            if (c.ne.1) then
!           ----------> top Hz is fixed
              kk=kk+1
              call n_hzijk(nhx,nhy,l,m,n,i,j,c,ii)
              sp(kk)=rijcx
              ija(kk)=ii
            end if
! 4 Hziec
            if (c.ne.1) then
!           ----------> top Hz is fixed
              kk=kk+1
              call n_hzijk(nhx,nhy,l,m,n,i,e,c,ii)
              sp(kk)=-rijcx
              ija(kk)=ii
            end if
! 5 Hzijk
            kk=kk+1
            call n_hzijk(nhx,nhy,l,m,n,i,j,k,ii)
            sp(kk)=-rijkx
            ija(kk)=ii
! 6 Hziek
            kk=kk+1
            call n_hzijk(nhx,nhy,l,m,n,i,e,k,ii)
            sp(kk)=rijkx
            ija(kk)=ii

            ija(jj+1)=kk+1
            ija1(jj+1)=kk1+1
20    continue


! Now for Hz...
!
!     -------------
      do 30 k=2,n
!     -------------
        c=k-1
        f=k+1
        call leng_zijk(k,z,zz)
        zm=(z(k)+z(k-1))/2.0d0
        zp=(z(k)+z(k+1))/2.0d0

!         ++++++++++++++++
          j=1
!         < N-pole >
!         ++++++++++++++++
          e=j+1
            yp=(y(j)+y(j+1))/2.0d0
            call area_sijk(j,k,y,z,sijk)
            skij2=2.0d0*pi*(zp**2)*(1.d0-dcos(yp))
!           < N-pole cap >
! 0 Hzijk
            jj=jj+1
!toh
            dp(jj)=ciuw*skij2/zz
!toh
!           -----------
            do i=1,l
!           -----------
              if (i == 1) then
                a=l
              else
                a=i-1
              end if
              if (i == l) then
                d=1
              else
                d=i+1
              end if
              xipp=zp*dsin(yp)*x(i)
              xapp=zp*dsin(yp)*x(a)
              rijkx=(resist(i,j,k)*xipp+resist(a,j,k)*xapp)/2.0d0/sijk
!toh
              dp(jj)=dp(jj)+rijkx
!toh
! 2 Hziek
              kk=kk+1
              kk1=kk1+1
              call n_hzijk(nhx,nhy,l,m,n,i,e,k,ii)
              sp(kk)=-rijkx
!toh
              sp1(kk1)=dcmplx(sp(kk),0.d0)
!toh
              ija(kk)=ii
              ija1(kk1)=ii

!           -----------
            end do
!           -----------
!toh
            sp1(jj)=dp(jj)
!toh
            ija(jj+1)=kk+1
            ija1(jj+1)=kk1+1

!     -------------
        do j=2,m
!     -------------
!         ++++++++++++++++
!         < neither N-pole nor S-pole >
!         ++++++++++++++++
          b=j-1
          e=j+1
            yp=(y(j)+y(j+1))/2.0d0
            ym=(y(j)+y(j-1))/2.0d0
            ybp=zp*( y(b+1)-y(b) )
            yjp=zp*( y(j+1)-y(j) )
            call area_sijk(j,k,y,z,sijk)
            call area_sijk(b,k,y,z,sibk)
!     ---------------
          do i=1,l
!     ---------------

            if (i == 1) then
               a=L
            else
               a=i-1
            end if
            if (i == l) then
               d=1
            else
               d=i+1
            end if
            ximp=zp*dsin(ym)*x(i)
            xamp=zp*dsin(ym)*x(a)
            xipp=zp*dsin(yp)*x(i)
            xapp=zp*dsin(yp)*x(a)
            xm=x(a)/2.0d0
            xp=x(i)/2.0d0
!   -----------------------------------------------------
!             set of line elements
!   -----------------------------------------------------
            call area_sjki(i,j,k,x,y,z,sjki)
            call area_sjki(a,j,k,x,y,z,sjka)
!   -----------------------------------------------------
!             set of surface elements 1
!   -----------------------------------------------------
            call area_skij2(a,i,b,k,x,y,z,skij2)
!           skij2=(zp**2)*(cos(ym)-cos(yp))*(xm+xp)
!           < S(k+1/2)(i-1/2)(j-1/2) >
!   -----------------------------------------------------
!             set of surface elements 2
!   -----------------------------------------------------

            rijky=(resist(i,j,k)*yjp+resist(i,b,k)*ybp)/2.0d0/sjki
            rijkx=(resist(i,j,k)*xipp+resist(a,j,k)*xapp)/2.0d0/sijk
            rajky=(resist(a,j,k)*yjp+resist(a,b,k)*ybp)/2.0d0/sjka
            ribkx=(resist(i,b,k)*ximp+resist(a,b,k)*xamp)/2.0d0/sibk

! 0 Hzijk
            jj=jj+1
!toh
            dp(jj)=ciuw*skij2/zz+rijky+rajky+ribkx+rijkx
            sp1(jj)=dp(jj)
!toh
! 1 Hzdjk
           if (d.ne.1) then
            kk=kk+1
            kk1=kk1+1
            call n_hzijk(nhx,nhy,l,m,n,d,j,k,ii)
            sp(kk)=-rijky
!toh
              sp1(kk1)=dcmplx(sp(kk),0.d0)
!toh
            ija(kk)=ii
            ija1(kk1)=ii
           end if
! 1.5 Hzajk
           if (a == l) then
            kk=kk+1
            kk1=kk1+1
            call n_hzijk(nhx,nhy,l,m,n,a,j,k,ii)
            sp(kk)=-rajky
!toh
              sp1(kk1)=dcmplx(sp(kk),0.d0)
!toh
            ija(kk)=ii
            ija1(kk1)=ii
           end if
! 2 Hziek
            kk=kk+1
            kk1=kk1+1
            call n_hzijk(nhx,nhy,l,m,n,i,e,k,ii)
            sp(kk)=-rijkx
!toh
              sp1(kk1)=dcmplx(sp(kk),0.d0)
!toh
            ija(kk)=ii
            ija1(kk1)=ii

            ija(jj+1)=kk+1
            ija1(jj+1)=kk1+1
!           -----------
            end do
!         -------------
          end do
!         ------
!         +++++++++++++++++++++++
          j=m+1
!         < S-pole >
!         +++++++++++++++++++++++
          b=j-1
            ym=(y(j)+y(j-1))/2.0d0
            call area_sijk(b,k,y,z,sibk)
            skij2=2.0d0*pi*(zp**2)*(dcos(ym)+1.0d0)
!           < S-pole cap >
! 0 Hzijk
            jj=jj+1
!toh
            dp(jj)=ciuw*skij2/zz
!toh
!           -----------
            do i=1,l
!           -----------
              if (i  == 1) then
                a=l
              else
                a=i-1
              end if
              if (i == l) then
                d=1
              else
                d=i+1
              end if
              ximp=zp*dsin(ym)*x(i)
              xamp=zp*dsin(ym)*x(a)
              ribkx=(resist(i,b,k)*ximp+resist(a,b,k)*xamp)/2.0d0/sibk
!toh
              dp(jj)=dp(jj)+ribkx
!toh
!           -----------
            end do
!           -----------
!toh
            sp1(jj)=dp(jj)
!toh
            ija(jj+1)=kk+1
            ija1(jj+1)=kk1+1

 30   continue

      ! AK: Compute the maximum "row length" = max(ija(jj+1) - ija(jj))
      ijalen = jj+1
      rowlen = 0
      do jj=1,ijalen-1
            if(ija(jj+1)-ija(jj) > rowlen) then
                rowlen = ija(jj+1)-ija(jj)
            end if
      end do

      return
      end subroutine set_op3u


! *****************************************************************************


      subroutine sum_divBC3(l,m,n,ar,ai,z,r)
!-------------------------------------------------------------
! subroutine to compute a summation of div
!     normalized by volume concerned
!-------------------------------------------------------------

	  integer					 :: l,m,n
      real(8),dimension(np1)     :: ar,ai
      real(8),dimension(n+1)     :: z
      real(8)                    :: r,vol
      integer                    :: i,ii,j,k

      r=0.0d0
      do k=2,n
!toh:29/FEB/2000
        j=1
            call n_dvpotijkBC3(l,m,n,1,j,k,ii)
            r=r+dsqrt( (ar(ii)**2)+(ai(ii)**2) )
!       end j=1
        do j=2,m
          do i=1,l
            call n_dvpotijkBC3(l,m,n,i,j,k,ii)
            r=r+dsqrt( (ar(ii)**2)+(ai(ii)**2) )
          end do
        end do
        j=m+1
            call n_dvpotijkBC3(l,m,n,1,j,k,ii)
            r=r+dsqrt( (ar(ii)**2)+(ai(ii)**2) )
!       end j=M+1
      end do

      call volume_sph(1,n+1,z,n+1,vol)
!toh:29/FEB/2000
      r=r/vol

      return
      end subroutine sum_divBC3



! *****************************************************************************




      subroutine updateh3BC3(l,m,n,pot,hx,hy,hz,x,y,z,job)
!-------------------------------------------------------------
! subroutine to make a correction to H: H=H0+grad(pot)
! the potential comes from requiring divH=0
!
!     job=0 update real part
!     job=1 update imaginary part
!
! modified to assume boundary values at top of air layers
! are fixed to their starting values
!
! modified for real space values
!toh:29/FEB/2000
!-------------------------------------------------------------


      integer                          :: i,j,k,l,m,n,job
      real(8),dimension(l)	           :: x
      real(8),dimension(m+1)           :: y
      real(8),dimension(n+1)           :: z
      real(8),dimension(np1)           :: pot
      real(8)                          :: grad
      complex(8),dimension(np1)        :: hx,hy,hz
      complex(8)                       :: grad1

      integer                          :: a,d
      integer                          :: nijf,nijk,niek,ndjk,ic
      real(8)                          :: xx,yy,zz

! update Hx

      do i=1,l
        if (i.ne.1) then
           a=i-1
        else
           a=l
        end if
        if (i.ne.l) then
           d=i+1
        else
           d=1
        end if
        do k=2,n
!toh:29/FEB/2000
          do j=2,m
            call leng_xijk(i,j,k,x,y,z,xx)
            call n_allpotijk(l,m,n,i,j,k,nijk)
            call n_allpotijk(l,m,n,d,j,k,ndjk)
            grad=pot(ndjk)-pot(nijk)
            grad=grad/xx
            if (job == 0) then
              grad1=dcmplx(grad,0.0d0)
            else
              grad1=dcmplx(0.0d0,grad)
            end if
            call n_allhxijk(l,m,n,i,j,k,ic)
            hx(ic)=hx(ic)+grad1
          end do
        end do
      end do

! update Hy

      do j=1,m
        do k=2,n
!toh:29/FEB/2000
          call leng_yijk(j,k,y,z,yy)
          do i=1,l
            call n_allpotijk(l,m,n,i,j,k,nijk)
            call n_allpotijk(l,m,n,i,j+1,k,niek)
            grad=pot(niek)-pot(nijk)
            grad=grad/yy
            if (job == 0) then
              grad1=dcmplx(grad,0.0d0)
            else
              grad1=dcmplx(0.0d0,grad)
            end if
            call n_allhyijk(l,m,n,i,j,k,ic)
            hy(ic)=hy(ic)+grad1
          end do
        end do
      end do

! update Hz

      do k=1,n
!toh:29/FEB/2000
        call leng_zijk(k,z,zz)
        j=1
            call n_allpotijk(l,m,n,0,j,k,nijk)
            call n_allpotijk(l,m,n,0,j,k+1,nijf)
            grad=pot(nijf)-pot(nijk)
            grad=grad/zz
            if (job == 0) then
              grad1=dcmplx(grad,0.0d0)
            else
              grad1=dcmplx(0.0d0,grad)
            end if
            call n_allhzijk(l,m,n,0,j,k,ic)
            hz(ic)=hz(ic)+grad1
        do j=2,m
          do i=1,l
            call n_allpotijk(l,m,n,i,j,k,nijk)
            call n_allpotijk(l,m,n,i,j,k+1,nijf)
            grad=pot(nijf)-pot(nijk)
            grad=grad/zz
            if (job == 0) then
              grad1=dcmplx(grad,0.0d0)
            else
              grad1=dcmplx(0.0d0,grad)
            end if
            call n_allhzijk(l,m,n,i,j,k,ic)
            hz(ic)=hz(ic)+grad1
          end do
        end do
        j=m+1
            call n_allpotijk(l,m,n,0,j,k,nijk)
            call n_allpotijk(l,m,n,0,j,k+1,nijf)
            grad=pot(nijf)-pot(nijk)
            grad=grad/zz
            if (job == 0) then
              grad1=dcmplx(grad,0.0d0)
            else
              grad1=dcmplx(0.0d0,grad)
            end if
            call n_allhzijk(l,m,n,0,j,k,ic)
            hz(ic)=hz(ic)+grad1
      end do

      return
      end subroutine updateh3BC3

end module coreFwd
