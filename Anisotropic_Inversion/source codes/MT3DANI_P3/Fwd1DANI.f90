module Fwd1DANI

  use math_constants
  use utilities

  type :: rscalar_1d


     ! Typical usage:  conductivity averaged on centers of
     ! staggered grid
     ! v: dimension Nx, Ny, Nz
     ! Note that the arrays are defined through dynamic memory allocation
     real(kind=prec), pointer, dimension(:)        :: v

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction:
     integer                                            ::  nz =0

     ! allocated:  .true.  v array has been allocated
     logical		                                :: allocated = .false.

     ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
     ! (probably will not be needed in the future when compilers will support
     ! ISO/IEC 15581 - the "allocatable array extension")
     logical										:: temporary = .false.

  end type rscalar_1d

  ! ***************************************************************************
  ! type rscalar defines scalar for either edge or face in a staggered grid as
  ! a real field
  type :: rscalar_2d


     ! Typical usage:  conductivity averaged on centers of
     ! staggered grid
     ! v: dimension Nx, Ny, Nz
     ! Note that the arrays are defined through dynamic memory allocation
     real(kind=prec), pointer, dimension(:,:)        :: v

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction:
     integer                                            :: ny = 0, nz =0

     ! allocated:  .true.  v array has been allocated
     logical		                                :: allocated = .false.

     ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
     ! (probably will not be needed in the future when compilers will support
     ! ISO/IEC 15581 - the "allocatable array extension")
     logical										:: temporary = .false.

  end type rscalar_2d

  type ::  grid2d_t
      integer	:: Nz,Ny,Nza
      real(kind=prec), pointer, dimension(:) :: Dy,Dz
  end type


Contains

  ! **********set boundary condition for 2d forward problem**********
  subroutine boufld(m,n,nlayr,nb,y,z,ihpol,period,cond,ex,ey,hx,hy,hz)
    implicit none
    integer m,n,nlayr,nb,ihpol
	real(kind=prec)    :: period,y(*),z(*)
    type(rscalar_2d)   :: cond(6)
	complex(kind=prec) :: ex(n,m),ey(n,m)
    complex(kind=prec) :: hx(n,m),hy(n,m),hz(n,m) 
    ! local variables
	real(kind=prec)    :: sy(m),sz(n)
    integer nlbl,nlbr,i,j,k
    real(kind=prec)    :: hbl(n),hbr(n),dire
	type(rscalar_1d)   :: sgbl(6),sgbr(6)
    complex(kind=prec) :: hx0,hy0
    complex(kind=prec) :: ex1d(n),ey1d(n),hx1d(n),hy1d(n)

	sy(1)=0.
    do j=2,m
      sy(j)=sy(j-1)+y(j-1)
    end do
!   in the air
    sz(nb)=0.
    do k=nb-1,1,-1
      sz(k)=sz(k+1)-z(k)
    end do
!   in the earth
    do k=nb+1,n
      sz(k)=sz(k-1)+z(k-1) 
    end do
    
	if(ihpol.eq.1) then
      hx0=C_ONE
      hy0=C_ZERO	
	elseif(ihpol.eq.2) then
      hx0=C_ZERO
      hy0=C_ONE
	else
	  call errStop('Polarization mode does not exist!')	
	end if 

    ! Construct layered models on the margins
	call boumedcell(m,n,nlayr,nb,z,cond,nlbl,nlbr,hbl,hbr,sgbl,sgbr)

    ! left boundary
    call em1abc(nlbl,hbl,sgbl,period,n,nb,sz,hx0,hy0,  &
	                                 ex1d,ey1d,hx1d,hy1d)
    do k=1,n
      ex(k,1)=ex1d(k)
      ey(k,1)=ey1d(k)
      hx(k,1)=hx1d(k)
	    hy(k,1)=hy1d(k)
      hz(k,1)=0.d0
	end do
	! right boundary
	call em1abc(nlbr,hbr,sgbr,period,n,nb,sz,hx0,hy0,  &
                                    ex1d,ey1d,hx1d,hy1d)
    do k=1,n
      ex(k,m)=ex1d(k)
      ey(k,m)=ey1d(k)
      hx(k,m)=hx1d(k)
	    hy(k,m)=hy1d(k)
      hz(k,m)=0.d0
    end do

    do k=1,n
      do j=2,m-1
        dire=sy(j)/sy(m)
        !
        ex(k,j)=ex(k,1)+dire*(ex(k,m)-ex(k,1))
        ey(k,j)=ey(k,1)+dire*(ey(k,m)-ey(k,1))
      enddo
    enddo

!   do j=1,m
!     do k=1,n
!       write(103,3000) ihpol,j,k,hy(k,j)
!       write(104,3000) ihpol,j,k,hz(k,j)
!     end do
!   end do
!
! 3000 FORMAT(I3,1X,I3,1X,I3,3X,2E11.3)

  end subroutine boufld

  ! *************Construct layered models on the margins***************
  subroutine boumedcell(m,n,nlayr,nb,z,cond,nlbl,nlbr,hbl,hbr,sgbl,sgbr)
    implicit none
    integer m,n,nlayr,nb,nlbl,nlbr
    real(kind=prec)     :: z(*)
    type(rscalar_2d)	:: cond(6)
    real(kind=prec)     :: hbl(n),hbr(n)
	type(rscalar_1d)	:: sgbl(6),sgbr(6)
    ! local variables
	integer ia,k  
    !
	do ia = 1,6

      call create_rscalar_1d(n, sgbl(ia))
      nlbl=1
	  do k = nb, n-1
	     hbl(nlbl) = z(k)
         sgbl(ia)%v(nlbl) = cond(ia)%v(k,1)
		 nlbl = nlbl + 1
      end do
	  hbl(nlbl)=0.0 
      sgbl(ia)%v(nlbl) = cond(ia)%v(n-1,1)

      call create_rscalar_1d(n, sgbr(ia))
      nlbr=1
	  do k = nb, n-1
	     hbr(nlbr) = z(k)
         sgbr(ia)%v(nlbr) = cond(ia)%v(k,m-1)
		 nlbr = nlbr + 1
      end do
	  hbr(nlbr)=0.0 
      sgbr(ia)%v(nlbr) = cond(ia)%v(n-1,m-1)

	end do

!   output for check
!    do ia=1,6
!	  write(101,*) nlbl
!      do k=1,nlbl
!	    write(101,*) hbl(k),sgbl(ia)%v(k)
!	  end do
!	  write(102,*) nlbr
!      do k=1,nlbr
!	    write(102,*) hbr(k),sgbr(ia)%v(k)
!	  end do
!	end do

  end subroutine boumedcell



  ! *******************1D boundary*******************
  subroutine em1abc(nl,h,sig,t,nzs,nzsd,zs,hx0,hy0,  & 
                                     ex1d,ey1d,hx1d,hy1d)
    implicit none
	integer nl,nzs,nzsd
    type(rscalar_1d)	:: sig(6)
	real(kind=prec)     :: t,h(nzs),zs(nzs)
    complex(kind=prec)  :: ic,hx0,hy0
    complex(kind=prec)  :: ex1d(nzs),ey1d(nzs),hx1d(nzs),hy1d(nzs)
	! local variables
    integer i,ii
	complex(kind=prec)  :: zpom(2,2,nzs),mgf(2,nzs),elf(2,nzs)
	real(kind=prec)     :: smax(nzs),smin(nzs),sazi(nzs),ommi0
    !
    call azim1damod(nzs,sig,smax,smin,sazi,nl)
    call z1adu(nzs,t,h,smax,smin,sazi,zpom,nl)
    mgf(1,1)=hx0
	mgf(2,1)=hy0
    call h1aud(nzs,t,h,smax,smin,sazi,zpom,mgf,elf,nl)

    ic=cmplx(0.0d0,1.0d0)
	ommi0=8.0d-7*PI*PI

	do i=nzsd,nzs
  	  ii=i-nzsd+1
	  hx1d(i)=mgf(1,ii)
	  hy1d(i)=mgf(2,ii)
	  ex1d(i)=elf(1,ii)
	  ey1d(i)=elf(2,ii)
    end do

	do i=1,nzsd
	  hx1d(i)=mgf(1,1)
	  hy1d(i)=mgf(2,1)
	  ex1d(i)=ic*ommi0*zs(i)*mgf(2,1)/t+elf(1,1)
	  ey1d(i)=-ic*ommi0*zs(i)*mgf(1,1)/t+elf(2,1)
    end do

  end subroutine em1abc

 
  !***********************************************
  subroutine azim1damod(nzs,sig,smax,smin,sazi,nl) 
    implicit none
	integer nl,nzs
    type(rscalar_1d)	:: sig(6) 
    real(kind=prec)     :: smax(nzs),smin(nzs),sazi(nzs)
    ! local variables
	integer l,ia
	real(kind=prec)     :: axx,axy,ayy
	real(kind=prec)     :: da0,da,sa0,db0,doubletiny
	!
	do l=1,nl
	  ! Axx = sigxx - sigxz * sigzx / sigzz  
      axx=sig(1)%v(l)-sig(5)%v(l)*sig(5)%v(l)/sig(3)%v(l)
	  ! Axy = sigxy - sigxz * sigzy / sigzz 
      axy=sig(4)%v(l)-sig(5)%v(l)*sig(6)%v(l)/sig(3)%v(l)
	  ! Ayy = sigyy - sigyz * sigzy / sigzz
	  ayy=sig(2)%v(l)-sig(6)%v(l)*sig(6)%v(l)/sig(3)%v(l)

	  sa0=axx+ayy
	  da0=axx-ayy
	  da=dsqrt(da0*da0+4.d0*axy*axy)
	  smax(l)=0.5d0*(sa0+da)
	  smin(l)=0.5d0*(sa0-da)
	  if(da.lt.tiny(doubletiny))then
	    sazi(l)=0.d0
	  else
	    db0=da0/da
	    sazi(l)=0.5d0*dacos(db0)
	    if(axy.lt.0.d0)sazi(l)=-sazi(l)
	  endif
	enddo
	do ia = 1,6
	  call deall_rscalar_1d(sig(ia))
	enddo
	
  end subroutine azim1damod


  !***********************************************
  subroutine z1adu(nzs,per,h,smax,smin,sazi,z,nl)
    implicit none
	integer nzs,nl
	real(kind=prec)     :: h(nzs),per
	real(kind=prec)     :: smax(nzs),smin(nzs),sazi(nzs)
	complex(kind=prec)  :: z(2,2,nzs)
	! local variables
	integer l,l1
	real(kind=prec)     :: omega,s2,c2,dsazi
	complex(kind=prec)  :: ic,kmax,kmin,commi,cimmo
	complex(kind=prec)  :: kapmax,kapmin
	complex(kind=prec)  :: sz1,sz2,dz1,dz2,dz,dz0,dsec,thmin,thmax
	!
    complex(kind=prec)  :: dth,cx
	dth(cx)=(1.d0-cdexp(-2.d0*cx))/(1.d0+cdexp(-2.d0*cx))
    !
	ic=cmplx(0.0d0,1.0d0)
    omega=2.d0*PI/per
    cimmo=ic*omega*MU_0
	commi=(1.d0-ic)*dsqrt(0.5d0*omega*MU_0)
	!
    do l=nl,1,-1
      kmax=commi*dsqrt(smax(l))
      kmin=commi*dsqrt(smin(l))
	  kapmax=kmax/cimmo
	  kapmin=kmin/cimmo
      if(l.eq.nl)then
        z(1,1,l)=0.d0
        z(1,2,l)=kmax/smax(l)
        z(2,1,l)=-kmin/smin(l)
        z(2,2,l)=0.d0
      else
        l1=l+1
        if(sazi(l).ne.sazi(l1))then
          dsazi=sazi(l)-sazi(l1)
          c2=dcos(2.d0*dsazi)
          s2=dsin(2.d0*dsazi)
          sz1=z(1,1,l1)+z(2,2,l1)
          sz2=z(1,2,l1)+z(2,1,l1)
          dz1=z(1,1,l1)-z(2,2,l1)
          dz2=z(1,2,l1)-z(2,1,l1)
          z(1,1,l1)=0.5d0*(sz1+dz1*c2+sz2*s2)
          z(1,2,l1)=0.5d0*(dz2+sz2*c2-dz1*s2)
          z(2,1,l1)=0.5d0*(-dz2+sz2*c2-dz1*s2)
          z(2,2,l1)=0.5d0*(sz1-dz1*c2-sz2*s2)
        endif
        thmax=dth(kmax*h(l))
        thmin=dth(kmin*h(l))
        dz=z(1,1,l1)*z(2,2,l1)-z(1,2,l1)*z(2,1,l1)
        dz0=1.d0-kapmax*thmax*z(1,2,l1)+kapmin*thmin*z(2,1,l1)+  &
                 kapmax*kapmin*thmax*thmin*dz
        dsec=4.d0*cdexp(-(kmax+kmin)*h(l))/   & 
             ((1.d0+cdexp(-2.d0*kmax*h(l)))*  & 
             (1.d0+cdexp(-2.d0*kmin*h(l))))
        z(1,1,l)=z(1,1,l1)*dsec/dz0
        z(1,2,l)=(z(1,2,l1)-thmax/kapmax-kapmin*thmin*dz-  & 
                 (kapmin/kapmax)*thmax*thmin*z(2,1,l1))/dz0
        z(2,1,l)=(z(2,1,l1)+thmin/kapmin+kapmax*thmax*dz-  & 
                 (kapmax/kapmin)*thmax*thmin*z(1,2,l1))/dz0
        z(2,2,l)=z(2,2,l1)*dsec/dz0
      endif
    enddo
    if(sazi(1).ne.0.d0)then
      dsazi=-sazi(1)
      c2=dcos(2.d0*dsazi)
      s2=dsin(2.d0*dsazi)
      sz1=z(1,1,1)+z(2,2,1)
      sz2=z(1,2,1)+z(2,1,1)
      dz1=z(1,1,1)-z(2,2,1)
      dz2=z(1,2,1)-z(2,1,1)
      z(1,1,1)=0.5d0*(sz1+dz1*c2+sz2*s2)
      z(1,2,1)=0.5d0*(dz2+sz2*c2-dz1*s2)
      z(2,1,1)=0.5d0*(-dz2+sz2*c2-dz1*s2)
      z(2,2,1)=0.5d0*(sz1-dz1*c2-sz2*s2)
    endif

  end subroutine z1adu


  !***********************************************
  subroutine h1aud(nzs,per,h,smax,smin,sazi,z,mgf,elf,nl)
    implicit none
	integer nzs,nl
	real(kind=prec)     :: h(nzs),per
 	real(kind=prec)     :: smax(nzs),smin(nzs),sazi(nzs) 
	complex(kind=prec)  :: z(2,2,nzs),mgf(2,nzs),elf(2,nzs)
	! local variables
	integer nl1,l,l1
	real(kind=prec)     :: omega,s1,c1,dsazi
	complex(kind=prec)  :: mgfu(2),mgfd(2),mgfpom,elfpom
	complex(kind=prec)  :: ic,kmax,kmin,commi,cimmo
	complex(kind=prec)  :: kapmax,kapmin
	complex(kind=prec)  :: dz,dz0,thmin,thmax,dzmax,dzmin
	!
    complex(kind=prec)  :: dth,cx
	dth(cx)=(1.d0-cdexp(-2.d0*cx))/(1.d0+cdexp(-2.d0*cx))
    !
	ic=cmplx(0.0d0,1.0d0)
    omega=2.d0*PI/per
    cimmo=ic*omega*MU_0
	commi=(1.d0-ic)*dsqrt(0.5d0*omega*MU_0)

	nl1=nl-1
	elf(1,1)=z(1,1,1)*mgf(1,1)+z(1,2,1)*mgf(2,1)
	elf(2,1)=z(2,1,1)*mgf(1,1)+z(2,2,1)*mgf(2,1)
	mgfu(1)=mgf(1,1)
	mgfu(2)=mgf(2,1)

	if(nl1.gt.0)then
      do l=1,nl1
	    l1=l+1
        kmax=commi*dsqrt(smax(l))
        kmin=commi*dsqrt(smin(l))
	    kapmax=kmax/cimmo
	    kapmin=kmin/cimmo

	    if(l.eq.1)then
	      if(sazi(l).ne.0.d0)then
	        c1=dcos(sazi(l))
	        s1=dsin(sazi(l))
	        mgfpom=mgfu(1)*c1+mgfu(2)*s1
	        mgfu(2)=-mgfu(1)*s1+mgfu(2)*c1
	        mgfu(1)=mgfpom
	      endif
	    else
	      mgfu(1)=mgfd(1)
		  mgfu(2)=mgfd(2)
		  if(sazi(l).ne.sazi(l-1))then
		    dsazi=sazi(l)-sazi(l-1)
	        c1=dcos(dsazi)
	        s1=dsin(dsazi)
		    mgfpom=mgfu(1)*c1+mgfu(2)*s1
		    mgfu(2)=-mgfu(1)*s1+mgfu(2)*c1
		    mgfu(1)=mgfpom
	      endif
	    endif

        thmax=dth(kmax*h(l))
        thmin=dth(kmin*h(l))
        dz=z(1,1,l1)*z(2,2,l1)-z(1,2,l1)*z(2,1,l1)
        dz0=1.d0-kapmax*thmax*z(1,2,l1)+kapmin*thmin*z(2,1,l1)+  &
               kapmax*kapmin*thmax*thmin*dz
	    dzmax=2.d0*cdexp(-kmax*h(l))/  &
     		  (dz0*(1.d0+cdexp(-2.d0*kmax*h(l))))
	    dzmin=2.d0*cdexp(-kmin*h(l))/  &
              (dz0*(1.d0+cdexp(-2.d0*kmin*h(l))))
	    mgfd(1)=(1.d0-kapmax*thmax*z(1,2,l1))*dzmin*mgfu(1)-  &
                 kapmin*thmin*dzmax*z(2,2,l1)*mgfu(2)
	    mgfd(2)=kapmax*thmax*z(1,1,l1)*dzmin*mgfu(1)+  &
                 (1.d0+kapmin*thmin*z(2,1,l1))*dzmax*mgfu(2)

	    elf(1,l1)=z(1,1,l1)*mgfd(1)+z(1,2,l1)*mgfd(2)
	    elf(2,l1)=z(2,1,l1)*mgfd(1)+z(2,2,l1)*mgfd(2)

        if(sazi(l).eq.0.d0)then
	      mgf(1,l1)=mgfd(1)
	      mgf(2,l1)=mgfd(2)
	    else
	      c1=dcos(sazi(l))
	      s1=dsin(sazi(l))
	      mgf(1,l1)=mgfd(1)*c1-mgfd(2)*s1
	      mgf(2,l1)=mgfd(1)*s1+mgfd(2)*c1
	      elfpom=elf(1,l1)*c1-elf(2,l1)*s1
	      elf(2,l1)=elf(1,l1)*s1+elf(2,l1)*c1
	      elf(1,l1)=elfpom
	    endif
	  enddo
	endif

  end subroutine h1aud


  !***********************************************
  subroutine create_rscalar_1d(n, E)

    implicit none
	integer n
    type (rscalar_1d), intent(inout)      :: E
    ! local variables
    integer                            :: status,nz


    if(E%allocated) then
       ! first deallocate memory for v
       deallocate(E%v, STAT=status)
    end if

    ! Grid dimensions
    nz = n

    E%nz = nz

    ! allocate memory for v
    ! E%allocated will be true if all allocations succeed

    allocate(E%v(nz), STAT=status)

    E%allocated = status .EQ. 0

    if (E%allocated) then
       E%v = R_ZERO
    end if

  end subroutine create_rscalar_1d


  !***********************************************
  subroutine deall_rscalar_1d(E)

    implicit none
	integer n
    type (rscalar_1d), intent(inout)      :: E
    ! local variables
    integer                            :: status


    ! first deallocate memory for v
    deallocate(E%v, STAT=status)

    E%nz = 0
    E%allocated = .false.

  end subroutine deall_rscalar_1d  

end module Fwd1DANI 
