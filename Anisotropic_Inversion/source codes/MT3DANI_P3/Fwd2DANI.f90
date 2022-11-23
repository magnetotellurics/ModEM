module Fwd2DANI
  
  use Fwd1DANI

Contains

  !*************************************************************************************
  subroutine d2atbp(np2_2d,np3_2d,m,n,nb,y,z,nlayr,period,ihpol,cond,ex,ey,ez,strike)
    implicit none
	integer np2_2d,np3_2d,m,n,nb,nlayr,ihpol
	real(kind=prec) :: period,y(*),z(*)
	! Array for 2D Conductivity
    type(rscalar_2d)	:: cond(6)
	complex(kind=prec) :: hx(n,m),hy(n,m),hz(n,m) 
	character*1 strike
    ! local variables
	integer j,k,node_e(m,n),node_h(m,n)
	integer neh,meh
    complex(kind=prec) :: ex(n,m),ey(n,m),ez(n,m)
    complex(kind=prec),allocatable :: apom(:), bpom(:)

! debug
!    do j=1,m-1
!	  write(1010,*) y(j)
!	end do
!    do k=1,n-1
!	  write(1010,*) z(k)
!	end do

    ! construct nodes number
	call classnodes(m,n,nb,node_e,node_h)

    neh=(m-2)*(2*n-nb-3)
    meh=2*n-nb

    ! set boundary condition for 2d forward problem
    call boufld(m,n,nlayr,nb,y,z,ihpol,period,cond,ex,ey,hx,hy,hz)
    do j = 1,m
      do k = 1,n
        ez(k,j) = C_ZERO
      enddo
    enddo

!    allocate(apom(np2_2d*np3_2d),bpom(np2_2d))
!
!!    call czero(np2_2d*np3_2d,apom)
!!	call czero(np2_2d,bpom)
!
!	if(strike.eq.'x') then
!
!      ! construct the coefficient matrix and r.h.s for X-strike problem
!	  call koef3_x(neh,meh,m,n,nlayr,nb,y,z,ihpol,period,  & 
!         cond,node_e,node_h,ex,hx,apom,bpom)
!
!	  WRITE(*,3000)PERIOD,IHPOL
!      
!	  call gaussr(m,n,nb,neh,meh,node_h,apom,bpom)
!    
!	  call getfld(m,n,nb,neh,ihpol,node_e,node_h,bpom,ex,hx)
!      ! compute hy and hz
!      call gethyhz(m,n,nlayr,ihpol,y,z,period,cond,ex,hx,hy,hz)
!                    
!
!	elseif(strike.eq.'y') then
!
!      ! construct the coefficient matrix and r.h.s for Y-strike problem
!      call koef3_y(neh,meh,m,n,nlayr,nb,y,z,ihpol,period,  & 
!             cond,node_e,node_h,ey,hy,apom,bpom)
!
!	  WRITE(*,3000)PERIOD,IHPOL
!      
!	  call gaussr(m,n,nb,neh,meh,node_h,apom,bpom)
!
!	  call getfld(m,n,nb,neh,ihpol,node_e,node_h,bpom,ey,hy)
!      ! compute hx and hz
!      call gethxhz(m,n,nlayr,ihpol,y,z,period,cond,ey,hy,hx,hz)                  
!
!    else
!      stop 'strike in d2atbp must be x or y'
!	end if
!
!3000 FORMAT(3X,'PERIOD',F12.4,' => PASS',I2,' OF DIRECT SOLUTION')
!
!    deallocate(apom,bpom)

  end subroutine d2atbp


  ! construct nodes number
  subroutine classnodes(m,n,nb,node_e,node_h)

	implicit none
	integer j,k,jk,m,n,nb
	integer node_e(m,n),node_h(m,n)

    jk=0
	do j=1,m
        do k=1,n
	    if((j.eq.1).or.(j.eq.m))then   !left and right boundary
	      node_e(j,k)=-2
	      if(k.lt.nb)then
	        node_h(j,k)=0
	      elseif(k.eq.nb)then
	        node_h(j,k)=-1
	      else
	        node_h(j,k)=-2
	      endif
	    else
	      if(k.eq.1)then
	        node_e(j,k)=-2	        
	        node_h(j,k)=0
	      elseif(k.lt.nb)then
	        jk=jk+1
	        node_e(j,k)=jk    
	        node_h(j,k)=0
	      elseif(k.eq.nb)then
	        jk=jk+1
	        node_e(j,k)=jk    
	        node_h(j,k)=-1
          elseif(k.lt.n)then
	        jk=jk+1
	        node_e(j,k)=jk
	        jk=jk+1    
	        node_h(j,k)=jk
		  else
	        node_e(j,k)=-2	        
	        node_h(j,k)=-2
		  endif
		endif
	  end do
	end do

  end subroutine classnodes


  ! construct the coefficient matrix and r.h.s for X-strike problem
  subroutine koef3_x(neh,meh,m,n,nlayr,nb,y,z,ihpol,period,  &
                   cond,node_e,node_h,ex,hx,apom,bpom)
    implicit none
    integer neh,meh,m,n,nlayr,nb,ihpol
	integer node_e(m,n),node_h(m,n)
    real(kind=prec)     :: period,y(*),z(*)
    type(rscalar_2d)	:: cond(6)
    complex(kind=prec)  :: ex(n,m),hx(n,m)
    complex(kind=prec)  :: apom(neh,meh), bpom(neh)
    ! local variables
	integer j,k,jk,jnode,knode,nodac,l,nodae,nodah
	integer md1max,md2max,mh1max,mh2max,mhmax,meh1max,meh3max
    real(kind=prec)     :: ommi0,ommi,dd,ommii,pe,qe,re,se,te,ve
	real(kind=prec)     :: slyem,slyep,slye,sryem,sryep,srye
	real(kind=prec)     :: suzem,suzep,suze,sdzem,sdzep,sdze,sce
	real(kind=prec)     :: pmh,pch,ph,pph,qch,qh,rch,rh
	real(kind=prec)     :: smh,sch,sh,sph,tch,th,vh
	real(kind=prec)     :: hym,hyp,hzm,hzp
	real(kind=prec)     :: sigxxmm,sigyymm,sigzzmm,sigxymm,sigxzmm,sigyzmm
	real(kind=prec)     :: sigxxpm,sigyypm,sigzzpm,sigxypm,sigxzpm,sigyzpm
	real(kind=prec)     :: sigxxmp,sigyymp,sigzzmp,sigxymp,sigxzmp,sigyzmp
	real(kind=prec)     :: sigxxpp,sigyypp,sigzzpp,sigxypp,sigxzpp,sigyzpp
	real(kind=prec)     :: sbymm,sbzmm,scymm,sczmm,saymm,sazmm
	real(kind=prec)     :: sbypm,sbzpm,scypm,sczpm,saypm,sazpm
	real(kind=prec)     :: sbymp,sbzmp,scymp,sczmp,saymp,sazmp
	real(kind=prec)     :: sbypp,sbzpp,scypp,sczpp,saypp,sazpp
	complex(kind=prec)  :: ic,ce,bouah
    !
	ommi0=8.0d-7*PI*PI

    md1max=nb-1
    md2max=nb-2
    mh2max=n-nb-1
    mh1max=mh2max+1   !n-nb
    mhmax=mh1max+1    !n-nb+1
    meh1max=n+mh2max-1      !meh=2*n-nb
    meh3max=meh1max+2

    ic=cmplx(0.0d0,1.0d0)
	ommi=ommi0/period          
    ommii=1./ommi
    ! magnetic field on the erath surface
    bouah=hx(nb,1)
	!
	jk=0
	do j=1,m-2
	  jnode=j+1
      hym=y(j)
      hyp=y(j+1)
      do k=1,n-2
        knode=k+1
        hzm=z(k)
        hzp=z(k+1)

		jk=jk+1
! -----------------------------------------------------------
! approximation of the quasi-e-mode equation at (jnode,knode)
! -----------------------------------------------------------		
!	    o o o
!	    o x o   e
!	    o o o
        nodac=node_e(jnode,knode)
        
		do l=1,meh3max
		  apom(jk,l)=0.
		end do
		bpom(jk)=0.
        
		pe=0.5*(hzm+hzp)/hym
        qe=0.5*(hym+hyp)/hzm
        re=0.5*(hym+hyp)/hzp
        se=0.5*(hzm+hzp)/hyp
        te=-(pe+qe+re+se)

        if(node_h(jnode,knode).eq.0)then
          ce=te
          apom(jk,1)=ce
!	      o o o
!	      o o o   e
!	      o x o
          nodae=node_e(jnode,knode+1)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=re
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-re*ex(n,j+1)
          endif        
!  	      o o o
!	      o o x   e
!	      o o o
          nodae=node_e(jnode+1,knode)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=se
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-se*ex(k+1,m)
          endif
!	      o o o
!	      x o o   e
!	      o o o
          nodae=node_e(jnode-1,knode)
          if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-pe*ex(k+1,1)
          endif
!	      o x o
!	      o o o   e
!	      o o o
          nodae=node_e(jnode,knode-1)
          if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-qe*ex(1,j+1)
          endif

		else if(node_h(jnode,knode).ne.0)then
          ! k,j
          sigxxmm=cond(1)%v(k,j)
		  sigyymm=cond(2)%v(k,j)
		  sigzzmm=cond(3)%v(k,j)
		  sigxymm=cond(4)%v(k,j)
		  sigxzmm=cond(5)%v(k,j)
		  sigyzmm=cond(6)%v(k,j)
		  if(sigxxmm.lt.1e-6.and.sigyymm.lt.1e-6.and.  &
		     sigzzmm.lt.1e-6) then
            sbymm=0.
            sbzmm=0.
            scymm=0.
            sczmm=0.
            saymm=0.
            sazmm=0.
		  else
		    dd=sigyymm*sigzzmm-sigyzmm*sigyzmm ! D
		    sbymm=sigyzmm/dd
		    sbzmm=sigyymm/dd
		    scymm=sigzzmm/dd
		    sczmm=sigyzmm/dd
            saymm=sigxzmm*sbymm-sigxymm*scymm  ! B
            sazmm=sigxzmm*sbzmm-sigxymm*sczmm  !-A
          end if
          ! k,j+1
          sigxxpm=cond(1)%v(k,j+1)
		  sigyypm=cond(2)%v(k,j+1)
		  sigzzpm=cond(3)%v(k,j+1)
		  sigxypm=cond(4)%v(k,j+1)
		  sigxzpm=cond(5)%v(k,j+1)
		  sigyzpm=cond(6)%v(k,j+1)
		  if(sigxxpm.lt.1e-6.and.sigyypm.lt.1e-6.and.  &
		     sigzzpm.lt.1e-6) then
            sbypm=0.
            sbzpm=0.
            scypm=0.
            sczpm=0.
            saypm=0.
            sazpm=0.
		  else
		    dd=sigyypm*sigzzpm-sigyzpm*sigyzpm
		    sbypm=sigyzpm/dd
		    sbzpm=sigyypm/dd
		    scypm=sigzzpm/dd
		    sczpm=sigyzpm/dd
            saypm=sigxzpm*sbypm-sigxypm*scypm
            sazpm=sigxzpm*sbzpm-sigxypm*sczpm
          end if
          ! k+1,j
          sigxxmp=cond(1)%v(k+1,j)
		  sigyymp=cond(2)%v(k+1,j)
		  sigzzmp=cond(3)%v(k+1,j)
		  sigxymp=cond(4)%v(k+1,j)
		  sigxzmp=cond(5)%v(k+1,j)
		  sigyzmp=cond(6)%v(k+1,j)
		  if(sigxxmp.lt.1e-6.and.sigyymp.lt.1e-6.and.  &
		     sigzzmp.lt.1e-6) then
            sbymp=0.
            sbzmp=0.
            scymp=0.
            sczmp=0.
            saymp=0.
            sazmp=0.
		  else
		    dd=sigyymp*sigzzmp-sigyzmp*sigyzmp
		    sbymp=sigyzmp/dd
		    sbzmp=sigyymp/dd
		    scymp=sigzzmp/dd
		    sczmp=sigyzmp/dd
            saymp=sigxzmp*sbymp-sigxymp*scymp
            sazmp=sigxzmp*sbzmp-sigxymp*sczmp
          end if
          ! k+1,j+1
          sigxxpp=cond(1)%v(k+1,j+1)
		  sigyypp=cond(2)%v(k+1,j+1)
		  sigzzpp=cond(3)%v(k+1,j+1)
		  sigxypp=cond(4)%v(k+1,j+1)
		  sigxzpp=cond(5)%v(k+1,j+1)
		  sigyzpp=cond(6)%v(k+1,j+1)
		  if(sigxxpp.lt.1e-6.and.sigyypp.lt.1e-6.and.  &
		     sigzzpp.lt.1e-6) then
            sbypp=0.
            sbzpp=0.
            scypp=0.
            sczpp=0.
            saypp=0.
            sazpp=0.
		  else
		    dd=sigyypp*sigzzpp-sigyzpp*sigyzpp
		    sbypp=sigyzpp/dd
		    sbzpp=sigyypp/dd
		    scypp=sigzzpp/dd
		    sczpp=sigyzpp/dd
            saypp=sigxzpp*sbypp-sigxypp*scypp
            sazpp=sigxzpp*sbzpp-sigxypp*sczpp
          end if

	      ve=0.25*ommi*((sigxxmm+sigxymm*saymm-sigxzmm*sazmm)*hym*hzm  &
		     +(sigxxpm+sigxypm*saypm-sigxzpm*sazpm)*hyp*hzm  &
		     +(sigxxmp+sigxymp*saymp-sigxzmp*sazmp)*hym*hzp  &
		     +(sigxxpp+sigxypp*saypp-sigxzpp*sazpp)*hyp*hzp)
          ce=cmplx(te,ve)

          slyem=-(sigxymm*sbymm-sigxzmm*sbzmm)*hzm
          slyep=-(sigxymp*sbymp-sigxzmp*sbzmp)*hzp
          slye=0.25*ommi*(slyem+slyep)

          sryem=(sigxypm*sbypm-sigxzpm*sbzpm)*hzm
          sryep=(sigxypp*sbypp-sigxzpp*sbzpp)*hzp
          srye=0.25*ommi*(sryem+sryep)

          suzem=-(sigxymm*scymm-sigxzmm*sczmm)*hym
		  suzep=-(sigxypm*scypm-sigxzpm*sczpm)*hyp
          suze=0.25*ommi*(suzem+suzep)

          sdzem=(sigxymp*scymp-sigxzmp*sczmp)*hym
		  sdzep=(sigxypp*scypp-sigxzpp*sczpp)*hyp
          sdze=0.25*ommi*(sdzem+sdzep)

          sce=-(slye+suze+srye+sdze)

          apom(jk,1)=ce
!	      o x o
!	      o o o   e
!	      o o o
          nodae=node_e(jnode,knode-1)
          if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-qe*ex(1,j+1)
          endif
!	      o o o
!	      x o o   e
!	      o o o
          nodae=node_e(jnode-1,knode)
          if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-pe*ex(k+1,1)
          endif
!	      o o o
!	      o o o   e
!	      o x o
          nodae=node_e(jnode,knode+1)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=re
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-re*ex(n,j+1)
          endif
!	      o o o
!	      o o x   e
!	      o o o
          nodae=node_e(jnode+1,knode)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=se
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-se*ex(k+1,m)
          endif
!	      o o o
!	      o x o   h
!	      o o o
          nodah=node_h(jnode,knode)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=cmplx(0.,sce)
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ic*sce*bouah
          endif
!	      o o o
!	      o o o   h
!	      o x o
          nodah=node_h(jnode,knode+1)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=cmplx(0.,sdze)
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ic*sdze*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-ic*sdze*hx(n,j+1)
          endif
!	      o o o
!	      o o x   h
!	      o o o
          nodah=node_h(jnode+1,knode)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=cmplx(0.,srye)
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ic*srye*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-ic*srye*hx(k+1,m)
          endif
!	      o x o
!	      o o o   h
!	      o o o
          nodah=node_h(jnode,knode-1)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ic*suze*bouah
          endif
!	      o o o
!	      x o o   h
!	      o o o
          nodah=node_h(jnode-1,knode)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ic*slye*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-ic*slye*hx(k+1,1)
          endif

		end if

        do l=1,meh3max
          apom(jk,l)=-ic*ommii*apom(jk,l)
        end do
        bpom(jk)=-ic*ommii*bpom(jk)

! -----------------------------------------------------------
! approximation of the quasi-h-mode equation at (jnode,knode)
! -----------------------------------------------------------
!	    o o o
!	    o x o   h
!	    o o o
        nodac=node_h(jnode,knode)

        if(nodac.gt.0)then

		  jk=jk+1

		  do l=1,meh3max
		    apom(jk,l)=0.
	  	  end do
		  bpom(jk)=0.
    
          pmh=0.25*(sczmm+sbymm)
          pch=0.5*(sbzmm*hzm+sbzmp*hzp)/hym
          ph=pch
          pph=-0.25*(sczmp+sbymp)
          qch=0.5*(scymm*hym+scypm*hyp)/hzm
          qh=qch
          rch=0.5*(scymp*hym+scypp*hyp)/hzp
          rh=rch
          smh=-0.25*(sczpm+sbypm)
          sch=0.5*(sbzpm*hzm+sbzpp*hzp)/hyp
          sh=sch
          sph=0.25*(sczpp+sbypp)
          tch=-(pch+qch+rch+sch)
          th=tch+0.25*(sczmp+sczpm-sczmm-sczpp+sbymp+sbypm-sbymm-sbypp)
          vh=0.25*ommi*(hym+hyp)*(hzm+hzp)

		  apom(jk,1)=cmplx(th,vh)

!	      o o o
!	      o x o   e
!	      o o o
          nodae=node_e(jnode,knode)
!	      o o o
!	      o o o   h
!	      o x o
          nodah=node_h(jnode,knode+1)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=rh
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-rh*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-rh*hx(n,j+1)
          endif
!	      o o x
!	      o o o   h
!	      o o o
          nodah=node_h(jnode+1,knode-1)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=smh
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-smh*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-smh*hx(k,m)
          endif
!	      o o o
!	      o o x   h
!	      o o o
          nodah=node_h(jnode+1,knode)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=sh
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-sh*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-sh*hx(k+1,m)
          endif
!	      o o o
!	      o o o   h
!	      o o x
          nodah=node_h(jnode+1,knode+1)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=sph
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-sph*bouah
          else if(nodah.eq.-2)then
            if(j.lt.m-2)then
              bpom(jk)=bpom(jk)-sph*hx(n,j+2)
            else
              bpom(jk)=bpom(jk)-sph*hx(k+2,m)
            endif
          endif
!	      o x o
!	      o o o   h
!	      o o o
          nodah=node_h(jnode,knode-1)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-qh*bouah
          endif
!	      x o o
!	      o o o   h
!	      o o o
          nodah=node_h(jnode-1,knode-1)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-pmh*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-pmh*hx(k,1)
          endif
!	      o o o
!	      x o o   h
!	      o o o
          nodah=node_h(jnode-1,knode)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ph*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-ph*hx(k+1,1)
          endif
!	      o o o
!	      o o o   h
!	      x o o
          nodah=node_h(jnode-1,knode+1)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-pph*bouah
          else if(nodah.eq.-2)then
            if(j.gt.1)then
              bpom(jk)=bpom(jk)-pph*hx(n,j)
            else
              bpom(jk)=bpom(jk)-pph*hx(k+2,1)
            endif
          endif
!	      o o o
!	      o o o   e
!	      o x o
          nodae=node_e(jnode,knode+1)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=-ommii*sdze
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)+ommii*sdze*ex(n,j+1)
          endif
!	      o o o
!	      o o x   e
!	      o o o
          nodae=node_e(jnode+1,knode)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=-ommii*srye
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)+ommii*srye*ex(k+1,m)
          endif
!	      o x o
!	      o o o   e
!	      o o o
          nodae=node_e(jnode,knode-1)
!	      o o o
!	      x o o   e
!	      o o o
          nodae=node_e(jnode-1,knode)
          if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)+ommii*slye*ex(k+1,1)
          endif
		    
		end if

      end do ! k
    end do ! j

  end subroutine koef3_x


  ! construct the coefficient matrix and r.h.s for X-strike problem
  subroutine koef3_y(neh,meh,m,n,nlayr,nb,y,z,ihpol,period,  &
                   cond,node_e,node_h,ex,hx,apom,bpom)
    implicit none
    integer neh,meh,m,n,nlayr,nb,ihpol
	integer node_e(m,n),node_h(m,n)
    real(kind=prec)     :: period,y(*),z(*)
    type(rscalar_2d)	:: cond(6)
    complex(kind=prec)  :: ex(n,m),hx(n,m)
    complex(kind=prec)  :: apom(neh,meh), bpom(neh)
    ! local variables
	integer j,k,jk,jnode,knode,nodac,l,nodae,nodah
	integer md1max,md2max,mh1max,mh2max,mhmax,meh1max,meh3max
    real(kind=prec)     :: ommi0,ommi,dd,ommii,pe,qe,re,se,te,ve
	real(kind=prec)     :: slyem,slyep,slye,sryem,sryep,srye
	real(kind=prec)     :: suzem,suzep,suze,sdzem,sdzep,sdze,sce
	real(kind=prec)     :: pmh,pch,ph,pph,qch,qh,rch,rh
	real(kind=prec)     :: smh,sch,sh,sph,tch,th,vh
	real(kind=prec)     :: hym,hyp,hzm,hzp
	real(kind=prec)     :: sigxxmm,sigyymm,sigzzmm,sigxymm,sigxzmm,sigyzmm
	real(kind=prec)     :: sigxxpm,sigyypm,sigzzpm,sigxypm,sigxzpm,sigyzpm
	real(kind=prec)     :: sigxxmp,sigyymp,sigzzmp,sigxymp,sigxzmp,sigyzmp
	real(kind=prec)     :: sigxxpp,sigyypp,sigzzpp,sigxypp,sigxzpp,sigyzpp
	real(kind=prec)     :: sbymm,sbzmm,scymm,sczmm,saymm,sazmm
	real(kind=prec)     :: sbypm,sbzpm,scypm,sczpm,saypm,sazpm
	real(kind=prec)     :: sbymp,sbzmp,scymp,sczmp,saymp,sazmp
	real(kind=prec)     :: sbypp,sbzpp,scypp,sczpp,saypp,sazpp
	complex(kind=prec)  :: ic,ce,bouah
    !
	ommi0=8.0d-7*PI*PI

    md1max=nb-1
    md2max=nb-2
    mh2max=n-nb-1
    mh1max=mh2max+1   !n-nb
    mhmax=mh1max+1    !n-nb+1
    meh1max=n+mh2max-1      !meh=2*n-nb
    meh3max=meh1max+2

    ic=cmplx(0.0d0,1.0d0)
	ommi=ommi0/period          
    ommii=1./ommi
    ! magnetic field on the erath surface
    bouah=hx(nb,1)
	!
	jk=0
	do j=1,m-2
	  jnode=j+1
      hym=y(j)
      hyp=y(j+1)
      do k=1,n-2
        knode=k+1
        hzm=z(k)
        hzp=z(k+1)

		jk=jk+1
! -----------------------------------------------------------
! approximation of the quasi-e-mode equation at (jnode,knode)
! -----------------------------------------------------------		
!	    o o o
!	    o x o   e
!	    o o o
        nodac=node_e(jnode,knode)
        
		do l=1,meh3max
		  apom(jk,l)=0.
		end do
		bpom(jk)=0.
        
		pe=0.5*(hzm+hzp)/hym
        qe=0.5*(hym+hyp)/hzm
        re=0.5*(hym+hyp)/hzp
        se=0.5*(hzm+hzp)/hyp
        te=-(pe+qe+re+se)

        if(node_h(jnode,knode).eq.0)then
          ce=te
          apom(jk,1)=ce
!	      o o o
!	      o o o   e
!	      o x o
          nodae=node_e(jnode,knode+1)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=re
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-re*ex(n,j+1)
          endif        
!  	      o o o
!	      o o x   e
!	      o o o
          nodae=node_e(jnode+1,knode)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=se
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-se*ex(k+1,m)
          endif
!	      o o o
!	      x o o   e
!	      o o o
          nodae=node_e(jnode-1,knode)
          if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-pe*ex(k+1,1)
          endif
!	      o x o
!	      o o o   e
!	      o o o
          nodae=node_e(jnode,knode-1)
          if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-qe*ex(1,j+1)
          endif

		else if(node_h(jnode,knode).ne.0)then
          ! k,j
          sigxxmm=cond(1)%v(k,j)
		  sigyymm=cond(2)%v(k,j)
		  sigzzmm=cond(3)%v(k,j)
		  sigxymm=cond(4)%v(k,j)
		  sigxzmm=cond(5)%v(k,j)
		  sigyzmm=cond(6)%v(k,j)
		  if(sigxxmm.lt.1e-6.and.sigyymm.lt.1e-6.and.  &
		     sigzzmm.lt.1e-6) then
            sbymm=0.
            sbzmm=0.
            scymm=0.
            sczmm=0.
            saymm=0.
            sazmm=0.
		  else
		    dd=sigxxmm*sigzzmm-sigxzmm*sigxzmm ! D
		    sbymm=sigyzmm/dd
		    sbzmm=sigxxmm/dd
		    scymm=sigzzmm/dd
		    sczmm=sigxzmm/dd
            saymm=sigxzmm*sbymm-sigxymm*scymm  ! B
            sazmm=sigyzmm*sbzmm-sigxymm*sczmm  !-A
          end if
          ! k,j+1
          sigxxpm=cond(1)%v(k,j+1)
		  sigyypm=cond(2)%v(k,j+1)
		  sigzzpm=cond(3)%v(k,j+1)
		  sigxypm=cond(4)%v(k,j+1)
		  sigxzpm=cond(5)%v(k,j+1)
		  sigyzpm=cond(6)%v(k,j+1)
		  if(sigxxpm.lt.1e-6.and.sigyypm.lt.1e-6.and.  &
		     sigzzpm.lt.1e-6) then
            sbypm=0.
            sbzpm=0.
            scypm=0.
            sczpm=0.
            saypm=0.
            sazpm=0.
		  else
		    dd=sigxxpm*sigzzpm-sigxzpm*sigxzpm
		    sbypm=sigyzpm/dd
		    sbzpm=sigxxpm/dd
		    scypm=sigzzpm/dd
		    sczpm=sigxzpm/dd
            saypm=sigxzpm*sbypm-sigxypm*scypm
            sazpm=sigyzpm*sbzpm-sigxypm*sczpm
          end if
          ! k+1,j
          sigxxmp=cond(1)%v(k+1,j)
		  sigyymp=cond(2)%v(k+1,j)
		  sigzzmp=cond(3)%v(k+1,j)
		  sigxymp=cond(4)%v(k+1,j)
		  sigxzmp=cond(5)%v(k+1,j)
		  sigyzmp=cond(6)%v(k+1,j)
		  if(sigxxmp.lt.1e-6.and.sigyymp.lt.1e-6.and.  &
		     sigzzmp.lt.1e-6) then
            sbymp=0.
            sbzmp=0.
            scymp=0.
            sczmp=0.
            saymp=0.
            sazmp=0.
		  else
		    dd=sigxxmp*sigzzmp-sigxzmp*sigxzmp
		    sbymp=sigyzmp/dd
		    sbzmp=sigxxmp/dd
		    scymp=sigzzmp/dd
		    sczmp=sigxzmp/dd
            saymp=sigxzmp*sbymp-sigxymp*scymp
            sazmp=sigyzmp*sbzmp-sigxymp*sczmp
          end if
          ! k+1,j+1
          sigxxpp=cond(1)%v(k+1,j+1)
		  sigyypp=cond(2)%v(k+1,j+1)
		  sigzzpp=cond(3)%v(k+1,j+1)
		  sigxypp=cond(4)%v(k+1,j+1)
		  sigxzpp=cond(5)%v(k+1,j+1)
		  sigyzpp=cond(6)%v(k+1,j+1)
		  if(sigxxpp.lt.1e-6.and.sigyypp.lt.1e-6.and.  &
		     sigzzpp.lt.1e-6) then
            sbypp=0.
            sbzpp=0.
            scypp=0.
            sczpp=0.
            saypp=0.
            sazpp=0.
		  else
		    dd=sigxxpp*sigzzpp-sigxzpp*sigxzpp
		    sbypp=sigyzpp/dd
		    sbzpp=sigxxpp/dd
		    scypp=sigzzpp/dd
		    sczpp=sigxzpp/dd
            saypp=sigxzpp*sbypp-sigxypp*scypp
            sazpp=sigyzpp*sbzpp-sigxypp*sczpp
          end if

	      ve=0.25*ommi*((sigyymm+sigxymm*saymm-sigyzmm*sazmm)*hym*hzm  &
		     +(sigyypm+sigxypm*saypm-sigyzpm*sazpm)*hyp*hzm  &
		     +(sigyymp+sigxymp*saymp-sigyzmp*sazmp)*hym*hzp  &
		     +(sigyypp+sigxypp*saypp-sigyzpp*sazpp)*hyp*hzp)
          ce=cmplx(te,ve)

          slyem=(sigxymm*sczmm-sigyzmm*sbzmm)*hzm
          slyep=(sigxymp*sczmp-sigyzmp*sbzmp)*hzp
          slye=0.25*ommi*(slyem+slyep)

          sryem=-(sigxypm*sczpm-sigyzpm*sbzpm)*hzm
          sryep=-(sigxypp*sczpp-sigyzpp*sbzpp)*hzp
          srye=0.25*ommi*(sryem+sryep)

          suzem=(sigxymm*scymm-sigxzmm*sbymm)*hym
		  suzep=(sigxypm*scypm-sigxzpm*sbypm)*hyp
          suze=0.25*ommi*(suzem+suzep)

          sdzem=-(sigxymp*scymp-sigxzmp*sbymp)*hym
		  sdzep=-(sigxypp*scypp-sigxzpp*sbypp)*hyp
          sdze=0.25*ommi*(sdzem+sdzep)

          sce=-(slye+suze+srye+sdze)

          apom(jk,1)=ce
!	      o x o
!	      o o o   e
!	      o o o
          nodae=node_e(jnode,knode-1)
          if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-qe*ex(1,j+1)
          endif
!	      o o o
!	      x o o   e
!	      o o o
          nodae=node_e(jnode-1,knode)
          if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-pe*ex(k+1,1)
          endif
!	      o o o
!	      o o o   e
!	      o x o
          nodae=node_e(jnode,knode+1)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=re
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-re*ex(n,j+1)
          endif
!	      o o o
!	      o o x   e
!	      o o o
          nodae=node_e(jnode+1,knode)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=se
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)-se*ex(k+1,m)
          endif
!	      o o o
!	      o x o   h
!	      o o o
          nodah=node_h(jnode,knode)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=cmplx(0.,sce)
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ic*sce*bouah
          endif
!	      o o o
!	      o o o   h
!	      o x o
          nodah=node_h(jnode,knode+1)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=cmplx(0.,sdze)
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ic*sdze*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-ic*sdze*hx(n,j+1)
          endif
!	      o o o
!	      o o x   h
!	      o o o
          nodah=node_h(jnode+1,knode)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=cmplx(0.,srye)
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ic*srye*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-ic*srye*hx(k+1,m)
          endif
!	      o x o
!	      o o o   h
!	      o o o
          nodah=node_h(jnode,knode-1)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ic*suze*bouah
          endif
!	      o o o
!	      x o o   h
!	      o o o
          nodah=node_h(jnode-1,knode)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ic*slye*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-ic*slye*hx(k+1,1)
          endif

		end if

        do l=1,meh3max
          apom(jk,l)=-ic*ommii*apom(jk,l)
        end do
        bpom(jk)=-ic*ommii*bpom(jk)

! -----------------------------------------------------------
! approximation of the quasi-h-mode equation at (jnode,knode)
! -----------------------------------------------------------
!	    o o o
!	    o x o   h
!	    o o o
        nodac=node_h(jnode,knode)

        if(nodac.gt.0)then

		  jk=jk+1

		  do l=1,meh3max
		    apom(jk,l)=0.
	  	  end do
		  bpom(jk)=0.
    
          pmh=0.5*sczmm
          pch=0.5*(sbzmm*hzm+sbzmp*hzp)/hym
          ph=pch
          pph=-0.5*sczmp
          qch=0.5*(scymm*hym+scypm*hyp)/hzm
          qh=qch
          rch=0.5*(scymp*hym+scypp*hyp)/hzp
          rh=rch
          smh=-0.5*sczpm
          sch=0.5*(sbzpm*hzm+sbzpp*hzp)/hyp
          sh=sch
          sph=0.5*sczpp
          tch=-(pch+qch+rch+sch)
          th=tch+0.5*(sczmp+sczpm-sczmm-sczpp)
          vh=0.25*ommi*(hym+hyp)*(hzm+hzp)

		  apom(jk,1)=cmplx(th,vh)

!	      o o o
!	      o x o   e
!	      o o o
          nodae=node_e(jnode,knode)
!	      o o o
!	      o o o   h
!	      o x o
          nodah=node_h(jnode,knode+1)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=rh
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-rh*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-rh*hx(n,j+1)
          endif
!	      o o x
!	      o o o   h
!	      o o o
          nodah=node_h(jnode+1,knode-1)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=smh
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-smh*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-smh*hx(k,m)
          endif
!	      o o o
!	      o o x   h
!	      o o o
          nodah=node_h(jnode+1,knode)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=sh
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-sh*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-sh*hx(k+1,m)
          endif
!	      o o o
!	      o o o   h
!	      o o x
          nodah=node_h(jnode+1,knode+1)
          if(nodah.gt.0)then
            apom(jk,1+nodah-nodac)=sph
          else if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-sph*bouah
          else if(nodah.eq.-2)then
            if(j.lt.m-2)then
              bpom(jk)=bpom(jk)-sph*hx(n,j+2)
            else
              bpom(jk)=bpom(jk)-sph*hx(k+2,m)
            endif
          endif
!	      o x o
!	      o o o   h
!	      o o o
          nodah=node_h(jnode,knode-1)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-qh*bouah
          endif
!	      x o o
!	      o o o   h
!	      o o o
          nodah=node_h(jnode-1,knode-1)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-pmh*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-pmh*hx(k,1)
          endif
!	      o o o
!	      x o o   h
!	      o o o
          nodah=node_h(jnode-1,knode)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-ph*bouah
          else if(nodah.eq.-2)then
            bpom(jk)=bpom(jk)-ph*hx(k+1,1)
          endif
!	      o o o
!	      o o o   h
!	      x o o
          nodah=node_h(jnode-1,knode+1)
          if(nodah.eq.-1)then
            bpom(jk)=bpom(jk)-pph*bouah
          else if(nodah.eq.-2)then
            if(j.gt.1)then
              bpom(jk)=bpom(jk)-pph*hx(n,j)
            else
              bpom(jk)=bpom(jk)-pph*hx(k+2,1)
            endif
          endif
!	      o o o
!	      o o o   e
!	      o x o
          nodae=node_e(jnode,knode+1)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=-ommii*sdze
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)+ommii*sdze*ex(n,j+1)
          endif
!	      o o o
!	      o o x   e
!	      o o o
          nodae=node_e(jnode+1,knode)
          if(nodae.gt.0)then
            apom(jk,1+nodae-nodac)=-ommii*srye
          else if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)+ommii*srye*ex(k+1,m)
          endif
!	      o x o
!	      o o o   e
!	      o o o
          nodae=node_e(jnode,knode-1)
!	      o o o
!	      x o o   e
!	      o o o
          nodae=node_e(jnode-1,knode)
          if(nodae.eq.-2)then
            bpom(jk)=bpom(jk)+ommii*slye*ex(k+1,1)
          endif
		    
		end if

      end do ! k
    end do ! j

  end subroutine koef3_y


  ! Gauss Eliminition Solver for 2D anisotropy problem
  subroutine gaussr(m,n,nb,neh,meh,node_h,apom,bpom)
    implicit none
	integer m,n,nb
	integer neh,meh,node_h(m,n)
	complex(kind=prec)  :: apom(neh,meh),bpom(neh)
	! local variables
	integer m1,nk,m1stor,mh2,nk1,nkm,mg,nkm1
	integer i,i1,me,k,ik,j,l,ne
	integer n9,j9,j9m
	real(kind=prec)     :: aik,ank
	complex(kind=prec)  :: c
    !
	nk=node_h(m-1,n-1)
    m1stor=n-1
    mh2=n-nb-1
    m1=meh
    
	nk1=nk-1
    nkm=nk-m1+2
    mg=m1
    nkm1=nkm-1

    do 210 i=1,nk1
      i1=i-1
      if(i.ge.nkm)mg=mg-1
      me=2
      do 220 k=2,mg
        aik=cdabs(apom(i,k))
        if(aik.eq.0.) goto 220
        c=apom(i,k)/apom(i,1)
        ik=i1+k
        j=0
        do 230 l=me,mg
          j=j+1
230     apom(ik,j)=apom(ik,j)-c*apom(i,l)
        bpom(ik)=bpom(ik)-c*bpom(i)
220   me=me+1
210 continue

    ne=nk+1
310 ne=ne-1
    if(ne)350,350,320
320 l=ne
    do 340 k=2,m1
      l=l+1
      if(l.gt.nk)goto 360
      ank=cdabs(apom(ne,k))
      if(ank)330,340,330
330   bpom(ne)=bpom(ne)-apom(ne,k)*bpom(l)
340 continue
360 bpom(ne)=bpom(ne)/apom(ne,1)
    goto 310
350 continue

  end subroutine gaussr


  ! get the electric and magnetic fields from the r.h.s
  subroutine getfld(m,n,nb,neh,ihpol,node_e,node_h,bpom,ex,hx)
    implicit none
    integer m,n,nb,neh,ihpol
	integer node_e(m,n),node_h(m,n)
    complex(kind=prec)  :: ex(n,m),hx(n,m),bpom(neh)
    ! local variables
	integer j,k
	!
    do j=2,m-1
	  do k=2,n-1
        if(k.le.nb)then
          ex(k,j)=bpom(node_e(j,k))
        else 
          ex(k,j)=bpom(node_e(j,k))
          hx(k,j)=bpom(node_h(j,k))
        endif
	  end do
    end do

! debug
!   do j=1,m
!     do k=1,n
!       write(107,3000) ihpol,j,k,ex(k,j)
!       write(108,3000) ihpol,j,k,hx(k,j)
!     end do
!	end do
!
! 3000 FORMAT(I3,1X,I3,1X,I3,3X,2E11.3)

  end subroutine getfld


  ! compute hy and hz
  subroutine gethyhz(m,n,nlayr,ihpol,y,z,period,cond,ex,hx,hy,hz)
    implicit none
	integer m,n,nlayr,ihpol
    real(kind=prec)    :: period,y(*),z(*)
    type(rscalar_2d)   :: cond(6)
	complex(kind=prec) :: ex(n,m),hx(n,m),hy(n,m),hz(n,m) 
    ! local variables
    integer j,k
    real(kind=prec)    :: hym,hyp,hzm,hzc
	complex(kind=prec) :: pem,pec,pep,pel,per,phm,phc,php,phl,phr
	complex(kind=prec) :: omdedzb,omdedzc,omdedz,omdedyb,omdedyc,omdedy
    !
	do j=2,m-1
	  do k=2,n-1
        pem=ex(k-1,j)   
        pec=ex(k,j)
        pep=ex(k+1,j)
        pel=ex(k,j-1)
        per=ex(k,j+1)
        phm=hx(k-1,j)     
        phc=hx(k,j)
        php=hx(k+1,j)
        phl=hx(k,j-1)
        phr=hx(k,j+1)
!
        hym=y(j-1)
        hyp=y(j)
        hzm=z(k-1)
        hzc=z(k)

	    call deredz_x(j,k,pem,pec,pep,phc,phm,php,phl,phr,hzm,  & 
             hzc,hym,hyp,nlayr,period,cond,omdedzb,omdedzc,omdedz)
	    hy(k,j)=omdedz

        call deredy_x(j,k,pel,pec,per,phc,phl,phr,phm,php,hym,  &
             hyp,hzm,hzc,nlayr,period,cond,omdedyb,omdedyc,omdedy)
		hz(k,j)=omdedy

	  end do
    end do

! debug
!   do j=1,m
!     do k=1,n
!       write(103,3000) ihpol,j,k,hy(k,j)
!       write(104,3000) ihpol,j,k,hz(k,j)
!     end do
!	end do
!
! 3000 FORMAT(I3,1X,I3,1X,I3,3X,2E11.3)

  end subroutine gethyhz


  ! *****************************************************************
  subroutine deredz_x(jnode,knode,pem,pec,pep,phc,phm,php,phl,phr,  &
          hzm,hzc,hym,hyp,nlayr,period,cond,omdedzb,omdedzc,omdedz)
    implicit none
    integer jnode,knode,nlayr
    real(kind=prec)    :: hym,hyp,hzm,hzc,period
    type(rscalar_2d)   :: cond(6)
	complex(kind=prec) :: pem,pec,pep,phm,phc,php,phl,phr
	complex(kind=prec) :: omdedzb,omdedzc,omdedz
	! local variables
	integer is
	real(kind=prec)    :: ommi0,ommi,dd,ommii,hzmc,hymp
	real(kind=prec)    :: sigxxmm,sigyymm,sigzzmm,sigxymm,sigxzmm,sigyzmm
	real(kind=prec)    :: sigxxpm,sigyypm,sigzzpm,sigxypm,sigxzpm,sigyzpm
	real(kind=prec)    :: sigxxmp,sigyymp,sigzzmp,sigxymp,sigxzmp,sigyzmp
	real(kind=prec)    :: sigxxpp,sigyypp,sigzzpp,sigxypp,sigxzpp,sigyzpp
	real(kind=prec)    :: sbymm,sbzmm,scymm,sczmm,saymm,sazmm,sigmm
	real(kind=prec)    :: sbypm,sbzpm,scypm,sczpm,saypm,sazpm,sigpm
	real(kind=prec)    :: sbymp,sbzmp,scymp,sczmp,saymp,sazmp,sigmp
	real(kind=prec)    :: sbypp,sbzpp,scypp,sczpp,saypp,sazpp,sigpp
	real(kind=prec)    :: sigm,sigp,saym,sayp
	complex(kind=prec) :: ic
	!
	ommi0=8.0d-7*PI*PI
    ic=cmplx(0.0d0,1.0d0)
	ommi=ommi0/period          
    ommii=1./ommi
	is=2

    call der3po(pem,pec,pep,hzm,hzc,is,omdedzb)
    omdedzb=-ic*ommii*omdedzb

    hzmc=hzm+hzc
	hymp=hym+hyp

    ! k,j
    sigxxmm=cond(1)%v(knode-1,jnode-1)
	sigyymm=cond(2)%v(knode-1,jnode-1)
	sigzzmm=cond(3)%v(knode-1,jnode-1)
	sigxymm=cond(4)%v(knode-1,jnode-1)
	sigxzmm=cond(5)%v(knode-1,jnode-1)
	sigyzmm=cond(6)%v(knode-1,jnode-1)
	if(sigxxmm.lt.1e-6.and.sigyymm.lt.1e-6.and.  &
		     sigzzmm.lt.1e-6) then
      sbymm=0.
      sbzmm=0.
      scymm=0.
      sczmm=0.
      saymm=0.
      sazmm=0.
	  sigmm=0.
	else
	  dd=sigyymm*sigzzmm-sigyzmm*sigyzmm ! D
	  sbymm=sigyzmm/dd
	  sbzmm=sigyymm/dd
	  scymm=sigzzmm/dd
	  sczmm=sigyzmm/dd
      saymm=sigxzmm*sbymm-sigxymm*scymm  ! B
      sazmm=sigxzmm*sbzmm-sigxymm*sczmm  !-A
	  sigmm=sigxxmm+sigxymm*saymm-sigxzmm*sazmm
    end if
    ! k,j+1
    sigxxpm=cond(1)%v(knode-1,jnode)
	sigyypm=cond(2)%v(knode-1,jnode)
	sigzzpm=cond(3)%v(knode-1,jnode)
	sigxypm=cond(4)%v(knode-1,jnode)
	sigxzpm=cond(5)%v(knode-1,jnode)
	sigyzpm=cond(6)%v(knode-1,jnode)
	if(sigxxpm.lt.1e-6.and.sigyypm.lt.1e-6.and.  &
		     sigzzpm.lt.1e-6) then
      sbypm=0.
      sbzpm=0.
      scypm=0.
      sczpm=0.
      saypm=0.
      sazpm=0.
	  sigpm=0.
	else
	  dd=sigyypm*sigzzpm-sigyzpm*sigyzpm
	  sbypm=sigyzpm/dd
	  sbzpm=sigyypm/dd
	  scypm=sigzzpm/dd
	  sczpm=sigyzpm/dd
      saypm=sigxzpm*sbypm-sigxypm*scypm
      sazpm=sigxzpm*sbzpm-sigxypm*sczpm
	  sigpm=sigxxpm+sigxypm*saypm-sigxzpm*sazpm
    end if
    ! k+1,j
    sigxxmp=cond(1)%v(knode,jnode-1)
	sigyymp=cond(2)%v(knode,jnode-1)
	sigzzmp=cond(3)%v(knode,jnode-1)
    sigxymp=cond(4)%v(knode,jnode-1)
	sigxzmp=cond(5)%v(knode,jnode-1)
	sigyzmp=cond(6)%v(knode,jnode-1)
	if(sigxxmp.lt.1e-6.and.sigyymp.lt.1e-6.and.  &
		     sigzzmp.lt.1e-6) then
      sbymp=0.
      sbzmp=0.
      scymp=0.
      sczmp=0.
      saymp=0.
      sazmp=0.
	  sigmp=0.
	else
	  dd=sigyymp*sigzzmp-sigyzmp*sigyzmp
	  sbymp=sigyzmp/dd
	  sbzmp=sigyymp/dd
	  scymp=sigzzmp/dd
	  sczmp=sigyzmp/dd
      saymp=sigxzmp*sbymp-sigxymp*scymp
      sazmp=sigxzmp*sbzmp-sigxymp*sczmp
	  sigmp=sigxxmp+sigxymp*saymp-sigxzmp*sazmp
    end if
    ! k+1,j+1
    sigxxpp=cond(1)%v(knode,jnode)
	sigyypp=cond(2)%v(knode,jnode)
	sigzzpp=cond(3)%v(knode,jnode)
	sigxypp=cond(4)%v(knode,jnode)
	sigxzpp=cond(5)%v(knode,jnode)
	sigyzpp=cond(6)%v(knode,jnode)
	if(sigxxpp.lt.1e-6.and.sigyypp.lt.1e-6.and.  &
		     sigzzpp.lt.1e-6) then
      sbypp=0.
      sbzpp=0.
      scypp=0.
      sczpp=0.
      saypp=0.
      sazpp=0.
	  sigpp=0.
	else
	  dd=sigyypp*sigzzpp-sigyzpp*sigyzpp
	  sbypp=sigyzpp/dd
	  sbzpp=sigyypp/dd
	  scypp=sigzzpp/dd
	  sczpp=sigyzpp/dd
      saypp=sigxzpp*sbypp-sigxypp*scypp
      sazpp=sigxzpp*sbzpp-sigxypp*sczpp
	  sigpp=sigxxpp+sigxypp*saypp-sigxzpp*sazpp
    end if

	sigm=(hym*sigmm+hyp*sigpm)/hymp
	sigp=(hym*sigmp+hyp*sigpp)/hymp
	omdedzc=(0.5*hzm*hzc*(sigp-sigm)/hzmc)*pec
	omdedzc=omdedzc+(0.5*hzm*hzc/(hzmc*hymp))*   &
          ((sazpm-sazpp)*(phr-phc)+(sazmm-sazmp)*(phc-phl))
	saym=(hym*saymm+hyp*saypm)/hymp
	sayp=(hym*saymp+hyp*saypp)/hymp
	omdedzc=omdedzc-0.5*(hzm*sayp*(php-phc)-hzc*saym*(phc-phm))/hzmc

	omdedz=omdedzb+omdedzc

  end subroutine deredz_x


  ! *****************************************************************
  subroutine deredy_x(jnode,knode,pel,pec,per,phc,phl,phr,phm,php,  & 
          hym,hyp,hzm,hzc,nlayr,period,cond,omdedyb,omdedyc,omdedy)
    implicit none
    integer jnode,knode,nlayr
    real(kind=prec)    :: hym,hyp,hzm,hzc,period
    type(rscalar_2d)   :: cond(6)
	complex(kind=prec) :: pel,pec,per,phc,phl,phr,phm,php
	complex(kind=prec) :: omdedyb,omdedyc,omdedy
	! local variables
	integer is
	real(kind=prec)    :: ommi0,ommi,dd,ommii,hzmc,hymp
	real(kind=prec)    :: sigxxmm,sigyymm,sigzzmm,sigxymm,sigxzmm,sigyzmm
	real(kind=prec)    :: sigxxpm,sigyypm,sigzzpm,sigxypm,sigxzpm,sigyzpm
	real(kind=prec)    :: sigxxmp,sigyymp,sigzzmp,sigxymp,sigxzmp,sigyzmp
	real(kind=prec)    :: sigxxpp,sigyypp,sigzzpp,sigxypp,sigxzpp,sigyzpp
	real(kind=prec)    :: sbymm,sbzmm,scymm,sczmm,saymm,sazmm,sigmm
	real(kind=prec)    :: sbypm,sbzpm,scypm,sczpm,saypm,sazpm,sigpm
	real(kind=prec)    :: sbymp,sbzmp,scymp,sczmp,saymp,sazmp,sigmp
	real(kind=prec)    :: sbypp,sbzpp,scypp,sczpp,saypp,sazpp,sigpp
	real(kind=prec)    :: sigl,sigr,sazl,sazr
	complex(kind=prec) :: ic
	!
	ommi0=8.0d-7*PI*PI
    ic=cmplx(0.0d0,1.0d0)
	ommi=ommi0/period          
    ommii=1./ommi
	is=2

	call der3po(pel,pec,per,hym,hyp,is,omdedyb)
	omdedyb=ic*ommii*omdedyb

    hzmc=hzm+hzc
	hymp=hym+hyp

    ! k,j
    sigxxmm=cond(1)%v(knode-1,jnode-1)
	sigyymm=cond(2)%v(knode-1,jnode-1)
	sigzzmm=cond(3)%v(knode-1,jnode-1)
	sigxymm=cond(4)%v(knode-1,jnode-1)
	sigxzmm=cond(5)%v(knode-1,jnode-1)
	sigyzmm=cond(6)%v(knode-1,jnode-1)
	if(sigxxmm.lt.1e-6.and.sigyymm.lt.1e-6.and.  &
		     sigzzmm.lt.1e-6) then
      sbymm=0.
      sbzmm=0.
      scymm=0.
      sczmm=0.
      saymm=0.
      sazmm=0.
	  sigmm=0.
	else
	  dd=sigyymm*sigzzmm-sigyzmm*sigyzmm ! D
	  sbymm=sigyzmm/dd
	  sbzmm=sigyymm/dd
	  scymm=sigzzmm/dd
	  sczmm=sigyzmm/dd
      saymm=sigxzmm*sbymm-sigxymm*scymm  ! B
      sazmm=sigxzmm*sbzmm-sigxymm*sczmm  !-A
	  sigmm=sigxxmm+sigxymm*saymm-sigxzmm*sazmm
    end if
    ! k,j+1
    sigxxpm=cond(1)%v(knode-1,jnode)
	sigyypm=cond(2)%v(knode-1,jnode)
	sigzzpm=cond(3)%v(knode-1,jnode)
	sigxypm=cond(4)%v(knode-1,jnode)
	sigxzpm=cond(5)%v(knode-1,jnode)
	sigyzpm=cond(6)%v(knode-1,jnode)
	if(sigxxpm.lt.1e-6.and.sigyypm.lt.1e-6.and.  &
		     sigzzpm.lt.1e-6) then
      sbypm=0.
      sbzpm=0.
      scypm=0.
      sczpm=0.
      saypm=0.
      sazpm=0.
	  sigpm=0.
	else
	  dd=sigyypm*sigzzpm-sigyzpm*sigyzpm
	  sbypm=sigyzpm/dd
	  sbzpm=sigyypm/dd
	  scypm=sigzzpm/dd
	  sczpm=sigyzpm/dd
      saypm=sigxzpm*sbypm-sigxypm*scypm
      sazpm=sigxzpm*sbzpm-sigxypm*sczpm
	  sigpm=sigxxpm+sigxypm*saypm-sigxzpm*sazpm
    end if
    ! k+1,j
    sigxxmp=cond(1)%v(knode,jnode-1)
	sigyymp=cond(2)%v(knode,jnode-1)
	sigzzmp=cond(3)%v(knode,jnode-1)
    sigxymp=cond(4)%v(knode,jnode-1)
	sigxzmp=cond(5)%v(knode,jnode-1)
	sigyzmp=cond(6)%v(knode,jnode-1)
	if(sigxxmp.lt.1e-6.and.sigyymp.lt.1e-6.and.  &
		     sigzzmp.lt.1e-6) then
      sbymp=0.
      sbzmp=0.
      scymp=0.
      sczmp=0.
      saymp=0.
      sazmp=0.
	  sigmp=0.
	else
	  dd=sigyymp*sigzzmp-sigyzmp*sigyzmp
	  sbymp=sigyzmp/dd
	  sbzmp=sigyymp/dd
	  scymp=sigzzmp/dd
	  sczmp=sigyzmp/dd
      saymp=sigxzmp*sbymp-sigxymp*scymp
      sazmp=sigxzmp*sbzmp-sigxymp*sczmp
	  sigmp=sigxxmp+sigxymp*saymp-sigxzmp*sazmp
    end if
    ! k+1,j+1
    sigxxpp=cond(1)%v(knode,jnode)
	sigyypp=cond(2)%v(knode,jnode)
	sigzzpp=cond(3)%v(knode,jnode)
	sigxypp=cond(4)%v(knode,jnode)
	sigxzpp=cond(5)%v(knode,jnode)
	sigyzpp=cond(6)%v(knode,jnode)
	if(sigxxpp.lt.1e-6.and.sigyypp.lt.1e-6.and.  &
		     sigzzpp.lt.1e-6) then
      sbypp=0.
      sbzpp=0.
      scypp=0.
      sczpp=0.
      saypp=0.
      sazpp=0.
	  sigpp=0.
	else
	  dd=sigyypp*sigzzpp-sigyzpp*sigyzpp
	  sbypp=sigyzpp/dd
	  sbzpp=sigyypp/dd
	  scypp=sigzzpp/dd
	  sczpp=sigyzpp/dd
      saypp=sigxzpp*sbypp-sigxypp*scypp
      sazpp=sigxzpp*sbzpp-sigxypp*sczpp
	  sigpp=sigxxpp+sigxypp*saypp-sigxzpp*sazpp
    end if
	
	sigl=(hzm*sigmm+hzc*sigmp)/hzmc
	sigr=(hzm*sigpm+hzc*sigpp)/hzmc
	omdedyc=-(0.5*hym*hyp*(sigr-sigl)/hymp)*pec
	omdedyc=omdedyc+(0.5*hym*hyp/(hzmc*hymp))*   &
           ((saypp-saymp)*(php-phc)+(saypm-saymm)*(phc-phm))
	sazl=-(hzm*sazmm+hzc*sazmp)/hzmc
	sazr=-(hzm*sazpm+hzc*sazpp)/hzmc
	omdedyc=omdedyc-0.5*(hym*sazr*(phr-phc)-hyp*sazl*(phc-phl))/hymp

	omdedy=omdedyb+omdedyc

!    if(jnode.eq.3) then
!	  write(107,*) omdedyb !debug
!	  write(107,*) sigl,sigr !debug
!      write(107,*) sazl,sazr !debug
!	  write(107,*) omdedyc !debug
!	  pause
!	end if
			  
  end subroutine deredy_x 


  ! compute hy and hz
  subroutine gethxhz(m,n,nlayr,ihpol,y,z,period,cond,ex,hx,hy,hz)
    implicit none
	integer m,n,nlayr,ihpol
    real(kind=prec)    :: period,y(*),z(*)
    type(rscalar_2d)   :: cond(6)
	complex(kind=prec) :: ex(n,m),hx(n,m),hy(n,m),hz(n,m) 
    ! local variables
    integer j,k
    real(kind=prec)    :: hym,hyp,hzm,hzc
	complex(kind=prec) :: pem,pec,pep,pel,per,phm,phc,php,phl,phr
	complex(kind=prec) :: omdedzb,omdedzc,omdedz,omdedyb,omdedyc,omdedy
    !
	do j=2,m-1
	  do k=2,n-1
        pem=ex(k-1,j)   
        pec=ex(k,j)
        pep=ex(k+1,j)
        pel=ex(k,j-1)
        per=ex(k,j+1)
        phm=hx(k-1,j)     
        phc=hx(k,j)
        php=hx(k+1,j)
        phl=hx(k,j-1)
        phr=hx(k,j+1)
!
        hym=y(j-1)
        hyp=y(j)
        hzm=z(k-1)
        hzc=z(k)

	    call deredz_y(j,k,pem,pec,pep,phc,phm,php,phl,phr,hzm,  & 
             hzc,hym,hyp,nlayr,period,cond,omdedzb,omdedzc,omdedz)
	    hy(k,j)=omdedz

        call deredy_y(j,k,pel,pec,per,phc,phl,phr,phm,php,hym,  &
             hyp,hzm,hzc,nlayr,period,cond,omdedyb,omdedyc,omdedy)
		hz(k,j)=omdedy

	  end do
    end do

! debug
!   do j=1,m
!     do k=1,n
!       write(103,3000) ihpol,j,k,hy(k,j)
!       write(104,3000) ihpol,j,k,hz(k,j)
!     end do
!	end do
!
! 3000 FORMAT(I3,1X,I3,1X,I3,3X,2E11.3)

  end subroutine gethxhz


  ! *****************************************************************
  subroutine deredz_y(jnode,knode,pem,pec,pep,phc,phm,php,phl,phr,  &
          hzm,hzc,hym,hyp,nlayr,period,cond,omdedzb,omdedzc,omdedz)
    implicit none
    integer jnode,knode,nlayr
    real(kind=prec)    :: hym,hyp,hzm,hzc,period
    type(rscalar_2d)   :: cond(6)
	complex(kind=prec) :: pem,pec,pep,phm,phc,php,phl,phr
	complex(kind=prec) :: omdedzb,omdedzc,omdedz
	! local variables
	integer is
	real(kind=prec)    :: ommi0,ommi,dd,ommii,hzmc,hymp
	real(kind=prec)    :: sigxxmm,sigyymm,sigzzmm,sigxymm,sigxzmm,sigyzmm
	real(kind=prec)    :: sigxxpm,sigyypm,sigzzpm,sigxypm,sigxzpm,sigyzpm
	real(kind=prec)    :: sigxxmp,sigyymp,sigzzmp,sigxymp,sigxzmp,sigyzmp
	real(kind=prec)    :: sigxxpp,sigyypp,sigzzpp,sigxypp,sigxzpp,sigyzpp
	real(kind=prec)    :: sbymm,sbzmm,scymm,sczmm,saymm,sazmm,sigmm
	real(kind=prec)    :: sbypm,sbzpm,scypm,sczpm,saypm,sazpm,sigpm
	real(kind=prec)    :: sbymp,sbzmp,scymp,sczmp,saymp,sazmp,sigmp
	real(kind=prec)    :: sbypp,sbzpp,scypp,sczpp,saypp,sazpp,sigpp
	real(kind=prec)    :: sigm,sigp,saym,sayp
	complex(kind=prec) :: ic
	!
	ommi0=8.0d-7*PI*PI
    ic=cmplx(0.0d0,1.0d0)
	ommi=ommi0/period          
    ommii=1./ommi
	is=2

    call der3po(pem,pec,pep,hzm,hzc,is,omdedzb)
    omdedzb=ic*ommii*omdedzb

    hzmc=hzm+hzc
	hymp=hym+hyp

    ! k,j
    sigxxmm=cond(1)%v(knode-1,jnode-1)
	sigyymm=cond(2)%v(knode-1,jnode-1)
	sigzzmm=cond(3)%v(knode-1,jnode-1)
	sigxymm=cond(4)%v(knode-1,jnode-1)
	sigxzmm=cond(5)%v(knode-1,jnode-1)
	sigyzmm=cond(6)%v(knode-1,jnode-1)
	if(sigxxmm.lt.1e-6.and.sigyymm.lt.1e-6.and.  &
		     sigzzmm.lt.1e-6) then
      sbymm=0.
      sbzmm=0.
      scymm=0.
      sczmm=0.
      saymm=0.
      sazmm=0.
	  sigmm=0.
	else
	  dd=sigxxmm*sigzzmm-sigxzmm*sigxzmm ! D
	  sbymm=sigyzmm/dd
	  sbzmm=sigxxmm/dd
	  scymm=sigzzmm/dd
	  sczmm=sigxzmm/dd
      saymm=sigxzmm*sbymm-sigxymm*scymm  ! B
      sazmm=sigyzmm*sbzmm-sigxymm*sczmm  !-A
	  sigmm=sigyymm+sigxymm*saymm-sigyzmm*sazmm
    end if
    ! k,j+1
    sigxxpm=cond(1)%v(knode-1,jnode)
	sigyypm=cond(2)%v(knode-1,jnode)
	sigzzpm=cond(3)%v(knode-1,jnode)
	sigxypm=cond(4)%v(knode-1,jnode)
	sigxzpm=cond(5)%v(knode-1,jnode)
	sigyzpm=cond(6)%v(knode-1,jnode)
	if(sigxxpm.lt.1e-6.and.sigyypm.lt.1e-6.and.  &
		     sigzzpm.lt.1e-6) then
      sbypm=0.
      sbzpm=0.
      scypm=0.
      sczpm=0.
      saypm=0.
      sazpm=0.
	  sigpm=0.
	else
	  dd=sigxxpm*sigzzpm-sigxzpm*sigxzpm
	  sbypm=sigyzpm/dd
	  sbzpm=sigxxpm/dd
	  scypm=sigzzpm/dd
	  sczpm=sigxzpm/dd
      saypm=sigxzpm*sbypm-sigxypm*scypm
      sazpm=sigyzpm*sbzpm-sigxypm*sczpm
	  sigpm=sigyypm+sigxypm*saypm-sigyzpm*sazpm
    end if
    ! k+1,j
    sigxxmp=cond(1)%v(knode,jnode-1)
	sigyymp=cond(2)%v(knode,jnode-1)
	sigzzmp=cond(3)%v(knode,jnode-1)
    sigxymp=cond(4)%v(knode,jnode-1)
	sigxzmp=cond(5)%v(knode,jnode-1)
	sigyzmp=cond(6)%v(knode,jnode-1)
	if(sigxxmp.lt.1e-6.and.sigyymp.lt.1e-6.and.  &
		     sigzzmp.lt.1e-6) then
      sbymp=0.
      sbzmp=0.
      scymp=0.
      sczmp=0.
      saymp=0.
      sazmp=0.
	  sigmp=0.
	else
	  dd=sigxxmp*sigzzmp-sigxzmp*sigxzmp
	  sbymp=sigyzmp/dd
	  sbzmp=sigxxmp/dd
	  scymp=sigzzmp/dd
	  sczmp=sigxzmp/dd
      saymp=sigxzmp*sbymp-sigxymp*scymp
      sazmp=sigyzmp*sbzmp-sigxymp*sczmp
	  sigmp=sigyymp+sigxymp*saymp-sigyzmp*sazmp
    end if
    ! k+1,j+1
    sigxxpp=cond(1)%v(knode,jnode)
	sigyypp=cond(2)%v(knode,jnode)
	sigzzpp=cond(3)%v(knode,jnode)
	sigxypp=cond(4)%v(knode,jnode)
	sigxzpp=cond(5)%v(knode,jnode)
	sigyzpp=cond(6)%v(knode,jnode)
	if(sigxxpp.lt.1e-6.and.sigyypp.lt.1e-6.and.  &
		     sigzzpp.lt.1e-6) then
      sbypp=0.
      sbzpp=0.
      scypp=0.
      sczpp=0.
      saypp=0.
      sazpp=0.
	  sigpp=0.
	else
	  dd=sigxxpp*sigzzpp-sigxzpp*sigxzpp
	  sbypp=sigyzpp/dd
	  sbzpp=sigxxpp/dd
	  scypp=sigzzpp/dd
	  sczpp=sigxzpp/dd
      saypp=sigxzpp*sbypp-sigxypp*scypp
      sazpp=sigyzpp*sbzpp-sigxypp*sczpp
	  sigpp=sigyypp+sigxypp*saypp-sigyzpp*sazpp
    end if

	sigm=(hym*sigmm+hyp*sigpm)/hymp
	sigp=(hym*sigmp+hyp*sigpp)/hymp
	omdedzc=-(0.5*hzm*hzc*(sigp-sigm)/hzmc)*pec
	omdedzc=omdedzc+(0.5*hzm*hzc/(hzmc*hymp))*   &
          ((sazpm-sazpp)*(phr-phc)+(sazmm-sazmp)*(phc-phl))
	saym=(hym*saymm+hyp*saypm)/hymp
	sayp=(hym*saymp+hyp*saypp)/hymp
	omdedzc=omdedzc-0.5*(hzm*sayp*(php-phc)-hzc*saym*(phc-phm))/hzmc

	omdedz=omdedzb+omdedzc

  end subroutine deredz_y


  ! *****************************************************************
  subroutine deredy_y(jnode,knode,pel,pec,per,phc,phl,phr,phm,php,  & 
          hym,hyp,hzm,hzc,nlayr,period,cond,omdedyb,omdedyc,omdedy)
    implicit none
    integer jnode,knode,nlayr
    real(kind=prec)    :: hym,hyp,hzm,hzc,period
    type(rscalar_2d)   :: cond(6)
	complex(kind=prec) :: pel,pec,per,phc,phl,phr,phm,php
	complex(kind=prec) :: omdedyb,omdedyc,omdedy
	! local variables
	integer is
	real(kind=prec)    :: ommi0,ommi,dd,ommii,hzmc,hymp
	real(kind=prec)    :: sigxxmm,sigyymm,sigzzmm,sigxymm,sigxzmm,sigyzmm
	real(kind=prec)    :: sigxxpm,sigyypm,sigzzpm,sigxypm,sigxzpm,sigyzpm
	real(kind=prec)    :: sigxxmp,sigyymp,sigzzmp,sigxymp,sigxzmp,sigyzmp
	real(kind=prec)    :: sigxxpp,sigyypp,sigzzpp,sigxypp,sigxzpp,sigyzpp
	real(kind=prec)    :: sbymm,sbzmm,scymm,sczmm,saymm,sazmm,sigmm
	real(kind=prec)    :: sbypm,sbzpm,scypm,sczpm,saypm,sazpm,sigpm
	real(kind=prec)    :: sbymp,sbzmp,scymp,sczmp,saymp,sazmp,sigmp
	real(kind=prec)    :: sbypp,sbzpp,scypp,sczpp,saypp,sazpp,sigpp
	real(kind=prec)    :: sigl,sigr,sazl,sazr
	complex(kind=prec) :: ic
	!
	ommi0=8.0d-7*PI*PI
    ic=cmplx(0.0d0,1.0d0)
	ommi=ommi0/period          
    ommii=1./ommi
	is=2

	call der3po(pel,pec,per,hym,hyp,is,omdedyb)
	omdedyb=-ic*ommii*omdedyb

    hzmc=hzm+hzc
	hymp=hym+hyp

    ! k,j
    sigxxmm=cond(1)%v(knode-1,jnode-1)
	sigyymm=cond(2)%v(knode-1,jnode-1)
	sigzzmm=cond(3)%v(knode-1,jnode-1)
	sigxymm=cond(4)%v(knode-1,jnode-1)
	sigxzmm=cond(5)%v(knode-1,jnode-1)
	sigyzmm=cond(6)%v(knode-1,jnode-1)
	if(sigxxmm.lt.1e-6.and.sigyymm.lt.1e-6.and.  &
		     sigzzmm.lt.1e-6) then
      sbymm=0.
      sbzmm=0.
      scymm=0.
      sczmm=0.
      saymm=0.
      sazmm=0.
	  sigmm=0.
	else
	  dd=sigxxmm*sigzzmm-sigxzmm*sigxzmm ! D
	  sbymm=sigyzmm/dd
	  sbzmm=sigxxmm/dd
	  scymm=sigzzmm/dd
	  sczmm=sigxzmm/dd
      saymm=sigxzmm*sbymm-sigxymm*scymm  ! B
      sazmm=sigyzmm*sbzmm-sigxymm*sczmm  !-A
	  sigmm=sigyymm+sigxymm*saymm-sigyzmm*sazmm
    end if
    ! k,j+1
    sigxxpm=cond(1)%v(knode-1,jnode)
	sigyypm=cond(2)%v(knode-1,jnode)
	sigzzpm=cond(3)%v(knode-1,jnode)
	sigxypm=cond(4)%v(knode-1,jnode)
	sigxzpm=cond(5)%v(knode-1,jnode)
	sigyzpm=cond(6)%v(knode-1,jnode)
	if(sigxxpm.lt.1e-6.and.sigyypm.lt.1e-6.and.  &
		     sigzzpm.lt.1e-6) then
      sbypm=0.
      sbzpm=0.
      scypm=0.
      sczpm=0.
      saypm=0.
      sazpm=0.
	  sigpm=0.
	else
	  dd=sigxxpm*sigzzpm-sigxzpm*sigxzpm
	  sbypm=sigyzpm/dd
	  sbzpm=sigxxpm/dd
	  scypm=sigzzpm/dd
	  sczpm=sigxzpm/dd
      saypm=sigxzpm*sbypm-sigxypm*scypm
      sazpm=sigyzpm*sbzpm-sigxypm*sczpm
	  sigpm=sigyypm+sigxypm*saypm-sigyzpm*sazpm
    end if
    ! k+1,j
    sigxxmp=cond(1)%v(knode,jnode-1)
	sigyymp=cond(2)%v(knode,jnode-1)
	sigzzmp=cond(3)%v(knode,jnode-1)
    sigxymp=cond(4)%v(knode,jnode-1)
	sigxzmp=cond(5)%v(knode,jnode-1)
	sigyzmp=cond(6)%v(knode,jnode-1)
	if(sigxxmp.lt.1e-6.and.sigyymp.lt.1e-6.and.  &
		     sigzzmp.lt.1e-6) then
      sbymp=0.
      sbzmp=0.
      scymp=0.
      sczmp=0.
      saymp=0.
      sazmp=0.
	  sigmp=0.
	else
	  dd=sigxxmp*sigzzmp-sigxzmp*sigxzmp
	  sbymp=sigyzmp/dd
	  sbzmp=sigxxmp/dd
	  scymp=sigzzmp/dd
	  sczmp=sigxzmp/dd
      saymp=sigxzmp*sbymp-sigxymp*scymp
      sazmp=sigyzmp*sbzmp-sigxymp*sczmp
	  sigmp=sigyymp+sigxymp*saymp-sigyzmp*sazmp
    end if
    ! k+1,j+1
    sigxxpp=cond(1)%v(knode,jnode)
	sigyypp=cond(2)%v(knode,jnode)
	sigzzpp=cond(3)%v(knode,jnode)
	sigxypp=cond(4)%v(knode,jnode)
	sigxzpp=cond(5)%v(knode,jnode)
	sigyzpp=cond(6)%v(knode,jnode)
	if(sigxxpp.lt.1e-6.and.sigyypp.lt.1e-6.and.  &
		     sigzzpp.lt.1e-6) then
      sbypp=0.
      sbzpp=0.
      scypp=0.
      sczpp=0.
      saypp=0.
      sazpp=0.
	  sigpp=0.
	else
	  dd=sigxxpp*sigzzpp-sigxzpp*sigxzpp
	  sbypp=sigyzpp/dd
	  sbzpp=sigxxpp/dd
	  scypp=sigzzpp/dd
	  sczpp=sigxzpp/dd
      saypp=sigxzpp*sbypp-sigxypp*scypp
      sazpp=sigyzpp*sbzpp-sigxypp*sczpp
	  sigpp=sigyypp+sigxypp*saypp-sigyzpp*sazpp
    end if
	
	sigl=(hzm*sigmm+hzc*sigmp)/hzmc
	sigr=(hzm*sigpm+hzc*sigpp)/hzmc
	omdedyc=(0.5*hym*hyp*(sigr-sigl)/hymp)*pec
	omdedyc=omdedyc+(0.5*hym*hyp/(hzmc*hymp))*   &
           ((saypp-saymp)*(php-phc)+(saypm-saymm)*(phc-phm))
	sazl=-(hzm*sazmm+hzc*sazmp)/hzmc
	sazr=-(hzm*sazpm+hzc*sazpp)/hzmc
	omdedyc=omdedyc-0.5*(hym*sazr*(phr-phc)-hyp*sazl*(phc-phl))/hymp

	omdedy=omdedyb+omdedyc
			  
  end subroutine deredy_y


  ! ***************************************************************** 
  subroutine der3po(y1,y2,y3,h12,h23,iy,der)
    implicit none
	integer iy
	real(kind=prec)    :: h12,h23
    complex(kind=prec) :: y1,y2,y3,der
    ! local variables
	real(kind=prec)    :: d12,d21
	complex(kind=prec) :: d1,d2
    !
    d1=(y2-y1)/h12
    d2=(y3-y2)/h23
    d12=h12/(h12+h23)
    d21=h23/(h12+h23)

	if(iy-2)1,2,3
1     der=(d12+d12+d21)*d1-d12*d2
      return
2     der=d21*d1+d12*d2
      return
3     der=-d21*d1+(d12+d21+d21)*d2
      return

  end subroutine der3po

  
!  ! **************************
!  subroutine czero(n,x)
!    implicit none
!	integer n
!	complex(kind=prec) :: x(n)
!    ! local variables
!	integer i
!
!	do i=1,n
!	  x(i)=C_ZERO
!	end do
!
!  end subroutine czero

end module Fwd2DANI