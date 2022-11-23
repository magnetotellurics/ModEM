!    routines which define mappings between the "natural" representation
!        of conductivity/resistivity on the model grid and the formal
!        model parameter structure: in general (allowing for a non-linear
!        mapping, as will occur when log condutivity is used for inversion;
!        in the implementation here,we combine the linear and non-linear
!        mappings in a single routine)
!        one needs the non-linear mapping, a linearized mapping, and the
!        adjoint of the linearized mapping.  If data functionals depend
!        on conductivity or resistivity additional routines are needed--
!        i.e., a real function (SigC here), with input arguments a conductivity
!        parameter and cell indicies, which returns the conductivity of
!        this cell (or resistivity), as well as what is essentially the
!        adjoint of the linearization of this function (sort of ...
!         see below; this is called QtoModelParam here.)
!       Note that the functionality is always needed, but names can be
!        flexible ... these routines are called by the interpolation
!        routines in EMfieldInterp* (for the latter two only), by the forward
!        modeling codes (the non-linear mapping only), and by SolverSens/DataSens
!        (the linearized mappings and adjoint) and thus names/interfaces
!        need to be kept consistent with what is used in these routines.
!        Note that for 2D MT different mapping routines are used for
!        TE and TM modes, as the natural representations for the two cases
!        are different (conductivity vs. resistivity)

   !**********************************************************************
   subroutine ModelParamToTensor(m,cCond,cRho,paramType,grid,AirRho)

     type(modelParam_t), intent(in)        :: m
     type(rscalar), intent(inout)	       :: cCond(6),cRho(6)
     character(80), intent(out), optional  :: paramType
     type(grid_t), intent(out), optional :: grid
     real(kind=prec), intent(out), optional :: AirRho
     ! local variables
     integer :: ia
     type(rscalar) :: tempRho(6)

     do ia = 1,6
     
       if(cRho(ia)%allocated) then
          if((cRho(ia)%Ny .ne. m%Ny).or. (cRho(ia)%Nx .ne. m%Nx) .or. &
		      (cRho(ia)%Nz .ne. m%grid%Nz)) then
             call deall_rscalar(cRho(ia))
             call create_rscalar(m%grid,cRho(ia),CENTER)
          endif
       else
          call create_rscalar(m%grid,cRho(ia),CENTER)
       endif
       
       if(cCond(ia)%allocated) then
          if((cCond(ia)%Ny .ne. m%Ny).or. (cCond(ia)%Nx .ne. m%Nx) .or. &
		      (cCond(ia)%Nz .ne. m%grid%Nz)) then
             call deall_rscalar(cCond(ia))
             call create_rscalar(m%grid,cCond(ia),CENTER)
          endif
       else
          call create_rscalar(m%grid,cCond(ia),CENTER)
       endif
       
       if(tempRho(ia)%allocated) then
          if((tempRho(ia)%Ny .ne. m%Ny).or. (tempRho(ia)%Nx .ne. m%Nx) .or. &
		      (tempRho(ia)%Nz .ne. m%grid%Nz)) then
             call deall_rscalar(tempRho(ia))
             call create_rscalar(m%grid,tempRho(ia),CELL_EARTH)
          endif
       else
          call create_rscalar(m%grid,tempRho(ia),CELL_EARTH)
       endif
       
     enddo            
     
     if(m%paramType .eq. LOGE) then
        do ia = 1,3
          cRho(ia)%v(:,:,1:m%grid%NzAir) = exp(m%AirRho)
          cCond(ia)%v(:,:,1:m%grid%NzAir) = exp(-m%AirRho)
          tempRho(ia)%v = exp(m%cellRho(ia)%v)
        enddo
     else
        do ia = 1,3
          cRho(ia)%v(:,:,1:m%grid%NzAir) = m%AirRho
          cCond(ia)%v(:,:,1:m%grid%NzAir) = ONE/m%AirRho
          tempRho(ia)%v = m%cellRho(ia)%v
        enddo
     endif
     
     do ia = 4,6
       cRho(ia)%v(:,:,1:m%grid%NzAir) = R_ZERO
       cCond(ia)%v(:,:,1:m%grid%NzAir) = R_ZERO
       tempRho(ia)%v = m%cellRho(ia)%v
     enddo

     call TensorCal(m%grid,tempRho,cCond,cRho)
     
     if (present(paramType)) then
        paramType = m%paramType
     end if

     if (present(grid)) then
        if (grid%allocated) then
           call deall_grid(grid)
        end if
        grid = m%grid
     end if

     if (present(AirRho)) then
        AirRho = m%AirRho
     end if
     
     do ia = 1,6
       call deall_rscalar(tempRho(ia)) 
     enddo

   end subroutine ModelParamToTensor
   
   
   !**********************************************************************
   subroutine SensModelParamToTensor(m,dxx,dxy,dxz,dyy,dyz,dzz)
   
     type(modelParam_t), intent(in) :: m
     type(rscalar), intent(inout) :: dxx(6),dxy(6),dxz(6),dyy(6),dyz(6),dzz(6)
     ! local variables
	   integer :: ia
     character(80) :: paramType = 'LINEAR'  !!!!!	   
	   type(rscalar) :: tempRho(6)
	   
     do ia=1,6

       if(dxx(ia)%allocated) then
          if((dxx(ia)%Ny .ne. m%Ny).or. (dxx(ia)%Nx .ne. m%Nx) .or. &
					(dxx(ia)%Nz .ne. m%grid%NzEarth)) then
             call deall_rscalar(dxx(ia))
             call create_rscalar(m%grid,dxx(ia),CELL_EARTH)
          endif
       else
          call create_rscalar(m%grid,dxx(ia),CELL_EARTH)
       endif

       if(dxy(ia)%allocated) then
          if((dxy(ia)%Ny .ne. m%Ny).or. (dxy(ia)%Nx .ne. m%Nx) .or. &
					(dxy(ia)%Nz .ne. m%grid%NzEarth)) then
             call deall_rscalar(dxy(ia))
             call create_rscalar(m%grid,dxy(ia),CELL_EARTH)
          endif
       else
          call create_rscalar(m%grid,dxy(ia),CELL_EARTH)
       endif

       if(dxz(ia)%allocated) then
          if((dxz(ia)%Ny .ne. m%Ny).or. (dxz(ia)%Nx .ne. m%Nx) .or. &
					(dxz(ia)%Nz .ne. m%grid%NzEarth)) then
             call deall_rscalar(dxz(ia))
             call create_rscalar(m%grid,dxz(ia),CELL_EARTH)
          endif
       else
          call create_rscalar(m%grid,dxz(ia),CELL_EARTH)
       endif

       if(dyy(ia)%allocated) then
          if((dyy(ia)%Ny .ne. m%Ny).or. (dyy(ia)%Nx .ne. m%Nx) .or. &
					(dyy(ia)%Nz .ne. m%grid%NzEarth)) then
             call deall_rscalar(dyy(ia))
             call create_rscalar(m%grid,dyy(ia),CELL_EARTH)
          endif
       else
          call create_rscalar(m%grid,dyy(ia),CELL_EARTH)
       endif

       if(dyz(ia)%allocated) then
          if((dyz(ia)%Ny .ne. m%Ny).or. (dyz(ia)%Nx .ne. m%Nx) .or. &
					(dyz(ia)%Nz .ne. m%grid%NzEarth)) then
             call deall_rscalar(dyz(ia))
             call create_rscalar(m%grid,dyz(ia),CELL_EARTH)
          endif
       else
          call create_rscalar(m%grid,dyz(ia),CELL_EARTH)
       endif
       
       if(dzz(ia)%allocated) then
          if((dzz(ia)%Ny .ne. m%Ny).or. (dzz(ia)%Nx .ne. m%Nx) .or. &
					(dzz(ia)%Nz .ne. m%grid%NzEarth)) then
             call deall_rscalar(dzz(ia))
             call create_rscalar(m%grid,dzz(ia),CELL_EARTH)
          endif
       else
          call create_rscalar(m%grid,dzz(ia),CELL_EARTH)
       endif

       if(tempRho(ia)%allocated) then
          if((tempRho(ia)%Ny .ne. m%Ny).or. (tempRho(ia)%Nx .ne. m%Nx) .or. &
					(tempRho(ia)%Nz .ne. m%grid%NzEarth)) then
             call deall_rscalar(tempRho(ia))
             call create_rscalar(m%grid,tempRho(ia),CELL_EARTH)
          endif
       else
          call create_rscalar(m%grid,tempRho(ia),CELL_EARTH)
       endif
              
		 end do
		 
		 call getValue_modelParam(m,paramType,tempRho)  
     
     call AniSensCal(m%grid,tempRho,dxx,dxy,dxz,dyy,dyz,dzz)
     
     do ia = 1,6
       call deall_rscalar(tempRho(ia)) 
     enddo
          
   end subroutine SensModelParamToTensor
   

   ! subrountine to calculate the tensor conductivity and resistivity in anisotropy
   subroutine TensorCal(grid,RhoIn,CondOut,RhoOut)
     ! modified from the code of Pek and Santos (2002)
     type(grid_t), intent(in)     :: grid
     type(rscalar), intent(inout)	 :: CondOut(6),RhoOut(6)
	   type(rscalar), intent(in)     :: RhoIn(6)

	   !local variables
	   integer i, j, k, km, nx, ny, nz, nzAir
	   real(kind=prec) :: ron1,ron2,ron3,rstr,rdip,r3
	   real(kind=prec) :: sps,cps,sth,cth,sfi,cfi,pom1,pom2,pom3
	   real(kind=prec) :: c2ps,s2ps,c2th,s2th,csps,csth 

     nx = grid%Nx
     ny = grid%Ny
     nz = grid%Nz
	   nzAir = grid%nzAir

!	 if(CondOut(1)%Nx.ne.RhoIn(1)%Nx) then
!       call errStop('Dimensions do not match!')
!	 end if

     ! In the Earth
     do k = nzAir+1,nz
       km = k - nzAir
	     do j = 1,ny
	       do i = 1,nx

           ! axial cond
           ron1 = ONE/RhoIn(1)%v(i,j,km) ! ONE/100.0 !
		       ron2 = ONE/RhoIn(2)%v(i,j,km)
		       ron3 = ONE/RhoIn(3)%v(i,j,km)
           ! anisotropy angles
		       rstr = RhoIn(4)%v(i,j,km)
		       rdip = RhoIn(5)%v(i,j,km)
		       r3 = RhoIn(6)%v(i,j,km)
           
           sps=dsin(rstr)
           cps=dcos(rstr)
           sth=dsin(rdip)
           cth=dcos(rdip)
           sfi=dsin(r3)
           cfi=dcos(r3)

           c2ps=cps*cps
           s2ps=sps*sps
           c2th=cth*cth
           s2th=sth*sth
           csps=cps*sps
           csth=cth*sth

           pom1=ron1*cfi*cfi+ron2*sfi*sfi
           pom2=ron1*sfi*sfi+ron2*cfi*cfi
           pom3=(ron1-ron2)*sfi*cfi

           CondOut(1)%v(i,j,k) = pom1*c2ps+pom2*s2ps*c2th-  &
                       2.*pom3*cth*csps+ron3*s2th*s2ps  !xx

           CondOut(2)%v(i,j,k) = pom1*s2ps+pom2*c2ps*c2th+  &
                       2.*pom3*cth*csps+ron3*s2th*c2ps  !yy

           CondOut(3)%v(i,j,k) = pom2*s2th+ron3*c2th  !zz

           CondOut(4)%v(i,j,k) = pom1*csps-pom2*c2th*csps+  &
                       pom3*cth*(c2ps-s2ps)-ron3*s2th*csps  !xy

           CondOut(5)%v(i,j,k) = -pom2*csth*sps+pom3*sth*cps+  &
                       ron3*csth*sps  !xz

           CondOut(6)%v(i,j,k) = pom2*csth*cps+pom3*sth*sps-  &
                       ron3*csth*cps  !yz

!       Debugging
!		write(1111,*) i,j,k
!		write(1111,*) CondOut(1)%v(i,j,k),CondOut(4)%v(i,j,k),CondOut(5)%v(i,j,k)
!		write(1111,*) CondOut(4)%v(i,j,k),CondOut(2)%v(i,j,k),CondOut(6)%v(i,j,k)
!		write(1111,*) CondOut(5)%v(i,j,k),CondOut(6)%v(i,j,k),CondOut(3)%v(i,j,k)

           ! axial rho
           ron1 = RhoIn(1)%v(i,j,km) ! 100.0d0 !
		       ron2 = RhoIn(2)%v(i,j,km)
		       ron3 = RhoIn(3)%v(i,j,km)

           pom1=ron1*cfi*cfi+ron2*sfi*sfi
           pom2=ron1*sfi*sfi+ron2*cfi*cfi
           pom3=(ron1-ron2)*sfi*cfi

           RhoOut(1)%v(i,j,k) = pom1*c2ps+pom2*s2ps*c2th-  &
                       2.*pom3*cth*csps+ron3*s2th*s2ps  !xx

           RhoOut(2)%v(i,j,k) = pom1*s2ps+pom2*c2ps*c2th+  &
                       2.*pom3*cth*csps+ron3*s2th*c2ps  !yy

           RhoOut(3)%v(i,j,k) = pom2*s2th+ron3*c2th  !zz

           RhoOut(4)%v(i,j,k) = pom1*csps-pom2*c2th*csps+  &
                       pom3*cth*(c2ps-s2ps)-ron3*s2th*csps  !xy

           RhoOut(5)%v(i,j,k) = -pom2*csth*sps+pom3*sth*cps+  &
                       ron3*csth*sps  !xz

           RhoOut(6)%v(i,j,k) = pom2*csth*cps+pom3*sth*sps-  &
                       ron3*csth*cps  !yz

!       Debugging
!		write(1111,*) i,j,k
!		write(1111,*) RhoOut(1)%v(i,j,k),RhoOut(4)%v(i,j,k),RhoOut(5)%v(i,j,k)
!		write(1111,*) RhoOut(4)%v(i,j,k),RhoOut(2)%v(i,j,k),RhoOut(6)%v(i,j,k)
!		write(1111,*) RhoOut(5)%v(i,j,k),RhoOut(6)%v(i,j,k),RhoOut(3)%v(i,j,k)	
	         
         end do
	   end do
	 end do


  end subroutine TensorCal
  
  
  ! subroutine to calculate the derivatives of resistivity tensor w.r.t axial res and angles 
   subroutine AniSensCal(grid,RhoIn,dxx,dxy,dxz,dyy,dyz,dzz)
     ! modified from the code of Pek and Santos (2002)   
     type(grid_t), intent(in)      :: grid
     type(rscalar), intent(in)     :: RhoIn(6)
     type(rscalar), intent(inout)  :: dxx(6),dxy(6),dxz(6),dyy(6),dyz(6),dzz(6)
     ! local variables
     integer i,j,k
     integer Nx,Ny,Nz
     real(kind=prec) :: sgpx,sgpy,sgpz,rstr,rdip,rsla
     real(kind=prec) :: asg,bsg,csg,dsg,dbsgpz,dsgpxy
     real(kind=prec) :: cstr,sstr,cdip,sdip,csla,ssla
     real(kind=prec) :: c2str,s2str,c2dip,s2dip,c2sla,s2sla
      	     
     Nx = grid%Nx
     Ny = grid%Ny
     Nz = grid%NzEarth
     
!	 if(CondOut(1)%Nx.ne.dxx(1)%Nx) then
!       call errStop('Dimensions do not match!') !debug
!	 end if     
 
 	 
     ! In the Earth 
     do k = 1, Nz
       do j = 1, Ny
         do i = 1, Nx
             
  
					 ! axial resistivities    
           sgpx=RhoIn(1)%v(i,j,k)
           sgpy=RhoIn(2)%v(i,j,k)
           sgpz=RhoIn(3)%v(i,j,k)
           ! anisotropy angles in radian
           rstr=RhoIn(4)%v(i,j,k)
           rdip=RhoIn(5)%v(i,j,k)
           rsla=RhoIn(6)%v(i,j,k)
!
           cstr=dcos(rstr)
           sstr=dsin(rstr)
           cdip=dcos(rdip)
           sdip=dsin(rdip)
           csla=dcos(rsla)
           ssla=dsin(rsla)
           c2str=dcos(2.d0*rstr)
           s2str=dsin(2.d0*rstr)
           c2dip=dcos(2.d0*rdip)
           s2dip=dsin(2.d0*rdip)
           c2sla=dcos(2.d0*rsla)
           s2sla=dsin(2.d0*rsla)
!
           dsgpxy=sgpx-sgpy
           asg=sgpx*csla**2.d0+sgpy*ssla**2.d0
           bsg=sgpx*ssla**2.d0+sgpy*csla**2.d0
           csg=bsg*cdip**2.d0+sgpz*sdip**2.d0
           dsg=bsg*sdip**2.d0+sgpz*cdip**2.d0
           dbsgpz=bsg-sgpz  


!
!> sigxx and its derivatives
!
					dxx(1)%v(i,j,k)=(csla*cstr-ssla*cdip*sstr)**2.d0
					dxx(2)%v(i,j,k)=(ssla*cstr+csla*cdip*sstr)**2.d0
					dxx(3)%v(i,j,k)=(sdip*sstr)**2.d0
					dxx(4)%v(i,j,k)=(-asg+csg)*s2str-dsgpxy*cdip*s2sla*c2str
					dxx(5)%v(i,j,k)=-dbsgpz*s2dip*sstr**2.d0+  &
                        0.5d0*dsgpxy*sdip*s2sla*s2str
					dxx(6)%v(i,j,k)=dsgpxy*(((sstr*cdip)**2.d0-cstr**2.d0)*s2sla-  &
                        cdip*s2str*c2sla)
!
!> sigxy and its derivatives
!
					dxy(1)%v(i,j,k)=0.5d0*((csla**2.d0-(ssla*cdip)**2.d0)*s2str+  &
                        cdip*c2str*s2sla)
					dxy(2)%v(i,j,k)=0.5d0*((ssla**2.d0-(csla*cdip)**2.d0)*s2str-  &
                        cdip*c2str*s2sla)
					dxy(3)%v(i,j,k)=-0.5d0*sdip**2.d0*s2str
					dxy(4)%v(i,j,k)=(asg-csg)*c2str-dsgpxy*cdip*s2str*s2sla
					dxy(5)%v(i,j,k)=0.5d0*dbsgpz*s2dip*s2str-  &
                        0.5d0*dsgpxy*sdip*c2str*s2sla
					dxy(6)%v(i,j,k)=0.5d0*dsgpxy*(2.d0*cdip*c2str*c2sla-  &
                       (1.d0+cdip**2.d0)*s2str*s2sla)           
!
!> sigxz and its derivatives
!
					dxz(1)%v(i,j,k)=(cstr*csla-sstr*cdip*ssla)*ssla*sdip
					dxz(2)%v(i,j,k)=-(cstr*ssla+sstr*cdip*csla)*csla*sdip
					dxz(3)%v(i,j,k)=0.5d0*s2dip*sstr
					dxz(4)%v(i,j,k)=-0.5d0*dsgpxy*sdip*sstr*s2sla-  &
                        0.5d0*dbsgpz*s2dip*cstr
					dxz(5)%v(i,j,k)=0.5d0*dsgpxy*cdip*cstr*s2sla-  &
                        dbsgpz*c2dip*sstr
					dxz(6)%v(i,j,k)=0.5d0*dsgpxy*(2.d0*sdip*cstr*c2sla-s2dip*sstr*s2sla)
!
!> sigyy and its derivatives
!
					dyy(1)%v(i,j,k)=(sstr*csla+cdip*cstr*ssla)**2.d0
					dyy(2)%v(i,j,k)=(sstr*ssla-cdip*cstr*csla)**2.d0
					dyy(3)%v(i,j,k)=(cstr*sdip)**2.d0
					dyy(4)%v(i,j,k)=(asg-csg)*s2str+dsgpxy*cdip*c2str*s2sla
					dyy(5)%v(i,j,k)=-dbsgpz*s2dip*cstr**2.d0-  &
                        0.5d0*dsgpxy*sdip*s2str*s2sla
					dyy(6)%v(i,j,k)=dsgpxy*(cdip*s2str*c2sla+((cdip*cstr)**2.d0-  &
                        sstr**2.d0)*s2sla)
!
!> sigyz and its derivatives
!
					dyz(1)%v(i,j,k)=(sstr*csla+cdip*cstr*ssla)*sdip*ssla
					dyz(2)%v(i,j,k)=-(sstr*ssla-cdip*cstr*csla)*sdip*csla
					dyz(3)%v(i,j,k)=-0.5d0*s2dip*cstr
					dyz(4)%v(i,j,k)=0.5d0*dsgpxy*sdip*cstr*s2sla-  &
                        0.5d0*dbsgpz*s2dip*sstr
					dyz(5)%v(i,j,k)=0.5d0*dsgpxy*cdip*sstr*s2sla+  &
                        dbsgpz*c2dip*cstr
					dyz(6)%v(i,j,k)=0.5d0*dsgpxy*(2.d0*sdip*sstr*c2sla+s2dip*cstr*s2sla)
!
!> sigzz and its derivatives
!
					dzz(1)%v(i,j,k)=(sdip*ssla)**2.d0
					dzz(2)%v(i,j,k)=(sdip*csla)**2.d0
					dzz(3)%v(i,j,k)=cdip**2.d0
					dzz(4)%v(i,j,k)=0.d0
					dzz(5)%v(i,j,k)=dbsgpz*s2dip
					dzz(6)%v(i,j,k)=dsgpxy*s2sla*sdip**2.d0         
                  
        end do         
      end do
    end do
     
     
  end subroutine AniSensCal   