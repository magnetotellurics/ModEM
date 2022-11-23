! *****************************************************************************
! This module implements boundary conditions and maybe first guess.

module boundary_Pek
  use sg_vector
  use sg_scalar
  use sg_boundary
  use Fwd2DANI

Contains


  !****************************************************************************
  ! Generates boundary conditions and initial guess for the iterative solution.
  ! X-strike version.
  subroutine BC_x0_Pek(imode,period,grid3D,Cond3D,E0,BC)
    implicit none
    !  Input mode, period
    integer, intent(in)		:: imode
    real(kind=prec)	:: period
    !  Input 3D grid
    type(grid_t), intent(in)	:: grid3D
    !  Input 3D conductivity in cells
    type(rscalar), intent(in)	:: Cond3D(3)

    ! Output electric field first guess (for iterative solver)
    type(cvector), intent(inout)	:: E0
    ! Output boundary conditions
    type(cboundary), intent(inout)	:: BC

    ! local variables
    ! 2D grid definitions
    type(grid2d_t)		        :: grid2D
	integer				:: nl,nh,nv,nl1,nh1,nv1
    integer				:: iSlice,ih,iz,nSlice,nzap1,ia,iof
	integer             :: ll_2d,ln,ln1,np1_2d,np2_2d,np3_2d
    integer				:: IER,istat
    complex(kind=prec),allocatable :: hx2d(:),hy2d(:),hz2d(:)

    ! Array for 2D Conductivity
    type(rscalar_2d)	:: Cond2D(6)

    nl = grid3D%nx
    ! nh is number of horizontal grid cells
    nh = grid3D%ny
	grid2D%Ny = nh
	nSlice = grid3D%nx
    allocate(grid2D%Dy(nh))
    do ih = 1,nh
       grid2D%Dy(ih) = grid3D%dy(ih)
!	   write(1010,*) grid2D%Dy(ih)  ! debug
    enddo
    !  nv is number of vertical levels in 2D
    nv = grid3D%nz  ! nz=nzEarth+nzAir
    grid2D%nz = nv
    grid2D%nza = grid3D%nzAir
	nzap1 = grid2D%nza+1
    allocate(grid2D%Dz(nv))
    do iz = 1,nv
       grid2D%Dz(iz) = grid3D%Dz(iz)
!	   write(1010,*) grid2D%Dz(iz)  ! debug
    enddo

    nl1=nl+1
	nh1=nh+1
	nv1=nv+1
    ll_2d=max(nl+1,nh+1)
    ln=ll_2d*(nv+2)
	ln1=ll_2d*(nv+3)
	np1_2d=ln1*3
    np2_2d=2*(ll_2d-2)*nv
	np3_2d=2*(ll_2d-2)+3

  
  do iSlice = 1,nSlice+1  
    call slice(nl,nh,nv,Cond3D,iSlice,'y',Cond2D)
!    !debug
!    do ia = 1,6
!      iof = 1088 + ia
!      call write_rscalar_2d(Cond2D(ia),iof)
!    enddo
!    pause
!    !debug    
	  allocate(hx2d(ln),hy2d(ln),hz2d(ln))
	  hx2d=C_ZERO
	  hy2d=C_ZERO
	  hz2d=C_ZERO
    call d2atbp(np2_2d,np3_2d,nh1,nv1,nzap1,grid2D%Dy,  &
         grid2D%Dz,nv,period,imode,Cond2D,hx2d,hy2d,hz2d,'x') 
    call copy_h2e0(nl1,nh1,nv1,hx2d,hy2d,hz2d,E0,iSlice,'x')							 
	  deallocate(hx2d,hy2d,hz2d)	
    do ia = 1,6
      call deall_rscalar_2d(Cond2D(ia))
    enddo
    
	enddo
	
	  Call getBC(E0,BC)

    deallocate(grid2D%Dy,grid2D%Dz)
    
    
  end subroutine BC_x0_Pek

  !****************************************************************************
  ! Generates boundary conditions and initial guess for the iterative solution.
  ! Y-strike version.
  subroutine BC_y0_Pek(imode,period,grid3D,Cond3D,E0,BC)
    implicit none
    !  Input mode, period
    integer, intent(in)		:: imode
    real(kind=prec)	:: period
    !  Input 3D grid
    type(grid_t), intent(in)	:: grid3D
    !  Input 3D conductivity in cells
    type(rscalar), intent(in)	:: Cond3D(3)

    ! Output electric field first guess (for iterative solver)
    type(cvector), intent(inout)	:: E0
    ! Output boundary conditions
    type(cboundary), intent(inout)	:: BC

    ! local variables
    ! 2D grid definitions
    type(grid2d_t)		        :: grid2D
	integer             :: nh,nm,nv,nh1,nm1,nv1
    integer				:: iSlice,ih,iz,nSlice,nzap1,ia
	integer             :: ll_2d,ln,ln1,np1_2d,np2_2d,np3_2d
    integer				:: IER,istat
    complex(kind=prec),allocatable :: hx2d(:),hy2d(:),hz2d(:)


    ! Array for 2D Conductivity
    type(rscalar_2d)	:: Cond2D(6)

    nm = grid3D%ny
    ! nh is number of horizontal grid cells
    nh = grid3D%nx
	grid2D%Ny = nh
	nSlice = grid3D%ny
    allocate(grid2D%Dy(nh))
    do ih = 1,nh
       grid2D%Dy(ih) = grid3D%dx(ih)
    enddo
    !  nv is number of vertical levels in 2D
    nv = grid3D%nz
    grid2D%nz = nv
    grid2D%nza = grid3D%nzAir
	nzap1 = grid2D%nza+1
    allocate(grid2D%Dz(nv))
    do iz = 1,nv
       grid2D%Dz(iz) = grid3D%Dz(iz)
    enddo

    nh1=nh+1
	nm1=nm+1
	nv1=nv+1
    ll_2d=max(nm+1,nh+1)
    ln=ll_2d*(nv+2)
	ln1=ll_2d*(nv+3)
	np1_2d=ln1*3
    np2_2d=2*(ll_2d-2)*nv
	np3_2d=2*(ll_2d-2)+3

  do iSlice = 1,nSlice+1
    call slice(nh,nm,nv,Cond3D,iSlice,'x',Cond2D) ! direction
    allocate(hx2d(ln),hy2d(ln),hz2d(ln))
	  hx2d=C_ZERO
	  hy2d=C_ZERO
	  hz2d=C_ZERO
    call d2atbp(np2_2d,np3_2d,nh1,nv1,nzap1,grid2D%Dy,  &
         grid2D%Dz,nv,period,imode,Cond2D,hx2d,hy2d,hz2d,'y') ! strike  
    call copy_h2e0(nh1,nm1,nv1,hx2d,hy2d,hz2d,E0,iSlice,'y') ! strike
	  deallocate(hx2d,hy2d,hz2d)
	  do ia = 1,6
      call deall_rscalar_2d(Cond2D(ia))
    enddo
	enddo
    
    Call getBC(E0,BC)

    deallocate(grid2D%Dy,grid2D%Dz)

  end subroutine BC_y0_Pek


  subroutine slice(l,m,nlayr,Cond3D,isl,direction,Cond2D)
!***************************************************************
    implicit none
    integer l,m,isl,nlayr
    type(rscalar), intent(in)	:: Cond3D(3)
	character*1 direction
	type(rscalar_2d), intent(inout) :: Cond2D(6)
	! local variables
    integer i,j,k,ia

    do ia = 1,3
      if(direction.eq.'x')then
        if(isl.gt.(m+1) .or. isl.lt.1)   &  
           stop 'x outside model in subroutine slice'

        call create_rscalar_2d(l, nlayr, Cond2D(ia)) 
        if(isl.eq.(m+1)) then
          do i=1,l
            do k=1,nlayr
              Cond2D(ia)%v(k,i)=Cond3D(ia)%v(i,isl-1,k)
            end do
          end do
        else
          do i=1,l
            do k=1,nlayr
              Cond2D(ia)%v(k,i)=Cond3D(ia)%v(i,isl,k)
            end do
          end do 
        endif 
      elseif(direction.eq.'y')then
        if(isl.gt.(l+1) .or. isl.lt.1)  &   
           stop 'y outside model in subroutine slice'

		    call create_rscalar_2d(m, nlayr, Cond2D(ia))
		    if(isl.eq.(l+1)) then
          do j=1,m
            do k=1,nlayr
              Cond2D(ia)%v(k,j)=Cond3D(ia)%v(isl-1,j,k)
            end do
          end do
        else
          do j=1,m
            do k=1,nlayr
              Cond2D(ia)%v(k,j)=Cond3D(ia)%v(isl,j,k)
            end do
          end do
        endif
      else
        stop 'ERROR in slice direction. Must be x or y'
      endif
    end do
 
    do ia = 4,6
      if(direction.eq.'x')then
        if(isl.gt.(m+1) .or. isl.lt.1)   &  
           stop 'x outside model in subroutine slice'

        call create_rscalar_2d(l, nlayr, Cond2D(ia))
        do i=1,l
          do k=1,nlayr
            Cond2D(ia)%v(k,i)=R_ZERO
          end do
        end do
      elseif(direction.eq.'y')then
        if(isl.gt.(l+1) .or. isl.lt.1)  &   
           stop 'y outside model in subroutine slice'

	    	call create_rscalar_2d(m, nlayr, Cond2D(ia))
        do j=1,m
          do k=1,nlayr
            Cond2D(ia)%v(k,j)=R_ZERO
!			write(1010,*) Cond2D(ia)%v(k,j)  ! debug
          end do
        end do
      else
        stop 'ERROR in slice direction. Must be x or y'
      endif
    end do
       
  end subroutine slice


  ! copy 2D boundary magnetic fields to cboundary
  subroutine copy_h2cb(nx,ny,nz,hx,hy,hz,BC,strike,flag)
    ! ------- Kong, W.X., 2017-4-23 -------
    implicit none
	integer nx,ny,nz
	complex(kind=prec), intent(in)  :: hx(nz,*),hy(nz,*),hz(nz,*)
	type(cboundary), intent(inout)	:: BC
	character*1 strike
	character*3 flag
	! local variables
	integer i,j,k
	!
	if(strike.eq.'x') then

      do j=1,ny ! actually, ny+1
	    do k=1,nz ! actually, nz+1
		  if(flag.eq.'Max') then
		    if(j.ne.ny) then
		      BC%yXMax(j,k)=(hy(k,j)+hy(k,j+1))/2.0d0 ! E%yXMax(ny,nz+1)
			  !write(101,*) BC%yXMax(j,k) !debug
			endif
			if(k.ne.nz) then
              BC%zXMax(j,k)=(hz(k,j)+hz(k+1,j))/2.0d0 ! E%zXMax(ny+1,nz)
			  !write(102,*) BC%zXMax(j,k) !debug
			endif
		  elseif(flag.eq.'Min') then 
		    if(j.ne.ny) then
		      BC%yXMin(j,k)=(hy(k,j)+hy(k,j+1))/2.0d0 ! E%yXMin(ny,nz+1) 
			  !write(103,*) BC%yXMin(j,k) !debug
			endif
			if(k.ne.nz) then
              BC%zXMin(j,k)=(hz(k,j)+hz(k+1,j))/2.0d0 ! E%zXMin(ny+1,nz)
			  !write(104,*) BC%zXMin(j,k) !debug
			endif			  
		  else
		    call errStop('Flag is not compatible')
		  endif   
		end do
	  end do

	elseif(strike.eq.'y') then

      do i=1,nx ! actually, nx+1
	    do k=1,nz ! actually, nz+1
		  if(flag.eq.'Max') then
		    if(i.ne.nx) then
		      BC%xYMax(i,k)=(hx(k,i)+hx(k,i+1))/2.0d0 ! E%xYMax(nx,nz+1)
			  !write(105,*) BC%xYMax(i,k) !debug 
			endif
			if(k.ne.nz) then
              BC%zYMax(i,k)=(hz(k,i)+hz(k+1,i))/2.0d0 ! E%zYMax(nx+1,nz)
			  !write(106,*) BC%zYMax(i,k) !debug
			endif
		  elseif(flag.eq.'Min') then 
		    if(i.ne.nx) then
		      BC%xYMin(i,k)=(hx(k,i)+hx(k,i+1))/2.0d0 ! E%xYMin(nx,nz+1)
			  !write(107,*) BC%xYMin(i,k) !debug
			endif
			if(k.ne.nz) then
              BC%zYMin(i,k)=(hz(k,i)+hz(k+1,i))/2.0d0 ! E%zYMin(nx+1,nz)
			  !write(108,*) BC%zYMin(i,k) !debug
			endif					  
		  else
		    call errStop('Flag is not compatible')
		  endif   
		end do
	  end do

	else

      call errStop('Strike must be x or y')

	endif

  end subroutine copy_h2cb


  ! copy 2D boundary electric fields to E0
  subroutine copy_h2e0(nx,ny,nz,hx,hy,hz,E0,iSlice,strike)
    ! ------- Kong, W.X., 2017-4-23 -------
    implicit none
	integer nx,ny,nz
	complex(kind=prec), intent(in)  :: hx(nz,*),hy(nz,*),hz(nz,*)
	type(cvector), intent(inout)	:: E0
	! local variables
	integer i,j,k,iSlice
	character*1 strike
	!
	if(strike.eq.'x') then
	  if(iSlice.ne.nx) then
      do j = 1,ny
        do k = 1,nz
          E0%x(iSlice,j,k) = hx(k,j)
        enddo
      enddo
    endif
    do j = 1,ny-1
      do k = 1,nz
        E0%y(iSlice,j,k) = (hy(k,j)+hy(k,j+1))/2.0d0
      enddo
    enddo
    do j = 1,ny
      do k = 1,nz-1
        E0%z(iSlice,j,k) = (hz(k,j)+hz(k+1,j))/2.0d0
      enddo
    enddo
  elseif(strike.eq.'y') then
  	  if(iSlice.ne.ny) then
      do i = 1,nx
        do k = 1,nz
          E0%y(i,iSlice,k) = hy(k,i)
        enddo
      enddo
    endif
    do i = 1,nx-1
      do k = 1,nz
        E0%x(i,iSlice,k) = (hx(k,i)+hx(k,i+1))/2.0d0
      enddo
    enddo
    do i = 1,nx
      do k = 1,nz-1
        E0%z(i,iSlice,k) = (hz(k,i)+hz(k+1,i))/2.0d0
      enddo
    enddo
  else

      call errStop('Strike must be x or y')

	endif
	
  end subroutine copy_h2e0
  
  !****************************************************************************
  ! create_rscalar creates variable of derived type rscalar,
  ! using grid definition in structure "grid" ;
  ! allocates memory in v component array
  ! gridType is a character string to describe intended usage
  subroutine create_rscalar_2d(m, n, E)

    implicit none
	integer m,n
    type (rscalar_2d), intent(inout)      :: E
    ! local variables
    integer                            :: status,ny,nz


    if(E%allocated) then
       ! first deallocate memory for v
       deallocate(E%v, STAT=status)
    end if

    ! Grid dimensions
    ny = m
    nz = n

    E%ny = ny
    E%nz = nz

    ! allocate memory for v
    ! E%allocated will be true if all allocations succeed

    allocate(E%v(nz,ny), STAT=status)

    E%allocated = status .EQ. 0

    if (E%allocated) then
       E%v = R_ZERO
    end if

  end subroutine create_rscalar_2d
  
  
  subroutine deall_rscalar_2d(E)
    implicit none
    type (rscalar_2d)  :: E
    integer	    :: status

    if(E%allocated) then
       ! deallocate memory for v
       deallocate(E%v,STAT=status)
    end if

    E%ny = 0
    E%nz = 0
    E%allocated = .false.
    
  end subroutine deall_rscalar_2d

  !*****************************************        
  complex(kind=prec) function dimped(rho,omega)
      implicit none

      real(kind=prec) :: rho,omega
      real(kind=prec) :: period
      real(kind=prec) :: c3
      complex(kind=prec) :: c1mI

      c3 = PI*MU_0
      c1mI = cmplx(1.d0,-1.d0)
      period = 2 * PI / omega
      dimped = c1mI*sqrt(c3*rho/period)
  end function dimped  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_rscalar_2d(E,iof)

    implicit none
	integer iof
    type (rscalar_2d), intent(in)      :: E
    ! local variables
    integer j,k,ny,nz,nair

    ! Grid dimensions
    ny = E%ny
    nz = E%nz
    nair = 10
    do j = 1,ny
      do k = nair,nz
        write(iof,*) j,k,E%v(k,j)
      end do
    end do   

  end subroutine write_rscalar_2d
  
  
end module boundary_Pek