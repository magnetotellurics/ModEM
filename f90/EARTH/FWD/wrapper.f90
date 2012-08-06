! *****************************************************************************
module wrapper
  ! Magnetic fields in this program are stored in three different data structures,
  ! each used for its own purpose. One is dimension (np2) vectorh, that is used
  ! by the iterative forward solver. Another is dimension (np1) Hx,Hy,Hz, that is
  ! used by the divergence calculations for the divergence correction only.
  ! The last one is type (cvector) H, which is a convenient data structure to
  ! use for the Jacobian calculations and everywhere outside the iterative forward
  ! solver.
  ! The boundary values are also stored separately in type (sparsevecc) bvH.
  ! Data structure vectorh does not contain the boundary values.
  ! The subroutines in this module convert all these data structures from one to the
  ! other. They also multiply or divide by length or perpendicular surface elements.
  ! For this latter purpose, separate routines have been devised.
  ! No pre-multiplication is any longer performed during conversion from one data
  ! structure to another. They all contain _pure_fields_ (unless one of the
  ! multiplication or division routines has been called on them).
  ! The date of this modification is Oct 05, 2005.
  ! The data of last modification is Nov 27, 2005.

  use dimensions
  use elements
  use griddef
  use sg_vector
  use sg_sparse_vector
  implicit none

  ! NB: length np2 vector contains no information whatsoever regarding the
  ! boundary values of magnetic field. They are ignored (saved in vector b instead).
  ! However, type (cvector) computed from it is aimed to contain full information
  ! about the field, so it reads H%x,y,z at k=1 and H%x,y at k=nz+1 from the saved
  ! array bvH and stores it.


Contains


! *****************************************************************************
! * vectorh -> Hx,Hy,Hz; boundary values are not copied but saved in Hx,Hy,Hz
! *
! * No division by edge lengths, contrary to the previous versions.
! * Last mod.: Oct 05, 2005

  subroutine copyd1_d3_b(l,m,n,rx,ry,rz,rvec,x,y,z)

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
            !rx(ii)=rvec(ic)/xlen
			!print *, 'd1->d3: x ',i,j,k,ii,ic,rvec(ic)
			rx(ii)=rvec(ic)
          end do
        end do
      end do

      do j=1,m
        do k=2,n
          call leng_yijk(j,k,y,z,ylen)
          do i=1,l
            ic=ic+1
            call n_allhyijk(l,m,n,i,j,k,ii)
            !ry(ii)=rvec(ic)/ylen
			!print *, 'd1->d3: y ',i,j,k,ii,ic,rvec(ic)
            ry(ii)=rvec(ic)
          end do
        end do
      end do

      do k=2,n
        ic=ic+1
        call leng_zijk(k,z,zlen)
        call n_allhzijk(l,m,n,1,1,k,ii)
        !rz(ii)=rvec(ic)/zlen
        rz(ii)=rvec(ic)
        do j=2,m
          do i=1,l
            ic=ic+1
            call n_allhzijk(l,m,n,i,j,k,ii)
            !rz(ii)=rvec(ic)/zlen
			!print *, 'd1->d3: z ',i,j,k,ii,ic,rvec(ic)
            rz(ii)=rvec(ic)
          end do
        end do
        ic=ic+1
        call n_allhzijk(l,m,n,1,m+1,k,ii)
        !rz(ii)=rvec(ic)/zlen
        rz(ii)=rvec(ic)
      end do

      return

  end subroutine copyd1_d3_b


! *****************************************************************************
! * Hx,Hy,Hz -> vectorh; boundary values are not saved anywhere
! *
! * No multiplication by edge lengths, contrary to the previous versions.
! * Last mod.: Oct 05, 2005

  subroutine copyd3_d1_b(l,m,n,rx,ry,rz,rvec,x,y,z)

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
            !rvec(ic)=rx(ii)*xlen
			!print *, 'd3->d1: x ',i,j,k,ii,ic,rvec(ic)
            rvec(ic)=rx(ii)
          end do
        end do
      end do

      do j=1,m
        do k=2,n
          call leng_yijk(j,k,y,z,ylen)
          do i=1,l
            ic=ic+1
            call n_allhyijk(l,m,n,i,j,k,ii)
            !rvec(ic)=ry(ii)*ylen
			!print *, 'd3->d1: y ',i,j,k,ii,ic,rvec(ic)
            rvec(ic)=ry(ii)
          end do
        end do
      end do

      do k=2,n
        ic=ic+1
        call leng_zijk(k,z,zlen)
        call n_allhzijk(l,m,n,1,1,k,ii)
        !rvec(ic)=rz(ii)*zlen
		!i=1;j=1
		!print *, 'd3->d1: z ',i,j,k,ii,ic,rvec(ic)
        rvec(ic)=rz(ii)
        do j=2,m
          do i=1,l
            ic=ic+1
            call n_allhzijk(l,m,n,i,j,k,ii)
            !rvec(ic)=rz(ii)*zlen
			!print *, 'd3->d1: z ',i,j,k,ii,ic,rvec(ic)
            rvec(ic)=rz(ii)
          end do
        end do
        ic=ic+1
        call n_allhzijk(l,m,n,1,m+1,k,ii)
        !rvec(ic)=rz(ii)*zlen
		!i=1;j=m+1
		!print *, 'd3->d1: z ',i,j,k,ii,ic,rvec(ic)
        rvec(ic)=rz(ii)
      end do

      return

  end subroutine copyd3_d1_b


! *****************************************************************************
! * divides vectorh by edge length elements < l(i,j,k) >
! *
! * Last mod.: Oct 05, 2005

  subroutine divide_vec_by_l(l,m,n,rvec,x,y,z)

	  integer					:: l,m,n
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
            rvec(ic)=rvec(ic)/xlen
          end do
        end do
      end do

      do j=1,m
        do k=2,n
          call leng_yijk(j,k,y,z,ylen)
          do i=1,l
            ic=ic+1
            rvec(ic)=rvec(ic)/ylen
          end do
        end do
      end do

      do k=2,n
        ic=ic+1
        call leng_zijk(k,z,zlen)
        rvec(ic)=rvec(ic)/zlen
        do j=2,m
          do i=1,l
            ic=ic+1
            rvec(ic)=rvec(ic)/zlen
          end do
        end do
        ic=ic+1
        rvec(ic)=rvec(ic)/zlen
      end do

      return
      end subroutine divide_vec_by_l


! *****************************************************************************
! * multiplies vectorh by edge length elements < l(i,j,k) >
! *
! * Last mod.: Oct 05, 2005

  subroutine mult_vec_by_l(l,m,n,rvec,x,y,z)

	  integer					:: l,m,n
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
            rvec(ic)=rvec(ic)*xlen
          end do
        end do
      end do

      do j=1,m
        do k=2,n
          call leng_yijk(j,k,y,z,ylen)
          do i=1,l
            ic=ic+1
            rvec(ic)=rvec(ic)*ylen
          end do
        end do
      end do

      do k=2,n
        ic=ic+1
        call leng_zijk(k,z,zlen)
        rvec(ic)=rvec(ic)*zlen
        do j=2,m
          do i=1,l
            ic=ic+1
            rvec(ic)=rvec(ic)*zlen
          end do
        end do
        ic=ic+1
        rvec(ic)=rvec(ic)*zlen
      end do

      return
      end subroutine mult_vec_by_l


! *****************************************************************************
! * divides vectorh by perp. surface area elements < S(i+1/2)(j-1/2)(k-1/2) >
! *
! * Last mod.: Oct 05, 2005

  subroutine divide_vec_by_S(l,m,n,rvec,x,y,z)

	  integer					:: l,m,n
      real(8),dimension(l)	    :: x
      real(8),dimension(m+1)    :: y
      real(8),dimension(n+1)    :: z
      complex(8),dimension(np2) :: rvec
      real(8)                   :: sijk2,sjki2,skij2
      real(8)                   :: ym,yp,zp
      integer                   :: ic,i,j,k,ii

      ic=0
      do i=1,l
        do k=2,n
          do j=2,m
            ic=ic+1
            call area_sijk2(j-1,k-1,y,z,sijk2)
            rvec(ic)=rvec(ic)/sijk2
          end do
        end do
      end do

      do j=1,m
        do k=2,n
		  ! Zero longitude
          ic=ic+1
          call area_sjki2(l,1,j,k-1,x,y,z,sjki2)
          rvec(ic)=rvec(ic)/sjki2
          do i=2,l
            ic=ic+1
            call area_sjki2(i-1,i,j,k-1,x,y,z,sjki2)
            rvec(ic)=rvec(ic)/sjki2
          end do
        end do
      end do

      do k=2,n ! BOUNDARY CONDITIONS ERROR: should loop from k=1
        ic=ic+1
		! North pole cap
		j=1
		zp=(z(k)+z(k+1))/2.0d0
        yp=(y(j)+y(j+1))/2.0d0
        skij2=2.0d0*pi*(zp**2)*(1.0d0-dcos(yp))
		rvec(ic)=rvec(ic)/skij2
        do j=2,m
		  ! Zero longitude
		  ic=ic+1
          call area_skij2(l,1,j-1,k,x,y,z,skij2)
		  rvec(ic)=rvec(ic)/skij2
          do i=2,l
            ic=ic+1
            call area_skij2(i-1,i,j-1,k,x,y,z,skij2)
			rvec(ic)=rvec(ic)/skij2
          end do
        end do
        ic=ic+1
		! South pole cap
		j=m+1
		zp=(z(k)+z(k+1))/2.0d0
        ym=(y(j)+y(j-1))/2.0d0
        skij2=2.0d0*pi*(zp**2)*(dcos(ym)+1.0d0)
		rvec(ic)=rvec(ic)/skij2
      end do

      return
      end subroutine divide_vec_by_S


! *****************************************************************************
! * multiplies vectorh by perp. surface area elements < S(i+1/2)(j-1/2)(k-1/2) >
! *
! * Last mod.: Oct 05, 2005

  subroutine mult_vec_by_S(l,m,n,rvec,x,y,z)

	  integer					:: l,m,n
      real(8),dimension(l)	    :: x
      real(8),dimension(m+1)    :: y
      real(8),dimension(n+1)    :: z
      complex(8),dimension(np2) :: rvec
      real(8)                   :: sijk2,sjki2,skij2
      real(8)                   :: ym,yp,zp
      integer                   :: ic,i,j,k,ii

      ic=0
      do i=1,l
        do k=2,n
          do j=2,m
            ic=ic+1
            call area_sijk2(j-1,k-1,y,z,sijk2)
            rvec(ic)=rvec(ic)*sijk2
          end do
        end do
      end do

      do j=1,m
        do k=2,n
		  ! Zero longitude
          ic=ic+1
          call area_sjki2(l,1,j,k-1,x,y,z,sjki2)
          rvec(ic)=rvec(ic)*sjki2
          do i=2,l
            ic=ic+1
            call area_sjki2(i-1,i,j,k-1,x,y,z,sjki2)
            rvec(ic)=rvec(ic)*sjki2
          end do
        end do
      end do

      do k=2,n
        ic=ic+1
		! North pole cap
		j=1
		zp=(z(k)+z(k+1))/2.0d0
        yp=(y(j)+y(j+1))/2.0d0
        skij2=2.0d0*pi*(zp**2)*(1.0d0-dcos(yp))
		rvec(ic)=rvec(ic)*skij2
        do j=2,m
		  ! Zero longitude
		  ic=ic+1
          call area_skij2(l,1,j-1,k,x,y,z,skij2)
		  rvec(ic)=rvec(ic)*skij2
          do i=2,l
            ic=ic+1
            call area_skij2(i-1,i,j-1,k,x,y,z,skij2)
			rvec(ic)=rvec(ic)*skij2
          end do
        end do
        ic=ic+1
		! South pole cap
		j=m+1
		zp=(z(k)+z(k+1))/2.0d0
        ym=(y(j)+y(j-1))/2.0d0
        skij2=2.0d0*pi*(zp**2)*(dcos(ym)+1.0d0)
		rvec(ic)=rvec(ic)*skij2
      end do

      return
      end subroutine mult_vec_by_S


! ***************************************************************************
! * copyd3_d1_d converts a type (cvector) H into dimension (np2) vectorh.
! * required by SolveMaxwells; rvec is assumed to have been allocated.
! * x- and y- components of the field on the upper and lower domain boundaries,
! * and z-component on the upper boundary only are not saved in rvec.
! * Instead, they are saved in the type (sparsevecc) bvH, and used to
! * reconstruct the original vecE by the complementary subroutine copyd1_d3_d.
! *
! * No multiplication by edge length elements, contrary to the previous versions.
! * Last mod.: Oct 05, 2005

  subroutine copyd3_d1_d(rvec,vecE,igrid)


	complex(8), dimension(:), intent(out)				  :: rvec !np2
	type (cvector), intent(in)							  :: vecE
	type (grid_t), intent(in)						  :: igrid
	integer												  :: l,m,n
    real(8),dimension(:),allocatable					  :: x
    real(8),dimension(:),allocatable					  :: y
    real(8),dimension(:),allocatable					  :: z
    integer												  :: ic,i,j,k,ii

	! Check whether input vector is defined on edges
	if (vecE%gridType /= EDGE) then
	  write(0, *) 'Error: (copyd3_d1_d) input vector not defined on edges'
	  stop
	end if

	! Check whether output vector is of the right size
	if (size(rvec) /= np2) then
	  write(0, *) 'Error: (copyd3_d1_d) output vector length is not np2'
	  stop
	end if

	l = vecE%nx
	m = vecE%ny
	n = vecE%nz

	allocate(x(l),y(m+1),z(n+1))

	x = igrid%x
	y = igrid%y
	z = igrid%z

!
! Hx
!
      ic=0
      do i=1,l
        do k=2,n
          do j=2,m
            ic=ic+1
			rvec(ic)=vecE%x(i,j,k)
			!print *, 'd3->d1: x ',i,j,k,ii,ic,rvec(ic)
          end do
        end do
      end do
!
! Hy
!
      do j=1,m
        do k=2,n
          do i=1,l
            ic=ic+1
			rvec(ic)=vecE%y(i,j,k)
			!print *, 'd3->d1: y ',i,j,k,ii,ic,rvec(ic)
          end do
        end do
      end do
!
! Hz
!
      do k=2,n
        ic=ic+1
		rvec(ic)=vecE%z(1,1,k)
        do j=2,m
          do i=1,l
            ic=ic+1
			rvec(ic)=vecE%z(i,j,k)
			!print *, 'd3->d1: z ',i,j,k,ii,ic,rvec(ic)
          end do
        end do
        ic=ic+1
		rvec(ic)=vecE%z(1,m+1,k)
      end do

      return

  end subroutine copyd3_d1_d  ! 3D -> 1D


! ***************************************************************************
! * copyd1_d3_d converts a dimension (np2) vectorh into a type (cvector) H,
! * required by Jacobian calculations; H is allocated and initialized inside.
! * x- and y- components of the field on the upper and lower domain boundaries,
! * and z-component on the upper boundary only are copied from bvH for H fields.
! *
! * No division by edge length elements, contrary to the previous versions.
! * Last mod.: Oct 13, 2005

  subroutine copyd1_d3_d(rvec,vecE,igrid,bc)

	complex(8), dimension(:), intent(in)				  :: rvec !np2
	type (grid_t), intent(in)						  :: igrid
	type (cvector), intent(inout)							  :: vecE
	type (sparsevecc), intent(in), optional				  :: bc
	integer												  :: i,j,k,ii,ic
	integer												  :: l,m,n
    real(8)												  :: xlen,ylen,zlen


	! Check whether input vector is of the right size
	if (size(rvec) /= np2) then
	  write(0, *) 'Error: (copyd1_d3_d) input vector length is not np2'
	  stop
	end if

	! Deallocate output vector, if needed
	call deall_cvector(vecE)

	! Create output vector
	call create_cvector(igrid, vecE, EDGE)
	if (present(bc)) then
	  ! Insert boundary conditions
	  call add_scvector(C_ONE,bc,vecE)
	end if

	l = vecE%nx
	m = vecE%ny
	n = vecE%nz

!
! Hx
!
    ic=0
	do i=1,l
	  do k=2,n
		do j=2,m
		  ic=ic+1
          !call n_allhijk(l,m,n,i,j,k,1,ic)
		  vecE%x(i,j,k) = rvec(ic)
        end do
      end do
	end do
!
! Hy
!
    do j=1,m
      do k=2,n
        do i=1,l
		  ic=ic+1
          !call n_allhijk(l,m,n,i,j,k,2,ic)
		  vecE%y(i,j,k) = rvec(ic)
        end do
      end do
    end do
!
! Hz
!
    do k=2,n
	  ic=ic+1
      !call n_allhijk(l,m,n,1,1,k,3,ic)
	  do i=1,l
		vecE%z(i,1,k) = rvec(ic)
	  end do
      do j=2,m
        do i=1,l
		  ic=ic+1
		  !call n_allhijk(l,m,n,i,j,k,3,ic)
		  vecE%z(i,j,k) = rvec(ic)
        end do
      end do
	  ic=ic+1
      !call n_allhijk(l,m,n,1,m+1,k,3,ic)
	  do i=1,l
		vecE%z(i,m+1,k) = rvec(ic)
	  end do
    end do

	! redefine fields at zero longitude to fill type (cvector) variable
	i=l+1
	  ! x-component: fields undefined
	  ! y-component:
	  do j=1,m
		do k=1,n+1
		  vecE%y(i,j,k) = vecE%y(1,j,k)
		end do
	  end do
	  ! z-component:
	  do j=1,m+1
		do k=1,n
		  vecE%z(i,j,k) = vecE%z(1,j,k)
		end do
	  end do

  end subroutine copyd1_d3_d	! 1D -> 3D

! *****************************************************************************
! Boundary conditions subroutines that serve a temporary purpose


  ! ***************************************************************************
  ! * initBoundaryValues is the routine to save the initial values of the
  ! * magnetic fields on boundaries as boundary conditions
  subroutine initBoundaryValues(bvH,Hx,Hy,Hz,grid)

	type (sparsevecc), intent(out)			:: bvH
	type (grid_t), intent(in)				:: grid
	complex(8), dimension(:), intent(in)	:: Hx,Hy,Hz
	integer									:: i,j,k,ii,ib
	integer									:: bv1,bv2,bv3
	integer, parameter						:: x=1,y=2,z=3
	integer									:: l,m,n

	bv1= 2*grid%nx*(grid%ny-1)
	bv2= 2*grid%nx*grid%ny
	bv3=   grid%nx*(grid%ny-1)+2*grid%nx !of which nx*(ny-1)+2 are stored in Hz
	l=grid%nx
	m=grid%ny
	n=grid%nz

	call create_sparsevecc(bv1+bv2+bv3,bvH,EDGE)

    ib=0
! Hx
	k=1
	  do i=1,l
        do j=2,m
          ib=ib+1
          call n_allhxijk(l,m,n,i,j,k,ii)
		  call newValueC_sparsevecc(bvH,ib,Hx(ii),i,j,k,x)
		  !print *, 'd3->d1: x ',i,j,k,ii,ib,Hx(ii)
        end do
      end do
	k=n+1
	  do i=1,l
        do j=2,m
          ib=ib+1
          call n_allhxijk(l,m,n,i,j,k,ii)
		  call newValueC_sparsevecc(bvH,ib,Hx(ii),i,j,k,x)
		  !print *, 'd3->d1: x ',i,j,k,ii,ib,Hx(ii)
        end do
      end do
! Hy
	k=1
	  do i=1,l
        do j=1,m
          ib=ib+1
          call n_allhyijk(l,m,n,i,j,k,ii)
		  call newValueC_sparsevecc(bvH,ib,Hy(ii),i,j,k,y)
		  !print *, 'd3->d1: y ',i,j,k,ii,ib,Hy(ii)
        end do
      end do
	k=n+1
	  do i=1,l
        do j=1,m
          ib=ib+1
          call n_allhyijk(l,m,n,i,j,k,ii)
		  call newValueC_sparsevecc(bvH,ib,Hy(ii),i,j,k,y)
		  !print *, 'd3->d1: y ',i,j,k,ii,ib,Hy(ii)
        end do
      end do
! Hz
	k=1
      call n_allhzijk(l,m,n,1,1,k,ii)
	  do i=1,l
		ib=ib+1
		call newValueC_sparsevecc(bvH,ib,Hz(ii),i,1,k,z)
	  end do
      do j=2,m
        do i=1,l
          ib=ib+1
          call n_allhzijk(l,m,n,i,j,k,ii)
		  call newValueC_sparsevecc(bvH,ib,Hz(ii),i,j,k,z)
        end do
      end do
      call n_allhzijk(l,m,n,1,m+1,k,ii)
	  do i=1,l
		ib=ib+1
		call newValueC_sparsevecc(bvH,ib,Hz(ii),i,m+1,k,z)
	  end do

  end subroutine initBoundaryValues


  ! ***************************************************************************
  ! * insertBoundaryValues is the routine to insert the boundary values of the
  ! * magnetic fields into the Hx,Hy,Hz representation of the fields

  subroutine insertBoundaryValues(bvH,Hx,Hy,Hz,grid)

	use elements
	type (sparsevecc), intent(in)			:: bvH
	complex(8), dimension(:), intent(out)	:: Hx,Hy,Hz
	type (grid_t), intent(in)               :: grid
	integer									:: i,j,k,ii,ib
	integer									:: bv1,bv2,bv3
	integer									:: l,m,n

	l=grid%nx
	m=grid%ny
	n=grid%nz
	!bv1= 2*l*(m-1)
	!bv2= 2*l*m
	!bv3=   l*(m-1)+2*l !of which nx*(ny-1)+2 are stored in Hz

	do ib=1,bvH%nCoeff

	  select case (bvH%xyz(ib))
	  case (1)
		call n_allhxijk(l,m,n,bvH%i(ib),bvH%j(ib),bvH%k(ib),ii)
		Hx(ii) = bvH%c(ib)
	  case (2)
		call n_allhyijk(l,m,n,bvH%i(ib),bvH%j(ib),bvH%k(ib),ii)
		Hy(ii) = bvH%c(ib)
	  case (3)
		if ((bvH%j(ib) == 1).or.(bvH%j(ib) == m+1)) then  !poles
		  if (bvH%i(ib) == 1) then
			call n_allhzijk(l,m,n,bvH%i(ib),bvH%j(ib),bvH%k(ib),ii)
			Hz(ii) = bvH%c(ib)
		  else
		  end if
		else
		  call n_allhzijk(l,m,n,bvH%i(ib),bvH%j(ib),bvH%k(ib),ii)
		  Hz(ii) = bvH%c(ib)
		end if
	  end select

	end do

  end subroutine insertBoundaryValues ! insertBoundaryValues


  ! ***************************************************************************
  ! * divide_bc_by_l multiplies the boundary values of the magnetic fields saved
  ! * in type (cvector) by the relevant edge lengths. This is an implementation
  ! * of the diagonal operator l_b.

  subroutine divide_bc_by_l(Hb,lHb,grid)

	use elements
	type (sparsevecc), intent(in)			:: Hb
	type (sparsevecc), intent(out)			:: lHb
    type (grid_t), intent(in)               :: grid
	integer									:: i,j,k,ib
	integer									:: bv1,bv2,bv3
	real(8)									:: xx,yy,zz

	!bv1= 2*l*(m-1)
	!bv2= 2*l*m
	!bv3=   l*(m-1)+2*l !of which nx*(ny-1)+2 are stored in Hz

	call copy_sparsevecc(lHb,Hb)
	lHb%c = C_ZERO

	do ib=1,Hb%nCoeff

	  select case (Hb%xyz(ib))
	  case (1)
        call leng_xijk(Hb%i(ib),Hb%j(ib),Hb%k(ib),grid%x,grid%y,grid%z,xx)
		lHb%c(ib) = Hb%c(ib)/xx
	  case (2)
        call leng_yijk(Hb%j(ib),Hb%k(ib),grid%y,grid%z,yy)
		lHb%c(ib) = Hb%c(ib)/yy
	  case (3)
        call leng_zijk(Hb%k(ib),grid%z,zz)
		lHb%c(ib) = Hb%c(ib)/zz
	  end select

	end do

  end subroutine divide_bc_by_l ! divide_bc_by_l


  ! ***************************************************************************
  ! * mult_bc_by_l multiplies the boundary values of the magnetic fields saved
  ! * in type (cvector) by the relevant edge lengths. This is an implementation
  ! * of the diagonal operator l_b.

  subroutine mult_bc_by_l(Hb,lHb,grid)

	use elements
	type (sparsevecc), intent(in)			:: Hb
	type (sparsevecc), intent(out)			:: lHb
    type (grid_t), intent(in)               :: grid
	integer									:: i,j,k,ib
	integer									:: bv1,bv2,bv3
	real(8)									:: xx,yy,zz

	!bv1= 2*l*(m-1)
	!bv2= 2*l*m
	!bv3=   l*(m-1)+2*l !of which nx*(ny-1)+2 are stored in Hz

	call copy_sparsevecc(lHb,Hb)
	lHb%c = C_ZERO

	do ib=1,Hb%nCoeff

	  select case (Hb%xyz(ib))
	  case (1)
        call leng_xijk(Hb%i(ib),Hb%j(ib),Hb%k(ib),grid%x,grid%y,grid%z,xx)
		lHb%c(ib) = xx*Hb%c(ib)
	  case (2)
        call leng_yijk(Hb%j(ib),Hb%k(ib),grid%y,grid%z,yy)
		lHb%c(ib) = yy*Hb%c(ib)
	  case (3)
        call leng_zijk(Hb%k(ib),grid%z,zz)
		lHb%c(ib) = zz*Hb%c(ib)
	  end select

	end do

  end subroutine mult_bc_by_l ! mult_bc_by_l


! *****************************************************************************
! * calculate vector dimension(np2) vectorb of boundary conditions on interior
! * nodes from the vector of fields on boundary nodes type (sparsevecc) bc
! *
! * Modified for spherical problem on 26/11/97
! * Modified for type (sparsevecc) boundary conditions on 05/09/05
! * Does *not* pre-multiply by edge lengths any longer
! * Last mod.: Oct 05, 2005

  subroutine calcb_from_bc(l,m,n,bc,bvec,resist,x,y,z,grid)


	  use griddef
      implicit none

      integer                     :: l,m,n,i,j,k,ii,jj
      integer                     :: a,b,c,d,e,f

      real(8),dimension(l,m,n)	  :: resist

      real(8),dimension(l)	      :: x
      real(8),dimension(m+1)      :: y
      real(8),dimension(n+1)      :: z
	  type (sparsevecc)			  :: bc,newbc
	  type (grid_t)               :: grid
      complex(8),dimension(:),allocatable   :: hx,hy,hz	!np1
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

	  integer					  :: istat

	  allocate(hx(np1),hy(np1),hz(np1),STAT=istat)
	  hx = C_ZERO
	  hy = C_ZERO
	  hz = C_ZERO
	  call mult_bc_by_l(bc,newbc,grid)
	  call insertBoundaryValues(newbc,hx,hy,hz,grid)
	  call deall_sparsevecc(newbc)

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
              !sum=sum+rijcy*xx*hx(ii)
              sum=sum+rijcy*hx(ii)
            end if
! 2 Hxijf
            if (k == n) then
!           ----------> bottom Hx
              call leng_xijk(i,j,f,x,y,z,xx)
              call n_allhxijk(L,M,N,i,j,f,ii)
              !sum=sum+rijky*xx*hx(ii)
              sum=sum+rijky*hx(ii)
            end if
! 3 Hzijc
           if (k == 2) then
!           ----------> top Hz
              call n_allhzijk(L,M,N,i,j,c,ii)
              !sum=sum-rijcy*zc*hz(ii)
              sum=sum-rijcy*hz(ii)
            end if
! 4 Hzdjc
           if (k == 2) then
!           ----------> top Hz
              call n_allhzijk(L,M,N,d,j,c,ii)
              !sum=sum+rijcy*zc*hz(ii)
              sum=sum+rijcy*hz(ii)
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
              !sum=sum+rijcx*yy*hy(ii)
              sum=sum+rijcx*hy(ii)
            end if
! 2 Hyijf
            if (k == n) then
!           ----------> bottom Hy
              call leng_yijk(j,f,y,z,yy)
              call n_allhyijk(l,m,n,i,j,f,ii)
              !sum=sum+rijkx*yy*hy(ii)
              sum=sum+rijkx*hy(ii)
            end if
! 3 Hzijc
            if (k == 2) then
!           ----------> top Hz
              call n_allhzijk(l,m,n,i,j,c,ii)
              !sum=sum-rijcx*zc*hz(ii)
              sum=sum-rijcx*hz(ii)
            end if
! 4 Hziec
            if (k == 2) then
!           ----------> top Hz
              call n_allhzijk(l,m,n,i,e,c,ii)
              !sum=sum+rijcx*zc*hz(ii)
              sum=sum+rijcx*hz(ii)
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
                !sum=sum-rijkx*yy*hy(ii)
                sum=sum-rijkx*hy(ii)
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
                !sum=sum+ribkx*yy*hy(ii)
                sum=sum+ribkx*hy(ii)
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
              !sum=sum-rijky*xx*hx(ii)
              sum=sum-rijky*hx(ii)
            end if
! 2 Hxajf
            if (k == n) then
!           ----------> bottom Hx
              call leng_xijk(a,j,f,x,y,z,xx)
              call n_allhxijk(l,m,n,a,j,f,ii)
              !sum=sum+rajky*xx*hx(ii)
              sum=sum-rijky*hx(ii)
            end if
! 3 Hyijf
            if (k == n) then
!           ----------> bottom Hy
              call leng_yijk(j,f,y,z,yy)
              call n_allhyijk(l,m,n,i,j,f,ii)
              !sum=sum-rijkx*yy*hy(ii)
              sum=sum-rijkx*hy(ii)
            end if
! 4 Hyibf
            if (k == n) then
!           ----------> bottom Hy
              call leng_yijk(b,f,y,z,yy)
              call n_allhyijk(l,m,n,i,b,f,ii)
              !sum=sum+ribkx*yy*hy(ii)
              sum=sum+ribkx*hy(ii)
            end if

            bvec(jj)=sum
!           -----------
            end do
!           -----------
!         ++++++++++++++++
          end if
!         ++++++++++++++++
30    continue

	  deallocate(hx,hy,hz)

      return

  end subroutine calcb_from_bc	! BC -> interior nodes


! *****************************************************************************

  ! ***************************************************************************
  ! * Utility subroutine to compare two type (cvector) vector fields and to
  ! * output any discrepancies to screen

  subroutine compare_fields(vecE1,vecE2)

	type (cvector), intent(in)			  :: vecE1,vecE2
	real									  :: eps
	integer									  :: i,j,k

	eps = 0.0000001

	do k=1,vecE1%nz+1
	  do j=1,vecE1%ny+1
		do i=1,vecE1%nx
		  ! Compare real parts
		  if (abs(dreal(vecE1%x(i,j,k)-vecE2%x(i,j,k)))>eps) then
			print *, 'x: ', i,j,k, dreal(vecE1%x(i,j,k)), dreal(vecE2%x(i,j,k))
		  end if
		  if (abs(dreal(vecE1%y(i,j,k)-vecE2%y(i,j,k)))>eps) then
			print *, 'y: ', i,j,k, dreal(vecE1%y(i,j,k)), dreal(vecE2%y(i,j,k))
		  end if
		  if (abs(dreal(vecE1%z(i,j,k)-vecE2%z(i,j,k)))>eps) then
			print *, 'z: ', i,j,k, dreal(vecE1%z(i,j,k)), dreal(vecE2%z(i,j,k))
		  end if
		  ! Compare imaginary parts
		  if (abs(dimag(vecE1%x(i,j,k)-vecE2%x(i,j,k)))>eps) then
			print *, 'x: ', i,j,k, dimag(vecE1%x(i,j,k)), dimag(vecE2%x(i,j,k))
		  end if
		  if (abs(dimag(vecE1%y(i,j,k)-vecE2%y(i,j,k)))>eps) then
			print *, 'y: ', i,j,k, dimag(vecE1%y(i,j,k)), dimag(vecE2%y(i,j,k))
		  end if
		  if (abs(dimag(vecE1%z(i,j,k)-vecE2%z(i,j,k)))>eps) then
			print *, 'z: ', i,j,k, dimag(vecE1%z(i,j,k)), dimag(vecE2%z(i,j,k))
		  end if
		end do
	  end do
	end do


  end subroutine compare_fields	! compare_fields


  ! ***************************************************************************
  ! * Utility subroutine to compare two type length np2 vector fields and to
  ! * output any discrepancies to screen

  subroutine compare_fields_np2(rvec1,rvec2)

	complex(8), dimension(:), intent(in)				  :: rvec1,rvec2 !np2
	real												  :: eps
	integer												  :: ii

	eps = 0.0000001

	do ii=1,np2
	  if (abs(dreal(rvec1(ii)-rvec2(ii)))>eps) then
		print *, ii, dreal(rvec1(ii)), dreal(rvec2(ii))
	  end if
	  if (abs(dimag(rvec1(ii)-rvec2(ii)))>eps) then
		print *, ii, dimag(rvec1(ii)), dimag(rvec2(ii))
	  end if
	end do

  end subroutine compare_fields_np2	! compare_fields


end module wrapper
