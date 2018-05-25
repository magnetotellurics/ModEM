! *****************************************************************************
module boundaries
  ! Module containing the subroutines that extract the boundary conditions from
  ! the full field vector, insert them back etc

  use sg_sparse_vector


Contains

  ! ***************************************************************************
  ! * createBC creates a sparse vector of zero boundary conditions

  subroutine createBC(Hb,grid)

	type (sparsevecc), intent(inout)		:: Hb
	type (grid_t), intent(in)				:: grid
	integer									:: i,j,k,ib
	integer									:: b1,b2,b3
	integer									:: nx,ny,nz
	integer, parameter						:: x=1,y=2,z=3

	nx = grid%nx
	ny = grid%ny
	nz = grid%nz
	b1 = 2*nx*(ny-1)  ! # of x-components in Hb
	b2 = 2*nx*ny	  ! # of y-components in Hb
	b3 = 1*nx*(ny+1)  ! of which unique nx*(ny-1)+2

	call create_sparsevecc(b1+b2+b3,Hb,EDGE)

    ib=0
! Hx
	k=1
	  do i=1,nx
        do j=2,ny
          ib=ib+1
		  call newValueC_sparsevecc(Hb,ib,C_ZERO,i,j,k,x)
        end do
      end do
	k=nz+1
	  do i=1,nx
        do j=2,ny
          ib=ib+1
		  call newValueC_sparsevecc(Hb,ib,C_ZERO,i,j,k,x)
        end do
      end do
! Hy
	k=1
	  do i=1,nx
        do j=1,ny
          ib=ib+1
		  call newValueC_sparsevecc(Hb,ib,C_ZERO,i,j,k,y)
        end do
      end do
	k=nz+1
	  do i=1,nx
        do j=1,ny
          ib=ib+1
		  call newValueC_sparsevecc(Hb,ib,C_ZERO,i,j,k,y)
        end do
      end do
! Hz
	k=1
	  do i=1,nx
        do j=1,ny+1
          ib=ib+1
		  call newValueC_sparsevecc(Hb,ib,C_ZERO,i,j,k,z)
        end do
      end do


  end subroutine createBC ! createBC


  ! ***************************************************************************
  ! * insertBC adds a sparse vector of boundary conditions to a full field vector

  subroutine insertBC(Hb,H)

	type (sparsevecc), intent(in)			:: Hb
	type (cvector), intent(inout)			:: H


	call add_scvector(C_ONE,Hb,H)

  end subroutine insertBC ! insertBC


  ! ***************************************************************************
  ! * extractBC extracts a sparse vector of boundary conditions from a full field
  ! * vector. In the future, b3 should be zero: no z-components should be stored
  ! * as boundary conditions. Currently (artifact of the original version of the
  ! * program) z-components for the upper boundary are stored in Hb.

  subroutine extractBC(Hb,H)

	type (sparsevecc), intent(out)			:: Hb
	type (cvector), intent(in)				:: H
	integer									:: i,j,k,ib
	integer									:: b1,b2,b3
	integer									:: nx,ny,nz
	integer, parameter						:: x=1,y=2,z=3

	nx = H%nx
	ny = H%ny
	nz = H%nz
	b1 = 2*nx*(ny-1)  ! # of x-components in Hb
	b2 = 2*nx*ny	  ! # of y-components in Hb
	b3 = 1*nx*(ny+1)  ! of which unique nx*(ny-1)+2

	call create_sparsevecc(b1+b2+b3,Hb,EDGE)

    ib=0
! Hx
	k=1
	  do i=1,nx
        do j=2,ny
          ib=ib+1
		  call copyValue_csvector(Hb,ib,H,i,j,k,x)
        end do
      end do
	k=nz+1
	  do i=1,nx
        do j=2,ny
          ib=ib+1
		  call copyValue_csvector(Hb,ib,H,i,j,k,x)
        end do
      end do
! Hy
	k=1
	  do i=1,nx
        do j=1,ny
          ib=ib+1
		  call copyValue_csvector(Hb,ib,H,i,j,k,y)
        end do
      end do
	k=nz+1
	  do i=1,nx
        do j=1,ny
          ib=ib+1
		  call copyValue_csvector(Hb,ib,H,i,j,k,y)
        end do
      end do
! Hz
	k=1
	  do i=1,nx
        do j=1,ny+1
          ib=ib+1
		  call copyValue_csvector(Hb,ib,H,i,j,k,z)
        end do
      end do

  end subroutine extractBC	! extractBC


end module
