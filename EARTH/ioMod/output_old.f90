! *****************************************************************************
module output_old
  ! Module containing the subroutines for output only

  ! NEED TO CLEAN UP THIS MODULE

  use iotypes
  use iospec
  use griddef
  use sg_vector
  use sg_sparse_vector
  use modeldef
  use global
  use responses
  use functionals
  implicit none


Contains

  ! ***************************************************************************
  ! * outputCellResp is a utility subroutine to compute and output C and D
  ! * responses at cells. Effectively, the subroutine builds an observatory
  ! * at every cell, such that the (i,j)'th observatory is located at the mid-pt
  ! * between (i,j,k)'th and (i,j+1,k)'th nodes. At each of these observatories,
  ! * the fields and then the responses are computed. Once that is done, we output
  ! * these values to a file. If a value is not scientifically reasonable
  ! * (ie, the observatory is located too close to one of the poles or equator,
  ! * in which case the C- and D- responses aren't stable), we output NaN instead.
  ! * We only need to compute these values for testing and cross-comparison with
  ! * the original version of the forward solver.

  subroutine outputCellResp(freq,fname,H)

	type (transmitter_t), intent(in)					:: freq
	character(len=*), intent(in)					:: fname
	type (cvector), intent(in)						:: H
	type (dataValue_t), dimension(:,:), allocatable :: C_ijk,D_ijk	
	integer											:: istat,i,j

	! compute and output C and D responses at cells
	allocate(C_ijk(grid%nx,grid%ny),D_ijk(grid%nx,grid%ny),STAT=istat)
	do i=1,grid%nx
	  do j=1,grid%ny
		C_ijk(i,j)%resp%value = dataResp_ijk(C_ratio,i,j,grid%nzAir+1,H)
		D_ijk(i,j)%resp%value = dataResp_ijk(D_ratio,i,j,grid%nzAir+1,H)
	  end do
	end do
	call outputResponsesAtCells(freq,trim(fname)//'.cresp',C_ijk%resp%value)
	call outputResponsesAtCells(freq,trim(fname)//'.dresp',D_ijk%resp%value)
	deallocate(C_ijk,D_ijk)

  end subroutine outputCellResp


  ! ***************************************************************************
  ! * outputResponses writes the chosen kind of responses calculated at every
  ! * observatory location to an output file
  subroutine outputResponsesAtObs(fn_response,Resp,ifreq,nfreq)

	character(80), intent(in)						:: fn_response
	complex(8), dimension(:), intent(in)			:: Resp
	integer, intent(in)								:: ifreq,nfreq
	integer											:: ios,iobs

	if (fn_response /= '') then
	  if (ifreq == 1) then
		open(ioResp,file=fn_response,status='unknown',form='formatted',iostat=ios)
	    write(ioResp,*) "#Output of earth3d (",nfreq,") frequency values";
		write(ioResp,*) "#Period Code GM_Lon GM_Lat Real(km) Imag(km) Error(km)";
	  else
		open(ioResp,file=fn_response,position='append', form='formatted',iostat=ios)
	  end if		
		!write(ioResp,*) "freq = ",freqList%info(ifreq)%value
		do iobs=1,size(Resp)
		  if (.not.obsList%info(iobs)%defined) then
			!write(0,*) 'Observatory ',trim(obsList%info(iobs)%code),' is not defined at output; ignore'
			cycle
		  end if
		  write(ioResp,'(f0.3,a12,4g15.7)') &
			  freqList%info(ifreq)%period,&
			  trim(obsList%info(iobs)%code),&
			  obsList%info(iobs)%lon,obsList%info(iobs)%lat,&
			  Resp(iobs) * m2km
		  !write(ioResp,*) trim(obsList%info(iobs)%code),Resp(iobs) * m2km
		end do
	  close(ioResp)
	else
	  do iobs=1,size(Resp)
		write(*,*) trim(obsList%info(iobs)%code),Resp(iobs) * m2km
	  end do
	end if

  end subroutine outputResponsesAtObs  ! outputResponses



  ! ***************************************************************************
  ! * outputResponsesAtNodes writes the chosen kind of responses calculated at
  ! * every surface grid cell (mid-point between (i,j,k)th and (i,j+1,k)th nodes)
  ! * to an output file
  subroutine outputResponsesAtCells(freq,fn_response,Resp)

	type (transmitter_t), intent(in)					:: freq
	character(len=*), intent(in)					:: fn_response
	complex(8), dimension(:,:), intent(in)			:: Resp
	integer											:: ios,i,j,k

	if (fn_response /= '') then
	  if (freq%i == 1) then
		open(ioResp,file=trim(fn_response),status='unknown',form='formatted',iostat=ios)
	  else
		open(ioResp,file=trim(fn_response),position='append', form='formatted',iostat=ios)
	  end if		
		write(ioResp,*) "freq = ",freq%value
		k=grid%nzAir+1
		  do i=1,grid%nx
			do j=1,grid%ny
			  if((j==1).or.(j==grid%ny)) then
			  ! if (.not.Resp(i,j)%obs%defined) then
			  	write(ioResp,'(2i5,2a10)') i,j,'NaN','NaN'
			  else
				write(ioResp,'(2i5,2g15.7)') i,j,Resp(i,j) * m2km
			  end if
			end do
		  end do
	  close(ioResp)
	else
	  k=grid%nzAir+1
		do i=1,grid%nx
		  do j=1,grid%ny
			write(*,*) '#f,i,j,response = ',freq%i,i,j,Resp(i,j) * m2km
		  end do
		end do
	end if

  end subroutine outputResponsesAtCells  ! outputResponsesAtCells


  subroutine output_res(res,ifreq,ifunc,iobs)

	type (dataValue_t), dimension(:,:,:), intent(in)	:: res
	integer, intent(in)									:: ifreq,ifunc,iobs

  	write(*,'(a40,4i3,5g15.7)') &
	'ifunc,iobs,i,j,real_res,imag_res,err,res/err = ',ifunc,iobs,&
	  res(ifreq,ifunc,iobs)%obs%i,res(ifreq,ifunc,iobs)%obs%j,&
	  dreal(res(ifreq,ifunc,iobs)%resp%value),&
	  dimag(res(ifreq,ifunc,iobs)%resp%value),res(ifreq,ifunc,iobs)%resp%err,&
	  dreal(res(ifreq,ifunc,iobs)%resp%value)/res(ifreq,ifunc,iobs)%resp%err,&
	  dimag(res(ifreq,ifunc,iobs)%resp%value)/res(ifreq,ifunc,iobs)%resp%err

  end subroutine output_res	! output_res


      subroutine outputH(fname,H,grid)

!-------------------------------------------------------------
!     to output all the magnetic fields with position
!-------------------------------------------------------------

	  implicit none

	  character(len=*), intent(in)			:: fname
	  type (cvector), intent(in)		:: H
	  type (grid_t), intent(in)		:: grid
	  integer						   	 :: l,m,n
      real(8),dimension(:),allocatable	             :: x
      real(8),dimension(:),allocatable             :: y
      real(8),dimension(:),allocatable             :: z
      integer				             :: nair
      real(8)                            :: xx0,xx1,xx,yy,zz
      integer                            :: k,j,i,ii,istat


	l=grid%nx
	m=grid%ny
	n=grid%nz
	nair=grid%nzAir
	allocate(x(l),y(m+1),z(n+1),STAT=istat)
	x=grid%x
	y=grid%y
	z=grid%z
	inquire(ioHx,OPENED=opened)
    if (.not.opened) then
  	  call initFileWrite(trim(trim(fname)//'.hx'),ioHx)
	end if
	inquire(ioHy,OPENED=opened)
    if (.not.opened) then
	  call initFileWrite(trim(trim(fname)//'.hy'),ioHy)
	end if
	inquire(ioHz,OPENED=opened)
	if (.not.opened) then
	  call initFileWrite(trim(trim(fname)//'.hz'),ioHz)
	end if

!
!   First Hx....
!
!      do k=nair+1,n+1
      do k=nair+1,nair+1 !temporarily only output surface fields
        zz=z(k)/1000.d0
        do j=2,m
          yy=y(j)*r2d
          xx0=0.d0
          do i=1,l
            xx1=(xx0+x(i)*r2d)
            xx=(xx0+xx1)/2.d0
            xx0=xx1
            write(ioHx,'(2i4,5e15.7)')i,j,xx,yy,zz,H%x(i,j,k) !!!!!
          end do
        end do
      end do
!
!   Now Hy....
!
!     do k=nair+1,n+1
      do k=nair+1,nair+1 !temporarily only output surface fields
        zz=z(k)/1000.d0
        do j=1,m
          yy=( y(j)+y(j+1) )*r2d/2.d0
          xx=0.d0
          do i=1,l
             write(ioHy,'(5e15.7)')xx,yy,zz,H%y(i,j,k)
             xx=xx+x(i)*r2d
          end do
        end do
      end do

!
!   Now Hz....
!
!      do k=nair+1,n
      do k=nair+1,nair+1 !temporarily only output surface fields
        zz=( z(k)+z(k+1) )/1000.d0/2.d0
        j=1
          yy=y(j)*r2d
          i=1
            xx=0.d0
            write(ioHz,'(5e15.7)')xx,yy,zz,H%z(i,j,k)

        do j=2,m
          yy=y(j)*r2d
          xx=0.d0
          do i=1,l
            write(ioHz,'(5e15.7)')xx,yy,zz,H%z(i,j,k)
            xx=xx+x(i)*r2d
          end do
        end do

        j=m+1
          yy=y(j)*r2d
          i=1
            xx=0.d0
            write(ioHz,'(5e15.7)')xx,yy,zz,H%z(i,j,k)
      end do
  
	  deallocate(x,y,z)

      return
      end subroutine outputH


end module output_old
