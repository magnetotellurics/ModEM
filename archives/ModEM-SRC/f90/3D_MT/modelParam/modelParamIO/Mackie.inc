! I/O routines for 3D_MT modelParam - extended Randie Mackie's format

  !******************************************************************
  subroutine read_modelParam_mackie(grid,m,cfile)

      !  open cfile on unit fid, writes out object of
      !   type modelParam in Randie Mackie's format, closes file
      ! we can update the grid here, but the grid as an input is critical
      ! for setting the pointer to the grid in the modelParam

      type(grid_t), target, intent(inout) :: grid
      type(modelParam_t), intent(out)	:: m
      character(*), intent(in)          :: cfile
      integer              :: ios
      ! local
      character(80)                     :: type
      type(rscalar)                     :: ccond

 	  ! Read input files and set up basic grid geometry, conductivities,
	  ! and frequencies (stored in the transmitter dictionary, txDictMT)
	  call readRMgridCond(ioPrm,cfile,ccond,grid,type,ios)

	  ! move cell conductivities read into rscalar object into a modelParam
	  ! object ... this dance needed to keep modelParam attributes private
	  !   First need to create model parameter
	  call create_modelParam(grid,type,m,ccond)

	  ! convert modelParam to the required paramType
	  ! call setType_modelParam(m,paramType)

	  ! now done with ccond, so deallocate
	  call deall_rscalar(ccond)

  end subroutine read_modelParam_mackie

  !******************************************************************
  subroutine write_modelParam_mackie(m,cfile)

      !  open cfile on unit fid, writes out object of
      !   type modelParam in Randie Mackie's format, closes file
      ! note that while this is in general a good interface for
      ! a write subroutine for modelParam, we do not have to write
      ! out the grid unless we want to do so

      character(*), intent(in)          :: cfile
      type(modelParam_t), intent(in)	:: m
      integer              :: ios
      ! local
      type(rscalar)                     :: ccond
      type(grid_t)                      :: grid
      character(80)                     :: type = ''

 	  ! Read input files and set up basic grid geometry, conductivities,
	  ! and frequencies (stored in the transmitter dictionary, txDictMT)

	  call getValue_modelParam(m,type,ccond)

	  grid = ccond%grid

	  call writeRMgridCond(ioPrm,cfile,ccond,grid,type,ios)

	  ! now done with ccond and grid, so deallocate
	  call deall_rscalar(ccond)
	  call deall_grid(grid)

  end subroutine write_modelParam_mackie

  ! ***************************************************************************
  subroutine readRMgridCond(fidRM,inputFile,Cond,grid,paramType,ioerr)
  ! this routine reads files in Mackie's 3D formats, returning the basic
  !   grid components, and optionally also conductivity
  !   If present, Cond is created during call

    implicit none

    integer,intent(in)                          :: fidRM
    character(*), intent(in)                   :: inputFile
    type (grid_t), intent(inout)             :: grid
    type (rscalar), intent(inout), optional     :: Cond
    character (len=80), intent(out), optional   :: paramType
    integer, intent(out), optional              :: ioerr

    real(kind=prec)                     :: origin(3)
    real(kind=prec),pointer,dimension(:)    :: res
    integer                                     :: whichLayer
    integer                                     :: ix,iy,iz,ip,i,j
    real(kind=prec)                     :: alpha = 3.
    integer                                     :: jj,Nx,Ny,Nz,NzEarth,NzAir
    integer                                     :: status

    character (len=80)                          :: ifValues = ''
    character (len=80)                          :: someChar = ''
    integer                                     :: jOne, jTwo
    logical                                     :: returnCond

    returnCond = present(Cond)
    paramType = ''

    ! Open file and read grid
    open(unit=fidRM,file=inputFile,status='old',ERR=9000,iostat=ioerr)

    read(fidRM,'(a80)') someChar
    read(someChar,*) Nx, Ny, NzEarth, nzAir, ifValues

    if (ifValues(1:6) /= 'VALUES') then
        write(0, *) 'Mapping not supported yet in:ReadGridInputRM'
        stop
    end if

	! By default assume 'LINEAR RHO' - Randie Mackie's linear resistivity format
    if (index(someChar,'LOGE')>0) then
       paramType = 'LOGE'
    else
       paramType = 'LINEAR'
    end if

    call create_grid(Nx,Ny,NzAir,NzEarth,grid)

    ! In Randy Mackie's format, dx is read forwwards, as is dy and dz
    read(fidRM,*) (grid%dx(ix),ix=1,grid%nx)
    read(fidRM,*) (grid%dy(iy),iy=1,grid%ny)
    read(fidRM,*) (grid%dz(iz),iz=grid%nzAir+1,grid%nzAir+grid%nzEarth)

    !   Following is Kush's approach to setting air layers:
    ! mirror imaging the dz values in the air layer with respect to
    ! earth layer as far as we can using the following formulation
    ! air layer(bottom:top) = (alpha)^(j-1) * earth layer(top:bottom)
    i = grid%nzAir+1
    j = 0
    do iz = grid%nzAir, 1, -1
        j = j + 1
        grid%dz(iz) = ((alpha)**(j-1))*grid%dz(i)
        i = i + 1
    end do

    ! the topmost air layer has to be atleast 30 km
    if (grid%dz(1).lt.30000) then
        grid%dz(1) = 30000
    end if

    if(returnCond) then
       call create_rscalar(grid,Cond,CELL_EARTH)
    endif

    allocate(res(grid%Nx))
    do iz = 1,grid%nzEarth
       read(fidRM, *) whichLayer
       do iy = 1,grid%ny
          ! in Randy Mackie's format, x varies the fastest
          read(fidRM,*) res
          if(returnCond) then
            if(index(paramType,'LOGE')>0) then
              Cond%v(:,iy,iz) = - res
            else
              Cond%v(:,iy,iz) = 1./res
            endif
          endif
       enddo     !iy
    enddo        !iz
    deallocate(res)

    ! skip the three lines: a) WINGLINK, b) site name, and c) block numbers
    ! read WINGLINK (it can also be a blank line).
    read(fidRM, *, IOSTAT = ioerr) someChar
    ! a) WINGLINK
    if (ioerr /= 0) then
        if (index(someChar,'WINGLINK')>0) then
           write(0, *) 'Model file created by Winglink'
        else if (index(someChar,'MATLAB')>0) then
           write(0, *) 'Model file created by Matlab'
        else if (index(someChar,'ModEM')>0) then
           write(0, *) 'Model file created by ModEM'
        end if
    end if

    someChar = ''
    ! b) site name
    read(fidRM, *, IOSTAT = ioerr) someChar
    ! c) the block numbers
    read(fidRM, *, IOSTAT = ioerr) jOne, jTwo

    read(fidRM, *, IOSTAT = ioerr) grid%ox, grid%oy, grid%oz
    ! the defualt from a file read through Randy Mackie's format
    ! in Randy Mackie's format, real coordinates are in kilometers
    ! defaults in case of missing data
    if (ioerr /= 0) then
        grid%ox = 0.0
        grid%oy = 0.0
        grid%oz = 0.0
    else
       grid%ox = grid%ox*1000.0
       grid%oy = grid%oy*1000.0
       grid%oz = grid%oz*1000.0
    endif

    read(fidRM, *, IOSTAT = ioerr) grid%rotdeg
    if (ioerr /= 0) then
        grid%rotdeg = 0.0
    end if

    CLOSE(fidRM)

    GOTO 9999

9000 CONTINUE
    WRITE(0,*) '!!! FILE CANNOT BE FOUND !!!'
    STOP

9999 CONTINUE

  end subroutine readRMgridCond

 ! ***************************************************************************
  subroutine writeRMgridCond(fidRM,cFile,Cond,grid,paramType,ioerr)
  ! this routine writes files in Mackie's 3D formats, including the basic
  !   grid components and conductivity

    implicit none

    integer,intent(in)                          :: fidRM
    character(*), intent(in)                   :: cFile
    type (grid_t), intent(inout)             :: grid
    type (rscalar), intent(inout)     			:: Cond
    character(*), intent(in)         			:: paramType
    integer, intent(out)                        :: ioerr
    real(kind=prec), allocatable       :: value(:,:,:)
    integer                                    :: i,j,k
    integer                                    :: Nx,Ny,NzEarth,Nza,Nz
    logical                                    :: newFile

      if (len_trim(cfile)>0) then
         newFile = .true.
         open(unit=fidRM,file=cfile,status='unknown',iostat=ioerr)
      end if

      ! write grid geometry definitions
      Nx=grid%nx
      Ny=grid%ny
      Nz = grid%nz
      NzEarth=grid%nz - grid%nzAir
      Nza = grid%nzAir
      allocate(value(Nx,Ny,NzEarth))

      write(fidRM,'(4i5,a8)',advance='no') Nx,Ny,NzEarth,Nza,'VALUES'
      write(fidRM,*) trim(paramType)

      ! write grid spacings
      do j=1,Nx
      	write(fidRM,'(f12.3)',advance='no') grid%Dx(j)
	  end do
	  write(fidRM, *)
	  do j=1,Ny
      	write(fidRM,'(f12.3)',advance='no') grid%Dy(j)
	  end do
	  write(fidRM, *)
	  do j=Nza+1,Nz
      	write(fidRM,'(f12.3)',advance='no') grid%Dz(j)
	  end do
	  write(fidRM, *)

      ! convert from conductivity to resistivity
      if (index(paramType,'LOGE')>0) then
      	value = - Cond%v
      else if (index(paramType,'LINEAR')>0) then
       	value = ONE/Cond%v
      else
        ! assume resistivity and do nothing
 	  end if

      ! write out resistivity values
      do k=1,nzEarth
      	write(fidRM,'(i5)') k
      	do j=1,Ny
        	write(fidRM,'(a2)',advance='no') '  '
        	do i=1,Nx
      			write(fidRM,'(es13.5)',advance='no') value(i,j,k)
      		end do
      		write(fidRM, *)
      	end do
      end do

      ! written by ModEM
      write(fidRM, *, iostat = ioerr) 'ModEM'

      ! some crap for Winglink: site name and block numbers
      write(fidRM, *, iostat = ioerr) 'site'
      write(fidRM, *, iostat = ioerr) '1 1'

      ! origin
      write(fidRM, '(3f12.3)', iostat = ioerr) grid%ox/1000.0, grid%oy/1000.0, grid%oz/1000.0

 	  ! rotation
      write(fidRM,  '(f12.3)', iostat = ioerr) grid%rotdeg

      if (newFile) then
         close(fidRM)
      end if
      deallocate(value)

  end subroutine writeRMgridCond
