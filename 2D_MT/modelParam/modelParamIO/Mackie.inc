! I/O routines for 2D_MT modelParam - extended Randie Mackie's format

      !******************************************************************
      subroutine write_modelParam_mackie(cond,cfile)

      !  open cfile on unit ioPrm, writes out object of
      !   type modelParam in Randie Mackie's format, closes file
      ! note that while this is in general a good interface for
      ! a write subroutine for modelParam, we do not have to write
      ! out the grid unless we want to do so

      type(modelParam_t), intent(in)	:: cond
      character(*), intent(in)          :: cfile
      integer                           :: ios
      ! local
      type(grid_t)              :: grid
      real (kind=prec), allocatable    :: value(:,:)
      Integer                     ::k,j,Nz,Nza,NzEarth,Ny
      logical                     :: newFile
      character(80)               :: paramType

      if (len_trim(cfile)>0) then
         newFile = .true.
         open(unit=ioPrm,file=cfile,status='unknown',iostat=ios)
      end if

      !  Read in grid geometry definitions, store in structure TEgrid
      !    first grid dimensions ...
      grid = cond%grid
      Ny=grid%ny
      Nz=grid%nz
      NzEarth=grid%nz - grid%nza
      Nza = grid%nza
      allocate(value(Ny,NzEarth))

      ! extract conductivity values from the modelParam structure
      call getValue_modelParam(cond,value,paramType)

      write(ioPrm,'(2i5)',advance='no',iostat=ios) Ny,NzEarth
      write(ioPrm,*,iostat=ios) trim(paramType)

      !    write grid spacings: NOTE that Randie Mackie inserts an empty
      !    line after every 100 values in the grid (i.e. 10 lines);
      !    we do not do that here to save extra coding

        write(ioPrm,'(10g11.4)',iostat=ios) (grid%Dy(j),j=1,Ny)

        write(ioPrm,'(10g11.4)',iostat=ios) (grid%Dz(j),j=Nza+1,Nz)

      !call write_grid(ioPrm,'',grid)
      !Ny = grid%Ny
      !NzEarth = grid%Nz-grid%Nza


      ! convert from conductivity to resistivity
      if (index(paramType,'LOGE')>0) then
      	value(:,:) = - value(:,:)
      else if (index(paramType,'LINEAR')>0) then
       	value(:,:) = ONE/value(:,:)
      else
        ! assume resistivity and do nothing
 	  end if

      ! write out resistivity values
      write(ioPrm,'(a1)') '0'
      do j=1,NzEarth
        write(ioPrm,'(a2)',advance='no',iostat=ios) '  '
      	write(ioPrm,'(19es13.5)',iostat=ios) (value(k,j), k=1,Ny)
      end do

      if (newFile) then
         close(ioPrm)
      end if
      deallocate(value)
      call deall_grid(grid)

      end subroutine write_modelParam_mackie

      !******************************************************************
      subroutine read_modelParam_mackie(grid,cond,cfile)

      !  open cfile on unit ioPrm, writes out object of
      !   type modelParam in Randie Mackie's format, closes file
      ! we can update the grid here, but the grid as an input is critical
      ! for setting the pointer to the grid in the modelParam

      type(grid_t), target, intent(inout)       :: grid
      type(modelParam_t), intent(out)	:: cond
      character(*), intent(in)          :: cfile
      integer                           :: ios
      ! local
      real (kind=prec), allocatable    :: value(:,:)
      Integer                     ::k,j,Nz,Nza,NzEarth,Ny
      logical                     :: newFile
      real (kind=prec)    :: airCond
      character(80)               :: line,paramType

      if (len_trim(cfile)>0) then
         newFile = .true.
         open(unit=ioPrm,file=cfile,status='old',iostat=ios)
      end if

           ! We are using Randie Mackie's format, which does not have information
      ! about the air layers. So we make it equal 10 in this routine.
      Nza = 10

      !  Read in grid geometry definitions, store in structure TEgrid
      !    first grid dimensions ...
      !  If paramType == 'LOGE', the file contains natural log resistivity,
      !  otherwise it's resistivity
      read(ioPrm,'(a80)',iostat=ios) line

      read(line,*,iostat=ios) Ny,NzEarth

      ! By default assume 'LINEAR RHO' - Randie Mackie's linear resistivity format
      if (index(line,'LOGE')>0) then
         paramType = 'LOGE'
      else
         paramType = 'LINEAR'
      end if

      Nz = NzEarth + Nza
      ! then allocate for grid
      call create_grid(Ny,Nz,Nza,grid)

      !    read in grid spacings

        read(ioPrm,*,iostat=ios) (grid%Dy(j),j=1,Ny)

        read(ioPrm,*,iostat=ios) (grid%Dz(j),j=Nza+1,Nz)

      ! set the air layers to something quasi-reasonable: in the future,
      ! will include them in the input ascii file
      !grid%Dz(Nza) = max(grid%Dz(Nza+1),10.0d0)
	  grid%Dz(Nza) = grid%Dz(Nza+1)
      do j = 1,Nza-1
      	grid%Dz(Nza-j) = 3.0 * grid%Dz(Nza-j+1)
      end do

      !call read_grid(ioPrm,'',grid)
            ! complete grid definition
      call setup_grid(grid)

      !Ny = grid%Ny
      !NzEarth = grid%Nz-grid%Nza
      allocate(value(Ny,NzEarth))

      ! read in resistivity values
      read(ioPrm,*)
      do j=1,NzEarth
      	read(ioPrm,*,iostat=ios) (value(k,j), k=1,Ny)
      end do

      ! convert to log conductivity
      if (trim(paramType) == 'LOGE') then
        value(:,:) = - value(:,:)
      else if (trim(paramType) == 'LINEAR') then
      	value(:,:) = ONE/value(:,:)
      else
        call errStop('Unable to understand the model parameter type in read_cond2d')
      end if

      ! set natural log of air conductivity
      airCond = log(SIGMA_AIR)

      ! write into the modelParam structure (according to paramType)
      call create_modelParam(grid,paramType,cond,value,airCond)

      if (newFile) then
         close(ioPrm)
      end if
      deallocate(value)

      end subroutine read_modelParam_mackie
