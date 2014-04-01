! Binary IO routines for 3D_MT modelParam (Gary Egbert's format)

  !******************************************************************
   subroutine write_modelParam_binary(m,cfile)

   !  open cfile on unit ioPrm, writes out object of
   !   type modelParam in standard *binary* format (comparable
   !   format written by write_Cond3D), then close file
   ! NOTE: uses grid for the model size; needs to be thought
   !       about if we decide to decouple grid and the model

      character(*), intent(in)  :: cfile
      type(modelParam_t), intent(in)              :: m
       integer               :: ios
      ! local
      type(rscalar)           :: ccond
      type(grid_t)          :: grid
      real(kind=prec) :: AirCond
      character(80)           :: paramType

      paramType = LOGE

      call getValue_modelParam(m,paramType,ccond,AirCond)

      grid = ccond%grid

      open(unit=ioPrm, file=cfile, form='unformatted', iostat=ios)
      write(ioPrm) grid%nx,grid%ny,grid%nzEarth,grid%nzAir
      write(ioPrm) grid%dx
      write(ioPrm) grid%dy
      write(ioPrm) grid%dz(grid%NzAir+1:grid%Nz)
      write(ioPrm) paramType
      write(ioPrm) AirCond
      write(ioPrm) ccond%v
      close(ioPrm)

      call deall_rscalar(ccond)
      call deall_grid(grid)

      end subroutine write_modelParam_binary

  !******************************************************************
   subroutine read_modelParam_binary(grid,m,cfile)

   !  open cfile on unit ioPrm, writes out object of
   !   type modelParam in standard *binary* format (comparable
   !   format written by write_Cond3D), then close file

      type(grid_t), target, intent(inout)  :: grid
      type(modelParam_t), intent(out)      :: m
      character(*), intent(in) 		       :: cfile
      integer                  :: ios

      integer		:: NzAir,Nz,Nx,Ny,NzEarth
      real(kind=prec) 	:: AirCond
      type(rscalar)             :: ccond
      character(80)             :: gridType,paramType

      open(unit=ioPrm, file=cfile, form='unformatted',status='old')

      read(ioPrm) Nx,Ny,NzEarth,NzAir

      call create_grid(Nx,Ny,NzAir,NzEarth,grid)
      read(ioPrm)  grid%dx
      read(ioPrm)  grid%dy
      read(ioPrm)  grid%dz(NzAir+1:NzAir+NzEarth)
      call setup_grid(grid)

      read(ioPrm) paramType

      gridType = CELL_EARTH
      call create_rscalar(grid,ccond,gridType)

      read(ioPrm) AirCond
      read(ioPrm) ccond%v

      close(ioPrm)

      ! move cell conductivities read into rscalar object into a modelParam
	  ! object ... this dance needed to keep modelParam attributes private
	  !   First need to create model parameter
	  call create_modelParam(grid,paramType,m,ccond,AirCond)

	  ! convert modelParam to the required paramType
	  call setType_modelParam(m,paramType)

	  ! now done with ccond, so deallocate
	  call deall_rscalar(ccond)

      end subroutine read_modelParam_binary
!******************************************************************
      subroutine write_JT_multi_vec_binary(m,header,cfile)

      !  open cfile on unit ioPrm, writes out one objects of
      !   type modelParam , closes file

      character(*), intent(in)          :: cfile, header
      type(modelParam_t), intent(in)      :: m
      integer                :: ios

      integer 		:: i,Nz,NzAir
      logical		:: opened


      inquire(file=cfile, opened=opened)
      if (.not. opened) then
      	open(unit=ioSens, file=cfile, form='unformatted', iostat=ios)
      endif

         Nz = m%grid%Nz
         NzAir = m%grid%Nz-m%NzEarth
         write(ioSens) m%paramType
         write(ioSens) m%Nx,m%Ny,m%NzEarth
         write(ioSens) m%grid%dx
         write(ioSens) m%grid%dy
         write(ioSens) m%grid%dz(NzAir+1:Nz)
         write(ioSens) m%AirCond
         write(ioSens) m%CellCond%v

      end subroutine write_JT_multi_vec_binary
     !******************************************************************
      subroutine writeVec_modelParam_binary(nSigma,m,header,cfile)

      !  open cfile on unit ioPrm, writes out nSigma objects of
      !   type modelParam , closes file

      integer, intent(in)               :: nSigma
      character(*), intent(in)          :: cfile, header
      type(modelParam_t), intent(in)      :: m(nSigma)
      integer                :: ios

      integer 		:: i,Nz,NzAir
      character(80) :: temp
      logical		:: opened

      temp = header

      inquire(file=cfile, opened=opened)
      if (.not. opened) then
      	open(unit=ioSens, file=cfile, form='unformatted', iostat=ios)
      endif
      !write(ioSens) temp
      !write(ioSens) nSigma

      do i = 1,nSigma
         Nz = m(i)%grid%Nz
         NzAir = m(i)%grid%Nz-m(i)%NzEarth
         write(ioSens) m(i)%paramType
         write(ioSens) m(i)%Nx,m(i)%Ny,m(i)%NzEarth
         write(ioSens) m(i)%grid%dx
         write(ioSens) m(i)%grid%dy
         write(ioSens) m(i)%grid%dz(NzAir+1:Nz)
         write(ioSens) m(i)%AirCond
         write(ioSens) m(i)%CellCond%v
      enddo
      ! for efficiency, will need to close the file in the calling program`
      !close(ioSens)
      end subroutine writeVec_modelParam_binary

     !******************************************************************
      subroutine readVec_modelParam_binary(grid,nSigma,m,header,cfile)

      !  open cfile on unit ioSens, writes out nSigma objects of
      !   type modelParam , closes file
      !  note - the code will not work unless the global grid already exists
      !         (it has to be initialized outside of this routine)

      type(grid_t), target, intent(inout)      :: grid
      integer, intent(in)              :: nSigma
      character(*), intent(in)          :: cfile
      character(80), intent(out)        :: header
      type(modelParam_t), intent(inout)   :: m(nSigma)
      integer                :: ios

      integer 		:: i,nS,Nx,Ny,Nz,NzAir,NzEarth,istat

      open(unit=ioSens, file=cfile, form='unformatted', iostat=ios)
      read(ioSens) header
      read(ioSens) nS

      if(nS .NE. nSigma) then
          call errStop('size of sigma does not agree with contents of file in readAll_Cond2D')
      endif

      do i = 1,nSigma
         read(ioSens) m(i)%paramType
         read(ioSens) Nx,Ny,NzEarth
         if((m(i)%Nx .NE. Nx) .OR. (m(i)%Ny .NE. Ny) .OR. (m(i)%NzEarth .NE. NzEarth)) then
            close(ioSens)
            call errStop('Size of cond does not agree with contents of file in readAll_Cond2D')
         else
            Nz = m(i)%grid%Nz
            NzAir = m(i)%grid%Nz-m(i)%NzEarth
            read(ioSens) m(i)%grid%dx
            read(ioSens) m(i)%grid%dy
            read(ioSens) m(i)%grid%dz(NzAir+1:Nz)
            read(ioSens) m(i)%AirCond
            read(ioSens) m(i)%CellCond%v
         endif
         m(i)%grid => grid
      enddo
      close(ioSens)
      end subroutine readVec_modelParam_binary
