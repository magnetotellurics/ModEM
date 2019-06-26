! Binary IO routines for 2D_MT modelParam (Gary Egbert's format)

      !******************************************************************
      subroutine write_modelParam_binary(cond,cfile)

      !  open cfile on unit ioPrm, writes out object of
      !   type modelParam in standard format (readable by matlab
      !   routine readCond_2D.m), closes file

      type(modelParam_t), intent(in)	:: cond
      character(*), intent(in)          :: cfile
      integer                           :: ios
      ! local
      type(grid_t)              :: grid
      real (kind=prec), allocatable :: v(:,:)
      Integer                     ::k,j,Nz,Nza,NzEarth,Ny,stat
      character*80                :: paramType
      real (kind=prec)    :: AirCond

      grid = cond%grid
      Ny = grid%Ny
      NzEarth = grid%Nz-grid%Nza
      allocate(v(Ny,NzEarth))

      ! extract conductivity values from the modelParam structure
      call getValue_modelParam(cond,v,paramType,AirCond)

      open(unit=ioPrm, file=cfile, form='unformatted',iostat=ios)
      ! write grid dimensions ...
      write(ioPrm) grid%Ny,grid%Nz,grid%Nza

      ! write grid spacings
      write(ioPrm) grid%Dy
      write(ioPrm) grid%Dz

      ! write model parameter
      write(ioPrm) paramType
      write(ioPrm) Ny,NzEarth
      write(ioPrm) v
      write(ioPrm) AirCond
      close(ioPrm)

      deallocate(v)
      call deall_grid(grid)

      end subroutine write_modelParam_binary

     !******************************************************************
      subroutine read_modelParam_binary(grid,cond,cfile)

      !  open cfile on unit ioPrm, read in object of
      !   type modelParam in standard format
      ! we can update the grid here, but the grid as an input is critical
      ! for setting the pointer to the grid in the modelParam

      type(grid_t), target, intent(inout)       :: grid
      type(modelParam_t), intent(out)	:: cond
      character(*), intent(in)          :: cfile
      integer                           :: ios
      ! local

      real (kind=prec), allocatable :: v(:,:)

      !  local variables
      integer 		:: Ny, Nz, Nza, NzEarth
      character*80	:: paramType
      real (kind=prec) :: AirCond

      open(unit=ioPrm, file=cfile, form='unformatted',status='old', iostat=ios)
      !  Read in grid geometry definitions, store in structure TEgrid
      !    first grid dimensions ...
      read(ioPrm) Ny,Nz,Nza

      ! then allocate for grid
      call create_grid(Ny,Nz,Nza,grid)

      !    read in grid spacings
      read(ioPrm) grid%Dy
      read(ioPrm) grid%Dz

      ! complete grid definition
      call setup_grid(grid)

      read(ioPrm) paramType
      read(ioPrm) Ny,NzEarth

      if((grid%Ny .NE. Ny).OR.(grid%Nz - grid%Nza .NE. NzEarth)) then
         close(ioPrm)
         call errStop('Grid size does not agree with model parameter size in readCond_2D')
      end if

      ! deallocate modelParam if it is allocated
      call deall_modelParam(cond)

      ! read info from file
      allocate(v(Ny,NzEarth))
      read(ioPrm) v
      read(ioPrm) AirCond
      close(ioPrm)

      ! write into the modelParam structure
      call create_modelParam(grid,paramType,cond,v,AirCond)

      ! clean up
      deallocate(v)

      end subroutine read_modelParam_binary
!******************************************************************
      subroutine write_JT_multi_vec_binary(sigma,header,cfile)

      !  open cfile on unit ioPrm, writes out nSigma objects of
      !   type modelParam , closes file

      character(*), intent(in)		:: header, cfile
      type(modelParam_t), intent(in)	:: sigma
      integer      :: ios

      integer i
      character(80) temp
      logical opened

      temp = header

      inquire(file=cfile, opened=opened)
      if (.not. opened) then
      	open(unit=ioSens, file=cfile, form='unformatted', iostat=ios)
      endif



         write(ioSens) sigma%paramType
         write(ioSens) sigma%Ny,sigma%NzEarth
         write(ioSens) sigma%v
         write(ioSens) sigma%AirCond

      end subroutine write_JT_multi_vec_binary
      !******************************************************************
      subroutine writeVec_modelParam_binary(nSigma,sigma,header,cfile)

      !  open cfile on unit ioPrm, writes out nSigma objects of
      !   type modelParam , closes file

      integer, intent(in)		:: nSigma
      character(*), intent(in)		:: header, cfile
      type(modelParam_t), intent(in)	:: sigma(nSigma)
      integer      :: ios

      integer i
      character(80) temp
      logical opened

      temp = header

      inquire(file=cfile, opened=opened)
      if (.not. opened) then
      	open(unit=ioSens, file=cfile, form='unformatted', iostat=ios)
      endif
      !write(ioSens) temp
      !write(ioSens) nSigma

      do i = 1,nSigma
         write(ioSens) sigma(i)%paramType
         write(ioSens) sigma(i)%Ny,sigma(i)%NzEarth
         write(ioSens) sigma(i)%v
         write(ioSens) sigma(i)%AirCond
      enddo
      ! for efficiency, will need to close the file in the calling program`
      ! close(ioSens)
      end subroutine writeVec_modelParam_binary

     !******************************************************************
      subroutine readVec_modelParam_binary(grid,nSigma,sigma,header,cfile)

      !  open cfile on unit ioPrm, read nSigma objects of
      !   type modelParam , closes file
      !  sigma(nsigma) must be allocated before calling
      !  note - the code will not work unless the global grid already exists
      !         (it has to be initialized outside of this routine)

      type(grid_t), target, intent(inout)   :: grid
      integer, intent(in)		:: nSigma
      type(modelParam_t), intent(inout) 	:: sigma(nSigma)
      character(*), intent(out)		:: header
      character(*), intent(in)		:: cfile
      integer		:: ios

      ! local variables
      integer i, nS, Ny,NzEarth

      open(unit=ioSens, file=cfile, form='unformatted', iostat=ios)
      read(ioSens) header
      read(ioSens) nS
      if(nS .NE. nSigma) then
          call errStop('size of sigma does not agree with contents of file in readAll_Cond2D')
      endif

      do i = 1,nSigma
         read(ioSens) sigma(i)%paramType
         read(ioSens) Ny,NzEarth
         if((sigma(i)%Ny .NE. Ny).OR. (sigma(i)%NzEarth .NE. NzEarth)) then
            close(ioSens)
            call errStop('Size of cond does not agree with contents of file in readAll_Cond2D')
         else
            read(ioSens) sigma(i)%v
            read(ioSens) sigma(i)%AirCond
         endif
         sigma(i)%grid => grid
      enddo
      close(ioSens)

      end subroutine readVec_modelParam_binary
