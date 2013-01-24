!   Derived from file_io.f90, with some additions to make
!   this more like Modular2D IO, but at the same time to allow
!   previously used 3D files formats to be handled.
!   Modified by Egbert, Nov. 2007, following previous modifications
!    to codes originally written by Kush Tandon.
!
!   Comments following are from file_io:
!   CHANGES FROM KUSH VERION:
!    1) Input is input, does not ever compare to
!   expected array sizes etc.  (as it was the input routine can only be
!    used when you already know sizes; it's useful to be able to just read the
!    file;  any consistency checking is kept separate.
!    2)  comments cleaned up
!    3) nMode is fixed at 2 ... if we get more general, this will all
!      have to be changed anyway.
!    4)  order of arguments changed most places
!    5)  got rid of routine for closing files: one executable stmt, make no sense
!    6)  BC IO only supported for BC data type: these routines don't
!        need to know about type excite
!    7) eliminate skipping and rewinding for output routines: in general
!         this makes no sense for writing a sequential binary file!
!      NOTE: we might want to at least allow the option to just read the next
!       records, assuming that the file is positioned correctly (as it will
!        be in all applications we now have)
!    8) format of impedance file changed: all site locations are in a second header
!          record, all impedances for one frequency are in one sequential record

! *****************************************************************************
module ioAscii
  ! This module contains io routines for boundary nodes,
  ! electric field solutions, impedances, and solver diagnostics

  use griddef
  use math_constants
  use emsolve3d
  use dataspace
  use datafunc ! for data types
  use ForwardSolver, only: solnVectorMTX_t
  use transmitters
  use receivers
  use datatypes
  use sg_vector_mg, only: mg2c, real_cvector_mg_f, imag_cvector_mg_f, deall_cvector_mg

  implicit none

  public                   :: FileReadInit, FileWriteInit
  public                   :: ZfileReadInit, ZfileWriteInit
  public                   :: BCfileRead, BCfileWrite
  public                   :: EfileRead, EfileWrite
  public                   :: DfileWrite
  public                   :: ZfileRead, ZfileWrite
  public                   :: write_solnVectorMTX

Contains

  ! ***************************************************************************
  ! open file, read header record; used for reading header of
  !     all output files except impedance ... this just reads, returns
  !     what is in header.  Check consistency elsewhere if desired.
  !    If you don't want to check consistency, or get the header info,
  !    just open the file directly, don't use this routine!
   subroutine FileReadInit(inFile,ioNum,fileGrid,filePer,fileMode,fileVersion,ios)

    implicit none
    !  file name and iouint
    character (len=80), intent(in)              :: inFile
    integer, intent(in)                 	:: ioNum

    !  file header contents
    character (len=20), intent(out)             :: fileVersion
    type (grid_t), intent(inout)             :: fileGrid
    integer, intent(out)                        :: filePer
    integer, intent(out)                        :: fileMode
    integer, intent(out)                        :: ios

    open (unit=ioNum,file=inFile,status='old', form ='unformatted',iostat=ios)

    IF( ios/=0) THEN
       WRITE(0,*) 'Error opening file: ', inFile
    else
       ! read the header
       read(ioNum) fileVersion, filePer, fileMode, &
	fileGrid%nx, fileGrid%ny, fileGrid%nz, fileGrid%nzAir, &
	fileGrid%ox, fileGrid%oy, fileGrid%oz, fileGrid%rotdeg

       ! the basic grid information is also stored in the header
       ! (don't we need to allocate first???)
       ! NM: Yes we do

       call create_grid(fileGrid%nx,fileGrid%ny,fileGrid%nzAir,fileGrid%nz-fileGrid%nzAir,fileGrid)
       
       read(ioNum) fileGrid%dx
       read(ioNum) fileGrid%dy
       read(ioNum) fileGrid%dz

    endif
    return
    end subroutine FileReadInit

  !**********************************************************************
  ! check consistency between grids, number of periods
  subroutine CheckgridConsist(runGrid, fileGrid, nPer, filePer, ier)
    type (grid_t), intent(in)         :: fileGrid, runGrid
    integer, intent(in)			:: nPer, filePer
    integer, intent(out)		:: ier
    logical				:: ok = .true.

    ier = 0
    ok = ok .and. fileGrid%nx .eq. runGrid%nx
    ok = ok .and. fileGrid%ny .eq. runGrid%ny
    ok = ok .and. fileGrid%nz .eq. runGrid%nz
    ok = ok .and. fileGrid%nzAir .eq. runGrid%nzAir
    if(.not. ok) then
       write(0,*) 'Mismatch in grid size '
       ier = -1
    endif
    ok = ok .and. fileGrid%ox .eq. runGrid%ox
    ok = ok .and. fileGrid%oy .eq. runGrid%oy
    ok = ok .and. fileGrid%oz .eq. runGrid%oz
    if(.not. ok) then
       write(0,*) 'Mismatch in grid origin'
       ier = -1
    endif

    if(fileGrid%rotdeg /= runGrid%rotdeg) then
       write(0,*) 'Mismatch in grid rotation'
       ier = -1
    endif

    if(nPer /= filePer) then
       write(0,*) 'Mismatch in number of periods'
       ier = -1
    endif

  end subroutine CheckGridConsist

  ! ***************************************************************************
  !  open MT impedance file and read header (no reason for consistency check!)
  subroutine ZfileReadInit(inFile, ioNum, nPer, nSites, sites, ifBzTF,ios)

    implicit none
    !  input file name and io unit number
    character (len=80), intent(in)		:: inFile
    integer, intent(in)				:: ioNum

    ! output info from header
    integer, intent(out)			:: nPer
    integer, intent(out)			:: nSites
    integer, intent(out)			:: ios
    logical, intent(out)			:: ifBzTF
    real(kind=prec), pointer, dimension(:,:) :: sites

    integer		:: k

    open (unit=ioNum,file=inFile,status='old', form ='unformatted',&
         iostat=ios)

    if( ios/=0) then
       write(0,*) 'Error opening file in ZfileReadInit: ', inFile
    endif

    ! read the header  ... why no version here?????
    read(ioNum) nPer, nSites, ifBzTF

    !  read in site locations
    allocate(sites(3,nSites))
    read(ioNum) (sites(:,k),k=1,nSites)

  end subroutine ZfileReadInit

  ! ***************************************************************************
  ! write header for MT_fwd output files (except impedance file)
  subroutine FileWriteInit(version, outFile, ioNum, inGrid, nPer, nMode,ios)

    implicit none
    character (len=20), intent(in)		:: version
    ! the version of the software in use
    character (len=80), intent(in)		:: outFile
    ! the filename
    type(grid_t), intent(in)			:: inGrid
    integer, intent(in)				:: nPer
    integer, intent(in)				:: nMode
    integer, intent(in)				:: ioNum
    integer, intent(out)			:: ios

    open (unit=ioNum,file=outFile,STATUS='unknown', form ='unformatted',&
         iostat=ios)

    if( ios/=0) then
       write(0,*) 'Error opening file in FileWriteInit: ', outFile
    else
       ! write the header (contains the basic information for the forward
       ! modeling). the header is 4 lines
       write(ioNum) version,nPer,nMode,inGrid%nx,inGrid%ny,inGrid%nz,inGrid%nzAir, &
       inGrid%ox, inGrid%oy, inGrid%oz, inGrid%rotdeg
       write(ioNum) inGrid%dx
       write(ioNum) inGrid%dy
       write(ioNum) inGrid%dz
    endif
  end subroutine FileWriteInit

  ! ***************************************************************************
  ! write header for MT impedance file ... why no version?????
  subroutine ZfileWriteInit(inFile, ioNum, nPer, nSites, sites, ifBzTF,ios)

    implicit none
    ! file name and io unit
    character (len=80), intent(in)		:: inFile
    integer, intent(in)				:: ioNum

    ! header contents
    integer, intent(in)				:: nPer,nSites
    logical, intent(in)				:: ifBzTF
    real(kind=prec), intent(inout)	:: sites(3,nSites)

    !  iostatus
    integer, intent(out)			:: ios

    ! local variables
    integer	:: k

    open (unit=ioNum,file=inFile,status='unknown', form ='unformatted',&
         iostat=ios)

    if( ios/=0) then
       write(0,*) 'Error opening file in ZfileWriteInit: ', inFile
    else
       ! write the header
       write(ioNum) nPer, nSites, ifBzTF
    endif

    write(ioNum) (sites(:,k),k=1,nSites)

  end subroutine ZfileWriteInit

  ! ***************************************************************************
  ! read BCs for one mode, frequency; open unit ioNum first ...
  subroutine BCfileRead(ioNum, ifreq, imode, fileOmega, inBC)

    implicit none
    integer, intent(in)				:: ioNum,ifreq,imode

    real (kind=prec), intent(out)			:: fileOmega
    type (cboundary), intent(inout)		:: inBC

    !  local variables
    integer					:: nRecSkip, iskip
!    integer					:: iskip

    !  following could be optional oututs
    integer				:: fileIfreq, fileMode
    character (len = 20)            	:: ModeName

    !  hard code number of modes for now
    integer		:: nMode = 2

    ! calculate number of records before the header for this frequency/mode
    nRecSkip = ((ifreq-1)*nMode+(imode-1))*13+4

    ! rewind file, and skip to header record
    rewind(ioNum)
    do iskip = 1,nRecSkip
       read(ioNum)
    enddo

    ! read the frequency header - 1 record
    read(ioNum) fileOmega, fileIfreq, fileMode, ModeName

    ! read BC data for one frequency/mode - 12 records
    read(ioNum) inBC%xYMax
    read(ioNum) inBC%zYMax
    read(ioNum) inBC%xYMin
    read(ioNum) inBC%zYMin
    read(ioNum) inBC%yXMax
    read(ioNum) inBC%zXMax
    read(ioNum) inBC%yXMin
    read(ioNum) inBC%zXMin
    read(ioNum) inBC%xZMin
    read(ioNum) inBC%yZMin
    read(ioNum) inBC%xZMax
    read(ioNum) inBC%yZMax

  end subroutine BCfileRead

  ! ***************************************************************************
  ! write BCs for one mode, frequency; open unit ioNum first ...
  !  Note: for a sequential binary file skipping around makes no sense in general
  !   Just write next record (only thing that works, as far as I know!)
  subroutine BCfileWrite(ioNum, Omega, iFreq, iMode, ModeName, outBC)

    implicit none
    real (kind=prec), intent(in)	:: Omega
    integer, intent(in)				:: ioNum
    integer, intent(in)				:: ifreq,imode
    character (len=20), intent(in)		:: ModeName
    type (cboundary), intent(inout)		:: outBC

    ! write the frequency header - 1 record
    write(ioNum) Omega, iFreq, iMode, ModeName

    ! write the BC data - 12 records
    write(ioNum) outBC%xYMax
    write(ioNum) outBC%zYMax
    write(ioNum) outBC%xYMin
    write(ioNum) outBC%zYMin
    write(ioNum) outBC%yXMax
    write(ioNum) outBC%zXMax
    write(ioNum) outBC%yXMin
    write(ioNum) outBC%zXMin
    write(ioNum) outBC%xZMin
    write(ioNum) outBC%yZMin
    write(ioNum) outBC%xZMax
    write(ioNum) outBC%yZMax

  end subroutine BCfileWrite

  ! ***************************************************************************
  ! read electric field solution for one mode, frequency; open unit ioNum first ...
  subroutine EfileRead(ioNum, ifreq, imode, fileOmega, inE)

    implicit none
    integer, intent(in)             :: ioNum,ifreq,imode
    type (cvector), intent(inout)       :: inE
    real (kind=prec), intent(out)   :: fileOmega

    !  local variables
    integer                 :: nRecSkip, iskip
!    integer                    :: iskip

    !  following could be optional oututs
    integer             :: fileIfreq, fileMode
    character (len = 20)                :: ModeName

    !  hard code number of modes for now
    integer     :: nMode = 2

    ! calculate number of records before the header for this frequency/mode
    nRecSkip = ((ifreq-1)*nMode+(imode-1))*4+4

    ! rewind the file, skip to header record
    rewind(ioNum)
    do iskip = 1,nRecSkip
       read(ioNum)
    enddo

    ! read frequency header - 1 record
    read(ioNum) fileOmega, fileIfreq, fileMode, ModeName

    ! read electrical field data - 3 records
    read(ioNum) inE%x
    read(ioNum) inE%y
    read(ioNum) inE%z

  end subroutine EfileRead


  ! ***************************************************************************
  ! write electrical field solution for one frequency/mode
  !    ... again no sense to skipping, rewinding for sequential binary write
  subroutine EfileWrite(ioNum,Omega, iFreq, iMode, ModeName, outE)

    implicit none
    real (kind=prec), intent(in)	:: Omega
    integer, intent(in)				:: ioNum,iFreq,iMode
    character (len=20), intent(in)		:: ModeName
    type (cvector), intent(in)			:: outE

    ! write the frequency header - 1 record
    write(ioNum) Omega, iFreq, iMode,ModeName

    ! write the electrical field data - 3 records
    write(ioNum) outE%x
    write(ioNum) outE%y
    write(ioNum) outE%z

  end subroutine EfileWrite

  ! ***************************************************************************
  ! read impedances for nSites for one frequency
  subroutine ZfileRead(ioNum,ifreq,ifBzTF,nSites,Z,fileOmega)

    implicit none
    integer, intent(in)				:: ioNum,iFreq
    logical, intent(in)				:: ifBzTF
    integer, intent(in)				:: nSites
    complex (kind=prec), intent(out)	:: Z(3,2,nSites)
    real (kind=prec), intent(out)	:: fileOmega

    integer					:: nRecSkip,fileIfreq
    integer					:: iskip,i,j,k

    ! rewind file, skip past records for other frequencies
    nRecSkip = ifreq*2
    rewind(ioNum)
    do iskip = 1,nRecSkip
       read(ioNum)
    enddo

    ! read the frequency header - 1 record
    read(ioNum) fileOmega, fileIfreq

    if (ifBzTF) then
        read(ioNum) Z
    else
        read(ioNum)(((Z(i,j,k),i=1,2),j=1,2),k=1,nSites)
    end if

  end subroutine ZfileRead

  ! ***************************************************************************
  ! writes impedances for nSites for one frequency
  subroutine ZfileWrite(ioNum,Omega,ifreq,nSites,Z,ifBzTF)
  !   NOTE: impedances are stored in real dataVec structure:
  !    first dimension is mode, second component (real/imag for
  !     two impedance components + BzTF  ... rows and
  !     columns are switched compared to usual impedance
  !     storage.  Output in the usual way


    implicit none
    integer, intent(in)				:: ioNum
    real (kind=prec), intent(in)	:: Omega
    integer, intent(in)				:: ifreq
    integer, intent(in)				:: nSites
    real (kind=prec), intent(in)	:: Z(2,6,nSites)
    logical, intent(in)				:: ifBzTF
    integer					:: i,j,k,nComp

    ! write the frequency header - 1 record
    write(ioNum) Omega, iFreq

    if (ifBzTF) then
       nComp = 6
    else
       nComp = 4
    endif

    write(ioNum) (((Z(i,j,k),j=1,nComp),i=1,2),k=1,nSites)

  end subroutine ZfileWrite

  ! ***************************************************************************
  ! writes solution diagnostics from unit opened by FileWriteInit. The prefix D
  ! is for diagnostics.
  subroutine DfileWrite(ioNum,Omega, iFreq, iMode, ModeName,solverDiag)

    implicit none
    real (kind=prec), intent(in)		:: Omega
    integer, intent(in)				:: ifreq,imode
    character (len=20), intent(in)		:: ModeName
    integer, intent(in)				:: ioNum
    type(emsolve_diag), intent(in)		:: solverDiag
    integer					:: nRecSkip, nMode
    integer					:: iskip
    integer					:: i, j

    ! write the frequency header - 1 record
    write(ioNum) Omega, iFreq, iMode, ModeName

    ! write the solution diagnostics - 3 records
    write(ioNum) solverDiag%nIterTotal, solverDiag%nDivCor
    write(ioNum) (solverDiag%EMrelErr(i),i = 1,solverDiag%nIterTotal)
    write(ioNum) ((solverDiag%divJ(i,j), i = 1,2), j = 1,solverDiag%nDivCor)

    !  why no output of divcor convergence?

  end subroutine DfileWrite ! DfileWrite

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!   Addtional routines, to allow greater compatability with Modular2D
!    code structure
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!******************************************************************
      subroutine read_solnVectorMTX(grid,eAll,cfile)

      ! reads an array of solution vectors for all transmitters & subgrids
      ! currently uses the old binary format; will switch to NetCDF when
      ! time allows

      character(*), intent(in)                    :: cfile
      type(solnVectorMTX_t), intent(inout)        :: eAll
      type(grid_t), intent(in), target            :: grid

      !   local variables
      integer           :: j,k,nMode = 2, ios,ig,cdot
      integer           :: iTx,nTx
      character (len=3)         :: igchar
      character (len=20)        :: version = '',ModeNames(2)
      character (len=200)       :: fn_input
      real (kind=prec)   :: omega

      ModeNames(1) = 'Ey'
      ModeNames(2) = 'Ex'

      nTx = size(txDict)

      call create_solnVectorMTX(nTx,eAll)
      do iTx=1,size(txDict)
         call create_solnVector(grid,iTx,eAll%solns(iTx))
      end do

      ! unformatted read for multigrid
      do ig = 1, grid%mgridSize

          write (igchar,'(i3.3)') ig
          cdot = index(cfile,'.')
          if (cdot > 0) then
            fn_input = cfile(1:cdot-1)//'_'//igchar//cfile(cdot:len_trim(cfile))
          else
            fn_input = trim(cfile)//'_'//igchar
          end if
          write(*,*) 'Reading E-fields for subgrid number ',ig,' from file: ',trim(fn_input)

          call FileReadInit(fn_input,ioE,eAll%solns(1)%grid%gridArray(ig),eAll%nTX,nMode,version,ios)
          do j = 1,nTx
             do k = 1,2
               omega = txDict(eAll%solns(j)%tx)%omega

               call EfileRead(ioE, j, k, omega, eAll%solns(j)%pol(k)%cvarray(ig))

             enddo
          enddo
          close(ioE)

      enddo

      end subroutine read_solnVectorMTX

!******************************************************************
      subroutine write_solnVectorMTX(eAll,cfile)

      !  open cfile on unit fid, writes out object of
      !   type cvector in standard format (readable by matlab
      !   routine readcvector.m), closes file
      !  NOT coded at present to specifically write out TE/TM
      !    solutions, periods, etc. (can get this infor from
      !    eAll%solns(j)%tx, but only with access to TXdict.

      character(*), intent(in)          :: cfile
      type(solnVectorMTX_t), intent(in)               :: eAll

      !   local variables
      integer           :: j,k,nMode = 2, ios,ig,cdot
      character (len=3)         :: igchar
      character (len=20) 		:: version = '',ModeNames(2)
      character (len=200)       :: fn_output
      real (kind=prec)   :: omega

      ModeNames(1) = 'Ey'
      ModeNames(2) = 'Ex'

      ! Maria C.: temporary output to be moved outside of this routine
      call plotFields(eAll)

      ! unformatted write for multigrid
      do ig = 1, eAll%solns(1)%grid%mgridSize

          write (igchar,'(i3.3)') ig

          cdot = index(cfile,'.')
          if (cdot > 0) then
            fn_output = cfile(1:cdot-1)//'_'//igchar//cfile(cdot:len_trim(cfile))
          else
            fn_output = trim(cfile)//'_'//igchar
          end if
          write(*,*) 'E-fields for subgrid number ',ig,' written to ',trim(fn_output)

          call FileWriteInit(version,fn_output,ioE,eAll%solns(1)%grid%gridArray(ig),eAll%nTX,nMode,ios)
          do j = 1,eAll%nTx
             do k = 1,2
               omega = txDict(eAll%solns(j)%tx)%omega

               call EfileWrite(ioE, omega, j,  k, ModeNames(k), eAll%solns(j)%pol(k)%cvarray(ig))

             enddo
          enddo
          close(ioE)

      enddo

      end subroutine write_solnVectorMTX

! *********************************************************************
      subroutine plotFields(eAll)

      ! created by Cherevatova (Nov,2012)
      ! to plot fields

      type(solnVectorMTX_t), intent(in)               :: eAll

      !   local variables
      integer  :: j, ig, iSurface,ix,iy,iz
      integer  :: nx,ny,Iox,Ioy

      real(kind=prec)    :: tempD
      type(cvector)    :: Etemp

     ! plot curve at the surface
     do ig = 1, eAll%solns(1)%grid%mgridSize
         Isurface = ig
       if (eAll%solns(1)%grid%gridarray(ig)%flag == 0)exit
     enddo
     ! iz =  eAll%solns(1)%grid%gridArray(Isurface)%nz-eAll%solns(1)%grid%gridArray(Isurface)%nzAir
      nx = eAll%solns(1)%grid%gridArray(Isurface)%nx
      ny = eAll%solns(1)%grid%gridArray(Isurface)%ny
      tempD= 0
      do ix = 1,nx-1
        tempD= tempD+eAll%solns(1)%grid%dx(ix)
          Iox = ix
        if (tempD == -(eAll%solns(1)%grid%ox)) exit
      enddo
      tempD= 0
      do iy = 1,ny-1
        tempD = tempD + eAll%solns(1)%grid%dy(iy)
          Ioy = iy
        if (tempD == -(eAll%solns(1)%grid%oy)) exit
      enddo
      Isurface = 2
      iz = 12
      open(01, file = 'Curve.dat')
        do j = 1,eAll%nTx   ! loop periods
          ! 1 - Ey
          ! 2 - Ex
          write(01,*)'Nmode, period, X, Y, Z, Ex, Ey, Ez'
                     !Mode Period  Orig Iy   Ex          Ey       Ez
          write(01,'(i2,es13.5,2x,i3,2x,i3,2x,es13.5,2x,es13.5,2x,es13.5,es13.5,2x,es13.5,2x,es13.5)') (1,txDict(eAll%solns(j)%tx)%period, Ioy, ix, &
                                 real(eAll%solns(j)%pol(1)%cvarray(Isurface)%x(ix,Ioy,iz)), imag(eAll%solns(j)%pol(1)%cvarray(Isurface)%x(ix,Ioy,iz)),&
                                 real(eAll%solns(j)%pol(1)%cvarray(Isurface)%y(ix,Ioy,iz)), imag(eAll%solns(j)%pol(1)%cvarray(Isurface)%y(ix,Ioy,iz)),&
                                 real(eAll%solns(j)%pol(1)%cvarray(Isurface)%z(ix,Ioy,iz)), imag(eAll%solns(j)%pol(1)%cvarray(Isurface)%z(ix,Ioy,iz)),ix=1,nx)
          write(01,'(i2,es13.5,2x,i3,2x,i3,2x,es13.5,2x,es13.5,2x,es13.5,es13.5,2x,es13.5,2x,es13.5)') (2,txDict(eAll%solns(j)%tx)%period,Ioy, ix, &
                                 real(eAll%solns(j)%pol(2)%cvarray(Isurface)%x(ix,Ioy,iz)), imag(eAll%solns(j)%pol(2)%cvarray(Isurface)%x(ix,Ioy,iz)), &
                                 real(eAll%solns(j)%pol(2)%cvarray(Isurface)%y(ix,Ioy,iz)), imag(eAll%solns(j)%pol(2)%cvarray(Isurface)%y(ix,Ioy,iz)), &
                                 real(eAll%solns(j)%pol(2)%cvarray(Isurface)%z(ix,Ioy,iz)), imag(eAll%solns(j)%pol(2)%cvarray(Isurface)%z(ix,Ioy,iz)),ix=1,nx)

        enddo
      close(01)

      call create(eAll%solns(1)%pol(1)%grid,Etemp,eAll%solns(1)%pol(1)%gridType)
      ! to plot a plane view of the fields
      ! XPer plane view, across the center (Y)
      open(02, file = 'XZYZ.dat')
        write(02,*)'Nmode, period, X, Y, Z, Ex, Ey, Ez'
        do j = 1,eAll%nTx   ! loop over periods
            Etemp = eAll%solns(j)%pol(1)
            write(02,'(i2,es13.5,2x,i3,2x,i3,2x,i3,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5)') ((1,txDict(eAll%solns(j)%tx)%period, Iox, iy,iz, &
                                 real(Etemp%x(Iox,iy,iz)), imag(Etemp%x(Iox,iy,iz)), abs(Etemp%x(Iox,iy,iz)),&
                                 real(Etemp%y(Iox,iy,iz)), imag(Etemp%y(Iox,iy,iz)), abs(Etemp%y(Iox,iy,iz)),&
                                 real(Etemp%z(Iox,iy,iz)), imag(Etemp%z(Iox,iy,iz)), abs(Etemp%z(Iox,iy,iz)),iy=1,Etemp%ny),iz=1,Etemp%nz)
            call zero(Etemp)
            Etemp = eAll%solns(j)%pol(2)
            write(02,'(i2,es13.5,2x,i3,2x,i3,2x,i3,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5)') ((2,txDict(eAll%solns(j)%tx)%period, ix, Ioy,iz, &
                                 real(Etemp%x(ix,Ioy,iz)),imag(Etemp%x(ix,Ioy,iz)), abs(Etemp%x(ix,Ioy,iz)), &
                                 real(Etemp%y(ix,Ioy,iz)),imag(Etemp%y(ix,Ioy,iz)), abs(Etemp%y(ix,Ioy,iz)),&
                                 real(Etemp%z(ix,Ioy,iz)),imag(Etemp%z(ix,Ioy,iz)), abs(Etemp%z(ix,Ioy,iz)),ix=1,Etemp%nx),iz=1,Etemp%nz)
            call zero(Etemp)
        enddo ! periods
        close(02)
      ! to plot a plane view of the fields
      ! XPer plane view, across the center (Y)
      open(03, file = 'XY.dat')
        ig = 2
        iz = 12
        write(03,*)'Nmode, period, X, Y, Z,Ex, Ey, Ez'
        do j = 1,eAll%nTx   ! loop over periods
            write(03,'(i2,es13.5,2x,i3,2x,i3,2x,i3,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5)') ((1,txDict(eAll%solns(j)%tx)%period, ix, iy,iz, &
                                 real(eAll%solns(j)%pol(1)%cvarray(ig)%x(ix,iy,iz)),imag(eAll%solns(j)%pol(1)%cvarray(ig)%x(ix,iy,iz)), abs(eAll%solns(j)%pol(1)%cvarray(ig)%x(ix,iy,iz)),&
                                 real(eAll%solns(j)%pol(1)%cvarray(ig)%y(ix,iy,iz)),imag(eAll%solns(j)%pol(1)%cvarray(ig)%y(ix,iy,iz)), abs(eAll%solns(j)%pol(1)%cvarray(ig)%y(ix,iy,iz)),&
                                 real(eAll%solns(j)%pol(1)%cvarray(ig)%z(ix,iy,iz)),imag(eAll%solns(j)%pol(1)%cvarray(ig)%z(ix,iy,iz)), abs(eAll%solns(j)%pol(1)%cvarray(ig)%z(ix,iy,iz)),&
                                 ix=1,eAll%solns(j)%pol(1)%cvarray(ig)%nx),iy=1,eAll%solns(j)%pol(1)%cvarray(ig)%ny)
            write(03,'(i2,es13.5,2x,i3,2x,i3,2x,i3,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5,2x,es13.5)') ((2,txDict(eAll%solns(j)%tx)%period, ix, iy,iz, &
                                 real(eAll%solns(j)%pol(2)%cvarray(ig)%x(ix,iy,iz)),imag(eAll%solns(j)%pol(2)%cvarray(ig)%x(ix,iy,iz)), abs(eAll%solns(j)%pol(2)%cvarray(ig)%x(ix,iy,iz)),&
                                 real(eAll%solns(j)%pol(2)%cvarray(ig)%y(ix,iy,iz)),imag(eAll%solns(j)%pol(2)%cvarray(ig)%y(ix,iy,iz)), abs(eAll%solns(j)%pol(2)%cvarray(ig)%y(ix,iy,iz)),&
                                 real(eAll%solns(j)%pol(2)%cvarray(ig)%z(ix,iy,iz)),imag(eAll%solns(j)%pol(2)%cvarray(ig)%z(ix,iy,iz)), abs(eAll%solns(j)%pol(2)%cvarray(ig)%z(ix,iy,iz)),&
                                 ix=1,eAll%solns(j)%pol(2)%cvarray(ig)%nx),iy=1,eAll%solns(j)%pol(2)%cvarray(ig)%ny)
        enddo ! periods
        close(03)

        call deall(Etemp)
      end subroutine plotFields


end module ioAscii
