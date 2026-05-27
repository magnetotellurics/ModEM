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
!
!   Modified by Kelbert, Apr. 2018 to add IO routines for RHS E-fields.
!   To achieve joint inversion generality in the future, replaced ModeNames
!   which were hardcoded to Ex, Ey by Pol_name(:), now stored in EMsoln & EMrhs.
!   And of course we can and should now use read_cvector/write_cvector instead
!   (this is not fixed yet). Ideally, convert everything to NetCDF or similar.
!
!   This is the HDF5 I/O version coded by Spencer Wilbur (USGS) under the direction
!   of Anna Kelbert, Jul. 2023. File intentionally misnamed; it is renamed at the
!   level of the configuration file.

! *****************************************************************************
module ioAscii
  ! This module contains io routines for boundary nodes,
  ! electric field solutions, impedances, and solver diagnostics
  use hdf5
  use griddef
  use math_constants
  use emsolve3d
  use dataspace
  use datafunc ! for data types
  use ForwardSolver, only: solnVectorMTX_t
  use transmitters
  use receivers
  use datatypes
  use modelParam_IO_HDF5

  implicit none

  public                   :: FileReadInit, FileWriteInit, FileWriteInit_hdf5
  public                   :: ZfileReadInit, ZfileWriteInit
  public                   :: BCfileRead, BCfileWrite
  public                   :: EfileRead, EfileWrite, EfileWrite_hdf5
  public                   :: DfileWrite
  public                   :: ZfileRead, ZfileWrite
  public                   ::  write_solnVectorMTX_hdf5
  public                   :: read_solnVectorMTX


  integer(HID_T), private, save   :: file_id, group_id, attr_id, dset_id, dspace_id, atype_id, aspace_id, dtype_id ! file, data set, and dataspace handles
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
    integer (kind=4)                            :: filePer4,fileMode4,nx,ny,nz,nzAir

    open (unit=ioNum,file=inFile,status='old', form ='unformatted',iostat=ios)

    IF( ios/=0) THEN
       WRITE(0,*) 'Error opening file: ', inFile
    else
       ! read the header
       read(ioNum) fileVersion, filePer4, fileMode4, &
	nx, ny, nz, nzAir, &
	fileGrid%ox, fileGrid%oy, fileGrid%oz, fileGrid%rotdeg

	   ! convert integer*4 to default integer precision
	   filePer = int(filePer4)
	   fileMode = int(fileMode4)
	   fileGrid%nx = int(nx)
	   fileGrid%ny = int(ny)
	   fileGrid%nz = int(nz)
	   fileGrid%nzAir = int(nzAir)

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
       write(ioNum) version,int(nPer,4),int(nMode,4),int(inGrid%nx,4),int(inGrid%ny,4),int(inGrid%nz,4),int(inGrid%nzAir,4), &
       inGrid%ox, inGrid%oy, inGrid%oz, inGrid%rotdeg
       write(ioNum) inGrid%dx
       write(ioNum) inGrid%dy
       write(ioNum) inGrid%dz
       write(0,*) inGrid%dx
    endif
  end subroutine FileWriteInit

  ! ***************************************************************************
  ! write header for MT_fwd output files (except impedance file)
  subroutine FileWriteInit_hdf5(version,iFreq, outFile, inGrid, nPer, nMode)

    implicit none
    character (len=20), intent(in)		:: version
    ! the version of the software in use
    character (len=80), intent(in)		:: outFile
    ! the filename
    type(grid_t), intent(in)			:: inGrid
    type(custom_att), target             :: attr_obj
    integer, intent(in)				    :: nPer
    integer, intent(in)				    :: iFreq
    integer, intent(in)				    :: nMode
    integer                      :: hdferr, i 
    INTEGER(8), DIMENSION(1)      :: dims
    INTEGER, DIMENSION(1), allocatable        :: grid(:)
    character(len = 5)                        :: gname
    
    ! character(len=13)                         :: grid_dim = "grid_geometry"
    ! character(len=16)                         :: grid_grp

      ! write(0,*) "The fileInit hdf5"
      ! grid_grp= trim(txDict(iFreq)%tx_type)//"/"//grid_dim

      CALL h5gcreate_f(file_id, 'Grid', group_id, hdferr)

      ! Assign attribute values 
      attr_obj%att_real = -inGrid%ox
      call write_hdf5_attr('dbl','model_origin_x', attr_obj, group_id)
      attr_obj%att_real = -inGrid%oy
      call write_hdf5_attr('dbl','model_origin_y', attr_obj, group_id)
      attr_obj%att_real = -inGrid%oz
      call write_hdf5_attr('dbl','model_origin_z', attr_obj, group_id)
      
      
      do i = 1, 3
        if (i == 1) then
      
          dims = size(inGrid%dx)
          CALL h5screate_simple_f(1, dims, dspace_id, hdferr)  !define the dataspace (not sure what the dimensions are yet) 
          Call h5dcreate_f(group_id, 'dx', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, inGrid%dx, dims, hdferr)
        else if (i == 2) then
          dims = size(inGrid%dy)
          CALL h5screate_simple_f(1, dims, dspace_id, hdferr)  !define the dataspace (not sure what the dimensions are yet) 
          Call h5dcreate_f(group_id, 'dy', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, inGrid%dy, dims, hdferr)
        else 
          dims = size(inGrid%dz)
          CALL h5screate_simple_f(1, dims, dspace_id, hdferr)  !define the dataspace (not sure what the dimensions are yet) 
          Call h5dcreate_f(group_id, 'dz', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, inGrid%dz, dims, hdferr)
        end if 
        
      end do 
   
     
  end subroutine FileWriteInit_hdf5

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
  subroutine EfileRead(ioNum, nRecSkip, fileOmega, inE)

    implicit none
    integer, intent(in)             :: ioNum,nRecSkip !ifreq,imode
    type (cvector), intent(inout)       :: inE
    real (kind=prec), intent(out)   :: fileOmega

    !  local variables
    !integer                 :: nRecSkip, iskip
    integer                    :: iskip

    !  following could be optional oututs
    integer (kind=4)        :: fileIfreq, fileMode
    character (len = 20)                :: ModeName

    !  hard code number of modes for now
    !integer     :: nMode = 2

    ! calculate number of records before the header for this frequency/mode
    !nRecSkip = ((ifreq-1)*nMode+(imode-1))*4+4

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
    write(ioNum) Omega, int(iFreq,4), int(iMode,4),ModeName

    ! write the electrical field data - 3 records
    write(ioNum) outE%x
    write(ioNum) outE%y
    write(ioNum) outE%z
    
  
  end subroutine EfileWrite

  ! ***************************************************************************
    ! write electrical field solution for one frequency/mode to hdf5 
  subroutine EfileWrite_hdf5(Omega, iFreq, iMode,ModeName, outE)

    implicit none
    real (kind=prec), intent(in)	  :: Omega
    integer, intent(in)				      :: iFreq,iMode
    integer                         :: hdferr, i
    character (len=20), intent(in)	:: ModeName
    INTEGER(8), DIMENSION(2)        :: dims
    type (cvector), intent(in)			:: outE
    REAL(8), DIMENSION(2), allocatable          :: efield(:,:,:)
    character(len= 4)               :: freq
    character(len = 6)              :: pol
    character (len = 1)             :: component 
    logical                         :: exists
    freq = trim("/"//txDict(iFreq)%tx_type)
    pol = trim(txDict(iFreq)%tx_type)//"/"//ModeName !this changes on every iteration of K in writeVectorSoln_hdf5
    

  !NEED to make if else to see if groups exist after first dataset is written 
  ! Read tipper definition, if exists in file
    call h5lexists_f(file_id, freq, exists, hdferr)
    if (exists) then
        CALL h5gopen_f(file_id, freq, group_id, hdferr)
    else
        call h5gcreate_f(file_id, freq, group_id, hdferr) ! create group for transmitter/frequency
    end if 

    CALL h5gcreate_f(file_id, pol, group_id, hdferr) ! creates the polerization groups for each dataset
  dims(1) = 3
  do i = 1, 3
    if (i == 1) then
      efield = outE%x
      dims(2) =size(efield)
      component = 'x'
    else if (i == 2) then
      efield = outE%y
      dims(2) =size(efield)
      component = 'y'
    else 
      efield = outE%z
      dims(2) =size(efield)
      component = 'z'
    end if 
    
    CALL h5screate_simple_f(2, dims, dspace_id, hdferr)  !define the dataspace (not sure what the dimensions are yet)
    Call h5dcreate_f(group_id, component, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr )
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, efield, dims, hdferr)
  end do
    write(0,*) "The Efilewrite hdf5"

  end subroutine EfileWrite_hdf5
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
!   Addtional routines, to allow greater compatability with Modular2D
!    code structure
!*******************************************************************************
!*******************************************************************************
      subroutine open_ef_hdf5(cfile)
        character(*), intent(in)                :: cfile
        integer                                 :: hdferr
        CHARACTER(LEN=7), PARAMETER  :: e_field = "E_Field"
        ! CHARACTER(LEN=7), PARAMETER  :: data_mt_group = "Data/MT"
        logical                      :: lexist
      
        inquire(file = cfile, exist = lexist)
        if (lexist) then 
            CALL h5open_f(hdferr)
            CALL h5fopen_f(cfile, H5F_ACC_RDWR_F, file_id, hdferr)
        else 
                CALL h5open_f(hdferr) ! throws error if cannot open 
                CALL h5fcreate_f(cfile, H5F_ACC_TRUNC_F, file_id, hdferr)!create the file using the variable cfile given to the terminal when running the script 
            
                ! !Create the groups for the dataset
                ! CALL h5gcreate_f(file_id, e_field, group_id, hdferr) 
                write(0,*) "The efield hdf5 file was created and opened"
            
        end if
        
      end subroutine open_ef_hdf5
  !*******************************************************************************
      subroutine close_ef_hdf5(cfile)
        character(*), intent(in)                :: cfile
        integer                                 :: hdferr
      
        CALL h5fclose_f(file_id, hdferr)
        CALL h5close_f(hdferr)
      
      end subroutine close_ef_hdf5

!******************************************************************
      subroutine read_solnVectorMTX(Larg_Grid,eAll,cfile)

      ! reads an array of solution vectors for all transmitters & subgrids
      ! currently uses the old binary format; will switch to NetCDF when
      ! time allows
      ! this SHOULD initialize the transmitter dictionary and the grids, as needed
      ! but currently it can only work if these are pre-allocated and the same
	  
	  !NM: "this SHOULD initialize the ...", That will not work for the nested modelling and in case of Joint inversion.
	  ! Also, the large_grid is initialized based on what we have in the E-file. Thus, the first thing to do is to read the meta data from the file,
	  ! initialize both the large_grid and corrsponding solution vectors based on the number of Tx present in the file and NOT in the txDic!

      character(*), intent(in)                    :: cfile
      type(solnVectorMTX_t), intent(inout)        :: eAll
      type(grid_t), intent(inout)            :: Larg_Grid

      !   local variables
      integer           :: j,k,nMode = 2, ios,ig,cdot,filePer
      integer           :: iTx,nTx,nRecSkip
      character (len=3)         :: igchar
      character (len=20)        :: version = ''
      character (len=200)       :: fn_input
      real (kind=prec)   :: omega

      fn_input = trim(cfile)

      write(*,*) 'Reading E-fields from file: ',trim(fn_input)
      call FileReadInit(fn_input,ioE,Larg_Grid,filePer,nMode,version,ios)
      call setup_grid(Larg_Grid)
      nTx = filePer   ! Get nTx from the file

       ! Test that the grid is read properly, at least on a basic level
      if (.not. Larg_Grid%allocated) then
          write(0,*) 'WARNING: grid is not allocated in read_solnVectorMTX'
      else if ((Larg_Grid%Nx <= 0) .or. (Larg_Grid%Ny <= 0) .or. (Larg_Grid%Nz <= 0)) then
          write(0,*) 'Grid information: Nx=',Larg_Grid%Nx,' Ny=',Larg_Grid%Ny,' Nz=',Larg_Grid%Nz
          write(0,*) 'WARNING: grid is not set up properly in read_solnVectorMTX'
      end if

      call create_solnVectorMTX(nTx,eAll)
      do iTx=1,nTx
          call create_solnVector(Larg_Grid,iTx,eAll%solns(iTx))
      end do		  
		  
		  nRecSkip = 0
		  
          do j = 1,nTx
       
             do k = 1,eAll%solns(j)%nPol

               nRecSkip = nRecSkip+4

               call EfileRead(ioE, nRecSkip, omega, eAll%solns(j)%pol(k))

               if (abs((omega - txDict(eAll%solns(j)%tx)%omega)/omega) > TOL6) then
                    write(0,*) 'Warning: frequencies don''t match on E-field input ',j
               endif

             enddo
          enddo
          close(ioE)

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
      character (len=20) 		:: version = ''
      character (len=200)       :: fn_output
      real (kind=prec)   :: omega

          fn_output = trim(cfile)

           ! Test that the output grid is well-defined, at least on a basic level
          if (.not. eAll%solns(1)%grid%allocated) then
              write(0,*) 'WARNING: E-field grid is not allocated in write_solnVectorMTX'
          else if ((eAll%solns(1)%grid%Nx <= 0) .or. (eAll%solns(1)%grid%Ny <= 0) .or. (eAll%solns(1)%grid%Nz <= 0)) then
              write(0,*) 'Grid information: Nx=',eAll%solns(1)%grid%Nx,' Ny=',eAll%solns(1)%grid%Ny,' Nz=',eAll%solns(1)%grid%Nz
              write(0,*) 'WARNING: E-field grid is not output properly in write_solnVectorMTX'
          end if

          write(*,*) 'E-fields written to ',trim(fn_output)

          call FileWriteInit(version,fn_output,ioE,eAll%solns(1)%grid,eAll%nTX,nMode,ios)
          do j = 1,eAll%nTx
             do k = 1,eAll%solns(j)%nPol
               omega = txDict(eAll%solns(j)%tx)%omega
              
               call EfileWrite(ioE, omega, j,  k, eAll%solns(j)%Pol_name(k), eAll%solns(j)%pol(k))

             enddo
          enddo
          close(ioE)

      end subroutine write_solnVectorMTX
!******************************************************************
      subroutine write_solnVectorMTX_hdf5(eAll,cfile)

        !  open cfile on unit fid, writes out object of
        !   type cvector in standard format (readable by matlab
        !   routine readcvector.m), closes file
        !  NOT coded at present to specifically write out TE/TM
        !    solutions, periods, etc. (can get this infor from
        !    eAll%solns(j)%tx, but only with access to TXdict.
  
        character(*), intent(in)                :: cfile
        type(solnVectorMTX_t), intent(in)       :: eAll
  
        !   local variables
        integer                         :: j,k,nMode = 2, ios,ig,cdot, hdferr
        character (len=3)               :: igchar
        character (len=20) 		          :: version = ''
        character (len=200)             :: fn_output
        real (kind=prec)                :: omega
        character(len= 4)               :: freq
        character(len = 4)              :: pol
        character(len=2)                :: padded_j, padded_k
        logical                         :: exists
        
          call open_ef_hdf5(cfile)
          ! CALL h5gopen_f(file_id, 'E_field', group_id, hdferr) ! open the main group
            
          fn_output = trim(cfile)

          ! Test that the output grid is well-defined, at least on a basic level
         if (.not. eAll%solns(1)%grid%allocated) then
             write(0,*) 'WARNING: E-field grid is not allocated in write_solnVectorMTX'
         else if ((eAll%solns(1)%grid%Nx <= 0) .or. (eAll%solns(1)%grid%Ny <= 0) .or. (eAll%solns(1)%grid%Nz <= 0)) then
             write(0,*) 'Grid information: Nx=',eAll%solns(1)%grid%Nx,' Ny=',eAll%solns(1)%grid%Ny,' Nz=',eAll%solns(1)%grid%Nz
             write(0,*) 'WARNING: E-field grid is not output properly in write_solnVectorMTX'
         end if

         write(*,*) 'E-fields written to HDF5 ',trim(fn_output)
            write(*,*) 'E-fields written to HDF5 ',trim(fn_output)
  
            
            call FileWriteInit_hdf5(version,j,fn_output,eAll%solns(1)%grid,eAll%nTX,nMode) !displays grid dimensions in its own group under transmitter group
            
            do j = 1, eAll%nTx
               do k = 1, eAll%solns(j)%nPol
                write(0,*) j, k
                 omega = txDict(eAll%solns(j)%tx)%omega
                if (j== 1) then
                 call EfileWrite_hdf5(omega, j, k,eAll%solns(j)%Pol_name(k), eAll%solns(j)%pol(k))
                else 
                  write(0,*) 'Not MT transmitter so no ModeNames available'
                end if 
               enddo
            enddo
            
            call close_ef_hdf5(cfile)
        end subroutine write_solnVectorMTX_hdf5
  

      !******************************************************************
      subroutine read_rhsVectorMTX(grid,bAll,cfile,format)

          ! reads an array of solution vectors for all transmitters & subgrids
          ! currently uses the old binary format; will switch to NetCDF when
          ! time allows
          ! will allocate the grid from a full format file. for sparse files,
          ! this will only work if the grid is preallocated and compatible

          character(*), intent(in)                    :: cfile
          type(rhsVectorMTX_t), intent(inout)         :: bAll
          type(grid_t), intent(inout)                 :: grid
          character(*), optional, intent(in)          :: format

          !   local variables
          integer           :: j,k,jj,kk,nMode = 2, ios,istat,ig,cdot,filePer, nPol
          integer           :: iTx,nTx,nRecSkip
          character (len=3)         :: igchar
          character (len=20)        :: version = '',source_type,tx_type,mode
          character (len=30)        :: str
          character (len=200)       :: fn_input
          logical                   :: sparse
          real (kind=prec)   :: omega, period

          sparse = .true.

          if (present(format)) then
            if (trim(format) .eq. 'full') then
                sparse = .false.
            end if
          end if

          if (.not. sparse) then

              fn_input = trim(cfile)

              write(*,*) 'Reading RHS E-fields from file: ',trim(fn_input)
              call FileReadInit(fn_input,ioE,grid,filePer,nMode,version,ios)
              call setup_grid(grid)
              nTx = filePer   ! Get nTx from the file

              call create_rhsVectorMTX(nTx,bAll)
              do iTx=1,nTx
                  bAll%combs(iTx)%nonzero_BC = .false.
                  bAll%combs(iTx)%nonzero_Source = .true.
                  bAll%combs(iTx)%sparse_Source = .false.
                  call create_rhsVector(grid,iTx,bAll%combs(iTx))
              end do


              nRecSkip = 0

              do j = 1,nTx

                  do k = 1,bAll%combs(j)%nPol

                      nRecSkip = nRecSkip+4

                      call EfileRead(ioE, nRecSkip, omega, bAll%combs(j)%b(k)%s)

                      if (abs(omega - txDict(bAll%combs(j)%tx)%omega) > R_TINY) then
                          write(0,*) 'Warning: frequencies don''t match on E-field input ',j
                      endif

                  enddo
              enddo
              close(ioE)

          else

              open (unit=ioE,file=cfile,STATUS='unknown', form ='formatted', iostat=ios)

              if( ios/=0) then
                  write(0,*) 'Error opening sparse vector output file in read_rhsVectorMTX: ', cfile
              else
                  read(ioE,'(a30,i5)',iostat=istat) str,nTx
                  call create_rhsVectorMTX(nTx,bAll)
                  do iTx=1,nTx
                      bAll%combs(iTx)%nonzero_BC = .true.
                      bAll%combs(iTx)%nonzero_Source = .true.
                      bAll%combs(iTx)%sparse_Source = .true.
                      call create_rhsVector(grid,iTx,bAll%combs(iTx))
                  end do


                  do j = 1,bAll%nTx

                      ! now that everything is allocated, set them both to false
                      ! because one of them could be missing from file...
                      bAll%combs(j)%nonzero_BC = .false.
                      bAll%combs(j)%nonzero_Source = .false.

                      do k = 1,bAll%combs(j)%nPol
                          ! now that everything is allocated, set them both to false
                          ! because one of them could be missing from file...
                          bAll%combs(j)%b(k)%nonzero_BC = .false.
                          bAll%combs(j)%b(k)%nonzero_Source = .false.

                          read(ioE,'(a30,a4,a10,i5)',iostat=istat) str, str, tx_type, nPol
                          if ((trim(txDict(bAll%combs(j)%tx)%tx_type) .ne. trim(tx_type)) &
                              .or. (bAll%combs(j)%nPol .ne. nPol)) then
                              write(0,*) 'Warning: transmitter types don''t match on RHS E-field input ',j
                          end if
                          read(ioE,'(a30,i5,es14.6)',iostat=istat) str, jj, period
                          write(*,*) str,jj,period
                          read(ioE,'(a30,i5,a5,a20)',iostat=istat) str, kk, str, mode
                          bAll%combs(j)%Pol_name(k) = mode
                          if (abs((period - txDict(bAll%combs(j)%tx)%period)/period) > TOL6) then
                              write(0,*) 'Warning: periods don''t match on RHS E-field input ',j
                              write(0,*) period, txDict(bAll%combs(j)%tx)%period
                          end if

                          read(ioE,*,iostat=istat) source_type
                          write(*,*) 'Reading... ',source_type

                          if (trim(adjustl(source_type)) .eq. 'BC') then
                              call read_cboundary(ioE,bAll%combs(j)%b(k)%bc,grid)
                              bAll%combs(j)%b(k)%nonzero_BC = .true.
                              bAll%combs(j)%nonzero_BC = .true.
                              write(*,'(a34,i5,a6,a20)') 'Completed reading BC for period #',jj,' mode ',mode
                          end if

                          if (trim(adjustl(source_type)) .eq. 'SOURCE') then
                              call read_sparsevecc(ioE,bAll%combs(j)%b(k)%sSparse,grid)
                              bAll%combs(j)%b(k)%nonzero_Source = .true.
                              bAll%combs(j)%nonzero_Source = .true.
                              write(*,'(a38,i5,a6,a20)') 'Completed reading SOURCE for period #',jj,' mode ',mode
                          end if

                      enddo
                  enddo
                  close(ioE)
              end if

          end if

      end subroutine read_rhsVectorMTX

      !******************************************************************
      subroutine write_rhsVectorMTX(bAll,cfile,format)

          !  open cfile on unit fid, writes out object of
          !   type cvector in standard format (readable by matlab
          !   routine readcvector.m), closes file
          !  NOT coded at present to specifically write out TE/TM
          !    solutions, periods, etc. (can get this infor from
          !    eAll%solns(j)%tx, but only with access to TXdict.

          character(*), intent(in)                    :: cfile
          type(rhsVectorMTX_t), intent(in)            :: bAll
          character(*), intent(in), optional          :: format

          !   local variables
          type (cvector)    :: rhsFull
          type (sparsevecc) :: rhsSparse
          integer           :: j,k,nMode = 2, ios,istat,ig,cdot
          character (len=3)         :: igchar
          character (len=20)        :: version = ''
          character (len=200)       :: fn_output
          logical            :: sparse
          real (kind=prec)   :: omega,period

          sparse = .true.

          if (present(format)) then
            if (trim(format) .eq. 'full') then
                sparse = .false.
            end if
          end if

          fn_output = trim(cfile)

          if (.not. sparse) then

              call FileWriteInit(version,fn_output,ioE,bAll%combs(1)%grid,bAll%nTX,nMode,ios)
              do j = 1,bAll%nTx
                  do k = 1,bAll%combs(j)%nPol
                      omega = txDict(bAll%combs(j)%tx)%omega

                      if (bAll%combs(j)%b(k)%nonzero_Source) then
                        call copy_cvector(rhsFull, bAll%combs(j)%b(k)%s)
                      elseif (bAll%combs(j)%b(k)%sparse_Source) then
                        call add_scvector(C_ONE, bAll%combs(j)%b(k)%sSparse, rhsFull)
                      elseif (bAll%combs(j)%b(k)%nonzero_BC) then
                        call copy_bcvector(bAll%combs(j)%b(k)%bc, rhsFull)
                      end if

                      call EfileWrite(ioE, omega, j,  k, bAll%combs(j)%Pol_name(k), rhsFull)

                  enddo
              enddo
              close(ioE)

          else

              open (unit=ioE,file=cfile,STATUS='unknown', form ='formatted', iostat=ios)

              if( ios/=0) then
                  write(0,*) 'Error opening sparse vector output file in write_solnVectorMTX: ', cfile
              else
                  write(ioE,'(a30,i5)',iostat=istat) 'Number of transmitters      : ',bAll%nTx
                  do j = 1,bAll%nTx
                      do k = 1,bAll%combs(j)%nPol
                          omega = txDict(bAll%combs(j)%tx)%omega
                          period = txDict(bAll%combs(j)%tx)%period
                          write(ioE,'(a30,a4,a10,i5)',iostat=istat) 'Transmitter type & nPol     : ','',txDict(bAll%combs(j)%tx)%tx_type, bAll%combs(j)%nPol
                          write(ioE,'(a30,i5,es14.6)',iostat=istat) 'Period number & value (secs): ',j,period
                          write(ioE,'(a30,i5,a5,a20)',iostat=istat) 'Polarization number & name  : ',k, '', bAll%combs(j)%Pol_name(k)

                          if (bAll%combs(j)%b(k)%nonzero_BC) then
                              write(ioE,'(a10)',iostat=istat) 'BC'
                              call write_cboundary(ioE,bAll%combs(j)%b(k)%bc)
                          end if

                          if (bAll%combs(j)%b(k)%nonzero_Source .and. bAll%combs(j)%b(k)%sparse_Source) then
                              write(ioE,'(a10)',iostat=istat) 'SOURCE'
                              call write_sparsevecc(ioE,bAll%combs(j)%b(k)%sSparse)
                          end if

                      enddo
                  enddo
                  close(ioE)
              end if

          end if

          write(*,*) 'RHS E-fields written to ',trim(fn_output)

      end subroutine write_rhsVectorMTX

end module ioAscii
