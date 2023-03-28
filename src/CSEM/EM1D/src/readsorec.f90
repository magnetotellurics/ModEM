!**********************************************************************
!  FD EM subroutine readsorec
!
!  Purpose:  read source or receiver specifications from ascii input file
!
!  Rita Streich 2009-2011
!**********************************************************************

subroutine readsorec(filename,src,comm)
 use, intrinsic :: iso_fortran_env
  implicit none

  !external variables
  character(len=namlen),intent(in)      :: filename
  type(sorec),dimension(:),pointer      :: src        !source or receiver structure
  integer(kind=int32),intent(in)        :: comm   !MPI communicator (not required for serial I/O)

  !internal variables
  integer(kind=int32)    :: lu           !file unit number
  integer(kind=int32)    :: ierr         !error index
  integer(kind=int32)    :: nsorec       !number of sources/receivers
  integer(kind=int32)    :: ishot        !counter
  character(len=6)       :: ishotstr     !string for error output
  integer(kind=int32)    :: isrc         !source/receiver index (not really used)


  !let only processs 0 access the file!
  !let all processes know all sources
  if (pid .eq. 0) then
    !open file
    lu = AvailableUnit()
    open(unit=lu,file=trim(adjustl(filename)),status='old',iostat=ierr)
    if (ierr.ne.0) call open_error(pid,'readsorec',filename,ierr)

    !read number of sources or receiver groups
    read(lu,*,iostat=ierr) nsorec
    if (ierr.ne.0) call readwrite_error(pid,'readsorec',filename,'r',ierr)
  endif

#ifdef USE_MPI
  call MPI_Bcast(nsorec,1,MPI_INTEGER,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'readsorec','MPI_Bcast nsorec',ierr)
#endif

  !allocate source / receiver structure
  allocate(src(nsorec), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'readsorec','source/receiver structure ',ierr)


  do ishot=1,nsorec

    if (pid .eq. 0) then
      !read line containing "shot" index, source type
      !RS 08.06.2011 new format: source name is in this line in sorpos files
      !receivers need individual names so there is no name needed here,
      !  but still for consistency with sorpos files we put a dummy name here,
      !  this could be used later to identify which receiver group to use for which source
      read(lu,*,iostat=ierr) isrc, src(ishot)%type, src(ishot)%srcname
      if (ierr.ne.0) call io_error(pid,filename,'reading source type and name from',ierr)
    endif

#ifdef USE_MPI
    call MPI_Bcast(src(ishot)%type,1,MPI_INTEGER,0,comm,ierr)
    if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'readsorec','MPI_Bcast source type',ierr)
    call MPI_Bcast(src(ishot)%srcname,namlen,MPI_CHARACTER,0,comm,ierr)
    if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'readsorec','MPI_Bcast source name',ierr)
#endif

    !read the source (or receiver) group
    select case (src(ishot)%type)

    !dipole sources
    case (dipole)
      call read_dipolesource(src(ishot),lu,comm)

    !(general) wire sources
    case (wire)
      call read_wiresource(src(ishot),lu,comm)

    !"multipod" wire sources with several (most typically 3) wires laid out from one center point,
    !  currents at center point adding up to zero --> fixed phase shifts between currents on each wire
    case (star)
      call read_starsource(src(ishot),lu,comm)


    !receivers
    case (receiver_type)
      call read_receivers(src(ishot),lu,comm)

    !invalid source type: write error message & exit
    case default

      write(ishotstr,'(i6)') ishot
      call invalid_error(pid,'readsorec',filename,'source/receiver type for "shot" '//trim(adjustl(ishotstr))//' ', &
        intnum=src(ishot)%type)

    end select

  enddo !sources or receiver groups


  if (pid .eq. 0) then
    !close the file
    close(lu, iostat=ierr)
    if (ierr.ne.0) call close_error(pid,'readgeometry',filename,ierr)
  endif

endsubroutine readsorec



!**********************************************************************
!  FD EM subroutine read_dipolesource
!
!  Purpose:  read specifications of a dipole source
!    the source can have many dipole elements
!
!  Rita Streich 2009
!**********************************************************************
subroutine read_dipolesource(src,lu,comm)

  implicit none

  !external variable
  type(sorec)    :: src  !the source
  integer(kind=int32),intent(in)   :: lu !file unit number
  integer(kind=int32),intent(in)        :: comm   !MPI communicator (not required for serial I/O)

  !internal variables
  integer(kind=int32)              :: ierr   !error index
  character(len=namlen)            :: filename
  integer(kind=int32)              :: nelem  !number of dipole elements
  integer(kind=int32)              :: ielem  !dipole element counter


  !1-element vector for number of dipole elements (for wire source this would be nr of elements for each wire)
  allocate(src%nelem(1), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'read_dipolesource','nelem',ierr)


  !get number of dipole elements
  if (pid .eq. 0) then
    !get file name
    inquire(lu, NAME=filename)

    read(lu,*,iostat=ierr) nelem
    if (ierr.ne.0) call io_error(pid,filename,'reading number of dipole elements from ',ierr)

    !check number
    if (nelem.le.0) call invalid_error(pid,'read_dipolesource',filename,'number of dipole elements ',intnum=nelem)
  endif

  !communicate nr of dipole elements
#ifdef USE_MPI
  call MPI_Bcast(nelem,1,MPI_INTEGER,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_dipolesource','MPI_Bcast nr of dipole elements',ierr)
#endif

  src%nelem(1) = nelem


  !array for dipole positions
  allocate(src%pos(3,nelem), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'read_dipolesource','array for dipole positions',ierr)

  !info on source components
  allocate(src%ljx(nelem),src%ljy(nelem),src%ljz(nelem),src%akx(nelem),src%aky(nelem),src%akz(nelem), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'read_dipolesource','dipole source component vectors',ierr)


  if (pid .eq. 0) then
    !read dipole element specifications
    do ielem = 1, nelem
      read(lu,*,iostat=ierr) src%pos(:,ielem), &
                             src%ljx(ielem),src%ljy(ielem),src%ljz(ielem), &
                             src%akx(ielem),src%aky(ielem),src%akz(ielem)
      if (ierr.ne.0) call io_error(pid,filename,'reading dipole components from ',ierr)
    enddo
  endif

#ifdef USE_MPI
  call MPI_Bcast(src%pos,3*nelem,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_dipolesource','MPI_Bcast dipole positions',ierr)
  call MPI_Bcast(src%ljx,nelem,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_dipolesource','MPI_Bcast ljx',ierr)
  call MPI_Bcast(src%ljy,nelem,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_dipolesource','MPI_Bcast ljy',ierr)
  call MPI_Bcast(src%ljz,nelem,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_dipolesource','MPI_Bcast ljz',ierr)
  call MPI_Bcast(src%akx,nelem,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_dipolesource','MPI_Bcast akx',ierr)
  call MPI_Bcast(src%aky,nelem,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_dipolesource','MPI_Bcast aky',ierr)
  call MPI_Bcast(src%akz,nelem,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_dipolesource','MPI_Bcast akz',ierr)
#endif

  !this is used when looping over sources and computing length of output files
  src%ncur = 1

  !check which source components are present
  call checksourcecomp(src)

  !offset of x coordinate from line x=0, for 2.5D
  src%shiftx = sum(src%pos(1,:)) / nelem

endsubroutine read_dipolesource


!**********************************************************************
!  FD EM subroutine read_wiresource
!
!  Purpose:  read specifications of a wire source
!    the source can consist of several wire segments
!
!  Rita Streich 2009-2011
!**********************************************************************
subroutine read_wiresource(src,lu,comm)

  implicit none

  !external variables
  type(sorec)           :: src    !the source
  integer(kind=int32),intent(in)   :: lu !file unit number
  integer(kind=int32),intent(in)        :: comm   !MPI communicator (not required for serial I/O)

  !internal variables
  integer(kind=int32)   :: ierr   !error index
  character(len=namlen) :: filename
  integer(kind=int32)   :: nwire  !number of wires
  integer(kind=int32)   :: iwire  !counter for wires
  real(kind=real64)     :: x0,y0,z0,x1,y1,z1 !wire end points
  real(kind=real64)     :: dlwx,dlwy    !x,y lengths of wire elements
  real(kind=real64)     :: lw,lx,ly     !wire length, total and in x,y directions
  real(kind=real64),dimension(:,:),allocatable  :: postmp  !temp positions and input lengths of all wire elements
  integer(kind=int32)   :: ielem  !wire element counter
  integer(kind=int32)   :: nelem  !tmp number of wire elements


  !get number of long wire segments composing this source
  if (pid .eq. 0) then
    !get file name
    inquire(lu, NAME=filename)

    read(lu,*,iostat=ierr) nwire
    if (ierr.ne.0) call io_error(pid,filename,'reading number of wires from ',ierr)

    !check number
    if (nwire.le.0) call invalid_error(pid,'read_wiresource',filename,'number of wires within source ',intnum=nwire)
  endif

  !communicate nr of long wire segments
#ifdef USE_MPI
  call MPI_Bcast(nwire,1,MPI_INTEGER,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_wiresource','MPI_Bcast nr of wires',ierr)
#endif

  src%nwire = nwire


  !vectors for number of wire elements and wire element lengths
  allocate(src%nelem(nwire),src%wire(nwire), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'read_wiresource','nelem for wire source',ierr)

  !temp positions and element lengths
  allocate(postmp(7,nwire), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'read_wiresource','postmp',ierr)

  !read wire specifications
  if (pid .eq. 0) then
    do iwire=1,nwire
      read(lu,*,iostat=ierr) postmp(:,iwire)
      if (ierr.ne.0) call io_error(pid,filename,'position of wire source',ierr)

      !for now, we can only handle horizontal wires! Relax this later!
      if (postmp(3,iwire) .ne. postmp(6,iwire)) then
        call invalid_error(pid,'read_wiresource',filename,'wirepos',realnum=postmp(6,iwire))
      endif

      !get number of wire elements
      x0 = postmp(1,iwire)
      y0 = postmp(2,iwire)
      z0 = postmp(3,iwire)
      x1 = postmp(4,iwire)
      y1 = postmp(5,iwire)
      z1 = postmp(6,iwire)
      lw = sqrt((x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2)

      !number of wire elements: exclude the possibility that the number of wire elements becomes zero
      src%nelem(iwire) = max(1,nint(lw/postmp(7,iwire)))

      !actual wire element length
      postmp(7,iwire) = lw / real(src%nelem(iwire))
    enddo
  endif  !pid is 0

#ifdef USE_MPI
  call MPI_Bcast(postmp,7*nwire,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_wiresource','MPI_Bcast wire positions',ierr)
  call MPI_Bcast(src%nelem,nwire,MPI_INTEGER,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_wiresource','MPI_Bcast nr of wire elements',ierr)
#endif


  do iwire = 1,nwire
    !put coordinates in src structure 
    src%wire(iwire)%endpos(:,1) = postmp(1:3,iwire)
    src%wire(iwire)%endpos(:,2) = postmp(4:6,iwire)
    src%wire(iwire)%dlw         = postmp(7,iwire)
    src%wire(iwire)%nelem       = src%nelem(iwire)
    nelem = src%wire(iwire)%nelem

    !array for positions of wire elements
    allocate(src%wire(iwire)%elempos(nelem,3), stat=ierr)
    if (ierr.ne.0) call alloc_error(pid,'read_wiresource','elempos',ierr)

    !x, y wire lengths
    lx = src%wire(iwire)%endpos(1,2) - src%wire(iwire)%endpos(1,1)
    ly = src%wire(iwire)%endpos(2,2) - src%wire(iwire)%endpos(2,1)

    !x, y lengths of wire elements
    dlwx = lx / real(nelem, kind=real64)
    dlwy = ly / real(nelem, kind=real64)

    !temp start coordinates of wire: 1/2 element length from actual start point to simplify the calculations
    x0 = src%wire(iwire)%endpos(1,1) - dlwx/2._real64
    y0 = src%wire(iwire)%endpos(2,1) - dlwy/2._real64

    !x,y, positions of all wire elements
    do ielem = 1,nelem
      src%wire(iwire)%elempos(ielem,1) = x0 + real(ielem)*dlwx
      src%wire(iwire)%elempos(ielem,2) = y0 + real(ielem)*dlwy
    enddo

    !equal depths for all wire elements
    src%wire(iwire)%elempos(:,3) = src%wire(iwire)%endpos(3,1)

  enddo

  !this is used when looping over sources and computing length of output files
  src%ncur = 1

  deallocate(postmp, stat=ierr)

  !offset of x coordinate from line x=0, for 2.5D modeling
  src%shiftx = 0._real64
  do iwire = 1,nwire
    src%shiftx = src%shiftx + sum(src%wire(iwire)%endpos(1,:))
  enddo
  src%shiftx = src%shiftx / (2*nwire)
      
endsubroutine read_wiresource


!**********************************************************************
!  FD EM subroutine read_starsource
!
!  Purpose:  read specifications of a "star" source:
!    source consists of several wire segments with a common center point
!    and currents phase-shifted by specific angles such that currents add up to zero
!    at the center point
!
!  Rita Streich 2009-2011
!**********************************************************************
subroutine read_starsource(src,lu,comm)

  implicit none

  !external variables
  type(sorec)           :: src    !the source
  integer(kind=int32)   :: lu     !file unit number
  integer(kind=int32),intent(in)        :: comm   !MPI communicator (not required for serial I/O)

  !internal variables
  integer(kind=int32)   :: ierr   !error index
  character(len=namlen) :: filename
  integer(kind=int32)   :: nwire  !number of wires
  integer(kind=int32)   :: iwire  !counter for wires
  real(kind=real64)     :: x0,y0,z0,x1,y1,z1 !wire start and end point
  real(kind=real64)     :: dlwx,dlwy    !x,y lengths of wire elements
  real(kind=real64)     :: lw,lx,ly     !wire length, total and x,y components
  integer(kind=int32)   :: ncur   !number of source currents
  integer(kind=int32)   :: icur,iicur   !counter for source currents
  integer(kind=int32)   :: wavidx       !switch for new wavelet: 1=yes, 0=no
  integer(kind=int32),parameter  :: isnew = 1, isold = 0 !indicators for new or old wavelet
  real(kind=real64),dimension(:,:),allocatable  :: postmp  !temp positions and input lengths of all wire elements
  integer(kind=int32)   :: ielem  !wire element counter
  integer(kind=int32)   :: nelem  !tmp number of wire elements


  !get number of wires composing this source
  if (pid .eq. 0) then
    !get file name
    inquire(lu, NAME=filename)

    read(lu,*,iostat=ierr) nwire
    if (ierr.ne.0) call io_error(pid,filename,'reading number of wires from ',ierr)

    !check number
    if (nwire.le.0) call invalid_error(pid,'read_starsource',filename,'number of wires within source ',intnum=nwire)
  endif

  !communicate nr of long wire segments
#ifdef USE_MPI
  call MPI_Bcast(nwire,1,MPI_INTEGER,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_starsource','MPI_Bcast nr of wires',ierr)
#endif

  src%nwire = nwire

  !vectors for number of wire elements and wire element lengths
  allocate(src%nelem(nwire),src%wire(nwire), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'read_starsource','nelem for star source',ierr)

  !temp positions and element lengths
  allocate(postmp(4,nwire), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'read_starsource','postmp',ierr)


  !get center point := end point of all wires
  if (pid .eq. 0) then
    read(lu,*,iostat=ierr) src%wire(1)%endpos(:,2)
    if (ierr.ne.0) call io_error(pid,filename,'center point of star source',ierr)
  endif

  !communicate center point
#ifdef USE_MPI
  call MPI_Bcast(src%wire(1)%endpos(:,2),3,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_starsource','MPI_Bcast wire center',ierr)
#endif

  !copy end point to all other wires
  do iwire=2,nwire
    src%wire(iwire)%endpos(:,2) = src%wire(1)%endpos(:,2)
  enddo


  !get wire start points and wire element lengths
  if (pid .eq. 0) then
    !remember end point
    x1 = src%wire(1)%endpos(1,2)
    y1 = src%wire(1)%endpos(2,2)
    z1 = src%wire(1)%endpos(3,2)

    do iwire=1,nwire
      !read start point and length of wire elements
      read(lu,*,iostat=ierr) postmp(:,iwire)
      if (ierr.ne.0) call io_error(pid,filename,'reading start point of wire from',ierr)

      !for now, we can only handle horizontal wires! Relax this later!
      if (postmp(3,iwire) .ne. src%wire(1)%endpos(3,2)) then
        call invalid_error(pid,'read_starsource',filename,'wirepos',realnum=postmp(3,iwire))
      endif

      !get number of wire elements
      x0 = postmp(1,iwire)
      y0 = postmp(2,iwire)
      z0 = postmp(3,iwire)
      lw = sqrt((x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2)

      !make sure that there is at least one element
      src%nelem(iwire) = max(1,nint(lw/postmp(4,iwire)))

      !actual wire element length, "as close as possible" to the length read from input file
      postmp(4,iwire) = lw / real(src%nelem(iwire))
    enddo !wires
  endif !pid is 0


#ifdef USE_MPI
  call MPI_Bcast(postmp,4*nwire,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_starsource','MPI_Bcast wire start positions and elem lengths',ierr)

  call MPI_Bcast(src%nelem,nwire,MPI_INTEGER,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_starsource','MPI_Bcast number of wire elements',ierr)
#endif


  do iwire = 1,nwire
    !put coordinates in src structure 
    src%wire(iwire)%endpos(:,1) = postmp(1:3,iwire)
    src%wire(iwire)%dlw         = postmp(4,iwire)
    src%wire(iwire)%nelem       = src%nelem(iwire)
    nelem = src%wire(iwire)%nelem

    !array for positions of wire elements
    allocate(src%wire(iwire)%elempos(nelem,3), stat=ierr)
    if (ierr.ne.0) call alloc_error(pid,'read_wiresource','elempos',ierr)

    !x, y wire lengths
    lx = src%wire(iwire)%endpos(1,2) - src%wire(iwire)%endpos(1,1)
    ly = src%wire(iwire)%endpos(2,2) - src%wire(iwire)%endpos(2,1)

    !x, y lengths of wire elements
    dlwx = lx / real(nelem, kind=real64)
    dlwy = ly / real(nelem, kind=real64)

    !temp start coordinates of wire: 1/2 element length from actual start point to simplify the calculations
    x0 = src%wire(iwire)%endpos(1,1) - dlwx/2._real64
    y0 = src%wire(iwire)%endpos(2,1) - dlwy/2._real64

    !x,y, positions of all wire elements
    do ielem = 1,nelem
      src%wire(iwire)%elempos(ielem,1) = x0 + real(ielem)*dlwx
      src%wire(iwire)%elempos(ielem,2) = y0 + real(ielem)*dlwy
    enddo

    !equal depths for all wire elements
    src%wire(iwire)%elempos(:,3) = src%wire(iwire)%endpos(3,1)

  enddo


  !offset of x coordinate from line x=0, for 2.5D modeling
  src%shiftx = 0._real64
  do iwire = 1,nwire
    src%shiftx = src%shiftx + sum(src%wire(iwire)%endpos(1,:))
  enddo
  src%shiftx = src%shiftx / (2*nwire)
      

  !get number of input currents
  if (pid .eq. 0) then
    read(lu,*,iostat=ierr) ncur
    if (ierr.ne.0) call io_error(pid,filename,'reading number of source currents from ',ierr)
  endif

#ifdef USE_MPI
  call MPI_Bcast(ncur,1,MPI_INTEGER,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_starsource','MPI_Bcast number of currents',ierr)
#endif

  src%ncur = ncur

  !rotation frequency and phase constant
  allocate(src%omegarot(ncur),src%fi0(ncur), src%wavnames(ncur), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'read_starsource','omegarot, fi0, wavnames', ierr)


  if (pid .eq. 0) then
    !read source current specifications
    do icur=1,ncur

      !wavelet at the modeling frequencies is not enough
      !instead, we need the full wavelet here to allow extracting I(omega+omegarot) and I(omega-omegarot)
      read(lu,*,iostat=ierr) iicur, wavidx, src%omegarot(icur), src%fi0(icur)
      if (ierr.ne.0) call io_error(pid,filename,'reading omegarot and fi0 from ',ierr)

      !in input file we define frequency f, now convert it to omega
      src%omegarot(icur) = src%omegarot(icur) * dtwopi
      !convert the phase angle to radians
      src%fi0(icur) = src%fi0(icur) * dpi / 180._real64


      !read wavelet file name if there is a new wavelet
      if (wavidx .eq. isnew) then
        read(lu,*,iostat=ierr) src%wavnames(icur)
        if (ierr.ne.0) call io_error(pid,filename,'reading name of wavelet file from ',ierr)
      else

        !force new wavelet for first current
        if (icur .eq. 1) call invalid_error(pid,'read_starsource',filename,'wavelet index for first current ',intnum=wavidx)

        !set name to none if there is no new wavelet for this current
        !  (i.e. only the rotation frequency or constant phase shift change)
        if (wavidx .eq. isold) then
          src%wavnames(icur) = 'none'

        !invalid value for wavelet index
        else
          call invalid_error(pid,'read_starsource',filename,'wavelet index ',intnum=wavidx)
        endif

      endif !new or old wavelet

    enddo !source currents
  endif !pid is 0


  !communicate omegarot and fi0
  !also communicate names of wavelet files since these are used to decide if readwavelet is entered or not
#ifdef USE_MPI
  call MPI_Bcast(src%omegarot,ncur,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_starsource','MPI_Bcast omegarot',ierr)
  call MPI_Bcast(src%fi0,ncur,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_starsource','MPI_Bcast fi0',ierr)

  call MPI_Bcast(src%wavnames,namlen*ncur,MPI_CHARACTER,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_starsource','MPI_Bcast wavnames',ierr)
#endif

  deallocate(postmp,stat=ierr)

endsubroutine read_starsource


!**********************************************************************
!  FD EM subroutine read_receivers
!
!  Purpose:  read receivers of one receiver group
!
!  Rita Streich 2009-2011
!**********************************************************************
subroutine read_receivers(rec,lu,comm)

  implicit none

  !external variables
  type(sorec)                     :: rec !the receiver group
  integer(kind=int32),intent(in)  :: lu  !file unit numbe
  integer(kind=int32),intent(in)        :: comm   !MPI communicator (not required for serial I/O)

  !internal variables
  integer(kind=int32)   :: ierr     !error index
  character(len=namlen) :: filename
  integer(kind=int32)   :: nelem    !number of receivers
  integer(kind=int32)   :: ielem    !counter for receivers


  !1-element vector for number of receivers
  allocate(rec%nelem(1), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'read_receivers','nelem for receivers',ierr)


  !get number of receivers
  if (pid .eq. 0) then
    !get file name
    inquire(lu, NAME=filename)

    read(lu,*,iostat=ierr) nelem
    if (ierr.ne.0) call io_error(pid,filename,'reading number of receivers from ',ierr)

    !check number
    if (nelem.le.0) call invalid_error(pid,'read_receivers',filename,'number of receivers ',intnum=nelem)
  endif

  !communicate nr of receivers
#ifdef USE_MPI
  call MPI_Bcast(nelem,1,MPI_INTEGER,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_receivers','MPI_Bcast nr of receivers',ierr)
#endif

  rec%nelem(1) = nelem


  !array for receiver positions
  allocate(rec%pos(3,nelem),rec%recnames(nelem), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'read_receivers','array for receiver positions',ierr)


  !read receiver positions
  if (pid .eq. 0) then
    !implicit loop is faster
    read(lu,*,iostat=ierr) (rec%pos(:,ielem),rec%recnames(ielem),ielem=1,nelem)
    if (ierr.ne.0) call io_error(pid,filename,'reading receiver positions from ',ierr)
  endif

  !communicate receiver positions and names
#ifdef USE_MPI
  call MPI_Bcast(rec%pos,3*nelem,MPI_DOUBLE_PRECISION,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_receivers','MPI_Bcast receiver positions',ierr)
  call MPI_Bcast(rec%recnames,namlen*nelem,MPI_CHARACTER,0,comm,ierr)
  if(ierr.ne.MPI_SUCCESS) call error_mpi(pid,'read_receivers','MPI_Bcast receiver names',ierr)
#endif

endsubroutine read_receivers



!**********************************************************************
!  FD EM subroutine checksourcecomp
!
!  Purpose:  check which source components are present in a dipole source
!
!  Rita Streich 2009
!**********************************************************************
subroutine checksourcecomp(src)

  implicit none

  !external variables
  type(sorec)          :: src    !the source (can have many dipole elements)

  !internal variables
  integer(kind=int32)  :: ipt    !counter
  integer(kind=int32)  :: npoint !total number of dipole elements


  npoint = src%nelem(1)


  !electric sources
  src%elsrc = .false.
  do ipt = 1,npoint
    !stop scanning source components as soon as there is one non-zero value
    if ((src%ljx(ipt).ne.0.) .or. (src%ljy(ipt).ne.0.) .or. (src%ljz(ipt).ne.0.)) then
      src%elsrc = .true.
      exit
    endif
  enddo

  !magnetic sources
  src%magsrc = .false.
  do ipt = 1,npoint
    if ((src%akx(ipt).ne.0.) .or. (src%aky(ipt).ne.0.) .or. (src%akz(ipt).ne.0.)) then
      src%magsrc = .true.
      exit
    endif
  enddo

endsubroutine checksourcecomp

