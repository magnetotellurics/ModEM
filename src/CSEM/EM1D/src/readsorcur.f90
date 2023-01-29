!**********************************************************************
!>  FD EM subroutine readsorcur
!
!>  Purpose:  read source currents
!>    can be simple complex numbers for the modeling frequencies
!>    or time or frequency domain wavelets
!
!>  Rita Streich 2009
!**********************************************************************

subroutine readsorcur(filename,freqdat,sources,comm)
 use, intrinsic :: iso_fortran_env
  implicit none

  !external variables
  character(len=namlen),intent( in ) :: filename  !file for definition of source currents
  type(freqdata),intent(inout) :: freqdat   !frequency dependent specifications
  type(sorec),dimension(:),pointer :: sources   !source specifications
  integer(kind=int32),intent( in ) :: comm      !MPI communicator

  !internal variables
  character(len=10) :: flenstr    !string for file size, used when error occurs
  character(len=10) :: freqstr    !string for number frequencies, used when error occurs
  character(len=10) :: nsrcstr    !string for number source points, used when error occurs
  integer(kind=int32) :: ierr       !error index
  integer(kind=int32) :: lu         !file unit number
  integer(kind=int32) :: ifreq,icur !counters
  integer(kind=int32) :: iicur    !counter for currents within one "star" source
  integer(kind=int32) :: nsrc       !number of sources
  integer(kind=int32) :: ncurgroup  !number of "source current groups", either one or nr of source elements / wires
  integer(kind=int32) :: isrc       !source counter
  complex(kind=real32),dimension(:,:),allocatable :: curtmp !temp currents from simple binary file
  integer(kind=int32) :: nrec       !number of records in source current file
  integer(kind=int32) :: nrecmax    !max number of records in source current file
  integer(kind=int32) :: irec       !file record counter
  integer(kind=int32) :: iwire,ielem       !wire/dipole element counter
  real(kind=real64) :: fi,dfi     !phase angles for "star" source
  integer(kind=int32) :: recinc     !record increment
  type(wavelet) :: wav
  complex(kind=real64) :: Iplus,Iminus      !current at frequencies (omega+-omegarot)
  real(kind=real64) :: Ire,Iim    !tmp real and imag parts of current
  real(kind=real64) :: omegarot   !temp copy of rotation frequency
  real(kind=real64) :: omegatmp   !temp angular frequency
  real(kind=real64) :: omegamin,omegamax !min and max angular frequency in input wavelet
  real(kind=real64) :: omegaperiod       !frequency period for input wavelets
  integer(kind=int64) :: filesize   !size of source current file
  integer(kind=int32) :: cplx_size  !size in bytes of complex numbers in input file
  logical :: haswav = .FALSE.  !indicates if source wavelet files have to be read
  integer(kind=int32) :: nfreq      !number of modeling frequencies (extract from model structure)
  type(sorec),pointer :: src        !pointer to one source
#ifdef USE_MPI
  integer(kind=int32) :: curtype    !mpi vector type for currents of 1 source element and all frequencies
  integer(kind=MPI_OFFSET_KIND) :: disp      !displacement into file
#else
  integer(kind=int64) :: iolen     !record length when not using MPI
#endif


#ifdef USE_MPI
  !get size of complex MPI datatype
  call MPI_Type_size(MPI_COMPLEX,cplx_size,ierr)
  if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'getinput','MPI_Type_size',ierr)
#else
  cplx_size = 8
#endif

  !number of sources
  nsrc = size(sources)

  !number of modeling frequencies
  nfreq = freqdat%nfreq

  !-----------------------------------------------------------
  !get maximum number of input source currents for all sources
  !-----------------------------------------------------------
  !> and max nuber of source current vectors for computations
  !for dipole element sources and "simple wire" sources: source currents are defined in sor_cur file 
  !> by specifying complex current values, either for each frequency or for each frequency and current element / wire
  !for "star sources" with rotating source current permitted, we need a full wavelet, 
  !> since the actual current at frequency omega is computed from input current at frequencies
  !> omega+omegarot and omega-omegarot
  !> we also need to translate the input current to currents differing by specific phase shifts for the different wires
  nrecmax = 0
  ncurgroup = 0
  nrec = 0
  do isrc=1,nsrc
    src => sources(isrc)

    select case(src%type)
    case(dipole)
      nrecmax = nrecmax + src%nelem(1)
      ncurgroup = ncurgroup + src%nelem(1)
    case(wire)
      nrecmax = nrecmax + src%nwire
      ncurgroup = ncurgroup + src%nwire
    case(star)
      !nothing in sor_cur file for this source, source current is read from separate wavelet file(s)
      !nrecmax = nrecmax + 0
      !there will be different currents for each input current and each wire
      ncurgroup = ncurgroup + src%nwire * src%ncur
      haswav = .TRUE.
    end select
  enddo


  !make sure we don't have wild pointers - gfortran is sensitive to this
  nullify(wav%re,wav%im,wav%spline_re,wav%spline_im,wav%omega)


  !read simple binary current file if there are any dipole or simple wire sources
  !nothing to do here if we only have "star" sources with separate wavelet files
  if(nrecmax .gt. 0) then

    !temp matrix for currents - this should be the only place where currents are not split by sources!
    allocate(curtmp(nrecmax,nfreq),stat=ierr)
    if(ierr.NE.0) call alloc_error(pid,'readsorcur','temp source current vector',ierr)


    !get file size
#ifdef USE_MPI
    call MPI_File_open(comm,trim(adjustl(filename)),MPI_MODE_RDONLY,MPI_INFO_NULL,lu,ierr)
    if(ierr .NE. MPI_SUCCESS) call error_mpi(pid,'readsorcur','MPI_File_open',ierr)

    call MPI_File_get_size(lu,filesize,ierr)
    if(ierr .NE. MPI_SUCCESS) call error_mpi(pid,'readsorcur','MPI_File_get_size',ierr)
#else
    filesize = getfilesize(filename)
#endif

    !number of complex numbers in file has to be either = number of frequencies or =nfreq*nrecmax
    filesize = filesize / cplx_size
    nrec = filesize/nfreq

    !check the file size:
    !we can either have the same source currents for all source points, in this case the 
    !>  file size has to be 8bytes times the number of frequencies
    !or we can have different source currents for each source point, in this case the
    !>  file size has to be 8bytes times the number of frequencies times the number of source elements / wires (nrecmax)
    if((filesize.NE.nfreq) .AND. (filesize.NE.(nfreq*nrecmax))) then
      write(unit=flenstr,fmt='(i9)') filesize
      write(unit=freqstr,fmt='(i6)') nfreq
      write(unit=nsrcstr,fmt='(i8)') nfreq*nrecmax
      write(*,'(a)') 'ERROR in readsorcur:'
      write(*,'(a)') 'Source current file '//trim(adjustl(filename))//' has '//trim(adjustl(flenstr))// &
        ' entries. This matches neither the number of frequencies ('// &
        trim(adjustl(freqstr))//') nor the number of frequencies times number of dipole elements / wires ('// &
        trim(adjustl(nsrcstr))//').'
      write(*,fmt='(a)') 'Exiting.'

#ifdef USE_MPI
      call MPI_Finalize(ierr)
#endif
      stop
    endif

#ifdef USE_MPI
    !if file is ok, read it (     while sorting into curtmp)
    call MPI_Type_vector(nfreq,1,nrecmax,MPI_COMPLEX,curtype,ierr)
    if(ierr .NE. MPI_SUCCESS) call error_mpi(pid,'readsorcur','MPI_Type_vector',ierr)
    call MPI_Type_commit(curtype,ierr)
    if(ierr .NE. MPI_SUCCESS) call error_mpi(pid,'readsorcur','MPI_Type_commit',ierr)

    disp = 0
    call MPI_File_set_view(lu,disp,MPI_COMPLEX,MPI_COMPLEX,'native',MPI_INFO_NULL,ierr)
    if(ierr .NE. MPI_SUCCESS) call error_mpi(pid,'readsorcur','MPI_File_set_view',ierr)

    do irec=1,nrec
      call MPI_File_read(lu,curtmp(irec,1),1,curtype,MPI_STATUS_IGNORE,ierr)
      if(ierr .NE. MPI_SUCCESS) call error_mpi(pid,'readsorcur','MPI_File_read',ierr)
    enddo

    call MPI_File_close(lu,ierr)
    if(ierr .NE. MPI_SUCCESS) call error_mpi(pid,'readsorcur','MPI_File_close',ierr)

    call MPI_Type_free(curtype,ierr)
    if(ierr .NE. MPI_SUCCESS) call error_mpi(pid,'readsorcur','MPI_Type_free',ierr)
#else    

    lu = AvailableUnit()
    !order of entries in source current file: fast index is frequencies, slow index (if present) is sources
    !order is transposed for the computations!
    inquire(iolength=iolen) curtmp(1,:)
    open(unit=lu,file=trim(adjustl(filename)),form='unformatted',access='direct', &
        status='old',recl=iolen,iostat=ierr)
    if(ierr.NE.0) call io_error(pid,filename,'opening ',ierr)

    !read data into temp array
    do irec=1,nrec
      read(unit=lu,rec=irec,iostat=ierr) (curtmp(irec,ifreq),ifreq=1,nfreq)
      if(ierr.NE.0) call io_error(pid,filename,'reading ',ierr)
    enddo

    !close source file
    close(unit=lu, iostat=ierr)
#endif

  endif !dipole or simple wire sources present, data needed from simple binary source current file


  !the same current used for all sources
  if(nrec .EQ. 1) then
    recinc = 0
    
  !different currents for each source element / wire
  else
    recinc = 1    
  endif


  !copy currents from temp matrix to src structure
  irec = 1
  sourceloop: do isrc=1,nsrc
    src => sources(isrc)

    select case(src%type)
    !dipole sources
    case(dipole)

      !allocate matrix for source currents
      allocate(src%cur(src%nelem(1),nfreq),stat=ierr)
      if(ierr.NE.0) call alloc_error(pid,'readsorcur','source current matrix',ierr)

      do ielem=1,src%nelem(1)
        !CAUTION: Loeseth's sign convention is opposite to the one I want to use!!!
        !but do NOT take complex conjugate here since we may need original currents too!!!
!!$        src%cur(ielem,:) = conjg(curtmp(irec,:))    !this was used in 1D code!!!
        src%cur(ielem,:) = curtmp(irec,:)   !this was used in 3D and 2,5D code
        irec = irec + recinc
      enddo        

    !wire sources
    case(wire)

      !allocate matrix for source currents
      allocate(src%cur(src%nwire,nfreq),stat=ierr)
      if(ierr.NE.0) call alloc_error(pid,'readsorcur','source current matrix',ierr)

      do iwire=1,src%nwire
        !CAUTION: Loeseth's sign convention is opposite to the one I want to use!!!
        !but do NOT take complex conjugate here since we may need original currents too!!!
!!$        src%cur(iwire,:) = conjg(curtmp(irec,:))    !this was used in 1D code!!!
        src%cur(iwire,:) = curtmp(irec,:)   !this was used in 3D and 2,5D code
        irec = irec + recinc
      enddo


    !special wire sources
    case(star)
      
      !allocate matrix for source currents
      allocate(src%cur(src%ncur*src%nwire,nfreq),stat=ierr)
      if(ierr.NE.0) call alloc_error(pid,'readsorcur','source current matrix',ierr)

      !phase angle separation between the wires
      dfi = (360._real64 / real(src%nwire,kind=real64)) * (dpi/180._real64)

      iicur = 1
      currents: do icur = 1,src%ncur

        if(src%wavnames(icur) .NE. 'none' ) then
          !output from readwavelet is a frequency domain wavelet,
          !>  vector length nf
          !>  angular frequencies specified in wav%omega
          !>  order from negative over zero to positive
          !frequency range:
          !>  for even nr of samples: (-nf/2)*df to (nf/2-1)*df
          !>  for odd nr of samples: (-nf+1)/2*df to (nf-1)/2*df
          call readwavelet(src%wavnames(icur),wav)
          !remember min and max frequency
          omegamin = wav%omega(lbound(wav%omega,1))
          omegamax = wav%omega(ubound(wav%omega,1))
          omegaperiod = wav%nf * wav%domega
        endif
        
        omegarot = src%omegarot(icur)
      
        wires: do iwire=1,src%nwire
        
          fi = real(iwire-1,kind=real64) * dfi + src%fi0(icur)
          

          !for non-zero rotation frequency, need a full input wavelet over a wider frequency range 
          !>  to allow extracting currents at frequencies omega+omegarot and omega-omegarot
          freqloop: do ifreq=1,nfreq

            omegatmp = freqdat%omega(ifreq) + omegarot

            !wrap omeagatmp into frequency interval [omegamin,omegamax]
            !> need to check if this approach is valid!!!
            call wrap_omega(omegatmp,omegamin,omegamax,omegaperiod)

            call splint(wav%omega,wav%re,wav%spline_re,wav%nf,omegatmp,Ire)
            call splint(wav%omega,wav%im,wav%spline_im,wav%nf,omegatmp,Iim)
            Iplus = complex(Ire,Iim)


            omegatmp = freqdat%omega(ifreq) - omegarot

            !wrap omeagatmp into frequency interval [omegamin,omegamax]
            !> need to check if this approach is valid!!!
            call wrap_omega(omegatmp,omegamin,omegamax,omegaperiod)

            call splint(wav%omega,wav%re,wav%spline_re,wav%nf,omegatmp,Ire)
            call splint(wav%omega,wav%im,wav%spline_im,wav%nf,omegatmp,Iim)
            Iminus = complex(Ire,Iim)

            !need complex conjugate to match Loeseth's sign convention
            !but: do NOT take complex conjugate here
            !plus and minus in exponent from fftw definition of Fourier transform
!!$            src%cur(iicur,ifreq) = conjg(Iplus * exp(-dci*fi) / 2._real64 + Iminus * exp(dci*fi) / 2._real64) !this was used in 1D code
            src%cur(iicur,ifreq) = Iplus * exp(-dci*fi) / 2._real64 + Iminus * exp(dci*fi) / 2._real64 !this was used in 3D and 2.5D code
          enddo freqloop
          
          !increment iicur here to get the current for each wire
          iicur = iicur + 1
        enddo wires

      enddo currents
      
      !from now we treat the source like a "normal" wire source???
      !> i.e. for both "star" and general wire sources we keep the contributions from different wires separate,
      !> only add them up at the very end!
      deallocate(src%omegarot,src%fi0,src%wavnames, stat=ierr)
      src%type = wire

    end select

  enddo sourceloop


  if(nrecmax.gt.0) deallocate(curtmp, stat=ierr)
  if(haswav) deallocate(wav%re,wav%im,wav%spline_re,wav%spline_im,wav%omega, stat=ierr)

endsubroutine readsorcur

