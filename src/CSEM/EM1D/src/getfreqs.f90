!**********************************************************************
!  FD EM subroutine getfreqs
!
!  Purpose:  get modeling frequencies
!    for regular frequency spacing (fdef option 1): from jobfile parameters
!    for (possibly) irregular frequency spacing (fdef option 2): read from separate ascii file
!
!  Rita Streich 2009
!**********************************************************************
!!$subroutine getfreqs(upars,freqdat)
subroutine getfreqs(fdef,nf,fstart,df,fname,freqdat)

  implicit none

  !external variables
!!$  type(userpars)       :: upars    !user input parameters
  integer(kind=int32)  :: fdef     !frequency definition mode: 1=regularly spaced, 2=read from file, may be irregular
  integer(kind=int32)  :: nf       !number of frequencies from general input file, may be replaced by nf from freq. file
  real(kind=real64)    :: fstart   !start frequency from general parameter file
  real(kind=real64)    :: df       !regular frequency spacing from general parameter file
  character(len=namlen)  :: fname  !frequency file name

  type(freqdata)       :: freqdat  !frequencies

  !internal variables
  integer(kind=int32)  :: ierr  !error index
  integer(kind=int32)  :: ifreq !frequency counter
  integer(kind=int32)  :: lu    !file unit number


  !equally spaced frequencies as defined in jobfile
  if (fdef.eq.1) then
    !use parameters from input file

    freqdat%nfreq = nf

    !allocate frequency vector
    allocate(freqdat%omega(nf), stat=ierr)
    if (ierr.ne.0) call alloc_error(pid,'getfreqs','frequency vector',ierr)

    !populate frequency vector
    do ifreq=1,freqdat%nfreq
      freqdat%omega(ifreq) = fstart+(ifreq-1)*df
    enddo

  elseif (fdef.eq.2) then
    !get all info from separate file
    lu = AvailableUnit()
    open(unit=lu,file=trim(adjustl(fname)),status='old',iostat=ierr)
    if (ierr.ne.0) call open_error(pid,'getfreqs',fname,ierr)

    !read number of frequencies
    read(lu,*,iostat=ierr) freqdat%nfreq
    if (ierr.ne.0) call readwrite_error(pid,'getfreqs',fname,'r',ierr,1_int64)
    nf = freqdat%nfreq

    !allocate frequency vector
    allocate(freqdat%omega(freqdat%nfreq), stat=ierr)
    if (ierr.ne.0) call alloc_error(pid,'getfreqs','frequency vector',ierr)

    !read frequencies from file
    do ifreq=1,freqdat%nfreq
      read(lu,*,iostat=ierr) freqdat%omega(ifreq)
      if (ierr.ne.0) call readwrite_error(pid,'getfreqs',fname,'r',ierr,ifreq+1_int64)
    enddo

    !close file
    close(lu,iostat=ierr)
    if (ierr.ne.0) call close_error(pid,'getfreqs',fname,ierr)

  endif

  !convert to angular frequency
  freqdat%omega = freqdat%omega * dtwopi

endsubroutine getfreqs
