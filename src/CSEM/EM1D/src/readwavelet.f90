!**********************************************************************
!>  FD EM subroutine readwavelet
!
!>  Purpose:  read a wavelet file
!>    wavelet can be in time or frequency domain
!>    output is always in frequency domain
!
!>  Rita Streich 2009
!**********************************************************************
subroutine readwavelet(filename,wav)

  implicit none

  !external variables
  character(len=namlen),intent( in ) :: filename  !name of wavelet file
  type(wavelet) :: wav       !wavelet

  !internal variables
  integer(kind=int32) :: ierr  !error index
  integer(kind=int32) :: lu    !file unit number
  integer(kind=int32) :: dom   !domain in which wavelet is defined
  integer(kind=int32) :: nt    !number of time steps in time domain wavelet
  real(kind=real64) :: dt    !time step
  real(kind=real64) :: df    !frequency spacing
  integer(kind=int32) :: if0,if1 !start and end frequency indices
  integer(kind=int32) :: isplit  !where to split freq domain wavelet when shifting it
  integer(kind=int32) :: nf_in !number of frequency compoennts in input wavelet file
  complex(kind=real64),dimension(:),allocatable :: twav   !time domain wavelet, define as complex for input to fftw
  complex(kind=real64),dimension(:),allocatable :: wavtmp !temp freq domain wavelet for fft
  integer(kind=int32) :: it,ifreq !time and frequency counters
  character(len=10) :: itstr !string for time step (used in error message only)
  real(kind=real64) :: samp  !one wavelet sample (need this since gfortran does not accept reading into wavelet vector directly)


  if(pid .EQ. 0) then
    !open wavelet file
    lu = AvailableUnit()
    open(unit=lu,file=trim(adjustl(filename)),status='old',iostat=ierr)
    if(ierr.NE.0) call io_error(pid,filename,'opening ',ierr)

    !read domain
    read(lu,*,iostat=ierr) dom
    if(ierr.NE.0) call io_error(pid,filename,'reading wavelet domain from ',ierr)
  endif

#ifdef USE_MPI
  call MPI_Bcast(dom,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'readwavelet','MPI_Bcast domain',ierr)
#endif

  !different wavelet formats for time and frequency domains
  select case(dom)
  !time domain wavelet: read wavelet to temp vector, then transform to frequency domain
  case(timedom)

    if(pid .EQ. 0) then
      !read number of samples and time step
      read(lu,*,iostat=ierr) nt,dt
      if(ierr.NE.0) call io_error(pid,filename,'reading nt, dt from ',ierr)
      if(nt.le.0) call invalid_error(pid,'readwavelet',filename,'number of time steps ',intnum=nt)
      if(dt.le.0.d0) call invalid_error(pid,'readwavelet',filename,'time step ',realnum=dt)
    
      !allocate according to C convention and such that frequency zero is at zero'th element
      allocate(twav(0:nt-1), stat=ierr)
      if(ierr.NE.0) call alloc_error(pid,'readwavelet','time domain wavelet',ierr)

      !read the wavelet
      do it=0,nt-1
        read(lu,*,iostat=ierr) samp
        if(ierr.NE.0) then
          write(itstr,'(i9)') it
          call io_error(pid,filename,'reading wavelet sample '//trim(adjustl(itstr))//' from ',ierr)
        endif
        twav(it) = samp
      enddo
    endif !pid is 0

#ifdef USE_MPI
    call MPI_Bcast(nt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'readwavelet','MPI_Bcast nt',ierr)
    call MPI_Bcast(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'readwavelet','MPI_Bcast dt',ierr)
#endif

    !reallocate wav if it is not allocated at the same length already
    if(associated(wav%re)) then
      wav%nf = size(wav%re)
    else !wavelet not pre-allocated
      wav%nf = 0
    endif


    if(wav%nf .NE. nt) then
      if(associated(wav%re)) deallocate(wav%re,wav%im,wav%spline_re,wav%spline_im,wav%omega, stat=ierr)
      !start from negative frequencies
      if(mod(nt,2) .EQ. 1) then
        isplit = 1
        if0 = -(nt-1)/2
        if1 = (nt-1)/2
      else
        isplit = 0
        if0 = -nt/2
        if1 = nt/2-1
      endif

      allocate(wav%re(if0:if1),wav%im(if0:if1),wav%spline_re(if0:if1),wav%spline_im(if0:if1), &
               wav%omega(if0:if1), stat=ierr)
      if(ierr.NE.0) call alloc_error(pid,'readwavelet','freq domain wavelet',ierr)
      wav%nf = nt

    endif

    !frequency spacing
    df = 1.d0 / (real(nt)*dt)
    wav%domega = df * dtwopi

    !transform the wavelet to frequency domain
    if(pid .EQ. 0) then
      allocate(wavtmp(if0:if1), stat=ierr)
      if(ierr.NE.0) call alloc_error(pid,'readwavelet','wavtmp',ierr)

      !transform wavelet to frequency domain
      call transform_wav(wav%nf, twav, wavtmp)
      
      deallocate(twav, stat=ierr)

      !shift such that zero freq. is in the middle
      wavtmp = (/wavtmp(isplit:),wavtmp(:isplit-1)/)
      wav%re = real(wavtmp)
      wav%im = aimag(wavtmp)

      deallocate(wavtmp, stat=ierr)  
    endif

    !communicate wavelet
#ifdef USE_MPI
    call MPI_Bcast(wav%re,nt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'readwavelet','MPI_Bcast wav real part',ierr)
    call MPI_Bcast(wav%im,nt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'readwavelet','MPI_Bcast wav imag part',ierr)
#endif
  
  !frequency domain wavelet: directly read wavelet into wavelet vector
  case(freqdom)

    if(pid .EQ. 0)   then
      !read number of frequency samples, start frequency and frequency step
      read(lu,*,iostat=ierr) nf_in,df
      if(ierr.NE.0) call io_error(pid,filename,'reading nf_in, df from ',ierr)
    endif
    wav%domega = df * dtwopi

#ifdef USE_MPI
    call MPI_Bcast(nf_in,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'readwavelet','MPI_Bcast nf_in',ierr)
    call MPI_Bcast(wav%domega,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'readwavelet','MPI_Bcast domega',ierr)
#endif

    !check if wav was previously allocated
    if(.NOT. associated(wav%re)) then !wavelet not pre-allocated
      wav%nf = 0
    !else
    !otherwise wav%nf already matches the size of the previous wavelet
    endif


    !reallocate wav if it is not allocated at the right length already
    !assume that:
    !>  input wavelet goes from zero to fmax=(nf/2)*df
    !>  output wavelet goes from -(nf/2)*df to (nf/2-1)*df
    if(wav%nf .NE. (2*(nf_in-1))) then
      if(wav%nf .ge. 1) deallocate(wav%re,wav%im,wav%spline_re,wav%spline_im,wav%omega, stat=ierr)
      if0 = -nf_in + 1
      if1 = nf_in - 2
      allocate(wav%re(if0:if1),wav%im(if0:if1),wav%spline_re(if0:if1),wav%spline_im(if0:if1),wav%omega(if0:if1), stat=ierr)
      if(ierr.NE.0) call alloc_error(pid,'readwavelet','freq domain wavelet',ierr)
      wav%nf = if1 - if0 + 1
    endif


    if(pid .EQ. 0) then
      !read wavelet
      !> frequency components 0 to nf/2-1
      do ifreq=0,nf_in-2
        read(lu,*,iostat=ierr) wav%re(ifreq),wav%im(ifreq)
        if(ierr.NE.0) call io_error(pid,filename,'reading a wavelet sample from ',ierr)
      enddo
      !> frequency component nf/2
      read(lu,*,iostat=ierr) wav%re(if0),wav%im(if0)
      if(ierr.NE.0) call io_error(pid,filename,'reading a wavelet sample from ',ierr)

    endif

#ifdef USE_MPI
    !communicate wavelet points 0 to nf/2-1
    call MPI_Bcast(wav%re(0:),nf_in-1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'readwavelet','MPI_Bcast wav real part',ierr)
    call MPI_Bcast(wav%im(0:),nf_in-1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'readwavelet','MPI_Bcast wav imag part',ierr)

    !communicate wavelet point -nf/2
    call MPI_Bcast(wav%re(if0),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'readwavelet','MPI_Bcast wav real part2',ierr)
    call MPI_Bcast(wav%im(if0),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'readwavelet','MPI_Bcast wav imag part2',ierr)
#endif

    !complex conjugate at negative frequencies
    wav%re(-nf_in+2:-1) = wav%re(nf_in-2:1:-1)
    wav%im(-nf_in+2:-1) = -wav%im(nf_in-2:1:-1)


  !invalid domain
  case default
    call invalid_error(pid,'readwavelet',filename,'domain ',intnum=dom)
  end select

  if(pid .EQ. 0) then
    !close the wavelet file
    close(lu, iostat=ierr)
    if(ierr.NE.0) call io_error(pid,filename,'closing',ierr)
  endif !pid is 0


  !populate frequency vector
  if0 = lbound(wav%omega, dim=1)
  if1 = ubound(wav%omega, dim=1)
  do ifreq=if0,if1
    wav%omega(ifreq) = real(ifreq)*wav%domega
  enddo


  !spline derivatives
  !CAUTION: we don't want spline interpolation for pseudo random binary sequences!
  call spline(wav%omega,wav%re,wav%nf,spl_endval,spl_endval,wav%spline_re)
  call spline(wav%omega,wav%im,wav%nf,spl_endval,spl_endval,wav%spline_im)


endsubroutine readwavelet


!**********************************************************************
!>  FD EM subroutine wrap_omega
!
!>  Purpose:  wrap a frequency into a specified interval
!
!>  Rita Streich 2009
!**********************************************************************
subroutine wrap_omega(omega,omegamin,omegamax,omegaperiod)

  implicit none

  !external variables
  real(kind=real64) :: omega  !the frequency to wrap
  real(kind=real64),intent( in ) :: omegamin,omegamax,omegaperiod  !interval to wrap frequency into


  !check upper bound of interval
  do
    if(omega .le. omegamax) exit
    omega = omega - omegaperiod
  enddo
  !check lower bound of interval
  !(always check both since omega may be negative)
  do
    if(omega .ge. omegamin) exit
    omega = omega + omegaperiod
  enddo

endsubroutine wrap_omega


!**********************************************************************
!>  FD EM subroutine transform_wav
!
!>  Purpose:  fft of input wavelet
!>  separate routine to allow for different calling conventions for MPI and fftw
!>   - not really necessary to separate this but it does no harm ...
!
!>  Rita Streich 2010
!**********************************************************************
subroutine transform_wav(nf,wavin,wavout)
  implicit none
  include 'fftw3.f'
#ifdef WINDOWS
#ifndef X64
!DEC$ ATTRIBUTES ALIAS:'_dfftw_plan_dft_1d' :: _dfftw_plan_dft_1d_
!DEC$ ATTRIBUTES ALIAS:'_dfftw_execute' :: _dfftw_execute_
!DEC$ ATTRIBUTES ALIAS:'_dfftw_destroy_plan' :: _dfftw_destroy_plan_
#else
!DEC$ ATTRIBUTES ALIAS:'dfftw_plan_dft_1d_' :: dfftw_plan_dft_1d
!DEC$ ATTRIBUTES ALIAS:'dfftw_execute_' :: dfftw_execute
!DEC$ ATTRIBUTES ALIAS:'dfftw_destroy_plan_' :: dfftw_destroy_plan
#endif
#endif

  !external variables
  integer(kind=int32) :: nf  !nr of freq. samples
  complex(kind=real64),dimension(:) :: wavin   !time domain wavelet, define as complex for input to fftw
  complex(kind=real64),dimension(:) :: wavout !temp freq domain wavelet for fft

  !internal variables
  integer(kind=int64) :: trafoplan !"plan" for fftw initialization
  
  
  call dfftw_plan_dft_1d(trafoplan,nf,wavin,wavout, FFTW_FORWARD, FFTW_ESTIMATE)
  call dfftw_execute(trafoplan, wavin, wavout)
  call dfftw_destroy_plan(trafoplan)

endsubroutine transform_wav
