!**********************************************************************
!  1D EM subroutine get_outputlength
!
!  Purpose:  get total number of receivers, added up over all sources
!
!  Rita Streich 2009
!**********************************************************************
integer(kind=int32) function get_outputlength(sources,receivers) result(nrecout)

  implicit none

  !external variables
  type(sorec),dimension(:),intent(in)    :: sources     !source definitions
  type(sorec),dimension(:),intent(in)    :: receivers   !receiver definitions

  !internal variables
  integer(kind=int32)       :: nsrc,nrec   !number of source and receiver groups
  integer(kind=int32)       :: recinc      !receiver increment
  integer(kind=int32)       :: isrc,irec   !source and receiver counters
  integer(kind=MPI_OFFSET_KIND)  :: nrec_out_tot  !total number of receivers, added up over all sources

  nsrc = size(sources)
  nrec = size(receivers)

  irec = 1
  if(nrec .eq. 1) then
    recinc = 0
  else
    recinc = 1
  endif

  nrec_out_tot = 0

  do isrc=1,nsrc
    !requiring that the number of source currents has been correctly set to 1 for "dipole" and "wire" sources"
    !for "star" sources, ncur has been read from source file
    nrec_out_tot = nrec_out_tot + sources(isrc)%ncur * receivers(irec)%nelem(1)
    irec = irec + recinc
  enddo

  nrecout = nrec_out_tot

endfunction get_outputlength


