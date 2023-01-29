!EM1d get timer, hide this because of #ifdef ....
!RS 2010
  
real(kind=real64) function get_time()
 use, intrinsic :: iso_fortran_env
  implicit none

  integer :: count, crate

#ifdef USE_MPI
  get_time = MPI_Wtime()
#else
!!!#ifndef X64
!!!>  call cpu_time(get_time)
!!!#else !cpu_time produced wrong results for win64 release version
  call SYSTEM_CLOCK(count, crate)
  get_time = real(count, kind=real64) / real(crate,kind=real64)
!!!#endif
#endif

endfunction get_time

