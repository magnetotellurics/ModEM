!****************************************************************
!
!> EM modeling for counting parallel frequencies, keeping track of
!>   how many have been done
!
!> Rita Streich 2011
!**************************************************************
module counter_mod

  use machine_dep_mod
  use util_mod
  implicit none

#ifdef USE_MPI
  public:: counter_service, counter_finish
#endif
  
  integer(kind=int32) :: countertag = 333    !tag for communicating and counting ifreq values
  integer(kind=int32) :: count_ifreqall      !overall frequency counter
  integer(kind=int32) :: count_ndone         !number of processes or process groups (except proc0) that are done with frequency loop
  integer(kind=int32) :: count_max           !maximum counter value, store here to make it accessible from places
                                              !where original max value is not present (e.g. nfreq in reflectivity routine)

#ifdef USE_MPI
contains
#include "counter.f90"
#endif

endmodule counter_mod
